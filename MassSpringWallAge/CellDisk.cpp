//
// This file is part of MorphoDynamX - http://www.MorphoDynamX.org
// Copyright (C) 2012-2020 Richard S. Smith and collaborators.
//
// If you use MorphoDynamX in your work, please cite:
//   http://dx.doi.org/10.7554/eLife.05864
//
// MorphoDynamX is free software, and is licensed under under the terms of the 
// GNU General (GPL) Public License version 2.0, http://www.gnu.org/licenses.
//
#include <CellDisk.hpp>

namespace CellDisk
{
  bool CellDisk::initialize(QWidget *parent)
  {
    // Get the current mesh
    mesh = currentMesh();
    if(!mesh)
      throw(QString("CellDisk::initialize No current mesh"));

    // Get all the processes required

    // Cell Division
    if(!parm("Divide Process").isEmpty()) {
      if(!getProcess(parm("Divide Process"), divideProcess))
        throw QString("%1::initialize Cannot make divide cell process").arg(name());
      divideProcess->initialize(parent);
    } else
      divideProcess = 0;

    // Growth Process
    if(!getProcess(parm("Growth Process"), growthProcess))
      throw QString("%1::initialize Cannot make growth process").arg(name());
    growthProcess->initialize(parent);

    // Tissue process, growth process initializes
    if(!getProcess(parm("Tissue Process"), tissueProcess))
      throw QString("%1::initialize Cannot make tissue process").arg(name());

    // Solver process
    if(!getProcess(parm("Solver Process"), solverProcess))
      throw QString("%1::initialize Cannot make mass-spring solver process").arg(name());
    solverProcess->initialize(parent);

    // Get the split edges process
    if(!parm("Split Edges Process").isEmpty()) {
      if(!getProcess(parm("Split Edges Process"), splitEdgesProcess))
        throw QString("%1::initialize Cannot make splitEdges process").arg(name());
    } else
      splitEdgesProcess = 0;

    // Get the quantify angles process
    if(!parm("Quantify Angles Process").isEmpty()) {
      if(!getProcess(parm("Quantify Angles Process"), quantifyAnglesProcess))
        throw QString("%1::initialize() Cannot make quantify angles process").arg(name());
    } else
      quantifyAnglesProcess = 0;

    return true;
  }

  bool CellDisk::step()
  {
    // Mass spring, return if not converged
    solverProcess->solve();
    double xNorm = solverProcess->calcXNorm();
    double dxNorm = solverProcess->calcDXNorm();
    mesh->updatePositions(tissueProcess->tissue().tissueName());
    if(dxNorm < parm("Converge Threshold").toDouble())
      mdxInfo << "Converged, X Norm:" << xNorm << " Dx Norm:" << dxNorm << endl;
    else
      return true;

    if(quantifyAnglesProcess)
      quantifyAnglesProcess->run();

    // Remove some cells
    double cellKill = parm("Cell Kill").toDouble();
    if(cellKill >  0) {
      auto &cs = tissueProcess->tissue().cellStructure();
      auto &indexAttr = mesh->indexAttr();
      CCIndexVec toKill;
      for(CCIndex c : cs.faces())
        if(norm(indexAttr[c].pos) > cellKill)
          toKill.push_back(c);
      for(CCIndex c : toKill)
        tissueProcess->tissue().removeCell(c);
    }

    // Split large-enough cells
    if(divideProcess)
      divideProcess->step();

    // Split edges
    if(splitEdgesProcess)
     splitEdgesProcess->run();

    // Grow
    growthProcess->step();

    // If the structure has changed, re-initialize the solver 
    solverProcess->initialize(0);

    mesh->updateAll(tissueProcess->tissue().tissueName());
    mesh->updateAll(tissueProcess->tissue().tissueDualName());

    steps++;
    std::cout << "Step " << steps << std::endl;

    if(steps > parm("Stop Steps").toInt()){
      return false;
    }

    return true;
  }

  bool CellDisk::rewind(QWidget *parent)
  {
    // To rewind, we'll reload the mesh
    mesh = currentMesh();
    if(!mesh or mesh->file().isEmpty())
      throw(QString("No current mesh, cannot rewind"));
    MeshLoad meshLoad(*this);
    meshLoad.setParm("File Name", mesh->file());
    return meshLoad.run();
  }

  bool CellDiskSolver::initialize(QWidget *parent)
  {
    // Mass spring process
    MassSpringDerivs *msDerivs;
    if(!getProcess(parm("Mass Spring Process"), msDerivs))
      throw QString("%1::initialize Cannot make mass-spring derivs process %2").arg(name()).arg(parm("Mass Spring Process"));
    msDerivs->initialize(parent);

    // Pressure process
    Pressure2Derivs *p2Derivs;
    if(!getProcess(parm("Pressure Process"), p2Derivs))
      throw QString("%1::initialize Cannot make pressure derivs process %2").arg(name()).arg(parm("Pressure Process"));
    p2Derivs->initialize(parent);

    // RSS FIXME
    Mesh *mesh = currentMesh();
    msDerivs->cs = &mesh->ccStructure("Tissue");
    msDerivs->indexAttr = &mesh->indexAttr();
    p2Derivs->cs = &mesh->ccStructure("Tissue");
    p2Derivs->indexAttr = &mesh->indexAttr();

    clearDerivs();
    addDerivs(msDerivs);
    addDerivs(p2Derivs);
    initSolver(msDerivs->cs);

    // Reset the Dt after growth
    Dt = parm("Dt").toDouble();

    return true;
  }

  void MassSpringDerivs::initDerivs(Solver<Point3d> &solver, Solver<Point3d>::VertexAttr &vAttr, Solver<Point3d>::EdgeAttr &eAttr)
  {
    Mesh *mesh = currentMesh();
    if(!mesh)
      throw QString("%1::initialize No current mesh").arg(name());

    // Get parameters from GUI
    k = parm("Spring Constant").toDouble();
    ageConstant = parm("Stiffness New Spring").toDouble();

    cs = &mesh->ccStructure("Tissue"); // RSS FIXME
    indexAttr  = &mesh->indexAttr();
    msEdgeAttr = &mesh->attributes().attrMap<CCIndex, MassSpringEdgeData>("Mass Spring Data");

    // Locate border and mark it, also grab cell for orientation
    auto edges = cs->edges();

    #pragma omp parallel for
    for(uint i = 0; i < edges.size(); i++) {
      CCIndex e = edges[i];
      MassSpringEdgeData &eMS = (*msEdgeAttr)[e];

      auto cb = cs->cobounds(e);
      if(cb.size() == 1) { // On the border
        eMS.border = true;
        eMS.f = *cb.begin();
        eMS.age = 1.0; // Border edges are automatically full grown
      } else {
        eMS.border = false;
        eMS.f = CCIndex::UNDEF;
      }
    }
  }

  // compute the ageFactor at age=0 the factor is ageConstant, at age=1 the factor is 1
  double calcAgeFactor(double edgeAge, double ageConstant)
  {
    double ageFactor = 1.0;
    if(edgeAge < 1){
      ageFactor = edgeAge * 1 + ageConstant * (1 - edgeAge);
    }
    return ageFactor;
  }

  void MassSpringDerivs::calcDerivatives(const SolverT &solver, CCIndex v, Point3d &values)
  {
    if(!cs)
      throw QString("%1::initDerivs Cell complex not set").arg(name());
    if(!indexAttr)
      throw QString("%1::initDerivs Index data attribute not set").arg(name());
    if(!msEdgeAttr)
      throw QString("%1::initDerivs Edge data attribute not set").arg(name());

    CCIndexData &vIdx = (*indexAttr)[v];
    for(const Flip &flip : cs->matchV(CCIndex::BOTTOM, v, CCIndex::Q, CCIndex::Q)) {
      CCIndex e = flip.interior;
      CCIndex n = flip.otherFacet(v);
      CCIndexData &nIdx = (*indexAttr)[n];

      MassSpringEdgeData &eMS = (*msEdgeAttr)[e];
      Point3d dir =  nIdx.pos - vIdx.pos;
      double length = norm(dir);
      if(length > 0)
        dir /= length;


      values += (length - eMS.restLen)/eMS.restLen * k * calcAgeFactor(eMS.age, ageConstant) * dir;
    }
  }

  void MassSpringDerivs::calcJacobian(const SolverT &solver, VertexAttr &vAttr, EdgeAttr &eAttr)
  {
    if(!cs)
      throw QString("%1::initDerivs Cell complex not set").arg(name());
    if(!indexAttr)
      throw QString("%1::initDerivs Index data attribute not set").arg(name());
    if(!msEdgeAttr)
      throw QString("%1::initDerivs Edge data attribute not set").arg(name());

    for(CCIndex edge: cs->edges()) {
      MassSpringEdgeData &eMS = (*msEdgeAttr)[edge];

      std::pair<CCIndex,CCIndex> eps = cs->edgeBounds(edge);
      SolverVertexData &vd0 = vAttr[eps.first], &vd1 = vAttr[eps.second];

      Point3d dir = vd1.x - vd0.x;
      double len = norm(dir);
      double hooke = k * calcAgeFactor(eMS.age, ageConstant) / eMS.restLen;
      double kl = hooke * (len - eMS.restLen) / len;
      double kl3 = hooke * eMS.restLen / (len * len * len);

      // J = kl * I + kl3 * (dir x dir)
      Matrix3d mat;
      for(uint x = 0 ; x < 3 ; x++) {
        mat[x][x] += kl;
        for(uint y = 0 ; y < 3 ; y++)
          mat[x][y] += kl3 * dir[x] * dir[y];
      }

      vd0.j -= mat;
      vd1.j -= mat;

      eAttr[+edge].j += mat;
      eAttr[-edge].j += mat;
    }
  }

  void Pressure2Derivs::initDerivs(Solver<Point3d> &solver, Solver<Point3d>::VertexAttr &vAttr, Solver<Point3d>::EdgeAttr &eAttr)
  {
    Mesh *mesh = currentMesh();
    if(!mesh)
      throw(QString("Pressure2Derivs::initialize No current mesh"));

    // Get parameters from GUI
    pressure = parm("Pressure").toDouble();

    cs = &mesh->ccStructure("Tissue"); // RSS FIXME
    indexAttr  = &mesh->indexAttr();
    msEdgeAttr = &mesh->attributes().attrMap<CCIndex, MassSpringEdgeData>("Mass Spring Data");

    // Locate border and mark it, also grab cell for orientation
    auto edges = cs->edges();

    #pragma omp parallel for
    for(uint i = 0; i < edges.size(); i++) {
      CCIndex e = edges[i];
      MassSpringEdgeData &eMS = (*msEdgeAttr)[e];

      auto cb = cs->cobounds(e);
      if(cb.size() == 1) { // On the border
        eMS.border = true;
        eMS.f = *cb.begin();
      } else {
        eMS.border = false;
        eMS.f = CCIndex::UNDEF;
      }
    }
  }

  void Pressure2Derivs::calcDerivatives(const SolverT &solver, CCIndex v, VectorT &values)
  {
    if(!cs)
      throw QString("%1::initDerivs Cell complex not set").arg(name());
    if(!indexAttr)
      throw QString("%1::initDerivs Index data attribute not set").arg(name());
    if(!msEdgeAttr)
      throw QString("%1::initDerivs Edge data attribute not set").arg(name());

    if(!(pressure > 0 && cs->onBorder(v))) return;

    CCIndexData &vIdx = (*indexAttr)[v];
    for(const Flip &flip : cs->matchV(CCIndex::BOTTOM, v, CCIndex::Q, CCIndex::Q)) {
      CCIndex e = flip.interior;
      CCIndex n = flip.otherFacet(v);
      CCIndexData &nIdx = (*indexAttr)[n];

      MassSpringEdgeData &eMS = (*msEdgeAttr)[e];

      // Pressure only on the outside edges
      if(eMS.border) {
       Point3d dir =  nIdx.pos - vIdx.pos;
       double length = norm(dir);
       if(length > 0)
         dir /= length;

        // Find orientation
        int ro = cs->ro(eMS.f, e) * cs->ro(e, v);
        // Apply pressure
        values += Point3d(dir.y(), -dir.x(), 0) * (ro == ccf::POS ? 1.0 : -1.0) * pressure * length * 0.25; // 1/2 pressure applied twice
      } 
    }
  }

  void Pressure2Derivs::calcJacobian(const SolverT &solver, CCIndex v, CCIndex nb, MatrixT &j) {
    if(!cs)
      throw QString("%1::initDerivs Cell complex not set").arg(name());
    if(!indexAttr)
      throw QString("%1::initDerivs Index data attribute not set").arg(name());
    if(!msEdgeAttr)
      throw QString("%1::initDerivs Edge data attribute not set").arg(name());

    if(!(pressure > 0)) return;

    CCSignedIndex oedge = cs->orientedEdge(v,nb);
    CCIndex edge = ~oedge;
    if(edge == CCIndex::UNDEF) return;

    MassSpringEdgeData &eMS = (*msEdgeAttr)[edge];
    if(eMS.border) {
      int ro = (cs->ro(eMS.f, edge) * oedge.orientation() == ccf::POS) ? 1 : -1;

      // only [x][y] and [y][x] are nonzero
      j[0][1] += 0.5 * pressure * ro;
      j[1][0] -= 0.5 * pressure * ro;
    }
  }

  bool CellDiskGrowth::initialize(QWidget *parent)
  {
    // Get the current mesh
    mesh = currentMesh();
    if(!mesh)
      throw(QString("CellDisk::initialize No current mesh"));

    // Initialize the tissue
    if(!getProcess(parm("Tissue Process"), tissueProcess))
      throw(QString("CellDiskGrowth::initialize Cannot make tissue process"));
    tissueProcess->initialize(parent);

    return true;
  }

  bool CellDiskGrowth::step()
  {
    if(!mesh)
      throw(QString("CellDiskGrowth::step No current mesh"));

    if(!tissueProcess)
      throw(QString("CellDiskGrowth::step No tissue process"));

    CCStructure &cs = tissueProcess->tissue().cellStructure();
    auto &msEdgeAttr = mesh->attributes().attrMap<CCIndex, MassSpringEdgeData>("Mass Spring Data");
    auto &indexAttr = mesh->indexAttr();

    double growthRate = parm("Growth Rate").toDouble();
    double dt = parm("Dt").toDouble();
    double wallAging = parm("Wall Aging").toDouble();
    double initalStiffness = parm("New Wall Stiffness Factor").toDouble();

    auto edges = cs.edges();
    for(uint i = 0; i < edges.size(); i++) {
      CCIndex e = edges[i];
      MassSpringEdgeData &eMS = msEdgeAttr[e];
      if(eMS.age < 1.0) {
        eMS.age += dt / wallAging;
        if(eMS.age > 1.0)
          eMS.age = 1.0;
      }
      CCIndexPair ep = cs.edgeBounds(e);
      CCIndexData &vIdx = indexAttr[ep.first];
      CCIndexData &wIdx = indexAttr[ep.second];
      double length = norm(vIdx.pos - wIdx.pos);
      // Strain-based growth
      if(length > eMS.restLen and growthRate > 0)
        eMS.restLen *= 1.0 + (length - eMS.restLen)/eMS.restLen * growthRate * dt;
    }

    double &growthTime = mesh->attributes().attrMap<QString, double>("Growth Time")["Growth Time"];
    growthTime += dt;
    mdxInfo << "Growth time: " << growthTime << endl;

    // Age cells and find average age
    auto &cellAge = mesh->attributes().attrMap<CCIndex, double>("Cell Age");
    double avgCellAge = 0;
    CCIndexSet cells;
    for(CCIndex c : cs.faces()) {
      avgCellAge += (cellAge[c] += dt);
      cells.insert(c);
    }
    avgCellAge /= cs.faces().size();

    mdxInfo << "Avg cell age: " << avgCellAge << endl;
      

    return true;
  }

  // Initialize to grab subdivider
  bool CellDiskDivide::initialize(QWidget *parent)
  {
    Mesh *mesh = currentMesh();
    if(!mesh)
      throw(QString("CellDiskDivide::step No current mesh"));

    // Call base initialize
    if(!CellTissueCell2dDivide::initialize(parent))
      return false;

    if(!getProcess(parm("Solver Process"), solverProcess))
      throw(QString("CellDiskDivide::initialize Cannot make solver process"));
    solverProcess->initialize(parent);

    if(!getProcess(parm("Tissue Process"), tissueProcess))
      throw(QString("CellDisk::initialize Cannot make tissue process"));

    // Setup subdivision object
    subdiv.mdx = *CellTissueCell2dDivide::subdivider();
    subdiv.ms = MassSpringSubdivide(&mesh->attributes().attrMap<CCIndex, MassSpringEdgeData>("Mass Spring Data"), &mesh->indexAttr());

    return true;
  }

  // Run a step of cell division
  bool CellDiskDivide::step() 
  { 
    // Pass our subdivide
    return CellTissueCell2dDivide::step(currentMesh(), &subdiv);
  }

  bool SetRestLength::run()
  {
    Mesh *mesh = currentMesh();
    if(!mesh)
      throw QString("%1::run() No current mesh").arg(name());

    auto &cs = mesh->ccStructure("Tissue"); // RSS FIXME
    auto &indexAttr  = mesh->indexAttr();
    auto &msEdgeAttr = mesh->attributes().attrMap<CCIndex, MassSpringEdgeData>("Mass Spring Data");

    auto edges = cs.edges();
    #pragma omp parallel for
    for(uint i = 0; i < edges.size(); i++) {
      CCIndex e = edges[i];
      auto &eMS = msEdgeAttr[e];
      auto cb = cs.edgeBounds(e);
      eMS.restLen = norm(indexAttr[cb.first].pos - indexAttr[cb.second].pos);
    }

    return true;
  }
  REGISTER_PROCESS(SetRestLength);

  bool QuantifyAngles::run()
  {
    Mesh *mesh = currentMesh();
    if(!mesh)
      throw QString("%1::initialize() No current mesh").arg(name());

    auto &cs = mesh->ccStructure("Tissue");
    auto &indexAttr = mesh->indexAttr();

    double avgDiff3 = 0, avgDiff4 = 0;
    uint count3 = 0, count4 = 0;
    for(CCIndex v : cs.vertices()) {
      uint cellCount = cs.incidentCells(v, 2).size();
      if(cellCount < 3) // Vertex on boundary
        continue;

      auto nbs = cs.neighbors(v);
      if(nbs.size() != cellCount) // Vertex on boundary
        continue;

      if(cellCount > 4) {
        mdxInfo << QString("Vertex found incident to %1 cells").arg(cellCount) << endl;
        continue;
      }
      
      Point3d pdir, dir;
      bool first = true;
      double avgDiff = 0;
      for(CCIndex n : nbs) {
        dir = normalized(indexAttr[n].pos - indexAttr[v].pos);
        if(first)
          first = false;
        else {
          double angle = acos(pdir * dir) * 180.0/M_PI;
          double diff = fabs(360.0/cellCount - angle);

          avgDiff += diff;
        }
        pdir = dir;
      }
      avgDiff /= nbs.size();

      if(cellCount == 3) {
        avgDiff3 += avgDiff;
        count3++;
      } else if(cellCount == 4) {
        avgDiff4 += avgDiff;
        count4++;
      } else
        throw QString("%1::run() Invalid cell count %2").arg(name()).arg(cellCount);
    }
    mdxInfo << endl;
    if(count3 > 0)
      mdxInfo << "Avg angle difference for 3-way junctions:" << avgDiff3 / count3 << endl;
    if(count4 > 0)
      mdxInfo << "Avg angle difference for 4-way junctions:" << avgDiff4 / count4 << endl;

    CCIndexDoubleAttr minWalls;

    // Find average minimum wall
    double avgMinWall = 0.0;
    uint countC = 0;
    for(CCIndex f : cs.faces()) {
      double minWall = HUGE_VAL;
      for(CCIndex e : cs.bounds(f)) {
        auto pr = cs.edgeBounds(e);
        double length = norm(indexAttr[pr.first].pos - indexAttr[pr.second].pos);
        if(minWall > length)
          minWall = length;
      }
      avgMinWall += minWall;
      countC++;
      minWalls[f] = minWall;
    }
    if(countC > 0)
      mdxInfo << "Avg min wall:" << avgMinWall / countC << endl;

    return true;
  }
  REGISTER_PROCESS(QuantifyAngles);

  bool QuantifyAnglesExport::initialize(QWidget *parent)
  {
    return true;
  }


  bool QuantifyAnglesExport::run()
  {


    // go though all edges, find the ones with age < 1

    // quantify their angles

    // save to a file
    filename = parm("Filename");
    QFile file(filename);
    if(!file.open(QIODevice::WriteOnly)){
      throw(QString("File cannot be opened"));
    }
    QTextStream out(&file);

    out << "Vertex,EdgeAge,Angle" << endl;

    Mesh *mesh = currentMesh();
    if(!mesh)
      throw QString("%1::initialize() No current mesh").arg(name());

    auto &msEdgeAttr = mesh->attributes().attrMap<CCIndex, MassSpringEdgeData>("Mass Spring Data");

    auto &cs = mesh->ccStructure("Tissue");
    auto &indexAttr = mesh->indexAttr();

    int counter = 1;

    for(CCIndex v : cs.vertices()) {
      uint cellCount = cs.incidentCells(v, 2).size();
      if(cellCount < 3) // Vertex on boundary
        continue;

      auto nbs = cs.neighbors(v);
      if(nbs.size() != cellCount) // Vertex on boundary
        continue;

      if(cellCount > 4) {
        mdxInfo << QString("Vertex found incident to %1 cells").arg(cellCount) << endl;
        continue;
      }
      
      Point3d pdir, dir;
      bool first = true;
      double avgDiff = 0;
      double age = 1;
      for(CCIndex n : nbs) {
        // get the edge and its age here
        auto edgeB = edgeBetween(cs, v, n);
        MassSpringEdgeData &eMS = msEdgeAttr[edgeB];
        if(eMS.age < age) age = eMS.age;

        dir = normalized(indexAttr[n].pos - indexAttr[v].pos);
        if(first)
          first = false;
        else {
          double angle = acos(pdir * dir) * 180.0/M_PI;
          double diff = fabs(360.0/cellCount - angle);

          avgDiff += diff;
        }
        pdir = dir;
      }
      avgDiff /= nbs.size();
      std::cout << "n1 " << age << "/" << avgDiff << std::endl;
      out << counter << "," << age << "," << avgDiff << endl;
      counter++;
      
    }

     


    return true;
  }
  REGISTER_PROCESS(QuantifyAnglesExport);


  bool MeasureMinWall::run(Mesh &mesh, CCStructure &cs, CCIndexDataAttr &indexAttr, IntDoubleAttr &heatAttr)
  {
    heatAttr.clear();
    for(CCIndex f : cs.faces()) {
      auto &fIdx = indexAttr[f];
      double minWall = HUGE_VAL;
      for(CCIndex e : cs.bounds(f)) {
        auto pr = cs.edgeBounds(e);
        double length = norm(indexAttr[pr.first].pos - indexAttr[pr.second].pos);
        if(minWall > length)
          minWall = length;
      }
      int label = mesh.getLabel(fIdx.label);
      if(label <= 0) 
        continue;
      heatAttr[label] = minWall;
    }
    return true;
  }
  REGISTER_PROCESS(MeasureMinWall);

  REGISTER_PROCESS(CellDisk);
  REGISTER_PROCESS(CellDiskDivide);
  REGISTER_PROCESS(CellDiskGrowth);
  REGISTER_PROCESS(CellDiskTissue);
  REGISTER_PROCESS(CellDiskSolver);
  REGISTER_PROCESS(CellDiskSplitEdges);
  REGISTER_PROCESS(MassSpringDerivs);
  REGISTER_PROCESS(Pressure2Derivs);
}
