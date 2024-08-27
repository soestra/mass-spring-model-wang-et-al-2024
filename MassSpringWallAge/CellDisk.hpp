//
// This file is part of MorphoDynamX - http://www.MorphoDynamX.org
// Copyright (C) 2012-2016 Richard S. Smith and collaborators.
//
// If you use MorphoDynamX in your work, please cite:
//   http://dx.doi.org/10.7554/eLife.05864
//
// MorphoDynamX is free software, and is licensed under under the terms of the 
// GNU General (GPL) Public License version 2.0, http://www.gnu.org/licenses.
//
#ifndef CELL_APEX_HPP
#define CELL_APEX_HPP

#include <Attributes.hpp>
#include <Function.hpp>

#include <Contour.hpp>

#include <Mesh.hpp>
#include <Process.hpp>
#include <Random.hpp>
#include <MeshProcessSystem.hpp>

#include <Solver.hpp>

#include <MDXProcessTissue.hpp>
#include <MDXProcessCellDivide.hpp>
#include <MDXProcessMorphogens.hpp>
#include <MeshProcessStructure.hpp>
#include <MeshProcessMeasures2D.hpp>

using namespace mdx;

namespace CellDisk
{
  class CellDiskDivide;
  class CellDiskTissue;
  class CellDiskGrowth;
  class CellDiskSolver;
  class CellDiskSplitEdges;
  class MassSpringDerivs;
  class Pressure2Derivs;
  class QuantifyAngles;

  // Store rest lengths in edges
  struct MassSpringEdgeData
  {
    double restLen; // Rest length
    CCIndex f; // Cell if edge is on border for quick lookup
    bool border; // Is edge on the border?
    double age; // Age of wall from 0.0 - 1.0

    MassSpringEdgeData() : 
      restLen(1.0), f(0), border(false), age(1.0) {}

		bool operator==(const MassSpringEdgeData &other)
    {
      if(restLen == other.restLen and f == other.f and border == other.border and age == other.age)
        return true;
      return false;
    }
  };
  typedef AttrMap<CCIndex,MassSpringEdgeData> MSEdgeDataAttr;

  // Class for subdivision of cells
  class MassSpringSubdivide : virtual public Subdivide
  {
  public:
    MassSpringSubdivide() {}
    MassSpringSubdivide(MSEdgeDataAttr *msEdgeAttr, CCIndexDataAttr *indexAttr) 
     : msEdgeData(msEdgeAttr), indexAttr(indexAttr) {}
    void splitCellUpdate(Dimension dim, const CCStructure &cs, const CCStructure::SplitStruct& ss,
                   CCIndex otherP = CCIndex(), CCIndex otherN = CCIndex(), double interpPos = 0.5) const
    {
      if(!msEdgeData) {
        mdxInfo << "MassSpringSubdivide::splitCellUpdate Edge data not defined" << endl;
        return;
      }
      if(!indexAttr) {
        mdxInfo << "MassSpringSubdivide::splitCellUpdate Index data not defined" << endl;
        return;
      }

      if(dim == 1) {
        // Propagate rest length
        MassSpringEdgeData &rMS = (*msEdgeData)[ss.parent];
        MassSpringEdgeData &pMS = (*msEdgeData)[ss.childP];
        MassSpringEdgeData &nMS = (*msEdgeData)[ss.childN];

        nMS.restLen = (1.0 - interpPos) * rMS.restLen;
        pMS.restLen = interpPos * rMS.restLen;
        pMS.age = nMS.age = rMS.age;
      } else if(dim  == 2) {
        // Assign rest length to new wall, just use actual length
        MassSpringEdgeData &eMS = (*msEdgeData)[ss.membrane];
        CCIndexPair ep = cs.edgeBounds(ss.membrane);
        eMS.restLen = norm((*indexAttr)[ep.first].pos - (*indexAttr)[ep.second].pos);
        eMS.age = 0.0;
      }
    }

    MSEdgeDataAttr *msEdgeData;
    CCIndexDataAttr *indexAttr;
  };

  class CellDiskSubdivide : virtual public Subdivide
  {
  public:
    CellDiskSubdivide() {}

    void splitCellUpdate(Dimension dim, const CCStructure &cs, const CCStructure::SplitStruct& ss, 
        CCIndex otherP = CCIndex(), CCIndex otherN = CCIndex(),  double interpPos = 0.5)
    {
      
      mdx.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);
      ms.splitCellUpdate(dim, cs, ss, otherP, otherN, interpPos);
    }

    MDXSubdivide mdx;
    MassSpringSubdivide ms;
  };
 
  class CellDisk : public Process
  {
  public:
    CellDisk(const Process &process) : Process(process) 
    {
      setName("Model/CCF/01 Cell Disk");
      setDesc("Growing cellular apex");
      setIcon(QIcon(":/images/CellDisk.png"));
      addParm("Converge Threshold", "Convergence threshold for mass-spring system", ".0001");
      addParm("Cell Kill", "Radius to remove cells, 0 prevents removal", "50.0");
      addParm("Solver Process", "Name of the process for the Meinhardt Solver", "Model/CCF/02 Mass Spring Solver");
      addParm("Tissue Process", "Name of process for Cell Tissue", "Model/CCF/10 Cell Tissue");
      addParm("Growth Process", "Name of the process for Growth", "Model/CCF/04 Mass Spring Growth");
      addParm("Divide Process", "Name of the process for Cell Division", "Model/CCF/05 Divide Cells");
      addParm("Split Edges Process", "Name of the process to split edges", "Model/CCF/06 Split Edges");
      addParm("Quantify Angles Process", "Name of the process to quantify angles", "Model/CCF/20 Quantify Angles");
      addParm("Stop Steps", "Stop after this number of steps", "1500");
    }
    bool initialize(QWidget *parent);
    bool step();
    bool rewind(QWidget *parent);

  private:
    Mesh *mesh = 0;

    CellDiskSolver *solverProcess = 0;
    CellDiskTissue *tissueProcess = 0;
    CellDiskGrowth *growthProcess = 0;
    CellDiskDivide *divideProcess = 0;
    CellDiskSplitEdges *splitEdgesProcess = 0;
    QuantifyAngles *quantifyAnglesProcess = 0;
    MSEdgeDataAttr *msEdgeAttr = 0;

    int steps = 0;
  };

  class CellDiskSolver :  public Process, public Solver<Point3d>
  {
  public:
    CellDiskSolver(const Process &process) : Process(process) 
    {
      setName("Model/CCF/02 Mass Spring Solver");
      addParm("Mass Spring Process", "Name of derivs process for mass-spring mechanics", "Model/CCF/03a Mass Spring Derivs");
      addParm("Pressure Process", "Name of derivs process for pressure mechanics", "Model/CCF/03b 2D Pressure Derivs");
      addParm("Tissue Process", "Name of process for Cell Tissue", "Model/CCF/10 Cell Tissue");
    }
    bool initialize(QWidget *parent);
  };

  class MassSpringDerivs : public Process, public SolverDerivs<Point3d>
  {
  public:
    MassSpringDerivs(const Process &process): Process(process) 
    {
      setName("Model/CCF/03a Mass Spring Derivs");
      insertParm("Spring Constant", "Stiffness of springs", "1.0", 1);
      insertParm("Stiffness New Spring", "Stiffness New Spring", "0.0", 1);
    }

    // Reimplemented solver methods
    void initDerivs(SolverT &solver, VertexAttr &vAttr, EdgeAttr &eAttr);
    void calcDerivatives(const SolverT &solver, CCIndex v, VectorT &values);
    void calcJacobian(const SolverT &solver, VertexAttr &vAttr, EdgeAttr &eAttr);

    void setValues(const SolverT &solver, CCIndex v, const VectorT &values)
    {
      (*indexAttr)[v].pos = values;
    }
  
    void getValues(const SolverT &solver, CCIndex v, VectorT &values)
    {
      values = (*indexAttr)[v].pos;
    }

    CCStructure const *cs = 0;
    CCIndexDataAttr *indexAttr = 0;
    MSEdgeDataAttr *msEdgeAttr = 0;

    double k = 1.0;

    double ageConstant = 0.0;
  };

  class Pressure2Derivs : public Process, public SolverDerivs<Point3d>
  {
  public:
    Pressure2Derivs(const Process &process): Process(process) 
    {
      setName("Model/CCF/03b 2D Pressure Derivs");
      insertParm("Pressure", "Pressure", "1.0", 1);
    }

    // Reimplemented solver methods
    void initDerivs(SolverT &solver, VertexAttr &vAttr, EdgeAttr &eAttr);
    void calcDerivatives(const SolverT &solver, CCIndex v, VectorT &values);
    void calcJacobian(const SolverT &solver, CCIndex v, MatrixT &j) {} // no diagonal Jacobian part
    void calcJacobian(const SolverT &solver, CCIndex v, CCIndex nb, MatrixT &j);

    void setValues(const SolverT &solver, CCIndex v, const VectorT &values)
    {
      (*indexAttr)[v].pos = values;
    }
  
    void getValues(const SolverT &solver, CCIndex v, VectorT &values)
    {
      values = (*indexAttr)[v].pos;
    }

    CCStructure const *cs = 0;
    CCIndexDataAttr *indexAttr = 0;
    MSEdgeDataAttr *msEdgeAttr = 0;

    double pressure = 0.5;
  };

  class CellDiskGrowth : public Process
  {
  public:
    CellDiskGrowth(const Process &process): Process(process) 
    {
      setName("Model/CCF/04 Mass Spring Growth");
      addParm("Growth Rate", "The rate of growth", "1.0");
      addParm("Dt", "Time step for growth", "0.1");
      addParm("Wall Aging", "Time until new walls are fully mature", "10.0");
      addParm("New Wall Stiffness Factor", "New Wall Stiffness Factor", "1.0");


      addParm("Tissue Process", "Name of the process for the Cell Tissue", "Model/CCF/10 Cell Tissue");
    }

    // Initialize to grab subdivider
    bool initialize(QWidget *parent);

    // Run a step of cell division
    bool step();

  private:
    Mesh *mesh = 0;
    CellDiskTissue *tissueProcess = 0;
  };

  class CellDiskDivide : public CellTissueCell2dDivide
  {
  public:
    CellDiskDivide(const Process &process) : CellTissueCell2dDivide(process)
    {
      setName("Model/CCF/05 Divide Cells");
      addParm("Solver Process", "Name of the process for the Meinhardt solver", "Model/CCF/02 Mass Spring Solver");
      addParm("Tissue Process", "Name of the process for the Cell Tissue", "Model/CCF/10 Cell Tissue");
    }

    // Initialize to grab subdivider
    bool initialize(QWidget *parent);

    // Run a step of cell division
    bool step();

    CellDiskSubdivide *subdivider() { return &subdiv; }

  private:
    CellDiskSubdivide subdiv;

    CellDiskSolver *solverProcess = 0;
    CellDiskTissue *tissueProcess = 0;
  };

  class CellDiskSplitEdges : public SplitEdges
  {
  public:
    CellDiskSplitEdges(const Process &process) : SplitEdges(process) 
    {
      setName("Model/CCF/06 Split Edges");
    }
  };

  class SetRestLength :  public Process
  {
  public:
    SetRestLength(const Process &process) : Process(process) 
    {
      setName("Model/CCF/07 Set Rest Length");
      setDesc("Set the rest length to the current configuration");
    }
    bool run();
  };

  class CellDiskTissue : public CellTissueProcess
  {
  public:
    CellDiskTissue(const Process &process) : CellTissueProcess(process) 
    {
      setName("Model/CCF/10 Cell Tissue");
    }
  };

  class QuantifyAngles :  public Process
  {
  public:
    QuantifyAngles(const Process &process) : Process(process) 
    {
      setName("Model/CCF/20 Quantify Angles");
      setDesc("Quantify the angles of vertices (excludes boundary vertices");
    }
    bool run();
  };

  class QuantifyAnglesExport :  public Process
  {
  public:
    QuantifyAnglesExport(const Process &process) : Process(process) 
    {
      setName("Model/CCF/20a Quantify Angles Export");
      setDesc("Quantify the angles of vertices (excludes boundary vertices");
      addParm("Filename", "Filename", "");
    }
    bool run();

    bool initialize(QWidget *parent);

    QString filename;

  };

  class MeasureMinWall : public Measure2D
  {
  public:
    MeasureMinWall(const Process& process) : Measure2D(process)
    {
      setName("Model/CCF/21 Min Wall");
      setDesc("Generate heat map based on minimum wall length. Label faces in the mesh before running the process.");
      setParmDefault("Heat Name", "Min Wall");
      heatUnit = UM;
    }
    bool run(Mesh &mesh, CCStructure &cs, CCIndexDataAttr &indexAttr, IntDoubleAttr &heatAttr);
  };
}
#endif 
