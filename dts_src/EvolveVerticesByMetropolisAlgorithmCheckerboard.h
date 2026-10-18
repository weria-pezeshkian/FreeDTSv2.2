#ifndef EVOLVE_VERTICES_BY_METROPOLIS_ALGORITHM_CHECKERBOARD_H
#define EVOLVE_VERTICES_BY_METROPOLIS_ALGORITHM_CHECKERBOARD_H

#include <memory>
#include <utility>
#include <vector>
#include "SimDef.h"
#include "Vec3D.h"
#include "RNG.h"
#include "AbstractVertexPositionIntegrator.h"
#include "CheckerboardSchedule.h"
#include "VertexMoveKernel.h"

class State;
class vertex;

// Moves vertices of one checkerboard colour in parallel, one serial Metropolis chain per cell;
// global totals are committed after each colour, so a move sees its own cell's earlier moves
// but not those of other cells in the colour.
class EvolveVerticesByMetropolisAlgorithmCheckerboard : public AbstractVertexPositionIntegrator {
public:
    EvolveVerticesByMetropolisAlgorithmCheckerboard(State* pState, double rate_surf, double rate_edge, double dr);
    ~EvolveVerticesByMetropolisAlgorithmCheckerboard();
    void Initialize();
    bool EvolveOneStep(int step);
    std::string CurrentState();
    inline std::string GetDerivedDefaultReadName() {return "MetropolisCheckerboard";}
    inline static std::string GetDefaultReadName() {return "MetropolisCheckerboard";}

private:
    // padded so threads filling neighbouring records do not share cache lines
    struct alignas(128) CellRecord {
        double dE, dV, dA, dC;
        int accepted, attempted;
        std::vector<std::pair<vertex*, Vec3D> > accepted_moves;
    };

    double MinimumCellWidth();
    void ResizeCells(int n_cells);
    void BuildCells();
    void EvolveCell(int cell, int sweep_seed, VertexMoveKernel& kernel);
    void MoveCellVertices(int cell, const std::vector<vertex*>& cell_vertices, double rate,
                          RNG& rng, VertexMoveKernel& kernel, CellRecord& record);
    void CommitCell(CellRecord& record);
    bool IsFrozen(vertex* p_vertex) const;

    State* m_pState;
    CheckerboardSchedule m_Schedule;
    std::vector<std::unique_ptr<VertexMoveKernel> > m_Kernels;
    std::vector<std::vector<vertex*> > m_CellSurfV;
    std::vector<std::vector<vertex*> > m_CellEdgeV;
    std::vector<CellRecord> m_CellRecords;
    std::vector<int> m_VertexCell;
    std::vector<int> m_OccupiedCells;
    std::vector<std::vector<int> > m_ColorCells;
    std::vector<std::vector<int> > m_ColorOrder;
    std::vector<char> m_Frozen;
};

#endif
