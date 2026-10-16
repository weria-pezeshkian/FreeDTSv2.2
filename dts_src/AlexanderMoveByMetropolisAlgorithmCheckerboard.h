#ifndef ALEXANDER_MOVE_BY_METROPOLIS_ALGORITHM_CHECKERBOARD_H
#define ALEXANDER_MOVE_BY_METROPOLIS_ALGORITHM_CHECKERBOARD_H

#include <memory>
#include <string>
#include <vector>
#include "SimDef.h"
#include "RNG.h"
#include "AbstractAlexanderMove.h"
#include "CheckerboardSchedule.h"
#include "EdgeFlipKernel.h"

class State;
class links;

// Flips edges of one checkerboard colour in parallel, one serial Metropolis chain per cell; an edge
// belongs to the cell of its midpoint, before and after the flip. Global totals are committed after
// each colour, as for the checkerboard vertex moves.
class AlexanderMoveByMetropolisAlgorithmCheckerboard : public AbstractAlexanderMove {
public:
    AlexanderMoveByMetropolisAlgorithmCheckerboard(State* pState, double rate);
    ~AlexanderMoveByMetropolisAlgorithmCheckerboard();
    bool Initialize();
    bool EvolveOneStep(int step);
    std::string CurrentState();
    inline std::string GetDerivedDefaultReadName() {return "MetropolisAlgorithmOpenMPCheckerboard";}
    inline static std::string GetDefaultReadName() {return "MetropolisAlgorithmOpenMPCheckerboard";}

private:
    struct CellRecord {
        double dE, dV, dA, dC;
        int accepted, attempted;
    };

    double MinimumCellWidth();
    void BuildCells();
    void EvolveCell(int cell, int sweep_seed, EdgeFlipKernel& kernel);
    int MidpointCell(vertex* a, vertex* b) const;
    void CommitCell(const CellRecord& record);

    State* m_pState;
    CheckerboardSchedule m_Schedule;
    std::vector<std::unique_ptr<EdgeFlipKernel> > m_Kernels;
    std::vector<std::vector<links*> > m_CellEdges;
    std::vector<CellRecord> m_CellRecords;
    std::vector<int> m_ActiveCells;
    std::vector<int> m_CellSizes;
};

#endif
