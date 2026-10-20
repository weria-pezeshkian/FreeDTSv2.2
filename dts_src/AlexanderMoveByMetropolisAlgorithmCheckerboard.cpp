#include <algorithm>
#include <cmath>
#ifdef _OPENMP
# include <omp.h>
#endif
#include "AlexanderMoveByMetropolisAlgorithmCheckerboard.h"
#include "State.h"
#include "./Registry/FactoryForRegistration.h"

AlexanderMoveByMetropolisAlgorithmCheckerboard::AlexanderMoveByMetropolisAlgorithmCheckerboard(State* pState, double rate)
    : m_pState(pState) {
    m_NumberOfMovePerStep = rate;
}
AlexanderMoveByMetropolisAlgorithmCheckerboard::~AlexanderMoveByMetropolisAlgorithmCheckerboard() {
}
bool AlexanderMoveByMetropolisAlgorithmCheckerboard::Initialize() {
    m_NumberOfAttemptedMoves = 0;
    m_AcceptedMoves = 0;

    int n_threads = std::max(1, m_pState->GetThreads_Number());
    m_Kernels.clear();
    for (int t = 0; t < n_threads; t++) {
        m_Kernels.emplace_back(new EdgeFlipKernel(m_pState));
        m_Kernels.back()->Initialize();
    }
    return true;
}
// a flip reads and writes within L_max of its quad, and the quad lies within L_max of the edge midpoint;
// the margin makes the bound strict
double AlexanderMoveByMetropolisAlgorithmCheckerboard::MinimumCellWidth() {
    return 3 * std::sqrt(m_pState->GetSimulation()->GetMaxL2()) * 1.01;
}
int AlexanderMoveByMetropolisAlgorithmCheckerboard::MidpointCell(vertex* a, vertex* b) const {
    Vec3D mid = CanonicalMidpoint(a, b, *(m_pState->GetMesh()->GetBox()));
    return m_Schedule.GetGrid().CellIndex(mid(0), mid(1), mid(2));
}
void AlexanderMoveByMetropolisAlgorithmCheckerboard::ResizeCells(int n_cells) {
    if ((int)m_CellEdges.size() == n_cells)
        return;
    m_CellEdges.assign(n_cells, std::vector<links*>());
    m_CellRecords.resize(n_cells);
    m_OccupiedCells.clear();
}
// m_EdgeCell holds the midpoint cell of every edge in GetRightL()
void AlexanderMoveByMetropolisAlgorithmCheckerboard::BuildCells() {
    for (std::vector<int>::const_iterator it = m_OccupiedCells.begin(); it != m_OccupiedCells.end(); ++it)
        m_CellEdges[*it].clear();
    m_OccupiedCells.clear();
    const std::vector<links*>& edges = m_pState->GetMesh()->GetRightL();
    for (int i = 0; i < (int)m_EdgeCell.size(); i++) {
        int cell = m_EdgeCell[i];
        if (m_CellEdges[cell].empty())
            m_OccupiedCells.push_back(cell);
        m_CellEdges[cell].push_back(edges[i]);
    }
    std::sort(m_OccupiedCells.begin(), m_OccupiedCells.end());
    m_Schedule.GroupByColor(m_OccupiedCells, m_ColorCells);
}
bool AlexanderMoveByMetropolisAlgorithmCheckerboard::EvolveOneStep(int) {

    m_Schedule.Draw(*m_pState->GetRandomNumberGenerator(), *(m_pState->GetMesh()->GetBox()), MinimumCellWidth());
    int sweep_seed = m_Schedule.GetSweepSeed();
    ResizeCells(m_Schedule.GetGrid().NumberOfCells());

    const std::vector<links*>& edges = m_pState->GetMesh()->GetRightL();
    int n_edges = (int)edges.size();
    m_EdgeCell.resize(n_edges);

    #pragma omp parallel num_threads((int)m_Kernels.size())
    {
        int thread_id = 0;
#ifdef _OPENMP
        thread_id = omp_get_thread_num();
#endif
        #pragma omp for schedule(static)
        for (int i = 0; i < n_edges; i++)
            m_EdgeCell[i] = MidpointCell(edges[i]->GetV1(), edges[i]->GetV2());
        #pragma omp single
        BuildCells();

        for (int c = 0; c < CheckerboardCellGrid::NumberOfColors; c++) {
            const std::vector<int>& cells = m_ColorCells[m_Schedule.GetColor(c)];
            #pragma omp for schedule(dynamic, 1)
            for (int i = 0; i < (int)cells.size(); i++)
                EvolveCell(cells[i], sweep_seed, *m_Kernels[thread_id]);
            // commit in cell order so the result does not depend on the thread count
            #pragma omp single
            {
                for (int i = 0; i < (int)cells.size(); i++)
                    CommitCell(m_CellRecords[cells[i]]);
            }
        }
    }
    return true;
}
void AlexanderMoveByMetropolisAlgorithmCheckerboard::EvolveCell(int cell, int sweep_seed, EdgeFlipKernel& kernel) {
    CellRecord& record = m_CellRecords[cell];
    record.dE = record.dV = record.dA = record.dC = 0;
    record.accepted = record.attempted = 0;

    const std::vector<links*>& cell_edges = m_CellEdges[cell];
    int n_e = (int)cell_edges.size();
    RNG rng(RNG::StreamSeed(sweep_seed, cell));
    double target = m_NumberOfMovePerStep * n_e;
    int n_moves = (int)std::floor(target);
    if (rng.UniformRNG(1.0) < target - n_moves)
        n_moves++;

    for (int m = 0; m < n_moves; m++) {
        links* p_edge = cell_edges[rng.IntRNG(n_e)];
        record.attempted++;
        // earlier flips in this cell may have moved the edge out
        if (MidpointCell(p_edge->GetV1(), p_edge->GetV2()) != cell)
            continue;
        // the edge after the flip; canonical order keeps the reverse flip's test bitwise equal
        if (p_edge->GetMirrorFlag() &&
            MidpointCell(p_edge->GetMirrorLink()->GetV3(), p_edge->GetV3()) != cell)
            continue;
        EdgeFlipOutcome outcome = kernel.TryFlip(p_edge, rng.UniformRNG(1.0), {record.dV, record.dA, record.dC});
        if (outcome.accepted) {
            record.accepted++;
            record.dE += outcome.dE_elastic;
            record.dV += outcome.dVolume;
            record.dA += outcome.dArea;
            record.dC += outcome.dCurvature;
        }
    }
}
void AlexanderMoveByMetropolisAlgorithmCheckerboard::CommitCell(const CellRecord& record) {
    m_pState->GetEnergyCalculator()->AddToTotalEnergy(record.dE);
    VAHGlobalMeshProperties* p_VAH = m_pState->GetVAHGlobalMeshProperties();
    if (p_VAH->GetCalculateVAH()) {
        p_VAH->Add2Volume(record.dV);
        p_VAH->Add2TotalArea(record.dA);
        p_VAH->Add2GlobalCurvature(record.dC);
    }
    m_AcceptedMoves += record.accepted;
    m_NumberOfAttemptedMoves += record.attempted;
}
std::string AlexanderMoveByMetropolisAlgorithmCheckerboard::CurrentState() {

    std::string state = GetBaseDefaultReadName() + " = " + GetDerivedDefaultReadName();
    state = state + " " + Nfunction::D2S(m_NumberOfMovePerStep);
    return state;
}
namespace
{
#ifdef _OPENMP

    AbstractAlexanderMove* Create_AM_OpenMPCheckerboard(
        std::istream& input,
        State* state)
    {
        double rate = 0.0;

        if (!(input >> rate)) {
            return nullptr;
        }

        std::string rest;
        std::getline(input, rest);

        return new AlexanderMoveByMetropolisAlgorithmCheckerboard(state, rate);
    }

#else

    AbstractAlexanderMove* Create_AM_OpenMPCheckerboard(
        std::istream&,
        State*)
    {
        std::cerr
            << "Error: '"
            << AlexanderMoveByMetropolisAlgorithmCheckerboard::GetDefaultReadName()
            << "' requires OpenMP, but this build does not support it.\n"
            << "Use MetropolisAlgorithm instead.\n";

        return nullptr;
    }

#endif

    const bool registered = []()
    {
        FactoryAlexanderMove::Instance().Register(
            AlexanderMoveByMetropolisAlgorithmCheckerboard::GetDefaultReadName(),
            &Create_AM_OpenMPCheckerboard
        );
        return true;
    }();
}
