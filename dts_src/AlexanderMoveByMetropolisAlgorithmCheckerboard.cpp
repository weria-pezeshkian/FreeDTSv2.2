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
void AlexanderMoveByMetropolisAlgorithmCheckerboard::BuildCells() {
    int n_cells = m_Schedule.GetGrid().NumberOfCells();
    m_CellEdges.resize(n_cells);
    m_CellRecords.resize(n_cells);
    m_CellSizes.resize(n_cells);
    for (int c = 0; c < n_cells; c++)
        m_CellEdges[c].clear();
    const std::vector<links*>& edges = m_pState->GetMesh()->GetRightL();
    for (std::vector<links*>::const_iterator it = edges.begin(); it != edges.end(); ++it)
        m_CellEdges[MidpointCell((*it)->GetV1(), (*it)->GetV2())].push_back(*it);
    for (int c = 0; c < n_cells; c++)
        m_CellSizes[c] = (int)m_CellEdges[c].size();
}
bool AlexanderMoveByMetropolisAlgorithmCheckerboard::EvolveOneStep(int) {

    m_Schedule.Draw(*m_pState->GetRandomNumberGenerator(), *(m_pState->GetMesh()->GetBox()), MinimumCellWidth());
    BuildCells();
    int sweep_seed = m_Schedule.GetSweepSeed();

    for (int c = 0; c < CheckerboardCellGrid::NumberOfColors; c++) {
        m_Schedule.ActiveCells(m_Schedule.GetColor(c), m_CellSizes, m_ActiveCells);
        int n_active = (int)m_ActiveCells.size();
        #pragma omp parallel for schedule(dynamic,1) num_threads((int)m_Kernels.size())
        for (int i = 0; i < n_active; i++) {
            int thread_id = 0;
#ifdef _OPENMP
            thread_id = omp_get_thread_num();
#endif
            EvolveCell(m_ActiveCells[i], sweep_seed, *m_Kernels[thread_id]);
        }
        // commit in cell order so the result does not depend on the thread count
        for (int i = 0; i < n_active; i++)
            CommitCell(m_CellRecords[m_ActiveCells[i]]);
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
        EdgeFlipOutcome outcome = kernel.TryFlip(p_edge, rng.UniformRNG(1.0));
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
