#include <algorithm>
#include <cmath>
#include <cstdlib>
#ifdef _OPENMP
# include <omp.h>
#endif
#include "EvolveVerticesByMetropolisAlgorithmCheckerboard.h"
#include "State.h"
#include "./Registry/FactoryVertexPositionIntegrator.h"

EvolveVerticesByMetropolisAlgorithmCheckerboard::EvolveVerticesByMetropolisAlgorithmCheckerboard(State* pState, double rate_surf, double rate_edge, double dr)
    : m_pState(pState) {
    m_NumberOfMovePerStep_Surf = rate_surf;
    m_NumberOfMovePerStep_Edge = rate_edge;
    m_DR = dr;
}
EvolveVerticesByMetropolisAlgorithmCheckerboard::~EvolveVerticesByMetropolisAlgorithmCheckerboard() {
}
void EvolveVerticesByMetropolisAlgorithmCheckerboard::Initialize() {
    // the kernel has no bond or vertex-vertex nonbonded terms, and bonds may span cells
    if (m_pState->GetBondedPotentialBetweenVertices()->GetDerivedDefaultReadName() != "No" ||
        m_pState->GetNonbondedInteractionBetweenVertices()->GetDerivedDefaultReadName() != "No") {
        std::cout << "---> error: " << GetDerivedDefaultReadName() << " does not support "
                  << AbstractBondedPotentialBetweenVertices::GetBaseDefaultReadName() << " or "
                  << AbstractNonbondedInteractionBetweenVertices::GetBaseDefaultReadName() << "; use MetropolisAlgorithm \n";
        exit(EXIT_FAILURE);
    }
    m_NumberOfAttemptedMoves = 0;
    m_AcceptedMoves = 0;

    int n_threads = std::max(1, m_pState->GetThreads_Number());
    m_Kernels.clear();
    for (int t = 0; t < n_threads; t++) {
        m_Kernels.emplace_back(new VertexMoveKernel(m_pState));
        m_Kernels.back()->Initialize();
    }
    const std::vector<vertex*>& active_v = m_pState->GetMesh()->GetActiveV();
    int max_id = -1;
    for (std::vector<vertex*>::const_iterator it = active_v.begin(); it != active_v.end(); ++it)
        max_id = std::max(max_id, (*it)->GetVID());
    m_Frozen.assign(max_id + 1, 0);
    for (std::vector<vertex*>::const_iterator it = active_v.begin(); it != active_v.end(); ++it) {
        if ((*it)->GetVID() >= 0 && m_FreezGroupName == (*it)->GetGroupName())
            m_Frozen[(*it)->GetVID()] = 1;
    }
}
bool EvolveVerticesByMetropolisAlgorithmCheckerboard::IsFrozen(vertex* p_vertex) const {
    int id = p_vertex->GetVID();
    if (id >= 0 && id < (int)m_Frozen.size())
        return m_Frozen[id];
    return m_FreezGroupName == p_vertex->GetGroupName();
}
// cells of one colour are at least one cell width apart; this width keeps their read/write sets
// (2-ring + links of the 1-ring, 27 voxels around the target voxel) disjoint
double EvolveVerticesByMetropolisAlgorithmCheckerboard::MinimumCellWidth() {
    Voxelization<vertex>* p_voxels = m_pState->GetVoxelization();
    double voxel_side = std::max(p_voxels->GetXSideVoxel(), std::max(p_voxels->GetYSideVoxel(), p_voxels->GetZSideVoxel()));
    double max_length = std::sqrt(m_pState->GetSimulation()->GetMaxL2());
    return std::max(3 * max_length, 4 * voxel_side) + 2 * m_DR;
}
void EvolveVerticesByMetropolisAlgorithmCheckerboard::ResizeCells(int n_cells) {
    if ((int)m_CellSurfV.size() == n_cells)
        return;
    m_CellSurfV.assign(n_cells, std::vector<vertex*>());
    m_CellEdgeV.assign(n_cells, std::vector<vertex*>());
    m_CellRecords.resize(n_cells);
    m_OccupiedCells.clear();
}
// m_VertexCell holds the cell of every surface vertex, then every edge vertex (-1 if frozen)
void EvolveVerticesByMetropolisAlgorithmCheckerboard::BuildCells() {
    for (std::vector<int>::const_iterator it = m_OccupiedCells.begin(); it != m_OccupiedCells.end(); ++it) {
        m_CellSurfV[*it].clear();
        m_CellEdgeV[*it].clear();
    }
    m_OccupiedCells.clear();
    const std::vector<vertex*>& surf_v = m_pState->GetMesh()->GetSurfV();
    const std::vector<vertex*>& edge_v = m_pState->GetMesh()->GetEdgeV();
    int n_surf = (int)surf_v.size();
    for (int i = 0; i < (int)m_VertexCell.size(); i++) {
        int cell = m_VertexCell[i];
        if (cell < 0)
            continue;
        if (m_CellSurfV[cell].empty() && m_CellEdgeV[cell].empty())
            m_OccupiedCells.push_back(cell);
        if (i < n_surf)
            m_CellSurfV[cell].push_back(surf_v[i]);
        else
            m_CellEdgeV[cell].push_back(edge_v[i - n_surf]);
    }
    std::sort(m_OccupiedCells.begin(), m_OccupiedCells.end());
    m_Schedule.GroupByColor(m_OccupiedCells, m_ColorCells);
    // largest cells first, so the last cells handed out are the short ones
    m_ColorOrder = m_ColorCells;
    for (int c = 0; c < CheckerboardCellGrid::NumberOfColors; c++) {
        std::stable_sort(m_ColorOrder[c].begin(), m_ColorOrder[c].end(), [this](int a, int b) {
            return m_CellSurfV[a].size() + m_CellEdgeV[a].size() > m_CellSurfV[b].size() + m_CellEdgeV[b].size();
        });
    }
}
bool EvolveVerticesByMetropolisAlgorithmCheckerboard::EvolveOneStep(int) {

    m_Schedule.Draw(*m_pState->GetRandomNumberGenerator(), *(m_pState->GetMesh()->GetBox()), MinimumCellWidth());
    const CheckerboardCellGrid& grid = m_Schedule.GetGrid();
    int sweep_seed = m_Schedule.GetSweepSeed();
    ResizeCells(grid.NumberOfCells());

    const std::vector<vertex*>& surf_v = m_pState->GetMesh()->GetSurfV();
    const std::vector<vertex*>& edge_v = m_pState->GetMesh()->GetEdgeV();
    int n_surf = (int)surf_v.size();
    int n_vertices = n_surf + (int)edge_v.size();
    m_VertexCell.resize(n_vertices);

    // one parallel region per sweep; forking per colour costs more than the work of a small colour
    #pragma omp parallel num_threads((int)m_Kernels.size())
    {
        int thread_id = 0;
#ifdef _OPENMP
        thread_id = omp_get_thread_num();
#endif
        #pragma omp for schedule(static)
        for (int i = 0; i < n_vertices; i++) {
            vertex* p_vertex = i < n_surf ? surf_v[i] : edge_v[i - n_surf];
            m_VertexCell[i] = IsFrozen(p_vertex) ? -1 : grid.CellIndex(p_vertex->GetXPos(), p_vertex->GetYPos(), p_vertex->GetZPos());
        }
        #pragma omp single
        BuildCells();

        for (int c = 0; c < CheckerboardCellGrid::NumberOfColors; c++) {
            const std::vector<int>& order = m_ColorOrder[m_Schedule.GetColor(c)];
            #pragma omp for schedule(dynamic, 1)
            for (int i = 0; i < (int)order.size(); i++)
                EvolveCell(order[i], sweep_seed, *m_Kernels[thread_id]);
            // commit in cell order so the result does not depend on the thread count
            #pragma omp single
            {
                const std::vector<int>& cells = m_ColorCells[m_Schedule.GetColor(c)];
                for (int i = 0; i < (int)cells.size(); i++)
                    CommitCell(m_CellRecords[cells[i]]);
            }
        }
    }
    return true;
}
void EvolveVerticesByMetropolisAlgorithmCheckerboard::EvolveCell(int cell, int sweep_seed, VertexMoveKernel& kernel) {
    CellRecord& record = m_CellRecords[cell];
    record.dE = record.dV = record.dA = record.dC = 0;
    record.accepted = record.attempted = 0;
    record.accepted_moves.clear();

    RNG rng(RNG::StreamSeed(sweep_seed, cell));
    MoveCellVertices(cell, m_CellSurfV[cell], m_NumberOfMovePerStep_Surf, rng, kernel, record);
    MoveCellVertices(cell, m_CellEdgeV[cell], m_NumberOfMovePerStep_Edge, rng, kernel, record);
}
void EvolveVerticesByMetropolisAlgorithmCheckerboard::MoveCellVertices(int cell, const std::vector<vertex*>& cell_vertices, double rate,
                                                                       RNG& rng, VertexMoveKernel& kernel, CellRecord& record) {
    int n_v = (int)cell_vertices.size();
    if (n_v == 0)
        return;
    double target = rate * n_v;
    int n_moves = (int)std::floor(target);
    if (rng.UniformRNG(1.0) < target - n_moves)
        n_moves++;

    for (int m = 0; m < n_moves; m++) {
        vertex* p_vertex = cell_vertices[rng.IntRNG(n_v)];
        double dx = 1 - 2 * rng.UniformRNG(1.0);
        double dy = 1 - 2 * rng.UniformRNG(1.0);
        double dz = 1 - 2 * rng.UniformRNG(1.0);
        Vec3D Dx(m_DR * dx, m_DR * dy, m_DR * dz);
        record.attempted++;
        if (!m_pState->GetBoundary()->MoveHappensWithinTheBoundary(Dx(0), Dx(1), Dx(2), p_vertex))
            continue;
        // leaving the cell would break the separation between concurrent cells; the rule is symmetric
        if (m_Schedule.GetGrid().CellIndex(p_vertex->GetXPos() + Dx(0), p_vertex->GetYPos() + Dx(1), p_vertex->GetZPos() + Dx(2)) != cell)
            continue;
        double thermal = rng.UniformRNG(1.0);
        VertexMoveOutcome outcome = kernel.TryMove(p_vertex, Dx, thermal);
        if (outcome.accepted) {
            record.accepted++;
            record.dE += outcome.dE_elastic;
            record.dV += outcome.dVolume;
            record.dA += outcome.dArea;
            record.dC += outcome.dCurvature;
            record.accepted_moves.push_back(std::make_pair(p_vertex, Dx));
        }
    }
}
void EvolveVerticesByMetropolisAlgorithmCheckerboard::CommitCell(CellRecord& record) {
    m_pState->GetEnergyCalculator()->AddToTotalEnergy(record.dE);
    VAHGlobalMeshProperties* p_VAH = m_pState->GetVAHGlobalMeshProperties();
    if (p_VAH->GetCalculateVAH()) {
        p_VAH->Add2Volume(record.dV);
        p_VAH->Add2TotalArea(record.dA);
        p_VAH->Add2GlobalCurvature(record.dC);
    }
    AbstractApplyConstraintBetweenGroups* p_constraint = m_pState->GetApplyConstraintBetweenGroups();
    for (std::vector<std::pair<vertex*, Vec3D> >::const_iterator it = record.accepted_moves.begin(); it != record.accepted_moves.end(); ++it)
        p_constraint->AcceptVertexMove(it->first, it->second);
    m_AcceptedMoves += record.accepted;
    m_NumberOfAttemptedMoves += record.attempted;
}
std::string EvolveVerticesByMetropolisAlgorithmCheckerboard::CurrentState() {

    std::string state = AbstractVertexPositionIntegrator::GetBaseDefaultReadName() + " = " + GetDerivedDefaultReadName();
    state = state + " " + Nfunction::D2S(m_NumberOfMovePerStep_Surf) + " " + Nfunction::D2S(m_NumberOfMovePerStep_Edge);
    state = state + " " + Nfunction::D2S(m_DR);
    return state;
}
namespace
{
#ifdef _OPENMP

    AbstractVertexPositionIntegrator* Create_VerPosMC_OpenMPCheckerboard(
        std::istream& input,
        State* state)
    {
        double rate_surf = 0.0;
        double rate_edge = 0.0;
        double dr = 0.0;

        if (!(input >> rate_surf >> rate_edge >> dr)) {
            return nullptr;
        }

        std::string rest;
        std::getline(input, rest);

        return new EvolveVerticesByMetropolisAlgorithmCheckerboard(
            state,
            rate_surf,
            rate_edge,
            dr
        );
    }

#else

    AbstractVertexPositionIntegrator* Create_VerPosMC_OpenMPCheckerboard(
        std::istream&,
        State*)
    {
        std::cerr
            << "Error: '"
            << EvolveVerticesByMetropolisAlgorithmCheckerboard::GetDefaultReadName()
            << "' requires OpenMP, but this build does not support it.\n"
            << "Use MetropolisAlgorithm instead.\n";

        return nullptr;
    }

#endif

    const bool registered = []()
    {
        FactoryVertexPositionIntegrator::Instance().Register(
            EvolveVerticesByMetropolisAlgorithmCheckerboard::GetDefaultReadName(),
            &Create_VerPosMC_OpenMPCheckerboard
        );
        return true;
    }();
}
