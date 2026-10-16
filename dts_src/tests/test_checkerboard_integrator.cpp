#include <algorithm>
#include "kernel_fixtures.h"
#include "../EvolveVerticesByMetropolisAlgorithmCheckerboard.h"
#include "../AlexanderMoveByMetropolisAlgorithmCheckerboard.h"
#include "../EdgeFlipKernel.h"

struct RunResult {
    std::vector<double> positions;
    double energy;
    double acceptance;
};

// replays the integrator's grid draws on a copy of the master rng, without consuming them
static CheckerboardCellGrid predict_grid(State* pState, double dr) {
    Voxelization<vertex>* p_voxels = pState->GetVoxelization();
    double voxel_side = std::max(p_voxels->GetXSideVoxel(), std::max(p_voxels->GetYSideVoxel(), p_voxels->GetZSideVoxel()));
    double min_width = std::max(3 * std::sqrt(pState->GetSimulation()->GetMaxL2()), 4 * voxel_side) + 2 * dr;
    const Vec3D& box = *(pState->GetMesh()->GetBox());
    RNG rng = *pState->GetRandomNumberGenerator();
    CheckerboardCellGrid grid;
    grid.Build(box, min_width, Vec3D(0, 0, 0));
    Vec3D shift;
    for (int a = 0; a < 3; a++)
        shift(a) = rng.UniformRNG(box(a) / grid.GetNCells(a));
    grid.Build(box, min_width, shift);
    return grid;
}

static std::vector<int> cell_indices(const CheckerboardCellGrid& grid, const std::vector<vertex*>& verts) {
    std::vector<int> cells;
    for (vertex* v : verts) cells.push_back(grid.CellIndex(v->GetXPos(), v->GetYPos(), v->GetZPos()));
    return cells;
}

static RunResult run(int nthreads, int sweeps) {
    const double dr = 0.05;
    State* pState = make_test_state("in.dts", "top.top", nthreads, "groups.inx");
    pState->GetVoxelization()->Voxelize(pState->GetMesh()->GetActiveV());
    pState->GetEnergyCalculator()->UpdateTotalEnergy(pState->GetEnergyCalculator()->CalculateAllLocalEnergy());
    EvolveVerticesByMetropolisAlgorithmCheckerboard integrator(pState, 1, 1, dr);
    integrator.UpdateFreezGroupName("G1");
    integrator.Initialize();

    const std::vector<vertex*>& active_v = pState->GetMesh()->GetActiveV();
    const std::vector<vertex*>& frozen = pState->GetMesh()->GetGroups().at("G1");
    std::vector<Vec3D> frozen_pos;
    for (vertex* v : frozen) frozen_pos.push_back(v->GetPos());
    int left_cell = 0;
    for (int s = 0; s < sweeps; s++) {
        CheckerboardCellGrid grid = predict_grid(pState, dr);
        std::vector<int> before = cell_indices(grid, active_v);
        integrator.EvolveOneStep(s);
        if (cell_indices(grid, active_v) != before) left_cell++;
    }
    CHECK(left_cell == 0);
    for (size_t i = 0; i < frozen.size(); i++) {
        CHECK(frozen[i]->GetXPos() == frozen_pos[i](0));
        CHECK(frozen[i]->GetYPos() == frozen_pos[i](1));
        CHECK(frozen[i]->GetZPos() == frozen_pos[i](2));
    }

    RunResult result;
    for (vertex* v : active_v) {
        result.positions.push_back(v->GetXPos());
        result.positions.push_back(v->GetYPos());
        result.positions.push_back(v->GetZPos());
    }
    result.energy = pState->GetEnergyCalculator()->GetEnergy();
    result.acceptance = integrator.GetAcceptanceRate(false);
    // stored total must match a full recomputation, i.e. no deltas were lost in the commit
    pState->GetCurvatureCalculator()->Initialize();
    double E_ref = pState->GetEnergyCalculator()->CalculateAllLocalEnergy();
    CHECK_NEAR(result.energy, E_ref, 1e-8 * std::max(1.0, std::fabs(E_ref)));
    std::cout << "nt " << nthreads << ": energy " << result.energy << " (recomputed " << E_ref
              << "), acceptance " << result.acceptance << "\n";
    delete pState;
    return result;
}

void test_thread_count_independent() {
    const int sweeps = 20;
    RunResult r1 = run(1, sweeps);
    RunResult r4 = run(4, sweeps);
    CHECK(r1.acceptance > 0.1);
    CHECK(r1.positions == r4.positions);
    CHECK(r1.energy == r4.energy);
    CHECK(r1.acceptance == r4.acceptance);
}

struct AlexanderProbe : public AlexanderMoveByMetropolisAlgorithmCheckerboard {
    using AlexanderMoveByMetropolisAlgorithmCheckerboard::AlexanderMoveByMetropolisAlgorithmCheckerboard;
    double Accepted() const { return m_AcceptedMoves; }
    double Attempted() const { return m_NumberOfAttemptedMoves; }
};

struct FlipRunResult {
    std::vector<double> positions;
    std::vector<std::vector<int> > neighbours;
    double energy, reference_energy;
    double accepted, attempted;
    int left_cell, changed;
};

// replays the flip schedule's draws on a copy of the master rng
static CheckerboardCellGrid predict_flip_grid(State* pState) {
    RNG rng = *pState->GetRandomNumberGenerator();
    CheckerboardSchedule schedule;
    schedule.Draw(rng, *pState->GetMesh()->GetBox(), 3 * std::sqrt(pState->GetSimulation()->GetMaxL2()) * 1.01);
    return schedule.GetGrid();
}

static FlipRunResult run_flips(int nthreads, int sweeps) {
    State* pState = make_test_state("in.dts", "top.top", nthreads);
    pState->GetVoxelization()->Voxelize(pState->GetMesh()->GetActiveV());
    pState->GetEnergyCalculator()->UpdateTotalEnergy(pState->GetEnergyCalculator()->CalculateAllLocalEnergy());
    EvolveVerticesByMetropolisAlgorithmCheckerboard vertex_moves(pState, 1, 1, 0.05);
    vertex_moves.Initialize();
    AlexanderProbe flips(pState, 1);
    flips.Initialize();

    const Vec3D& box = *pState->GetMesh()->GetBox();
    const std::vector<links*>& edges = pState->GetMesh()->GetRightL();
    FlipRunResult result;
    result.left_cell = result.changed = 0;
    for (int s = 0; s < sweeps; s++) {
        vertex_moves.EvolveOneStep(s);
        CheckerboardCellGrid grid = predict_flip_grid(pState);
        std::vector<std::pair<vertex*, vertex*> > ends;
        std::vector<int> cells;
        for (links* e : edges) {
            ends.push_back(std::make_pair(e->GetV1(), e->GetV2()));
            Vec3D m = CanonicalMidpoint(e->GetV1(), e->GetV2(), box);
            cells.push_back(grid.CellIndex(m(0), m(1), m(2)));
        }
        flips.EvolveOneStep(s);
        for (size_t i = 0; i < edges.size(); i++) {
            if (ends[i] == std::make_pair(edges[i]->GetV1(), edges[i]->GetV2())) continue;
            result.changed++;
            Vec3D m = CanonicalMidpoint(edges[i]->GetV1(), edges[i]->GetV2(), box);
            if (grid.CellIndex(m(0), m(1), m(2)) != cells[i]) result.left_cell++;
        }
    }
    for (vertex* v : pState->GetMesh()->GetActiveV()) {
        result.positions.push_back(v->GetXPos());
        result.positions.push_back(v->GetYPos());
        result.positions.push_back(v->GetZPos());
        std::vector<int> ids;
        for (vertex* n : v->GetVNeighbourVertex()) ids.push_back(n->GetVID());
        std::sort(ids.begin(), ids.end());
        result.neighbours.push_back(ids);
    }
    result.energy = pState->GetEnergyCalculator()->GetEnergy();
    result.accepted = flips.Accepted();
    result.attempted = flips.Attempted();
    result.reference_energy = reference_energy(pState);
    std::cout << "flips nt " << nthreads << ": energy " << result.energy << " (recomputed " << result.reference_energy
              << "), flips accepted " << result.accepted << " of " << result.attempted
              << ", edges changed " << result.changed << "\n";
    delete pState;
    return result;
}

static const int g_flip_sweeps = 20;
static FlipRunResult g_flip_nt1, g_flip_nt4;

void test_flip_determinism() {
    // the lattice fixture makes most flips length-inadmissible, so require a count, not a rate
    CHECK(g_flip_nt1.accepted >= 20);
    CHECK(g_flip_nt1.positions == g_flip_nt4.positions);
    CHECK(g_flip_nt1.neighbours == g_flip_nt4.neighbours);
    CHECK(g_flip_nt1.energy == g_flip_nt4.energy);
    CHECK(g_flip_nt1.accepted == g_flip_nt4.accepted);
    CHECK(g_flip_nt1.attempted == g_flip_nt4.attempted);
}

void test_flip_confinement() {
    CHECK(g_flip_nt1.changed > 0);
    CHECK(g_flip_nt1.left_cell == 0);
    CHECK(g_flip_nt4.left_cell == 0);
}

void test_flip_energy_consistent() {
    for (const FlipRunResult* r : {&g_flip_nt1, &g_flip_nt4})
        CHECK_NEAR(r->energy, r->reference_energy, 1e-8 * std::max(1.0, std::fabs(r->reference_energy)));
}

int main() {
    test_thread_count_independent();
    g_flip_nt1 = run_flips(1, g_flip_sweeps);
    g_flip_nt4 = run_flips(4, g_flip_sweeps);
    test_flip_determinism();
    test_flip_confinement();
    test_flip_energy_consistent();
    return check_report();
}
