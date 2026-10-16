#include <algorithm>
#include "check.h"
#include "../EvolveVerticesByMetropolisAlgorithmCheckerboard.h"

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

int main() {
    test_thread_count_independent();
    return check_report();
}
