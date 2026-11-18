#include <algorithm>
#include "kernel_fixtures.h"
#include "../EvolveVerticesByMetropolisAlgorithmCheckerboard.h"
#include "../AlexanderMoveByMetropolisAlgorithmCheckerboard.h"

// With symmetric proposals, detailed balance implies p(-U) = exp(-U) p(U) for the exponents of
// proposed moves in equilibrium. Per shell B > 0 the expected accepted uphill flux sum_{U in B} exp(-U)
// therefore equals the downhill count N(-B); unlike a ratio of histograms this has no binning bias.
struct ShellBalance {
    double max_z;
    long n_shells;
};

static ShellBalance shell_balance(const std::vector<double>& U, const char* label) {
    const double width = 0.5;
    const int n_bins = 8;
    std::vector<double> up(n_bins, 0), up2(n_bins, 0), down(n_bins, 0);
    for (double u : U) {
        int b = (int)std::floor(std::fabs(u) / width);
        if (u == 0 || b >= n_bins)
            continue;
        if (u > 0) {
            up[b] += std::exp(-u);
            up2[b] += std::exp(-2 * u);
        } else {
            down[b] += 1;
        }
    }
    ShellBalance result = {0, 0};
    std::cout << label << ": " << U.size() << " proposals\n   shell      uphill flux  downhill   z\n";
    for (int b = 0; b < n_bins; b++) {
        if (down[b] < 200)
            continue;
        double z = (up[b] - down[b]) / std::sqrt(up2[b] + down[b]);
        std::cout << "  [" << b * width << "," << (b + 1) * width << ")  " << up[b] << "  " << down[b] << "  " << z << "\n";
        result.max_z = std::max(result.max_z, std::fabs(z));
        result.n_shells++;
    }
    return result;
}

struct Logs {
    std::vector<double> vertex_U, flip_U;
};

static Logs run(const char* input, int nthreads, int burn_in, int sweeps) {
    State* pState = make_test_state(input, "top.top", nthreads);
    pState->GetVoxelization()->Voxelize(pState->GetMesh()->GetActiveV());
    pState->GetEnergyCalculator()->UpdateTotalEnergy(pState->GetEnergyCalculator()->CalculateAllLocalEnergy());
    EvolveVerticesByMetropolisAlgorithmCheckerboard vertex_moves(pState, 1, 1, 0.05);
    vertex_moves.Initialize();
    AlexanderMoveByMetropolisAlgorithmCheckerboard flips(pState, 1);
    flips.Initialize();

    Logs logs;
    for (int s = 0; s < burn_in + sweeps; s++) {
        if (s == burn_in) {
            vertex_moves.SetProposalLog(&logs.vertex_U);
            flips.SetProposalLog(&logs.flip_U);
        }
        vertex_moves.EvolveOneStep(s);
        flips.EvolveOneStep(s);
    }
    delete pState;
    return logs;
}

void test_local_energy_balance() {
    Logs logs = run("in.dts", 4, 200, 2000);
    ShellBalance v = shell_balance(logs.vertex_U, "vertex moves, in.dts");
    ShellBalance f = shell_balance(logs.flip_U, "edge flips, in.dts");
    CHECK(v.n_shells >= 3);
    CHECK(f.n_shells >= 3);
    // z assumes independent proposals; correlations inflate it, hence the loose bound
    CHECK(v.max_z < 4);
    CHECK(f.max_z < 4);
}

int main() {
    test_local_energy_balance();
    return check_report();
}
