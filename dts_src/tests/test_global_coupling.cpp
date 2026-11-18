#include "check.h"
#include "../GlobalCouplingChange.h"

static double plain_change(State* pState, double dV, double dA, double dC) {
    return GlobalCouplingEnergyChange(pState, PendingGlobalChange(), dV, dA, dC);
}

int main() {
    State* pState = make_test_state("coupling.dts", "top.top", 1);
    VAHGlobalMeshProperties* vah = pState->GetVAHGlobalMeshProperties();
    CHECK(vah->GetCalculateVAH());
    const double changes[][3] = {{0.3, -0.2, 0.05}, {-1.1, 0.4, -0.3}, {0.0, 0.0, 0.0}};
    const PendingGlobalChange pendings[] = {{0, 0, 0}, {2.5, -1.5, 0.7}, {-4.0, 3.0, -1.2}};
    for (const PendingGlobalChange& p : pendings) {
        for (const auto& c : changes) {
            double with_pending = GlobalCouplingEnergyChange(pState, p, c[0], c[1], c[2]);
            vah->Add2Volume(p.dV);
            vah->Add2TotalArea(p.dA);
            vah->Add2GlobalCurvature(p.dC);
            double committed = plain_change(pState, c[0], c[1], c[2]);
            vah->Add2Volume(-p.dV);
            vah->Add2TotalArea(-p.dA);
            vah->Add2GlobalCurvature(-p.dC);
            CHECK_NEAR(with_pending, committed, 1e-8 * (1 + std::fabs(committed)));
        }
    }
    // the pending change must matter, or the check above is vacuous
    CHECK(std::fabs(GlobalCouplingEnergyChange(pState, pendings[1], 0.3, -0.2, 0.05) - plain_change(pState, 0.3, -0.2, 0.05)) > 1e-6);
    return check_report();
}
