#include "GlobalCouplingChange.h"
#include "State.h"

static double CouplingEnergyChange(State* pState, double dV, double dA, double dC) {
    double dE = pState->GetVolumeCoupling()->GetEnergyChange(0, 0, dA, dV);
    dE += pState->GetTotalAreaCoupling()->CalculateEnergyChange(0, dA);
    dE += pState->GetGlobalCurvature()->CalculateEnergyChange(dA, dC);
    return dE;
}
// each coupling returns E(total + d) - E(total), so the difference of two calls is the change at total + pending
double GlobalCouplingEnergyChange(State* pState, const PendingGlobalChange& pending, double dV, double dA, double dC) {
    double dE = CouplingEnergyChange(pState, pending.dV + dV, pending.dA + dA, pending.dC + dC);
    if (pending.dV != 0 || pending.dA != 0 || pending.dC != 0)
        dE -= CouplingEnergyChange(pState, pending.dV, pending.dA, pending.dC);
    return dE;
}
