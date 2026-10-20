#ifndef GLOBAL_COUPLING_CHANGE_H
#define GLOBAL_COUPLING_CHANGE_H

class State;

// changes to the global totals that are accepted but not yet committed
struct PendingGlobalChange {
    double dV, dA, dC;
};

// energy change of the global couplings for a change on top of the pending one
double GlobalCouplingEnergyChange(State* pState, const PendingGlobalChange& pending, double dV, double dA, double dC);

#endif
