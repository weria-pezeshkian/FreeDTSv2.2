#include "check.h"

static vertex* find_vertex_outside_groups(State* pState) {
    for (vertex* pv : pState->GetMesh()->GetActiveV())
        if (pv->GetGroupName() != "G1" && pv->GetGroupName() != "G2") return pv;
    return nullptr;
}

void test_stale_de_not_reapplied() {
    State* pState = make_test_state("constraint.dts", "top.top", 1, "groups.inx");
    AbstractApplyConstraintBetweenGroups* c = pState->GetApplyConstraintBetweenGroups();
    vertex* g1 = pState->GetMesh()->GetGroups().at("G1")[0];
    vertex* other = find_vertex_outside_groups(pState);
    CHECK(other != nullptr);
    double E0 = c->GetEnergy();
    c->CalculateEnergyChange(g1, Vec3D(0, 0, 0.1));
    c->AcceptVertexMove(other, Vec3D(0, 0, 0.1));
    CHECK_NEAR(c->GetEnergy(), E0, 1e-12);
    delete pState;
}

void test_accept_matches_prediction() {
    State* pState = make_test_state("constraint.dts", "top.top", 1, "groups.inx");
    AbstractApplyConstraintBetweenGroups* c = pState->GetApplyConstraintBetweenGroups();
    vertex* g1 = pState->GetMesh()->GetGroups().at("G1")[0];
    Vec3D Dx(0.1, -0.2, 0.3);
    double E0 = c->GetEnergy();
    double dE = c->CalculateEnergyChange(g1, Dx);
    CHECK(std::fabs(dE) > 0);
    c->AcceptVertexMove(g1, Dx);
    CHECK_NEAR(c->GetEnergy(), E0 + dE, 1e-10);
    delete pState;
}

void test_calculate_is_pure() {
    State* pState = make_test_state("constraint.dts", "top.top", 1, "groups.inx");
    AbstractApplyConstraintBetweenGroups* c = pState->GetApplyConstraintBetweenGroups();
    std::vector<vertex*>& g2 = pState->GetMesh()->GetGroups().at("G2");
    Vec3D Dx1(0, 0, 0.3), Dx2(0, 0, -0.5);
    double first = c->CalculateEnergyChange(g2[0], Dx1);
    c->CalculateEnergyChange(g2[1], Dx2);
    CHECK_NEAR(c->CalculateEnergyChange(g2[0], Dx1), first, 1e-14);
    delete pState;
}

int main() {
    test_stale_de_not_reapplied();
    test_accept_matches_prediction();
    test_calculate_is_pure();
    return check_report();
}
