#include <algorithm>
#include <random>
#include "kernel_fixtures.h"
#include "../VertexMoveKernel.h"

static std::vector<vertex*> pick_vertices(State* pState, int n) {
    std::vector<vertex*> all = pState->GetMesh()->GetSurfV();
    std::mt19937 gen(7);
    std::shuffle(all.begin(), all.end(), gen);
    all.resize(n);
    return all;
}

static Vec3D random_dx(std::mt19937& gen, double size) {
    std::normal_distribution<double> n(0.0, 1.0);
    Vec3D d(n(gen), n(gen), n(gen));
    return d * (size / d.norm());
}

struct VertexSnapshot {
    double x, y, z, energy, c1, c2, area;
};

static VertexSnapshot snap(vertex* v) {
    return {v->GetXPos(), v->GetYPos(), v->GetZPos(), v->GetEnergy(),
            v->GetP1Curvature(), v->GetP2Curvature(), v->GetArea()};
}

static void check_same(const VertexSnapshot& a, const VertexSnapshot& b) {
    CHECK_NEAR(a.x, b.x, 1e-12);
    CHECK_NEAR(a.y, b.y, 1e-12);
    CHECK_NEAR(a.z, b.z, 1e-12);
    CHECK_NEAR(a.energy, b.energy, 1e-12);
    CHECK_NEAR(a.c1, b.c1, 1e-12);
    CHECK_NEAR(a.c2, b.c2, 1e-12);
    CHECK_NEAR(a.area, b.area, 1e-12);
}

// links whose interaction energy the move can change are the links of inclusion-owning neighbours
static bool has_live_interaction(vertex* v) {
    for (vertex* n : v->GetVNeighbourVertex()) {
        if (!n->VertexOwnInclusion()) continue;
        for (links* l : n->GetVLinkList())
            if (std::fabs(l->GetIntEnergy()) > 1e-12) return true;
    }
    return false;
}

static int check_rejected_move_restores_state(const char* input) {
    State* pState = make_voxelized_state(input);
    VAHGlobalMeshProperties* vah = pState->GetVAHGlobalMeshProperties();
    VertexMoveKernel kernel(pState);
    kernel.Initialize();
    std::vector<links*> links_all = all_links(pState);
    std::mt19937 gen(11);
    int rejected = 0, metropolis_rejected = 0, live_interaction = 0;
    for (vertex* v : pick_vertices(pState, 20)) {
        std::vector<vertex*> ring = v->GetVNeighbourVertex();
        bool done = false;
        Vec3D dx;
        double E_ref = 0;
        std::vector<VertexSnapshot> before;
        std::vector<double> int_before;
        for (int attempt = 0; attempt < 50 && !done; attempt++) {
            E_ref = reference_energy(pState);
            before = {snap(v)};
            for (vertex* n : ring) before.push_back(snap(n));
            int_before.clear();
            for (links* l : links_all) int_before.push_back(l->GetIntEnergy());
            dx = random_dx(gen, 0.05);
            done = !kernel.TryMove(v, dx, 1.0 - 1e-15).accepted;
        }
        CHECK(done);
        if (!done) continue;
        rejected++;
        if (has_live_interaction(v)) live_interaction++;
        double V = vah->GetTotalVolume(), A = vah->GetTotalArea(), C = vah->GetTotalMeanCurvature();
        check_same(snap(v), before[0]);
        for (size_t i = 0; i < ring.size(); i++) check_same(snap(ring[i]), before[i + 1]);
        for (size_t i = 0; i < links_all.size(); i++) CHECK_NEAR(links_all[i]->GetIntEnergy(), int_before[i], 1e-12);
        CHECK_NEAR(vah->GetTotalVolume(), V, 1e-12);
        CHECK_NEAR(vah->GetTotalArea(), A, 1e-12);
        CHECK_NEAR(vah->GetTotalMeanCurvature(), C, 1e-12);
        CHECK_NEAR(reference_energy(pState), E_ref, 1e-9 * std::max(1.0, std::fabs(E_ref)));
        // accepted at thermal 0 means the rejection above was energetic, not geometric
        if (kernel.TryMove(v, dx, 0.0).accepted) metropolis_rejected++;
    }
    std::cout << input << ": metropolis rejections " << metropolis_rejected << " of " << rejected
              << ", with live interaction " << live_interaction << "\n";
    CHECK(rejected == 20);
    CHECK(metropolis_rejected >= 10);
    delete pState;
    return live_interaction;
}

static int check_accepted_dE_matches_total_energy(const char* input) {
    State* pState = make_voxelized_state(input);
    VertexMoveKernel kernel(pState);
    kernel.Initialize();
    std::mt19937 gen(13);
    int checked = 0, live_interaction = 0;
    double max_diff = 0;
    for (vertex* v : pick_vertices(pState, 20)) {
        double E0 = reference_energy(pState);
        bool live = has_live_interaction(v);
        VertexMoveOutcome out = kernel.TryMove(v, random_dx(gen, 0.05), 0.0);
        if (!out.accepted) continue;
        double E1 = reference_energy(pState);
        CHECK_NEAR(out.dE_elastic, E1 - E0, 1e-8 * std::max(1.0, std::fabs(E0)));
        max_diff = std::max(max_diff, std::fabs(out.dE_elastic - (E1 - E0)));
        checked++;
        if (live) live_interaction++;
    }
    std::cout << input << ": dE checked " << checked << " vertices, max |diff| " << max_diff
              << ", with live interaction " << live_interaction << "\n";
    CHECK(checked >= 10);
    delete pState;
    return live_interaction;
}

void test_rejected_move_restores_state() {
    check_rejected_move_restores_state("in.dts");
    CHECK(check_rejected_move_restores_state("interaction.dts") > 0);
}

void test_accepted_dE_matches_total_energy() {
    check_accepted_dE_matches_total_energy("in.dts");
    CHECK(check_accepted_dE_matches_total_energy("interaction.dts") > 0);
}

void test_kernel_leaves_globals_untouched() {
    State* pState = make_voxelized_state("constraint.dts", "groups.inx");
    VAHGlobalMeshProperties* vah = pState->GetVAHGlobalMeshProperties();
    AbstractApplyConstraintBetweenGroups* c = pState->GetApplyConstraintBetweenGroups();
    // volume stays off: the flat periodic fixture has no enclosed volume
    vah->MakeAreaActive();
    vah->MakeGlobalCurvatureActive();
    VertexMoveKernel kernel(pState);
    kernel.Initialize();
    vertex* v = pState->GetMesh()->GetGroups().at("G1")[0];
    double V = vah->GetTotalVolume(), A = vah->GetTotalArea(), C = vah->GetTotalMeanCurvature();
    double Ec = c->GetEnergy();
    double V0, A0, C0, V1, A1, C1;
    vah->CalculateGlobalVariables(V0, A0, C0);

    std::mt19937 gen(17);
    VertexMoveOutcome out{};
    Vec3D dx;
    for (int attempt = 0; attempt < 50 && !out.accepted; attempt++) {
        dx = random_dx(gen, 0.05);
        out = kernel.TryMove(v, dx, 0.0);
    }
    CHECK(out.accepted);
    CHECK_NEAR(vah->GetTotalVolume(), V, 1e-14);
    CHECK_NEAR(vah->GetTotalArea(), A, 1e-14);
    CHECK_NEAR(vah->GetTotalMeanCurvature(), C, 1e-14);
    CHECK_NEAR(c->GetEnergy(), Ec, 1e-14);

    vah->CalculateGlobalVariables(V1, A1, C1);
    CHECK(std::fabs(out.dArea) > 0);
    CHECK(std::fabs(out.dCurvature) > 0);
    CHECK(out.dVolume == 0);
    CHECK_NEAR(out.dArea, A1 - A0, 1e-8);
    CHECK_NEAR(out.dCurvature, C1 - C0, 1e-8);
    delete pState;
}

int main() {
    test_rejected_move_restores_state();
    test_accepted_dE_matches_total_energy();
    test_kernel_leaves_globals_untouched();
    return check_report();
}
