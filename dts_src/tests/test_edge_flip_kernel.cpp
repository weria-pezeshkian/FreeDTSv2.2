#include <algorithm>
#include <random>
#include "kernel_fixtures.h"
#include "../EdgeFlipKernel.h"

static std::vector<links*> shuffled_edges(State* pState) {
    std::vector<links*> all = pState->GetMesh()->GetRightL();
    std::mt19937 gen(7);
    std::shuffle(all.begin(), all.end(), gen);
    return all;
}

static std::vector<links*> pick_edges(State* pState, int n) {
    std::vector<links*> all = shuffled_edges(pState);
    all.resize(n);
    return all;
}

// the fixtures are near-regular lattices where few quad diagonals fit [MinL, MaxL],
// so samples are counted over edges whose flipped length is admissible
static bool flipped_length_admissible(State* pState, links* e) {
    const Vec3D& L = *pState->GetMesh()->GetBox();
    Vec3D a = e->GetV3()->GetPos(), b = e->GetMirrorLink()->GetV3()->GetPos();
    double d2 = 0;
    for (int k = 0; k < 3; k++) {
        double d = b(k) - a(k);
        if (std::fabs(d) > L(k) / 2) d -= (d > 0 ? 1 : -1) * L(k);
        d2 += d * d;
    }
    return d2 >= pState->GetSimulation()->GetMinL2() && d2 <= pState->GetSimulation()->GetMaxL2();
}

struct VertexSnapshot {
    double x, y, z, energy, c1, c2, area;
    std::vector<int> neighbours;
    std::vector<links*> link_list;
    std::vector<triangle*> triangle_list;
};

static VertexSnapshot snap(vertex* v) {
    VertexSnapshot s = {v->GetXPos(), v->GetYPos(), v->GetZPos(), v->GetEnergy(),
                        v->GetP1Curvature(), v->GetP2Curvature(), v->GetArea(), {}, {}, {}};
    for (vertex* n : v->GetVNeighbourVertex()) s.neighbours.push_back(n->GetVID());
    s.link_list = v->GetVLinkList();
    s.triangle_list = v->GetVTraingleList();
    std::sort(s.neighbours.begin(), s.neighbours.end());
    std::sort(s.link_list.begin(), s.link_list.end());
    std::sort(s.triangle_list.begin(), s.triangle_list.end());
    return s;
}

static void check_same(const VertexSnapshot& a, const VertexSnapshot& b) {
    CHECK_NEAR(a.x, b.x, 1e-12);
    CHECK_NEAR(a.y, b.y, 1e-12);
    CHECK_NEAR(a.z, b.z, 1e-12);
    CHECK_NEAR(a.energy, b.energy, 1e-12);
    CHECK_NEAR(a.c1, b.c1, 1e-12);
    CHECK_NEAR(a.c2, b.c2, 1e-12);
    CHECK_NEAR(a.area, b.area, 1e-12);
    CHECK(a.neighbours == b.neighbours);
    CHECK(a.link_list == b.link_list);
    CHECK(a.triangle_list == b.triangle_list);
}

struct LinkSnapshot {
    vertex *v1, *v2, *v3;
    triangle* t;
    links *n1, *n2;
};

static LinkSnapshot snap(links* l) {
    return {l->GetV1(), l->GetV2(), l->GetV3(), l->GetTriangle(), l->GetNeighborLink1(), l->GetNeighborLink2()};
}

static void check_same(const LinkSnapshot& a, const LinkSnapshot& b) {
    CHECK(a.v1 == b.v1);
    CHECK(a.v2 == b.v2);
    CHECK(a.v3 == b.v3);
    CHECK(a.t == b.t);
    CHECK(a.n1 == b.n1);
    CHECK(a.n2 == b.n2);
}

static std::vector<vertex*> quad(links* e) {
    return {e->GetV1(), e->GetV2(), e->GetV3(), e->GetMirrorLink()->GetV3()};
}

// links whose interaction energy a flip can change are the links of inclusion-owning quad vertices
static bool has_live_interaction(links* e) {
    for (vertex* v : quad(e)) {
        if (!v->VertexOwnInclusion()) continue;
        for (links* l : v->GetVLinkList())
            if (std::fabs(l->GetIntEnergy()) > 1e-12) return true;
    }
    return false;
}

void test_midpoint_min_image() {
    Vec3D box(10, 10, 10);
    Vec3D m = MinImageMidpoint(Vec3D(9.9, 5, 5), Vec3D(0.1, 5, 5), box);
    CHECK(m(0) >= 0 && m(0) < 10);
    CHECK(std::min(m(0), 10 - m(0)) < 1e-12);
    CHECK_NEAR(m(1), 5, 1e-12);
    CHECK_NEAR(m(2), 5, 1e-12);
    m = MinImageMidpoint(Vec3D(0.1, 5, 5), Vec3D(9.9, 5, 5), box);
    CHECK(m(0) >= 0 && m(0) < 10);
    CHECK(std::min(m(0), 10 - m(0)) < 1e-12);
    m = MinImageMidpoint(Vec3D(2, 2, 2), Vec3D(4, 4, 4), box);
    CHECK_NEAR(m(0), 3, 1e-12);
    CHECK_NEAR(m(1), 3, 1e-12);
    CHECK_NEAR(m(2), 3, 1e-12);

    std::mt19937 gen(5);
    std::uniform_real_distribution<double> u(0.0, 10.0);
    for (int i = 0; i < 1000; i++) {
        Vec3D a(u(gen), u(gen), u(gen)), b(u(gen), u(gen), u(gen));
        Vec3D c = MinImageMidpoint(a, b, box);
        for (int k = 0; k < 3; k++) CHECK(c(k) >= 0 && c(k) < 10);
    }
    // the edge's own midpoint, wrapped in the box
    State* pState = make_voxelized_state("in.dts");
    Vec3D L = *pState->GetMesh()->GetBox();
    for (links* e : pick_edges(pState, 40)) {
        Vec3D c = EdgeMidpoint(e, L);
        Vec3D d = MinImageMidpoint(e->GetV1()->GetPos(), e->GetV2()->GetPos(), L);
        for (int k = 0; k < 3; k++) {
            CHECK(c(k) >= 0 && c(k) < L(k));
            CHECK(c(k) == d(k));
        }
    }
    delete pState;
}

static int check_rejected_flip_restores_state(const char* input) {
    State* pState = make_voxelized_state(input);
    VAHGlobalMeshProperties* vah = pState->GetVAHGlobalMeshProperties();
    EdgeFlipKernel kernel(pState);
    kernel.Initialize();
    std::vector<links*> links_all = all_links(pState);
    int rejected = 0, metropolis_rejected = 0, live_interaction = 0, sampled = 0, visited = 0;
    for (links* e : shuffled_edges(pState)) {
        if (sampled == 40) break;
        visited++;
        if (flipped_length_admissible(pState, e)) sampled++;
        bool done = false;
        double E_ref = 0;
        std::vector<vertex*> q;
        std::vector<VertexSnapshot> before;
        LinkSnapshot e_before, m_before;
        std::vector<double> int_before;
        bool live = false;
        // a flip and its reverse cannot both be accepted at thermal ~1 unless dE ~ 0
        for (int attempt = 0; attempt < 3 && !done; attempt++) {
            E_ref = reference_energy(pState);
            q = quad(e);
            before.clear();
            for (vertex* v : q) before.push_back(snap(v));
            e_before = snap(e);
            m_before = snap(e->GetMirrorLink());
            int_before.clear();
            for (links* l : links_all) int_before.push_back(l->GetIntEnergy());
            live = has_live_interaction(e);
            done = !kernel.TryFlip(e, 1.0 - 1e-15).accepted;
        }
        CHECK(done);
        if (!done) continue;
        rejected++;
        double V = vah->GetTotalVolume(), A = vah->GetTotalArea(), C = vah->GetTotalMeanCurvature();
        for (size_t i = 0; i < q.size(); i++) check_same(snap(q[i]), before[i]);
        check_same(snap(e), e_before);
        check_same(snap(e->GetMirrorLink()), m_before);
        for (size_t i = 0; i < links_all.size(); i++) CHECK_NEAR(links_all[i]->GetIntEnergy(), int_before[i], 1e-12);
        CHECK_NEAR(vah->GetTotalVolume(), V, 1e-12);
        CHECK_NEAR(vah->GetTotalArea(), A, 1e-12);
        CHECK_NEAR(vah->GetTotalMeanCurvature(), C, 1e-12);
        CHECK_NEAR(reference_energy(pState), E_ref, 1e-9 * std::max(1.0, std::fabs(E_ref)));
        // accepted at thermal 0 means the rejection above was energetic, not geometric
        if (kernel.TryFlip(e, 0.0).accepted) {
            metropolis_rejected++;
            if (live) live_interaction++;
        }
    }
    std::cout << input << ": metropolis rejections " << metropolis_rejected << " of " << sampled
              << " admissible (" << rejected << " edges checked), with live interaction " << live_interaction << "\n";
    CHECK(sampled == 40);
    CHECK(rejected == visited);
    CHECK(metropolis_rejected >= 10);
    delete pState;
    return live_interaction;
}

static int check_accepted_dE_matches_total_energy(const char* input) {
    State* pState = make_voxelized_state(input);
    EdgeFlipKernel kernel(pState);
    kernel.Initialize();
    int checked = 0, live_interaction = 0, sampled = 0;
    double max_diff = 0;
    for (links* e : shuffled_edges(pState)) {
        if (sampled == 40) break;
        if (!flipped_length_admissible(pState, e)) continue;
        sampled++;
        double E0 = reference_energy(pState);
        bool live = has_live_interaction(e);
        EdgeFlipOutcome out = kernel.TryFlip(e, 0.0);
        if (!out.accepted) continue;
        double E1 = reference_energy(pState);
        CHECK_NEAR(out.dE_elastic, E1 - E0, 1e-8 * std::max(1.0, std::fabs(E0)));
        max_diff = std::max(max_diff, std::fabs(out.dE_elastic - (E1 - E0)));
        checked++;
        if (live) live_interaction++;
    }
    CHECK(sampled == 40);
    std::cout << input << ": dE checked " << checked << " of " << sampled << " admissible edges, max |diff| " << max_diff
              << ", with live interaction " << live_interaction << "\n";
    CHECK(checked >= 10);
    delete pState;
    return live_interaction;
}

void test_rejected_flip_restores_state() {
    check_rejected_flip_restores_state("in.dts");
    CHECK(check_rejected_flip_restores_state("interaction.dts") > 0);
}

void test_accepted_dE_matches_total_energy() {
    check_accepted_dE_matches_total_energy("in.dts");
    CHECK(check_accepted_dE_matches_total_energy("interaction.dts") > 0);
}

void test_kernel_leaves_globals_untouched() {
    State* pState = make_voxelized_state("constraint.dts", "groups.inx");
    VAHGlobalMeshProperties* vah = pState->GetVAHGlobalMeshProperties();
    AbstractEnergy* energy = pState->GetEnergyCalculator();
    // volume stays off: the flat periodic fixture has no enclosed volume
    vah->MakeAreaActive();
    vah->MakeGlobalCurvatureActive();
    energy->UpdateTotalEnergy(energy->CalculateAllLocalEnergy());
    EdgeFlipKernel kernel(pState);
    kernel.Initialize();
    double V = vah->GetTotalVolume(), A = vah->GetTotalArea(), C = vah->GetTotalMeanCurvature();
    double E = energy->GetEnergy();
    double V0, A0, C0, V1, A1, C1;
    vah->CalculateGlobalVariables(V0, A0, C0);

    EdgeFlipOutcome out{};
    for (links* e : shuffled_edges(pState)) {
        out = kernel.TryFlip(e, 0.0);
        if (out.accepted) break;
    }
    CHECK(out.accepted);
    CHECK_NEAR(vah->GetTotalVolume(), V, 1e-14);
    CHECK_NEAR(vah->GetTotalArea(), A, 1e-14);
    CHECK_NEAR(vah->GetTotalMeanCurvature(), C, 1e-14);
    CHECK_NEAR(energy->GetEnergy(), E, 1e-14);

    vah->CalculateGlobalVariables(V1, A1, C1);
    CHECK(std::fabs(out.dArea) > 0);
    CHECK(std::fabs(out.dCurvature) > 0);
    CHECK(out.dVolume == 0);
    CHECK_NEAR(out.dArea, A1 - A0, 1e-8);
    CHECK_NEAR(out.dCurvature, C1 - C0, 1e-8);
    delete pState;
}

int main() {
    test_midpoint_min_image();
    test_rejected_flip_restores_state();
    test_accepted_dE_matches_total_energy();
    test_kernel_leaves_globals_untouched();
    return check_report();
}
