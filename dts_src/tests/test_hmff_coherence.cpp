#include <algorithm>
#include <cstring>
#include <fstream>
#include <random>
#include "kernel_fixtures.h"
#include "../VertexMoveKernel.h"
#include "../EdgeFlipKernel.h"

// a density slab around the membrane (z ~ 50) with in-plane structure and voxel noise, so neighbours see different energies
static void write_map(const char* path) {
    const int nx = 34, ny = 34, nz = 16;
    char header[1024];
    std::memset(header, 0, sizeof(header));
    int32_t dims[4] = {nx, ny, nz, 2};
    float cell[3] = {nx, ny, nz};
    std::memcpy(header, dims, sizeof(dims));
    std::memcpy(header + 40, cell, sizeof(cell));
    std::memcpy(header + 208, "MAP ", 4);
    std::vector<float> data(nx * ny * nz);
    std::mt19937 gen(3);
    std::uniform_real_distribution<float> noise(0.0f, 0.3f);
    for (int z = 0; z < nz; z++)
        for (int y = 0; y < ny; y++)
            for (int x = 0; x < nx; x++) {
                double dz = z + 42 - 50.0;
                double plane = 1 + 0.5 * std::sin(2 * M_PI * x / 8) * std::cos(2 * M_PI * y / 11);
                data[x + nx * (y + ny * z)] = std::exp(-dz * dz / 8) * plane + noise(gen);
            }
    std::ofstream out(path, std::ios::binary);
    out.write(header, sizeof(header));
    out.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(float));
}

static std::vector<vertex*> pick_vertices(State* pState, int n) {
    std::vector<vertex*> all = pState->GetMesh()->GetSurfV();
    std::mt19937 gen(7);
    std::shuffle(all.begin(), all.end(), gen);
    all.resize(n);
    return all;
}

static std::vector<links*> shuffled_edges(State* pState) {
    std::vector<links*> all = pState->GetMesh()->GetRightL();
    std::mt19937 gen(7);
    std::shuffle(all.begin(), all.end(), gen);
    return all;
}

void test_coherence_contributes() {
    State* with = make_voxelized_state("hmff.dts");
    State* without = make_voxelized_state("hmff_nocoh.dts");
    CHECK(std::fabs(reference_energy(with) - reference_energy(without)) > 1.0);
    delete with;
    delete without;
}

// the coherence term couples each vertex to its 1-ring, so a move changes the neighbours' energies too
void test_vertex_move_dE_with_coherence() {
    State* pState = make_voxelized_state("hmff.dts");
    VertexMoveKernel kernel(pState);
    kernel.Initialize();
    std::mt19937 gen(13);
    std::normal_distribution<double> n(0.0, 1.0);
    int checked = 0;
    for (vertex* v : pick_vertices(pState, 40)) {
        Vec3D dx(n(gen), n(gen), n(gen));
        dx = dx * (0.1 / dx.norm());
        double E0 = reference_energy(pState);
        VertexMoveOutcome out = kernel.TryMove(v, dx, 0.0);
        if (!out.accepted) continue;
        CHECK_NEAR(out.dE_elastic, reference_energy(pState) - E0, 1e-8 * std::max(1.0, std::fabs(E0)));
        checked++;
    }
    CHECK(checked >= 10);
    delete pState;
}

// a flip changes the 1-rings of its four vertices
void test_edge_flip_dE_with_coherence() {
    State* pState = make_voxelized_state("hmff.dts");
    EdgeFlipKernel kernel(pState);
    kernel.Initialize();
    int checked = 0;
    for (links* e : shuffled_edges(pState)) {
        if (checked == 20) break;
        double E0 = reference_energy(pState);
        EdgeFlipOutcome out = kernel.TryFlip(e, 0.0);
        if (!out.accepted) continue;
        CHECK_NEAR(out.dE_elastic, reference_energy(pState) - E0, 1e-8 * std::max(1.0, std::fabs(E0)));
        checked++;
    }
    CHECK(checked >= 10);
    delete pState;
}

int main() {
    write_map("hmff.mrc");
    test_coherence_contributes();
    test_vertex_move_dE_with_coherence();
    test_edge_flip_dE_with_coherence();
    return check_report();
}
