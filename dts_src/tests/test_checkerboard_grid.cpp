#include <random>
#include "check.h"
#include "../CheckerboardCellGrid.h"

void test_grid_even_cells() {
    CheckerboardCellGrid grid;
    grid.Build(Vec3D(33.79, 33.79, 100), 5.3, Vec3D(0, 0, 0));
    CHECK(grid.GetNCells(0) == 6);
    CHECK(grid.GetNCells(1) == 6);
    CHECK(grid.GetNCells(2) == 18);
    CHECK(grid.NumberOfCells() == 6 * 6 * 18);
}

void test_grid_odd_rounds_down() {
    CheckerboardCellGrid grid;
    grid.Build(Vec3D(35, 35, 35), 5, Vec3D(0, 0, 0));
    for (int a = 0; a < 3; a++) CHECK(grid.GetNCells(a) == 6);
}

void test_grid_single_cell_axis() {
    CheckerboardCellGrid grid;
    grid.Build(Vec3D(30, 30, 4), 5.3, Vec3D(0, 0, 0));
    CHECK(grid.GetNCells(2) == 1);
    for (int c = 0; c < grid.NumberOfCells(); c++) CHECK((grid.Color(c) & 4) == 0);
}

void test_grid_wraps_pbc() {
    CheckerboardCellGrid grid;
    grid.Build(Vec3D(30, 30, 30), 5, Vec3D(0, 0, 0));
    CHECK(grid.CellIndex(30.0, 0, 0) == grid.CellIndex(0, 0, 0));
    CHECK(grid.CellIndex(-1e-12, 0, 0) == grid.CellIndex(29.999, 0, 0));
    std::mt19937 gen(1);
    std::uniform_real_distribution<double> dist(-30.0, 60.0);
    for (int i = 0; i < 10000; i++) {
        int id = grid.CellIndex(dist(gen), dist(gen), dist(gen));
        CHECK(id >= 0 && id < grid.NumberOfCells());
    }
}

static void check_same_color_separated(const Vec3D& box, double min_width) {
    CheckerboardCellGrid grid;
    grid.Build(box, min_width, Vec3D(0, 0, 0));
    int n[3] = {grid.GetNCells(0), grid.GetNCells(1), grid.GetNCells(2)};
    int N = grid.NumberOfCells();
    for (int a = 0; a < N; a++) {
        for (int b = a + 1; b < N; b++) {
            if (grid.Color(a) != grid.Color(b)) continue;
            int ia[3] = {a % n[0], (a / n[0]) % n[1], a / (n[0] * n[1])};
            int ib[3] = {b % n[0], (b / n[0]) % n[1], b / (n[0] * n[1])};
            bool separated = false;
            for (int k = 0; k < 3; k++) {
                if (n[k] < 2) continue;
                int d = std::abs(ia[k] - ib[k]);
                if (std::min(d, n[k] - d) >= 2) separated = true;
            }
            CHECK(separated);
        }
    }
}

void test_same_color_separated() {
    check_same_color_separated(Vec3D(30, 30, 30), 5);
    check_same_color_separated(Vec3D(30, 30, 4), 5.3);
}

void test_shift_moves_boundaries() {
    CheckerboardCellGrid grid;
    grid.Build(Vec3D(30, 30, 30), 5, Vec3D(1, 0, 0));
    CHECK(grid.CellIndex(0.5, 0, 0) == grid.CellIndex(29.5, 0, 0));
}

int main() {
    test_grid_even_cells();
    test_grid_odd_rounds_down();
    test_grid_single_cell_axis();
    test_grid_wraps_pbc();
    test_same_color_separated();
    test_shift_moves_boundaries();
    return check_report();
}
