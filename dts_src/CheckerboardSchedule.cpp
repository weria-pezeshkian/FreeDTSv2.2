#include <climits>
#include <utility>
#include "CheckerboardSchedule.h"

void CheckerboardSchedule::Draw(RNG& master, const Vec3D& box, double min_width) {
    m_Grid.Build(box, min_width, Vec3D(0, 0, 0));
    Vec3D shift;
    for (int a = 0; a < 3; a++)
        shift(a) = master.UniformRNG(box(a) / m_Grid.GetNCells(a));
    m_Grid.Build(box, min_width, shift);

    m_SweepSeed = master.IntRNG(INT_MAX);
    for (int i = 0; i < CheckerboardCellGrid::NumberOfColors; i++)
        m_Colors[i] = i;
    for (int i = CheckerboardCellGrid::NumberOfColors - 1; i > 0; i--)
        std::swap(m_Colors[i], m_Colors[master.IntRNG(i + 1)]);
}
void CheckerboardSchedule::ActiveCells(int color, const std::vector<int>& cell_sizes, std::vector<int>& active) const {
    active.clear();
    int n_cells = m_Grid.NumberOfCells();
    for (int cell = 0; cell < n_cells; cell++) {
        if (m_Grid.Color(cell) == color && cell_sizes[cell] > 0)
            active.push_back(cell);
    }
}
