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
void CheckerboardSchedule::GroupByColor(const std::vector<int>& cells, std::vector<std::vector<int> >& by_color) const {
    by_color.resize(CheckerboardCellGrid::NumberOfColors);
    for (int c = 0; c < CheckerboardCellGrid::NumberOfColors; c++)
        by_color[c].clear();
    for (std::vector<int>::const_iterator it = cells.begin(); it != cells.end(); ++it)
        by_color[m_Grid.Color(*it)].push_back(*it);
}
