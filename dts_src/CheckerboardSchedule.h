#pragma once
#include <vector>
#include "Vec3D.h"
#include "RNG.h"
#include "CheckerboardCellGrid.h"

// Per-sweep scheduling shared by the checkerboard integrators.
class CheckerboardSchedule {
public:
    // Draws, from the master RNG and in a fixed order: grid shift, sweep seed, colour permutation.
    void Draw(RNG& master, const Vec3D& box, double min_width);
    const CheckerboardCellGrid& GetGrid() const { return m_Grid; }
    int GetSweepSeed() const { return m_SweepSeed; }
    int GetColor(int i) const { return m_Colors[i]; }
    // Splits ascending cells by colour, keeping the order.
    void GroupByColor(const std::vector<int>& cells, std::vector<std::vector<int> >& by_color) const;

private:
    CheckerboardCellGrid m_Grid;
    int m_SweepSeed = 0;
    int m_Colors[CheckerboardCellGrid::NumberOfColors] = {0, 1, 2, 3, 4, 5, 6, 7};
};
