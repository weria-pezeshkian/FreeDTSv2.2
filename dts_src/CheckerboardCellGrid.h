#pragma once
#include "Vec3D.h"

class CheckerboardCellGrid {
public:
    void Build(const Vec3D& box, double min_width, const Vec3D& shift);
    int GetNCells(int axis) const { return m_NCells[axis]; }
    int NumberOfCells() const { return m_NCells[0] * m_NCells[1] * m_NCells[2]; }
    int CellIndex(double x, double y, double z) const;
    int Color(int cell) const;
    static constexpr int NumberOfColors = 8;

private:
    int AxisIndex(double x, int axis) const;

    int m_NCells[3] = {1, 1, 1};
    double m_Box[3] = {1.0, 1.0, 1.0};
    double m_Width[3] = {1.0, 1.0, 1.0};
    double m_Shift[3] = {0.0, 0.0, 0.0};
};
