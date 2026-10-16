#include <cmath>
#include "CheckerboardCellGrid.h"

void CheckerboardCellGrid::Build(const Vec3D& box, double min_width, const Vec3D& shift) {
    for (int a = 0; a < 3; a++) {
        int n = (int)std::floor(box(a) / min_width);
        if (n >= 2) n -= n % 2;  // even count keeps colouring consistent across the periodic wrap
        if (n < 1) n = 1;
        m_NCells[a] = n;
        m_Box[a] = box(a);
        m_Width[a] = box(a) / n;
        m_Shift[a] = shift(a);
    }
}

int CheckerboardCellGrid::AxisIndex(double x, int axis) const {
    double L = m_Box[axis];
    double w = std::fmod(x - m_Shift[axis], L);
    if (w < 0) w += L;
    int i = (int)std::floor(w / m_Width[axis]);
    if (i >= m_NCells[axis]) i = m_NCells[axis] - 1;
    if (i < 0) i = 0;
    return i;
}

int CheckerboardCellGrid::CellIndex(double x, double y, double z) const {
    int ix = AxisIndex(x, 0);
    int iy = AxisIndex(y, 1);
    int iz = AxisIndex(z, 2);
    return ix + m_NCells[0] * (iy + m_NCells[1] * iz);
}

int CheckerboardCellGrid::Color(int cell) const {
    int ix = cell % m_NCells[0];
    int iy = (cell / m_NCells[0]) % m_NCells[1];
    int iz = cell / (m_NCells[0] * m_NCells[1]);
    int color = 0;
    if (m_NCells[0] > 1) color |= (ix & 1);
    if (m_NCells[1] > 1) color |= (iy & 1) << 1;
    if (m_NCells[2] > 1) color |= (iz & 1) << 2;
    return color;
}
