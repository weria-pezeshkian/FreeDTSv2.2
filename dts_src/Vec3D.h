#ifndef VEC3D_H
#define VEC3D_H

#include <iostream>
#include <cmath>

class Vec3D {
private:
    double m_data[3];  // Array-based storage for fast indexing

public:
    Vec3D(double x = 0.0, double y = 0.0, double z = 0.0);
    
    inline double& operator()(int n) { return m_data[n]; }
    inline const double& operator()(int n) const { return m_data[n]; }

    Vec3D operator+(const Vec3D&) const;
    Vec3D operator-(const Vec3D&) const;
    Vec3D operator*(const Vec3D&) const;  // Cross Product
    Vec3D operator*(double) const;
    Vec3D& operator=(const Vec3D& other);

    void normalize();
    double norm() const;
    
    static double dot(const Vec3D&, const Vec3D&);

    bool isbad() const;
    bool isgood() const;
    void print() const;

    friend std::ostream& operator<<(std::ostream&, const Vec3D&);
};

#endif
