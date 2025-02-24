#ifndef TENSOR2_H
#define TENSOR2_H

#include "Vec3D.h"
#include <iostream>

class Tensor2 {
public:
    // Constructors
    Tensor2(const Vec3D& x, const Vec3D& y, const Vec3D& z);
    Tensor2(); // Default constructor
    Tensor2(char t); // Constructor for special matrices (Identity or Zero)
    ~Tensor2(); // Destructor

    // Access functions
    double at(int n, int m) const; // Get element at (n, m)
    void put(int n, int m, double s); // Set element at (n, m)
    double& operator()(int n, int m); // Overloaded operator for access

    // Matrix operations
    Vec3D operator*(const Vec3D& A) const; // Matrix-vector multiplication
    Tensor2 operator+(const Tensor2& M) const; // Matrix addition
    Tensor2 operator-(const Tensor2& M) const; // Matrix subtraction
    Tensor2 operator*(const Tensor2& M) const; // Matrix multiplication
    Tensor2 operator*(double x) const; // Scalar multiplication

    // Transpose function
    Tensor2 Transpose() const; // Transpose the matrix

    // Function to create a tensor from a vector
    Tensor2 makeTen(Vec3D X); // Create tensor from a vector

private:
    double m_matrix[3][3]; // 3x3 matrix representation
};

#endif // TENSOR2_H
