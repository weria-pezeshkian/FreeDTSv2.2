#ifndef TriangularPrismBuilder_H
#define TriangularPrismBuilder_H
#include <array>
#include "SimDef.h"
#include "triangle.h"
#include "Vec3D.h"
#include "Tensor2.h"


/*
 * -----------------------------------------------------------------------------
 * Author: Weria Pezeshkian (weria.pezeshkian@gmail.com)
 * Date: June 2026
 *
 * Class: TriangularPrismBuilder
 *
 * Description:
 * Constructs and analyzes all possible triangular-prism-like connections
 * between two triangular faces.
 *
 * Given two triangles, the builder generates the possible vertex
 * correspondences and side-face triangulations that connect them into a
 * closed surface. The generated configurations can be used for geometric,
 * topological, or mesh-construction purposes.
 *
 * The class stores references to the two input triangles and the simulation
 * box required for geometric operations involving periodic boundary
 * conditions.
 *
 * Notes:
 * - Multiple valid prism configurations may exist for the same triangle pair.
 * - The class does not own the input triangles or simulation box.
 *
 * -----------------------------------------------------------------------------
 * License:
 * Copyright (c) Weria Pezeshkian, 2026. All rights reserved.
 * -----------------------------------------------------------------------------
 */
using triple = std::array<int, 3>;
inline triple make_triple(int a, int b, int c) {
    return {a,b,c};
}
struct TriangularPrism {
    std::array<triple, 6> tp;
};
class TriangularPrismBuilder {
public:

    /*
     * Construct a prism builder from two triangles.
     *
     * Parameters:
     *   id   - Unique identifier.
     *   t1   - First triangle.
     *   t2   - Second triangle.
     *   pbox - Simulation box used for geometric calculations.
     *   max distance 
     */
    TriangularPrismBuilder(int id,triangle* t1, triangle* t2,Vec3D* pbox, double &maxl2
    );

    ~TriangularPrismBuilder();

    //--------------------------------------------------------------------------
    // Accessors
    //--------------------------------------------------------------------------
    int GetID() const { return m_ID; }
    int GetNoTopology() const { return m_NoTopology; }

private:
void GenerateDistanceMap();  
bool CheckMap(TriangularPrism & tp);
bool CheckDistanceEdge(int i, int j);    
bool MakePrismMaps();  
private:

    //--------------------------------------------------------------------------
    // Data members
    //--------------------------------------------------------------------------
    int m_ID;
     triangle* m_T1; // Unique identifier.
     triangle* m_T2;
     Vec3D* m_pBox; // Ownership is not required.
    double m_Dist[3][3];
    int m_NoTopology;
    std::vector<TriangularPrism> m_TopologyMaps;
    double &m_MaxLength2;
    // just to keep as global
    std::array<Vec3D, 6> m_NormalVector;
    std::array<Vec3D, 6> m_Vertices;

    //std::vector<> m_ValidTopology;
    

};

#endif
