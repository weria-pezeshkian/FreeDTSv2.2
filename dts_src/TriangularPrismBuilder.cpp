#include "TriangularPrismBuilder.h"
#include "vertex.h"
#include "MESH.h"

// Constructor with initialization of TriangularPrismBuilder parameters
TriangularPrismBuilder::TriangularPrismBuilder(int id,  triangle* t1,  triangle* t2,  Vec3D *pbox, double &maxl2)
    : m_ID(id), m_T1(t1), m_T2(t2), m_pBox(pbox), m_MaxLength2(maxl2) {

    }


// Destructor
TriangularPrismBuilder::~TriangularPrismBuilder() = default;
bool TriangularPrismBuilder::CheckMap(TriangularPrism &tprism) {
    
    // check distances
    for (const auto& tri : tprism.tp)
    {
        for (int k = 0; k < 3; k++)
        {
            int i = tri[k];
            int j = tri[(k + 1) % 3];

            if (!CheckDistanceEdge(i, j))
                return false;
        }
    }
    
    // check faces
    
    

    return true;
}
bool TriangularPrismBuilder::CheckDistanceEdge(int i, int j)
{
    // same-triangle edges (optional skip)
    if ((i < 3 && j < 3) || (i >= 3 && j >= 3))
        return true;

    int ia = (i < 3) ? i : i - 3;
    int jb = (j < 3) ? j : j - 3;

    return (m_Dist[ia][jb] <= m_MaxLength2);
}
void TriangularPrismBuilder::GenerateDistanceMap(){
    
        // ---- vertex sharing rejection ----
     vertex * v1 = m_T1->GetV1();
     vertex * v2 = m_T1->GetV2();
     vertex * v3 = m_T1->GetV3();

     vertex * u1 = m_T2->GetV1();
     vertex * u2 = m_T2->GetV2();
     vertex * u3 = m_T2->GetV3();


    double dist_11 = MESH::SquareDistanceBetweenTwoVertices(v1, u1, *m_pBox);
    double dist_12 = MESH::SquareDistanceBetweenTwoVertices(v1, u2, *m_pBox);
    double dist_13 = MESH::SquareDistanceBetweenTwoVertices(v1, u3, *m_pBox);
    double dist_21 = MESH::SquareDistanceBetweenTwoVertices(v2, u1, *m_pBox);
    double dist_22 = MESH::SquareDistanceBetweenTwoVertices(v2, u2, *m_pBox);
    double dist_23 = MESH::SquareDistanceBetweenTwoVertices(v2, u3, *m_pBox);
    double dist_31 = MESH::SquareDistanceBetweenTwoVertices(v3, u1, *m_pBox);
    double dist_32 = MESH::SquareDistanceBetweenTwoVertices(v3, u2, *m_pBox);
    double dist_33 = MESH::SquareDistanceBetweenTwoVertices(v3, u3, *m_pBox);

    m_Dist[0][0] = dist_11;
    m_Dist[0][1] = dist_12;
    m_Dist[0][2] = dist_13;
    m_Dist[1][0] = dist_21;
    m_Dist[1][1] = dist_22;
    m_Dist[1][2] = dist_23;
    m_Dist[2][0] = dist_31;
    m_Dist[2][1] = dist_32;
    m_Dist[2][2] = dist_33;
    
    
    return;
}
bool TriangularPrismBuilder::MakePrismMaps()
{
    //==== map 1
    TriangularPrism map1 = {{
        make_triple(0, 3, 4),
        make_triple(0, 4, 1),
        make_triple(2, 1, 4),
        make_triple(2, 4, 5),
        make_triple(0, 2, 3),
        make_triple(2, 5, 3)
    }};

    m_TopologyMaps.push_back(map1);

    //==== map 2
    TriangularPrism map2 = {{
        make_triple(0, 3, 4),
        make_triple(0, 4, 1),
        make_triple(2, 1, 4),
        make_triple(2, 4, 5),
        make_triple(0, 5, 3),
        make_triple(0, 2, 5)
    }};

    m_TopologyMaps.push_back(map2);

    return true;
}