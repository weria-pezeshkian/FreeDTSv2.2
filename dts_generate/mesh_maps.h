#include "Vec3D.h"

struct Vertex_Map {    // data structure for vertex map (not vertex object)
    double x,y,z;
    int id,domain;
};
struct Triangle_Map {    // data structure for triangle map (not triangle object)
    int id,v1,v2,v3;
};
struct Inclusion_Map {    // data structure for inclusion map (not triangle object)
    double x,y;
    int id,vid,tid;
};
struct MeshBluePrint {    // data structure for the mesh blue print
    std::vector<Vertex_Map> bvertex;       // a vector of all vertices (only the blueprint not the object) in the mesh
    std::vector<Triangle_Map> btriangle;   // a vector of all triangles (only the blueprint not the object) in the mesh
    std::vector<Inclusion_Map> binclusion; // a vector of all inclusions (only the blueprint not the object) in the mesh
    Vec3D simbox;
};
