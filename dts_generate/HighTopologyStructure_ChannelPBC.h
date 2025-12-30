#ifndef HighTopologyStructure_ChannelPBC_H
#define HighTopologyStructure_ChannelPBC_H

#include "Vec3D.h"

struct Triangle_Map;
class HighTopologyStructure_ChannelPBC {
public:
    // Constructors and Destructor
    HighTopologyStructure_ChannelPBC();
    ~HighTopologyStructure_ChannelPBC();
    // Public Member Functions
    void Generate_Shape(int genus, int Ny, Vec3D box, std::string OutputFilename);


private:
    int idfromij(int s, int i, int j, int Nx, int Ny);
    bool MakeHole(int s,int i, int j, int Nx, int Ny, int genus);
    std::vector<Triangle_Map>  MakeTrianglesAroundHole(int id, int i, int j, int Nx, int Ny);

    int m_noLowerCreatedHole;
    int m_noUpperCreatedHole;

};

#endif
