#if !defined(AFX_Generate_H_INCLUDED_)
#define AFX_Generate_H_INCLUDED_
#include "Def.h"
#include "Nfunction.h"
#include "Vec3D.h"
#include "mesh_maps.h"

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 This class makes initial TS files for DTS simulations. The output could be both q or tsi file format.
 currently 3 types of morphologies are generated.
 */

class Generate
{
public:
    
	Generate(std::vector <std::string> argument);
    Generate();
    ~Generate();
    
    
    
    
    // Functions to generate specific morphology. They will be extended
private:
    void Tetrahedron();               // to make a close surface
    void Cylinder();                  // to make a 1D channel (it is pbc in one dimention)
    void FlatBilayer();              // making flat bilayers
    int idfromij(int, int i, int j);
    std::vector<Triangle_Map>  MakeTrianglesAroundHole(int, int,int);
    bool MakeHole(int,int,int);
    void HighTopologyStructure();
    void HighTopologyStructurePBC();

public:
    std::vector <std::string> m_Argument;
    int m_Seed;       // seed for random number generator
    double m_MinFaceAngle;          //  minimum angle between the face (smaller will results in error), this is the value of the cos
    double m_MinVerticesDistanceSquare; //  minimum distance allowed between two vertices  (smaller will results in error)
    double m_MaxLinkLengthSquare;       //  maximum distance allowed between two nighbouring vertices  (larger will results in error)
    std::string m_OutputFilename; //  output TS file
    Vec3D m_Box;
    std::string m_Type;
    int m_N;      // number of the meshes for different morphology
    int m_Ny;      // number of v in Y direction for channel pbc high gen

    bool m_Healthy;
    std::string m_GeneralOutputFilename; //  a general file flag for specific run
    
    int m_genus;   // topological genus of the surface
    int m_noUpperCreatedHole;
    int m_noLowerCreatedHole;
    void WriteQFile(std::string filename , MeshBluePrint blueprint);
    void WriteTSI( std::string filename,  MeshBluePrint blueprint);

private:
    void ExploreArguments();         // updates variables based on the command line arguments
    void HelpMessage();              // writes a help message
    std::string m_tsiPrecision;
    
    
   

    
    
    
    // Some temporary function to be used in the morphology gen functions. They mostly return a proper vertex id based on some index numbers.
    // Nothing special about this functions and they may be replaced with a better one.
private:
    int findid(int faceno,int i, int j);
    int findidface4(int i, int j);
    int findindex(int M,int N,int i,int j);
    int CylinderIndex(int M,int N,int i,int j);

};

#endif
