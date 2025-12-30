#include <iostream>
#include "HighTopologyStructure_ChannelPBC.h"
//#include "mesh_maps.h"
#include "Def.h"
#include "Nfunction.h"
#include "IO_Files.h"


// Constructor with initialization of bond parameters
HighTopologyStructure_ChannelPBC::HighTopologyStructure_ChannelPBC(){
    
    m_noLowerCreatedHole = 0;
    m_noUpperCreatedHole = 0;
}


// Destructor
HighTopologyStructure_ChannelPBC::~HighTopologyStructure_ChannelPBC() = default;

// Update vertices associated with the bond
void HighTopologyStructure_ChannelPBC::Generate_Shape(int genus, int Ny, Vec3D box, std::string OutputFilename) {
    
    srand (40072);

    int TopDegree = genus;
    std::vector<Vertex_Map> allV;
    std::vector<Triangle_Map> allT;
    std::vector<Inclusion_Map> allI;
    
    
    double edge_l = 1.01;
    int Nx = int (box(0)/edge_l);
    box(0) = double(Nx) * edge_l;
    if (edge_l * double(Ny) > box(1) ){
        std::cout<<"---> error: box should be larger than the Ny value \n";
    }
    double x_min = 0.5 * edge_l;
    double y_min = 0.5 * (box(1) - edge_l * double(Ny)) + 0.5 * edge_l;
    double z_min = 0.5 * box(2);
    double noise = 0.001;    // value of tinly ruoghness added to the surface

   
    int id = 0;
    
    // upper layer
    for (int j = 0; j < Ny; j++) {
    for (int i = 0; i < Nx; i++) {

            double x = i * edge_l + x_min;
            double y = j * edge_l + y_min;
            double z = z_min + edge_l/2;
            
            
            // build a vertex
            Vertex_Map v;
            v.id = id;
            v.domain = 0;
            double dx = double(rand()%1000)/1000.0;
            v.x = x + dx * noise;
            dx = double(rand()%1000)/1000.0;
            v.y = y + dx * noise;
            dx = double(rand()%1000)/1000.0;
            v.z = z + dx * noise;
            allV.push_back(v);
            
            // id for the next vertex
            id++;
            
        }
        
    }
    
    // lower layer
    for (int j = 0; j < Ny; j++) {
    for (int i = 0; i < Nx; i++) {
            
            double x = i * edge_l + x_min;
            double y = j * edge_l + y_min;
            double z = z_min - edge_l/2;
            // build a vertex
            Vertex_Map v;
            v.id = id;
            v.domain = 0;
            double dx = double(rand()%1000)/1000.0;
            v.x = x + dx * noise;
            dx = double(rand()%1000)/1000.0;
            v.y = y + dx * noise;
            dx = double(rand()%1000)/1000.0;
            v.z = z + dx * noise;
            allV.push_back(v);
            
            // id for the next vertex
            id++;
        }
        
    }
    
    
    //============================
    // Connect triangles
    //==============================
       int Tid = 0;
        bool makehole = false;
        
    for (int j = 0; j < Ny-1; j++)
    {
        for (int i = 0; i < Nx; i++)
        {

                

                makehole = MakeHole(1,i,j, Nx, Ny, genus);
                
                if(makehole == true)
                {
                      std::vector<Triangle_Map> temT = MakeTrianglesAroundHole(id, i, j, Nx, Ny);
                    for (std::vector<Triangle_Map>::iterator it = temT.begin() ; it != temT.end(); ++it)
                       allT.push_back(*it);

                    id = id + 8;
                }
                else
                {
                    int v1id = idfromij(1, i, j, Nx , Ny);
                    int v2id = idfromij(1, i, j+1, Nx , Ny);
                    int v3id = idfromij(1, i+1, j+1, Nx , Ny);
                    
                    
                    Triangle_Map T1;
                    T1.id = Tid;
                    T1.v1 = v1id;
                    T1.v2 = v2id;
                    T1.v3 = v3id;
                    allT.push_back(T1);
                    Tid++;
                    
                    v1id = idfromij(1,i,j, Nx , Ny);
                    v2id = idfromij(1,i+1,j+1, Nx , Ny);
                    v3id = idfromij(1,i+1,j, Nx , Ny);
                    Triangle_Map T2;
                    T2.id = Tid;
                    T2.v1 = v1id;
                    T2.v2 = v2id;
                    T2.v3 = v3id;
                    allT.push_back(T2);
                    Tid++;
                }


                
            }
        }
    
    // lower layer
    for (int j = 0; j<Ny-1; j++)
    {
    for (int i = 0 ; i < Nx; i++)
    {
            makehole = MakeHole(2,i,j, Nx, Ny, genus);

            if(!makehole)
            {
                
                int v1id=idfromij(2 ,i, j, Nx , Ny);
                int v3id=idfromij(2 ,i, j+1, Nx , Ny);
                int v2id=idfromij(2 ,i+1, j+1, Nx , Ny);
                
                Triangle_Map T1;
                T1.id = id;
                T1.v1 = v1id;
                T1.v2 = v2id;
                T1.v3 = v3id;
                allT.push_back(T1);
                id++;
                v1id=idfromij(2,i,j, Nx , Ny);
                v3id=idfromij(2,i+1,j+1, Nx , Ny);
                v2id=idfromij(2,i+1,j, Nx , Ny);
                Triangle_Map T2;
                T2.id = id;
                T2.v1 = v1id;
                T2.v2 = v2id;
                T2.v3 = v3id;
                allT.push_back(T2);
                id++;
            }
        }
    }
    
    // Edge Y
    for (int i = 0 ; i < Nx; i++)
    {
        
        int J = Ny - 1;
        int v1id=idfromij(1,i,J, Nx , Ny);
        int v2id=idfromij(2,i,J, Nx , Ny);
        int v3id=idfromij(2,i+1,J, Nx , Ny);
                
                Triangle_Map T1;
                T1.id = id;
                T1.v1 = v1id;
                T1.v2 = v2id;
                T1.v3 = v3id;
                allT.push_back(T1);
                id++;
        v1id=idfromij(1,i,J, Nx , Ny);
        v2id=idfromij(2,i+1,J, Nx , Ny);
        v3id=idfromij(1,i+1,J, Nx , Ny);
                Triangle_Map T2;
                T2.id = id;
                T2.v1 = v1id;
                T2.v2 = v2id;
                T2.v3 = v3id;
                allT.push_back(T2);
                id++;
    }
    for (int i = 0 ; i < Nx; i++)
    {
        
        int v1id=idfromij(1,i+1,0, Nx , Ny);
        int v2id=idfromij(2,i+1,0, Nx , Ny);
        int v3id=idfromij(2,i,0, Nx , Ny);
        Triangle_Map T1;
        T1.id = id;
        T1.v1 = v1id;
        T1.v2 = v2id;
        T1.v3 = v3id;
        allT.push_back(T1);
        id++;
        
        
        v1id=idfromij(1,i+1,0, Nx , Ny);
        v2id=idfromij(2,i,0, Nx , Ny);
        v3id=idfromij(1,i,0, Nx , Ny);
        Triangle_Map T2;
        T2.id = id;
        T2.v1 = v1id;
        T2.v2 = v2id;
        T2.v3 = v3id;
        allT.push_back(T2);
        id++;
        
    }


    
    // Mash making is finish, make a blue print and make the files
    std::cout<<" number of the vertex "<<allV.size()<<" total number of triangles "<< allT.size()<<"\n";
    MeshBluePrint BluePrint;
    BluePrint.bvertex = allV;
    BluePrint.btriangle = allT;
    BluePrint.binclusion = allI;
    BluePrint.simbox = box;
    
    std::string ext = OutputFilename.substr(OutputFilename.find_last_of(".") + 1);
    IO_Files io_files;
    if(ext == TSExt)
    {
        
        io_files.WriteQFile(OutputFilename , BluePrint);
    }
    else if(ext == TSIExt)
    {
        io_files.WriteTSI(OutputFilename , BluePrint);
    }
    else
    {
        std::cout<<"---> Error: output file with "<<ext<<" extension is not recognized.  It should have either "<<TSExt<<" or "<<TSIExt<<" extension. "<<std::endl;
    }

}
int HighTopologyStructure_ChannelPBC::idfromij(int s, int i, int j, int Nx, int Ny)
{
    
    if(i > Nx ||  j > Ny)
    {
        std::cout<<"error, this should happen \n";
    }
    if( j == Ny)
    {
        j = j%Ny;
    }
    if( i == Nx)
    {
        i = i%Nx;
    }

    int id = 0;
    
    if(s == 1) {
        
        id = Nx * j + i;
    }
    else if(s == 2){
        id = Nx * Ny + Nx * j + i;
    }

    return id;
}
bool HighTopologyStructure_ChannelPBC::MakeHole(int s,
                                                int i,
                                                int j,
                                                int Nx,
                                                int Ny,
                                                int genus)
{
    // No holes if genus is zero or negative
    
    if (genus <= 0)
        return false;

    // Avoid division by zero or meaningless grids
    if (Nx <= 2 || Ny <= 2)
        return false;

    int space = sqrt((Nx - 2) * (Ny - 2) / (genus + 1))-1;
    if (space <= 0)
        return false;

    // Bounds check (shared)
    if (i >= Nx - 1 || j >= Ny - 1)
        return false;

    bool makeHole = false;

    // Upper surface
    
    //std::cout<<space <<" "<<i<<" here 6 \n";

    if (s == 1 && m_noUpperCreatedHole < genus)
    {
        if (((i + 2) % space == 0) && ((j + 2) % space == 0)){
            makeHole = true;
        }
    }

    // Lower surface
    else if (s == 2 && m_noLowerCreatedHole < genus)
    {
        if (((i + 2) % space == 0) && ((j + 2) % space == 0))
            makeHole = true;
    }

    // Update counters
    if (makeHole)
    {
        if (s == 1){
            m_noUpperCreatedHole++;
            std::cout<< "--> neck number "<<m_noUpperCreatedHole<<" is found \n";
        }
        
        else if (s == 2)
            m_noLowerCreatedHole++;
    }

    return makeHole;
}
std::vector<Triangle_Map>  HighTopologyStructure_ChannelPBC::MakeTrianglesAroundHole(int id, int i, int j, int Nx, int Ny)
{
    
    std::vector<Triangle_Map> allT;

    int v1id=idfromij(1,i,j,Nx,Ny);
    int v2id=idfromij(1,i,j+1,Nx,Ny);
    int v3id=idfromij(2,i,j+1,Nx,Ny);
    Triangle_Map T1;
    T1.id = id;
    T1.v1 = v1id;
    T1.v2 = v2id;
    T1.v3 = v3id;
    allT.push_back(T1);
    id++;


    v1id=idfromij(1,i,j,Nx,Ny);
    v2id=idfromij(2,i,j+1,Nx,Ny);
    v3id=idfromij(2,i,j,Nx,Ny);
    T1.id = id;
    T1.v1 = v1id;
    T1.v2 = v2id;
    T1.v3 = v3id;
    allT.push_back(T1);
    id++;

    
    v1id=idfromij(1,i,j+1,Nx,Ny);
    v2id=idfromij(2,i+1,j+1,Nx,Ny);
    v3id=idfromij(2,i,j+1,Nx,Ny);
    T1.id = id;
    T1.v1 = v1id;
    T1.v2 = v2id;
    T1.v3 = v3id;
    allT.push_back(T1);
    id++;

    v1id=idfromij(1,i,j+1,Nx,Ny);
    v2id=idfromij(1,i+1,j+1,Nx,Ny);
    v3id=idfromij(2,i+1,j+1,Nx,Ny);
    T1.id = id;
    T1.v1 = v1id;
    T1.v2 = v2id;
    T1.v3 = v3id;
    allT.push_back(T1);
    id++;

    v1id=idfromij(1,i,j,Nx,Ny);
    v2id=idfromij(2,i,j,Nx,Ny);
    v3id=idfromij(2,i+1,j,Nx,Ny);
    T1.id = id;
    T1.v1 = v1id;
    T1.v2 = v2id;
    T1.v3 = v3id;
    allT.push_back(T1);
    id++;

    v1id=idfromij(1,i,j,Nx,Ny);
    v2id=idfromij(2,i+1,j,Nx,Ny);
    v3id=idfromij(1,i+1,j,Nx,Ny);
    T1.id = id;
    T1.v1 = v1id;
    T1.v2 = v2id;
    T1.v3 = v3id;
    allT.push_back(T1);
    id++;
    
    
    v1id=idfromij(1,i+1,j,Nx,Ny);
    v2id=idfromij(2,i+1,j+1,Nx,Ny);
    v3id=idfromij(1,i+1,j+1,Nx,Ny);
    T1.id = id;
    T1.v1 = v1id;
    T1.v2 = v2id;
    T1.v3 = v3id;
    allT.push_back(T1);
    id++;
    
    
    v1id=idfromij(1,i+1,j,Nx,Ny);
    v2id=idfromij(2,i+1,j,Nx,Ny);
    v3id=idfromij(2,i+1,j+1,Nx,Ny);
    T1.id = id;
    T1.v1 = v1id;
    T1.v2 = v2id;
    T1.v3 = v3id;
    allT.push_back(T1);

    return allT;
}
