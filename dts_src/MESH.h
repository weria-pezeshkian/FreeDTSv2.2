#if !defined(AFX_MESH_H_INCLUDED_)
#define AFX_MESH_H_INCLUDED_
#include "inclusion.h"
#include "triangle.h"
#include "vertex.h"
#include "links.h"
#include "SimDef.h"
#include "CreateMashBluePrint.h"

class MESH
{
public:

    MESH();
    ~MESH();

    inline double GetMinLength()                        const      {return m_MinLength;}
    inline double GetMaxLength()                        const      {return m_MaxLength;}
    inline double GetMinAngle()                         const      {return m_MinAngle;}
    inline std::vector<vertex*>     GetActiveV()        const      {return m_pActiveV;}
    inline std::vector<triangle*>   GetActiveT()        const      {return m_pActiveT;}
    inline std::vector<links*>      GetActiveL()        const      {return m_pActiveL;}

    

public:
    // this has not been completed yet
    bool MakeALinkGhost(links* l);
    bool MakeAVertexGhost(links* l);
    bool MakeATriangleGhost(links* l);
    

public:
    void GenerateMesh(MeshBluePrint meshblueprint);
    MeshBluePrint Convert_Mesh_2_BluePrint(MESH *mesh);
    void  CenterMesh();

    
    
    std::vector <InclusionType> m_InclusionType;
    std::vector <InclusionType*> m_pInclusionType;
    

    
    std::vector<vertex*>        m_pActiveV; // all the active vertices edge + surf
    std::vector<vertex*>        m_pSurfV; // all the active vertices  surf
    std::vector<vertex*>        m_pEdgeV;  // edge
    
    
    std::vector<links*>         m_pActiveL;   // all the links
    std::vector<links*>         m_pHL;
    std::vector<links*>         m_pMHL;
    std::vector<links*>         m_pEdgeL;
    
    
    std::vector<triangle*>      m_pActiveT;

    std::vector<inclusion*>     m_pInclusion;
    Vec3D                       *m_pBox;
    

    
    
private:
    std::vector<vertex>         m_Vertex;
    std::vector<triangle>       m_Triangle;
    std::vector<links>          m_Links;
    std::vector<inclusion>      m_Inclusion;
    Vec3D                       m_Box;
    std::vector<triangle>       m_GhostT; // Some trinagles for initial storing
    std::vector<links>          m_GhostL;
    std::vector<vertex>         m_GhostV;
    std::vector<triangle*>       m_pGhostT; // Some trinagles for initial storing
    std::vector<links*>          m_pGhostL;
    std::vector<vertex*>         m_pGhostV;
    double m_MinLength;
    double m_MaxLength;
    double m_MinAngle;
    
};



#endif
