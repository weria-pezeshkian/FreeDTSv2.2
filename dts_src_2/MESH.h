#if !defined(AFX_MESH_H_INCLUDED_)
#define AFX_MESH_H_INCLUDED_
#include <map>
#include "SimDef.h"
#include "inclusion.h"
#include "triangle.h"
#include "vertex.h"
#include "links.h"
#include "CreateMashBluePrint.h"
/*
 by Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 the mesh class has all the functions for manipulating the mesh stracture
*/
class MESH
{
public:

    MESH();
    ~MESH();

    inline std::vector<vertex*>     GetActiveV()        const       {return m_pActiveV;}
    inline std::vector<triangle*>   GetActiveT()        const       {return m_pActiveT;}
    inline std::vector<links*>      GetActiveL()        const       {return m_pActiveL;}
    inline std::vector<vertex*>&           GetSurfV()                 {return m_pSurfV;}
    inline std::vector<vertex*>&           GetEdgeV()                 {return m_pEdgeV;}
    inline std::vector<links*>      GetEdgeL()          const       {return m_pEdgeL;}
    inline std::vector<links*>      GetRightL()         const       {return m_pHL;}
    inline std::vector<links*>      GetLeftL()          const       {return m_pMHL;}
    inline std::vector<inclusion*>  GetInclusion()      const       {return m_pInclusion;}
    std::map<std::string, std::vector<vertex*> > GetGroups()  const  {return m_Groups;}


    inline Vec3D                   *GetBox()                        {return m_pBox;}
    inline const bool               GetHasCrossedPBC()  const       {return m_MeshCrossedPBC;}
    std::vector <InclusionType*>    GetInclusionType()     const    {return m_pInclusionType;}
    
    
    
    
    inline void UpdateCrossedPBC(bool newValue){
        if(!m_MeshCrossedPBC)
            m_MeshCrossedPBC = newValue;
        return;
    }


public:
    void RemoveFromLinkList(links* z, std::vector<links*> &vect);
    void RemoveFromTriangleList(triangle* z, std::vector<triangle*> &vect);
    void RemoveFromVertexList(vertex* z, std::vector<vertex*> &vect);
// -- it could be optimised
    void MakeALinkGhost(links* l);
    void MakeATriangleGhost(triangle* tri);
    // this has not been completed yet
    bool Flip(links* l);
    bool MakeAVertexGhost(vertex* v);
    bool EdgeV2SurfV(vertex* v);
//----
    bool UpdateGroupFromIndexFile(std::string &indexfilename);
    void CenterMesh();    // this function centers the mesh inside the box. For broken systems it may not work nicely
    void GenerateMesh(MeshBluePrint meshblueprint);
    MeshBluePrint Convert_Mesh_2_BluePrint(MESH *mesh);
    std::vector <InclusionType*> m_pInclusionType;
    std::vector <InclusionType> m_InclusionType;
    
protected:
    std::vector<vertex*>        m_pActiveV; // all the active vertices edge + surf
    std::vector<vertex*>        m_pSurfV; // all the active vertices  surf
    std::vector<vertex*>        m_pEdgeV;  // edge
    std::vector<links*>         m_pActiveL;   // all the links
    std::vector<links*>         m_pHL;
    std::vector<links*>         m_pMHL;
    std::vector<links*>         m_pEdgeL;
    std::vector<triangle*>      m_pActiveT;
    std::vector<inclusion*>     m_pInclusion;
    std::vector<triangle*>       m_pGhostT; // Some trinagles for initial storing
    std::vector<links*>          m_pGhostL;
    std::vector<vertex*>         m_pGhostV;
    Vec3D                       *m_pBox;
    std::map<std::string, std::vector<vertex*> > m_Groups;


    
private:
    
    bool m_MeshCrossedPBC;
//========== this variables should be fully hidden from anything =======================
    std::vector<vertex>         m_Vertex;                           //                ||
    std::vector<triangle>       m_Triangle;                           //              ||
    std::vector<links>          m_Links;                           //                 ||
    std::vector<inclusion>      m_Inclusion;                           //             ||
    Vec3D                       m_Box;                           //                   ||
    std::vector<triangle>       m_GhostT; // Some trinagles for initial storing       ||
    std::vector<links>          m_GhostL;                           //                ||
    std::vector<vertex>         m_GhostV;                           //                ||
//======================================================================================

};



#endif
