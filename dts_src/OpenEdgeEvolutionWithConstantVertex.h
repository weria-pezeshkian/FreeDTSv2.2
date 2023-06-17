#if !defined(AFX_OpenEdgeEvolutionWithConstantVertex_H_INCLUDED_)
#define AFX_OpenEdgeEvolutionWithConstantVertex_H_INCLUDED_
#include "inclusion.h"
#include "triangle.h"
#include "vertex.h"
#include "links.h"
#include "SimDef.h"
#include "MESH.h"
class OpenEdgeEvolutionWithConstantVertex
{
public:
    OpenEdgeEvolutionWithConstantVertex();
    OpenEdgeEvolutionWithConstantVertex(bool state,MESH* pmesh);    
    ~OpenEdgeEvolutionWithConstantVertex();

private:
    std::vector<triangle>       m_GhostT; // Some trinagles for initial storing
    std::vector<links>          m_GhostL;
    MESH* m_pMESH;
public:

    std::vector<triangle*>      m_pGhostT; // Some trinagles ....
    std::vector<links*>         m_pGhostL;  // some edges for  ...
    
    
    void MC_Move();
    

};
#endif
