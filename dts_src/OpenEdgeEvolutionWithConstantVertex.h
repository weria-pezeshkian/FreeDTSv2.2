#if !defined(AFX_OpenEdgeEvolutionWithConstantVertex_H_INCLUDED_)
#define AFX_OpenEdgeEvolutionWithConstantVertex_H_INCLUDED_
#include "inclusion.h"
#include "triangle.h"
#include "vertex.h"
#include "links.h"
#include "SimDef.h"
#include "MESH.h"
#include "RNG.h"

class OpenEdgeEvolutionWithConstantVertex
{
public:
    OpenEdgeEvolutionWithConstantVertex();
    OpenEdgeEvolutionWithConstantVertex(bool state, MESH* pmesh,int rate, double lambda, double k1,double k2);
    ~OpenEdgeEvolutionWithConstantVertex();

    
public:
    bool m_State;
    double m_WholeSize;
private:
    std::vector<triangle>       m_GhostT; // Some trinagles for initial storing
    std::vector<links>          m_GhostL;
    MESH* m_pMESH;
    int m_Rate;
    double m_Lambda;
    double m_K1;
    double m_k2;
    Vec3D * m_pBox;
    
    void RemoveFromLinkList(links* z, std::vector<links*> &vect);
    void RemoveFromVertexList(vertex* z, std::vector<vertex*> &vect);
    void RemoveFromTriangleList(triangle* z, std::vector<triangle*> &vect);
    
    void AddtoLinkList(links* z, std::vector<links*> &vect);
    void AddtoVertexList(vertex* z, std::vector<vertex*> &vect);
    void AddtoTriangleList(triangle* z, std::vector<triangle*> &vect);


public:

    std::vector<triangle*>      m_pGhostT; // Some trinagles ....
    std::vector<links*>         m_pGhostL;  // some edges for  ...
    void Initialize();
    void MC_Move(RNG* rng, double lmin, double lmax, double maxangle);
    
private:
    links* CreateALink(vertex *);
    void KillALink(links *);
    bool Linkisvalid(vertex *, double lmin, double lmax, double maxangle);
    



    
    
};
#endif
