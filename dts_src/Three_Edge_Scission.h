#if !defined(AFX_Three_Edge_Scission_H_INCLUDED_)
#define AFX_Three_Edge_Scission_H_INCLUDED_
#include "inclusion.h"
#include "triangle.h"
#include "vertex.h"
#include "links.h"
#include "SimDef.h"
#include "MESH.h"
#include "RNG.h"
#include "Energy.h"
#include "DynamicTopology.h"

class State;
class Three_Edge_Scission : public DynamicTopology { // to use for polymorphism
public:
    Three_Edge_Scission();
    Three_Edge_Scission(int period, State *pState);
    ~Three_Edge_Scission();

    
public:
    void initialize();
    bool MCMove(double * TotalEnergy, double temp, GenerateCNTCells *pGenCNT );

private:
    std::vector<triangle>       m_GhostT; // Some trinagles for initial storing
    std::vector<links>          m_GhostL;
    MESH* m_pMESH;
    Vec3D * m_pBox;
    State *m_pState;
    Energy *m_pEnergyCalculator;
    void RemoveFromLinkList(links* z, std::vector<links*> &vect);
    void RemoveFromVertexList(vertex* z, std::vector<vertex*> &vect);
    void RemoveFromTriangleList(triangle* z, std::vector<triangle*> &vect);
    
    void AddtoLinkList(links* z, std::vector<links*> &vect);
    void AddtoVertexList(vertex* z, std::vector<vertex*> &vect);
    void AddtoTriangleList(triangle* z, std::vector<triangle*> &vect);
    bool Anglevalid4Vhole(vertex *v1, double minangle);

    std::vector<triangle> FindThreeEdgedLoop(MESH* mesh);
    double m_Beta;


public:
    std::vector<triangle*>      m_pGhostT; // Some trinagles ....
    std::vector<links*>         m_pGhostL;  // some edges for  ...





    
    
};
#endif
