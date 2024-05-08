#ifndef EVOLVE_VERTICES_BY_METROPOLIS_ALGORITHM_H
#define EVOLVE_VERTICES_BY_METROPOLIS_ALGORITHM_H


#include "SimDef.h"
#include "MESH.h"
#include "Vec3D.h"
#include "AbstractVertexPositionIntegrator.h"
#include "AbstractSimulation.h"

class State;
class EvolveVerticesByMetropolisAlgorithm : public AbstractVertexPositionIntegrator, public MESH, public AbstractSimulation {
public:
    EvolveVerticesByMetropolisAlgorithm(State *pState);
    ~EvolveVerticesByMetropolisAlgorithm();
    void Initialize();
    bool EvolveOneStep(int step);
    std::string CurrentState();

private:
    bool EvolveOneVertex(int step, vertex *pvertex, double dx, double dy, double dz,double temp);
    bool VertexMoveIsFine(vertex* pvertex, double dx,double dy, double dz,  double mindist2, double maxdist2);
    bool CheckFacesAfterAVertexMove(double &minangle, vertex* p_vertex);
    inline  std::string GetDerivedDefaultReadName() {return "MetropolisAlgorithm";}
    inline static std::string GetDefaultReadName() {return "MetropolisAlgorithm";}

    double  SystemEnergy();  // it is for bug finding only; slow function, this is for development time, should be deleted
    State *m_pState;
    

    
private:
    bool do_Simulation(){
        std::cout<<" ---> error, 999o1o this should have been called \n";
        return false;
    }

};


#endif
