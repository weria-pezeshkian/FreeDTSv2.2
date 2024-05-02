#ifndef EVOLVE_VERTICES_BY_METROPOLIS_ALGORITHM_H
#define EVOLVE_VERTICES_BY_METROPOLIS_ALGORITHM_H


#include "SimDef.h"
#include "vertex.h"
#include "triangle.h"
#include "links.h"
#include "Vec3D.h"
#include "AbstractVertexPositionIntegrator.h"
class State;
class EvolveVerticesByMetropolisAlgorithm : public AbstractVertexPositionIntegrator {
public:
    EvolveVerticesByMetropolisAlgorithm();
    ~EvolveVerticesByMetropolisAlgorithm();
    bool Initialize(State *pState);
    bool EvolveOneStep(int step);

private:
    int EvolveOneVertex(int step, vertex *pvertex, double dx, double dy, double dz,double temp);
    bool VertexMoveIsFine(vertex* pvertex, double dx,double dy, double dz,  double mindist2, double maxdist2);
    bool CheckFacesAfterAVertexMove(double &minangle, vertex* p_vertex);

    double  SystemEnergy();  // it is for bug finding only; slow function, this is for development time, should be deleted
private:
    double *m_pLmin2;
    double *m_pLmax2;
    double *m_pminAngle;
    double *m_pBeta;
    double *m_RateOfVMovePerStep;
    State *m_pState;
    Vec3D *m_pBox;


};


#endif
