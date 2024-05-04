#ifndef ALEXANDER_MOVE_BY_METROPOLIS_ALGORITHM_H
#define ALEXANDER_MOVE_BY_METROPOLIS_ALGORITHM_H

#include "SimDef.h"
#include "vertex.h"
#include "triangle.h"
#include "links.h"
#include "Vec3D.h"
#include "AbstractAlexanderMove.h"

class State;

class AlexanderMoveByMetropolisAlgorithm : public AbstractAlexanderMove {
public:
    AlexanderMoveByMetropolisAlgorithm();
    ~AlexanderMoveByMetropolisAlgorithm();

    bool Initialize(State *pState);
    bool EvolveOneStep(int step);

private:
    int EvolveOneVertex(int step, vertex *pvertex, double dx, double dy, double dz, double temp);
    bool VertexMoveIsFine(vertex* pvertex, double dx, double dy, double dz, double mindist2, double maxdist2);
    bool CheckFacesAfterAVertexMove(double &minangle, vertex* p_vertex);
    double SystemEnergy(); // For bug finding only; slow function (should be deleted in production code)

    inline  std::string GetDerivedDefaultReadName() {return "MetropolisAlgorithm";}
    inline static std::string GetDefaultReadName() {return "MetropolisAlgorithm";}
    
private:
    double *m_pLmin2;
    double *m_pLmax2;
    double *m_pminAngle;
    double *m_pBeta;
    double *m_RateOfVMovePerStep;
    State *m_pState;
    Vec3D *m_pBox;
};

#endif // ALEXANDER_MOVE_BY_METROPOLIS_ALGORITHM_H
