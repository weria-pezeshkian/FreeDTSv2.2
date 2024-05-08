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
    inline  std::string GetDerivedDefaultReadName() {return "MetropolisAlgorithm";}
    inline static std::string GetDefaultReadName() {return "MetropolisAlgorithm";}
    
private:
    bool FlipOneEdge(int step, links *pedge, double temp);
    bool EdgeCanBeFliped(links *pedge, double mindist2, double maxdist2);
    bool CheckFacesAfterFlip(double &minangle, links* edge);
    double SystemEnergy(); // For bug finding only; slow function (should be deleted in production code)
    std::string CurrentState();

    
private:
    double *m_pLmin2;
    double *m_pLmax2;
    double *m_pminAngle;
    double *m_pBeta;
    State *m_pState;
    Vec3D *m_pBox;
};

#endif // ALEXANDER_MOVE_BY_METROPOLIS_ALGORITHM_H
