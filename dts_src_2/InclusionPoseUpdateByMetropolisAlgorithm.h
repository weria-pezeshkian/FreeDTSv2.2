
#ifndef AFX_InclusionPoseUpdateByMetropolisAlgorithm_H_INCLUDED_
#define AFX_InclusionPoseUpdateByMetropolisAlgorithm_H_INCLUDED_

#include "AbstractInclusionPoseIntegrator.h"
#include "SimDef.h"

class State;
class inclusion;
class links;

/*
 * @file InclusionPoseUpdateByMetropolisAlgorithm.h
 * @brief Declaration of the InclusionPoseUpdateByMetropolisAlgorithm class.
 *
 * This file contains the declaration of the InclusionPoseUpdateByMetropolisAlgorithm class,
 * which implements an inclusion pose update algorithm based on the Metropolis algorithm.
 *
 * @author Weria
 */


class InclusionPoseUpdateByMetropolisAlgorithm : public AbstractInclusionPoseIntegrator {
public:

    InclusionPoseUpdateByMetropolisAlgorithm();
    ~InclusionPoseUpdateByMetropolisAlgorithm();


    bool Initialize(State *pState);
    bool EvolveOneStep(int step);
    inline std::string GetDerivedDefaultReadName() { return "MetropolisAlgorithm"; }
    inline static std::string GetDefaultReadName() { return "MetropolisAlgorithm"; }
    std::string CurrentState();

private:
    State *m_pState; ///< Pointer to the simulation state.
    double m_Beta;   ///< Inverse temperature factor.

private:
    /**
     * @brief Perform a Kawasaki move for inclusion pose update.
     *
     * Performs a Kawasaki move to update the pose of an inclusion.
     *
     * @param step Current step of the simulation.
     * @param p_inc Pointer to the inclusion whose pose is to be updated.
     * @param d_links Pointer to the link associated with the inclusion.
     * @param thermal Thermal fluctuation factor.
     * @return True if the move is accepted, false otherwise.
     */
    bool KawasakiMove(int step, inclusion* p_inc, links* d_links, double thermal);
    bool RotationMove(int step, inclusion* p_inc, double dx, double dy, double thermal);
};

#endif
