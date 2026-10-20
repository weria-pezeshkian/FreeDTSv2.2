#ifndef VERTEX_MOVE_KERNEL_H
#define VERTEX_MOVE_KERNEL_H

#include <vector>
#include "SimDef.h"
#include "Vec3D.h"
#include "GlobalCouplingChange.h"

class State;
class vertex;
class links;

struct VertexMoveOutcome {
    bool   accepted;
    double dE_elastic;
    double dVolume, dArea, dCurvature;
};

// Type1's vertex move without global bookkeeping; holds scratch buffers, so use one per thread.
// one per thread; padded so their scratch buffers do not share cache lines
class alignas(128) VertexMoveKernel {
public:
    explicit VertexMoveKernel(State* pState);
    void Initialize();
    VertexMoveOutcome TryMove(vertex* p_vertex, const Vec3D& dx, double thermal, const PendingGlobalChange& pending = PendingGlobalChange());

private:
    bool VertexMoveIsFine(vertex* p_vertex, const Vec3D& dx);
    bool CheckFacesAfterAVertexMove(vertex* p_vertex);
    void GetEdgesWithInteractionChange(vertex* p_vertex, std::vector<links*>& edges);

    State* m_pState;
    Vec3D* m_pBox;
    double* m_pBeta;
    double* m_pDBeta;
    double* m_pMinLength2;
    double* m_pMaxLength2;
    double* m_pMinAngle;
    std::vector<links*> m_AffectedLinks;
    std::vector<links*> m_CandidateLinks;
};

#endif
