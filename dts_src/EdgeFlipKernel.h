#ifndef EDGE_FLIP_KERNEL_H
#define EDGE_FLIP_KERNEL_H

#include <vector>
#include "SimDef.h"
#include "Vec3D.h"

class State;
class links;
class vertex;

struct EdgeFlipOutcome {
    bool   accepted;
    double dE_elastic;
    double dVolume, dArea, dCurvature;
};

// the serial Alexander move without global bookkeeping; holds scratch buffers, so use one per thread.
// one per thread; padded so their scratch buffers do not share cache lines
class alignas(128) EdgeFlipKernel {
public:
    explicit EdgeFlipKernel(State* pState);
    void Initialize();
    EdgeFlipOutcome TryFlip(links* p_edge, double thermal);

private:
    bool EdgeCanBeFliped(links* p_edge);
    bool CheckFacesAfterFlip(links* p_edge);
    void GetEdgesWithInteractionChange(links* p_edge, std::vector<links*>& edges);

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

// a, b inside the box; the result is wrapped into [0, L)
Vec3D MinImageMidpoint(const Vec3D& a, const Vec3D& b, const Vec3D& box);
Vec3D EdgeMidpoint(links* p_edge, const Vec3D& box);
// lower vertex id first, so an edge and its reverse get bitwise-equal midpoints
Vec3D CanonicalMidpoint(vertex* a, vertex* b, const Vec3D& box);

#endif
