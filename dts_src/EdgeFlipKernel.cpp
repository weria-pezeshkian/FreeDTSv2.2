#include <cmath>
#include <utility>
#include "EdgeFlipKernel.h"
#include "State.h"

Vec3D MinImageMidpoint(const Vec3D& a, const Vec3D& b, const Vec3D& box) {
    Vec3D mid;
    for (int k = 0; k < 3; k++) {
        double L = box(k);
        double d = b(k) - a(k);
        if (std::fabs(d) > L / 2.0) {
            if (d < 0)
                d += L;
            else
                d -= L;
        }
        double m = a(k) + 0.5 * d;
        if (m >= L) {
            m -= L;
        } else if (m < 0) {
            m += L;
        }
        // -tiny + L rounds to L
        if (m >= L)
            m = 0;
        mid(k) = m;
    }
    return mid;
}
Vec3D EdgeMidpoint(links* p_edge, const Vec3D& box) {
    return MinImageMidpoint(p_edge->GetV1()->GetPos(), p_edge->GetV2()->GetPos(), box);
}
Vec3D CanonicalMidpoint(vertex* a, vertex* b, const Vec3D& box) {
    if (b->GetVID() < a->GetVID())
        std::swap(a, b);
    return MinImageMidpoint(a->GetPos(), b->GetPos(), box);
}
EdgeFlipKernel::EdgeFlipKernel(State* pState)
    : m_pState(pState), m_pBox(nullptr), m_pBeta(nullptr), m_pDBeta(nullptr),
      m_pMinLength2(nullptr), m_pMaxLength2(nullptr), m_pMinAngle(nullptr) {
}
void EdgeFlipKernel::Initialize() {
    m_pBox = m_pState->GetMesh()->GetBox();
    AbstractSimulation* p_sim = m_pState->GetSimulation();
    m_pBeta = &p_sim->GetBeta();
    m_pDBeta = &p_sim->GetDBeta();
    m_pMinLength2 = &p_sim->GetMinL2();
    m_pMaxLength2 = &p_sim->GetMaxL2();
    m_pMinAngle = &p_sim->GetMinAngle();
}
EdgeFlipOutcome EdgeFlipKernel::TryFlip(links* p_edge, double thermal, const PendingGlobalChange& pending) {

    EdgeFlipOutcome outcome = {false, 0, 0, 0, 0};
    if (!EdgeCanBeFliped(p_edge))
        return outcome;

    vertex* v1 = p_edge->GetV1();
    vertex* v2 = p_edge->GetV2();
    vertex* v3 = p_edge->GetV3();
    vertex* v4 = p_edge->GetMirrorLink()->GetV3();

    double old_energy = v1->GetEnergy() + v2->GetEnergy() + v3->GetEnergy() + v4->GetEnergy();
    old_energy += v1->GetBindingEnergy() + v2->GetBindingEnergy() + v3->GetBindingEnergy() + v4->GetBindingEnergy();

    GetEdgesWithInteractionChange(p_edge, m_AffectedLinks);
    for (std::vector<links*>::iterator it = m_AffectedLinks.begin(); it != m_AffectedLinks.end(); ++it) {
        (*it)->Copy_InteractionEnergy();
        (*it)->Copy_VFInteractionEnergy();
        old_energy += 2 * (*it)->GetIntEnergy();
        old_energy += 2 * (*it)->GetVFIntEnergy();
    }
    //--- link-triangle contributions to the global variables (not the totals)
    double old_Tvolume = 0, old_Tarea = 0, old_Tcurvature = 0;
    double new_Tvolume = 0, new_Tarea = 0, new_Tcurvature = 0;
    VAHGlobalMeshProperties* p_VAH = m_pState->GetVAHGlobalMeshProperties();
    if (p_VAH->GetCalculateVAH()) {
        p_VAH->CalculateALinkTrianglesContributionToGlobalVariables(p_edge, old_Tvolume, old_Tarea, old_Tcurvature);
    }
    links* ml = p_edge->GetMirrorLink();
    links* l1 = p_edge->GetNeighborLink1();
    links* l2 = p_edge->GetNeighborLink2();
    links* l3 = ml->GetNeighborLink1();
    links* l4 = ml->GetNeighborLink2();
    triangle* t1 = p_edge->GetTriangle();
    triangle* t2 = ml->GetTriangle();

    t1->ConstantMesh_Copy();
    t2->ConstantMesh_Copy();
    if (!p_edge->Flip(p_edge))
        return outcome;
    t1->UpdateNormal_Area(m_pBox);
    t2->UpdateNormal_Area(m_pBox);
    if (!CheckFacesAfterFlip(p_edge)) {
        p_edge->Reverse_Flip(p_edge);
        t1->ReverseConstantMesh_Copy();
        t2->ReverseConstantMesh_Copy();
        return outcome;
    }
    // Reverse_Flip restores the lists, so only constant-mesh copies are needed
    p_edge->ConstantMesh_Copy();
    l1->ConstantMesh_Copy();
    l2->ConstantMesh_Copy();
    l3->ConstantMesh_Copy();
    l4->ConstantMesh_Copy();
    v1->ConstantMesh_Copy();
    v2->ConstantMesh_Copy();
    v3->ConstantMesh_Copy();
    v4->ConstantMesh_Copy();
    v1->Copy_VFsBindingEnergy();
    v2->Copy_VFsBindingEnergy();
    v3->Copy_VFsBindingEnergy();
    v4->Copy_VFsBindingEnergy();

    p_edge->UpdateShapeOperator(m_pBox);
    l1->UpdateShapeOperator(m_pBox);
    l2->UpdateShapeOperator(m_pBox);
    l3->UpdateShapeOperator(m_pBox);
    l4->UpdateShapeOperator(m_pBox);
    AbstractCurvature* p_curvature = m_pState->GetCurvatureCalculator();
    p_curvature->UpdateVertexCurvature(v1);
    p_curvature->UpdateVertexCurvature(v2);
    p_curvature->UpdateVertexCurvature(v3);
    p_curvature->UpdateVertexCurvature(v4);

    AbstractEnergy* p_energy = m_pState->GetEnergyCalculator();
    double new_energy = p_energy->SingleVertexEnergy(v1);
    new_energy += p_energy->SingleVertexEnergy(v2);
    new_energy += p_energy->SingleVertexEnergy(v3);
    new_energy += p_energy->SingleVertexEnergy(v4);
    new_energy += p_energy->CalculateVectorFieldMembraneBindingEnergy(v1);
    new_energy += p_energy->CalculateVectorFieldMembraneBindingEnergy(v2);
    new_energy += p_energy->CalculateVectorFieldMembraneBindingEnergy(v3);
    new_energy += p_energy->CalculateVectorFieldMembraneBindingEnergy(v4);
    for (std::vector<links*>::iterator it = m_AffectedLinks.begin(); it != m_AffectedLinks.end(); ++it) {
        new_energy += p_energy->TwoInclusionsInteractionEnergy(*it);
        for (int vf_layer = 0; vf_layer < m_pState->GetMesh()->GetNoVFPerVertex(); vf_layer++) {
            new_energy += p_energy->TwoVectorFieldInteractionEnergy(vf_layer, *it);
        }
    }
    if (p_VAH->GetCalculateVAH()) {
        p_VAH->CalculateALinkTrianglesContributionToGlobalVariables(p_edge, new_Tvolume, new_Tarea, new_Tcurvature);
    }
    double dE_global = GlobalCouplingEnergyChange(m_pState, pending, new_Tvolume - old_Tvolume, new_Tarea - old_Tarea, new_Tcurvature - old_Tcurvature);

    double diff_energy = new_energy - old_energy;
    double tot_diff_energy = diff_energy + dE_global;
    double U = (*m_pBeta) * tot_diff_energy - (*m_pDBeta);
    if (U <= 0 || exp(-U) > thermal) {
        outcome.accepted = true;
        outcome.dE_elastic = diff_energy;
        outcome.dVolume = new_Tvolume - old_Tvolume;
        outcome.dArea = new_Tarea - old_Tarea;
        outcome.dCurvature = new_Tcurvature - old_Tcurvature;
        return outcome;
    }
    //---> reject: reverse the changes
    p_edge->Reverse_Flip(p_edge);
    t1->ReverseConstantMesh_Copy();
    t2->ReverseConstantMesh_Copy();
    p_edge->ReverseConstantMesh_Copy();
    l1->ReverseConstantMesh_Copy();
    l2->ReverseConstantMesh_Copy();
    l3->ReverseConstantMesh_Copy();
    l4->ReverseConstantMesh_Copy();
    v1->ReverseConstantMesh_Copy();
    v2->ReverseConstantMesh_Copy();
    v3->ReverseConstantMesh_Copy();
    v4->ReverseConstantMesh_Copy();
    v1->Reverse_VFsBindingEnergy();
    v2->Reverse_VFsBindingEnergy();
    v3->Reverse_VFsBindingEnergy();
    v4->Reverse_VFsBindingEnergy();
    for (std::vector<links*>::iterator it = m_AffectedLinks.begin(); it != m_AffectedLinks.end(); ++it) {
        (*it)->Reverse_InteractionEnergy();
        (*it)->Reverse_VFInteractionEnergy();
    }
    return outcome;
}
bool EdgeFlipKernel::CheckFacesAfterFlip(links* p_edge) {
    links* p_medge = p_edge->GetMirrorLink();
    const double min_angle = *m_pMinAngle;
    const Vec3D& n1 = p_edge->GetTriangle()->GetNormalVector();
    const Vec3D& n2 = p_medge->GetTriangle()->GetNormalVector();
    if (Vec3D::dot(n1, n2) < min_angle)
        return false;
    if (p_edge->GetNeighborLink1()->m_LinkType == 0 &&
        Vec3D::dot(n1, p_edge->GetNeighborLink1()->GetMirrorLink()->GetTriangle()->GetNormalVector()) < min_angle)
        return false;
    if (p_edge->GetNeighborLink2()->m_LinkType == 0 &&
        Vec3D::dot(n1, p_edge->GetNeighborLink2()->GetMirrorLink()->GetTriangle()->GetNormalVector()) < min_angle)
        return false;
    if (p_medge->GetNeighborLink1()->m_LinkType == 0 &&
        Vec3D::dot(n2, p_medge->GetNeighborLink1()->GetMirrorLink()->GetTriangle()->GetNormalVector()) < min_angle)
        return false;
    if (p_medge->GetNeighborLink2()->m_LinkType == 0 &&
        Vec3D::dot(n2, p_medge->GetNeighborLink2()->GetMirrorLink()->GetTriangle()->GetNormalVector()) < min_angle)
        return false;
    return true;
}
bool EdgeFlipKernel::EdgeCanBeFliped(links* p_edge) {
    vertex* v1 = p_edge->GetV1();
    vertex* v2 = p_edge->GetV2();
    vertex* v3 = p_edge->GetV3();
    vertex* v4 = p_edge->GetMirrorLink()->GetV3();
    if (v1->GetVNeighbourVertex().size() <= 3 || v2->GetVNeighbourVertex().size() <= 3)
        return false;
    const std::vector<vertex*>& v3_neighbours = v3->GetVNeighbourVertex();
    for (std::vector<vertex*>::const_iterator it = v3_neighbours.begin(); it != v3_neighbours.end(); ++it) {
        if (*it == v4)
            return false;
    }
    double d[3] = {v4->GetVXPos() - v3->GetVXPos(), v4->GetVYPos() - v3->GetVYPos(), v4->GetVZPos() - v3->GetVZPos()};
    double dist2 = 0;
    for (int k = 0; k < 3; k++) {
        if (std::fabs(d[k]) > (*m_pBox)(k) / 2.0) {
            if (d[k] < 0)
                d[k] += (*m_pBox)(k);
            else
                d[k] -= (*m_pBox)(k);
        }
        dist2 += d[k] * d[k];
    }
    if (dist2 > *m_pMaxLength2 || dist2 < *m_pMinLength2)
        return false;
    return true;
}
void EdgeFlipKernel::GetEdgesWithInteractionChange(links* p_edge, std::vector<links*>& edges) {
    edges.clear();
    m_CandidateLinks.clear();
    vertex* quad[4] = {p_edge->GetV1(), p_edge->GetV2(), p_edge->GetV3(), p_edge->GetMirrorLink()->GetV3()};
    for (int i = 0; i < 4; i++) {
        vertex* v = quad[i];
        if (v->VertexOwnInclusion() || v->GetNumberOfVF() != 0) {
            const std::vector<links*>& v_links = v->GetVLinkList();
            m_CandidateLinks.insert(m_CandidateLinks.end(), v_links.begin(), v_links.end());
            if (v->m_VertexType == 1) {
                m_CandidateLinks.push_back(v->m_pPrecedingEdgeLink);
            }
        }
    }
    //--- drop repeated links and mirrors of links already kept
    for (std::vector<links*>::iterator it = m_CandidateLinks.begin(); it != m_CandidateLinks.end(); ++it) {
        bool should_be_added = true;
        for (std::vector<links*>::iterator it2 = edges.begin(); it2 != edges.end(); ++it2) {
            if (*it2 == *it || ((*it2)->GetMirrorFlag() && (*it2)->GetMirrorLink() == *it)) {
                should_be_added = false;
            }
        }
        if (should_be_added)
            edges.push_back(*it);
    }
}
