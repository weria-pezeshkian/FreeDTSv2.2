#include <cmath>
#include "VertexMoveKernel.h"
#include "State.h"

VertexMoveKernel::VertexMoveKernel(State* pState)
    : m_pState(pState), m_pBox(nullptr), m_pBeta(nullptr), m_pDBeta(nullptr),
      m_pMinLength2(nullptr), m_pMaxLength2(nullptr), m_pMinAngle(nullptr) {
}
void VertexMoveKernel::Initialize() {
    m_pBox = m_pState->GetMesh()->GetBox();
    AbstractSimulation* p_sim = m_pState->GetSimulation();
    m_pBeta = &p_sim->GetBeta();
    m_pDBeta = &p_sim->GetDBeta();
    m_pMinLength2 = &p_sim->GetMinL2();
    m_pMaxLength2 = &p_sim->GetMaxL2();
    m_pMinAngle = &p_sim->GetMinAngle();
}
VertexMoveOutcome VertexMoveKernel::TryMove(vertex* p_vertex, const Vec3D& dx, double thermal) {

    VertexMoveOutcome outcome = {false, 0, 0, 0, 0};
    if (!VertexMoveIsFine(p_vertex, dx))
        return outcome;

    //--- obtain vertices energy terms and make copies
    double old_energy = p_vertex->GetEnergy();
    old_energy += p_vertex->GetBindingEnergy();
    p_vertex->ConstantMesh_Copy();
    p_vertex->Copy_VFsBindingEnergy();
    const std::vector<vertex*>& v_neighbours = p_vertex->GetVNeighbourVertex();
    for (std::vector<vertex*>::const_iterator it = v_neighbours.begin(); it != v_neighbours.end(); ++it) {
        (*it)->ConstantMesh_Copy();
        old_energy += (*it)->GetEnergy();
        old_energy += (*it)->GetBindingEnergy();
        (*it)->Copy_VFsBindingEnergy();
    }
    const std::vector<triangle*>& v_triangles = p_vertex->GetVTraingleList();
    for (std::vector<triangle*>::const_iterator it = v_triangles.begin(); it != v_triangles.end(); ++it) {
        (*it)->ConstantMesh_Copy();
    }
    const std::vector<links*>& v_links = p_vertex->GetVLinkList();
    for (std::vector<links*>::const_iterator it = v_links.begin(); it != v_links.end(); ++it) {
        (*it)->ConstantMesh_Copy();
        (*it)->GetNeighborLink1()->ConstantMesh_Copy();
    }
    if (p_vertex->GetVertexType() == 1) {
        p_vertex->GetPrecedingEdgeLink()->ConstantMesh_Copy();
    }
    GetEdgesWithInteractionChange(p_vertex, m_AffectedLinks);
    for (std::vector<links*>::iterator it = m_AffectedLinks.begin(); it != m_AffectedLinks.end(); ++it) {
        (*it)->Copy_InteractionEnergy();
        (*it)->Copy_VFInteractionEnergy();
        old_energy += 2 * (*it)->GetIntEnergy();
        old_energy += 2 * (*it)->GetVFIntEnergy();
    }
    //--- ring contributions to the global variables (not the totals)
    double old_Tvolume = 0, old_Tarea = 0, old_Tcurvature = 0;
    double new_Tvolume = 0, new_Tarea = 0, new_Tcurvature = 0;
    VAHGlobalMeshProperties* p_VAH = m_pState->GetVAHGlobalMeshProperties();
    if (p_VAH->GetCalculateVAH()) {
        p_VAH->CalculateAVertexRingContributionToGlobalVariables(p_vertex, old_Tvolume, old_Tarea, old_Tcurvature);
    }
    double dE_force_from_inc = m_pState->GetForceonVerticesfromInclusions()->Energy_of_Force(p_vertex, dx);
    double dE_force_from_vector_fields = m_pState->GetForceonVerticesfromVectorFields()->Energy_of_Force(p_vertex, dx);
    double dE_force_on_vertex = m_pState->GetForceonVertices()->Energy_of_Force(p_vertex, dx);

    //---> move the vertex
    p_vertex->PositionPlus(dx(0), dx(1), dx(2));
    for (std::vector<triangle*>::const_iterator it = v_triangles.begin(); it != v_triangles.end(); ++it) {
        (*it)->UpdateNormal_Area(m_pBox);
    }
    if (!CheckFacesAfterAVertexMove(p_vertex)) {
        for (std::vector<triangle*>::const_iterator it = v_triangles.begin(); it != v_triangles.end(); ++it) {
            (*it)->ReverseConstantMesh_Copy();
        }
        p_vertex->PositionPlus(-dx(0), -dx(1), -dx(2));
        return outcome;
    }
    for (std::vector<links*>::const_iterator it = v_links.begin(); it != v_links.end(); ++it) {
        (*it)->UpdateShapeOperator(m_pBox);
        (*it)->GetNeighborLink1()->UpdateShapeOperator(m_pBox);
    }
    if (p_vertex->GetVertexType() == 1) {
        p_vertex->GetPrecedingEdgeLink()->UpdateEdgeVector(m_pBox);
    }
    AbstractCurvature* p_curvature = m_pState->GetCurvatureCalculator();
    p_curvature->UpdateVertexCurvature(p_vertex);
    for (std::vector<vertex*>::const_iterator it = v_neighbours.begin(); it != v_neighbours.end(); ++it) {
        p_curvature->UpdateVertexCurvature(*it);
    }
    //---> new energies
    AbstractEnergy* p_energy = m_pState->GetEnergyCalculator();
    double new_energy = p_energy->SingleVertexEnergy(p_vertex);
    new_energy += p_energy->CalculateVectorFieldMembraneBindingEnergy(p_vertex);
    for (std::vector<vertex*>::const_iterator it = v_neighbours.begin(); it != v_neighbours.end(); ++it) {
        new_energy += p_energy->SingleVertexEnergy(*it);
        new_energy += p_energy->CalculateVectorFieldMembraneBindingEnergy(*it);
    }
    for (std::vector<links*>::iterator it = m_AffectedLinks.begin(); it != m_AffectedLinks.end(); ++it) {
        new_energy += p_energy->TwoInclusionsInteractionEnergy(*it);
        if (p_vertex->GetNumberOfVF() != 0) {
            for (int vf_layer = 0; vf_layer < m_pState->GetMesh()->GetNoVFPerVertex(); vf_layer++) {
                new_energy += p_energy->TwoVectorFieldInteractionEnergy(vf_layer, *it);
            }
        }
    }
    double dE_Cgroup = m_pState->GetApplyConstraintBetweenGroups()->CalculateEnergyChange(p_vertex, dx);

    if (p_VAH->GetCalculateVAH()) {
        p_VAH->CalculateAVertexRingContributionToGlobalVariables(p_vertex, new_Tvolume, new_Tarea, new_Tcurvature);
    }
    double dE_volume = m_pState->GetVolumeCoupling()->GetEnergyChange(old_Tarea, old_Tvolume, new_Tarea, new_Tvolume);
    double dE_t_area = m_pState->GetTotalAreaCoupling()->CalculateEnergyChange(old_Tarea, new_Tarea);
    double dE_g_curv = m_pState->GetGlobalCurvature()->CalculateEnergyChange(new_Tarea - old_Tarea, new_Tcurvature - old_Tcurvature);

    double diff_energy = new_energy - old_energy;
    double tot_diff_energy = diff_energy + dE_Cgroup + dE_force_on_vertex + dE_force_from_inc
                           + dE_force_from_vector_fields + dE_volume + dE_t_area + dE_g_curv;
    double U = (*m_pBeta) * tot_diff_energy - (*m_pDBeta);
    if (U <= 0 || exp(-U) > thermal) {
        if (!p_vertex->CheckVoxel()) {
            p_vertex->UpdateVoxelAfterAVertexMove();
        }
        outcome.accepted = true;
        outcome.dE_elastic = diff_energy;
        outcome.dVolume = new_Tvolume - old_Tvolume;
        outcome.dArea = new_Tarea - old_Tarea;
        outcome.dCurvature = new_Tcurvature - old_Tcurvature;
        return outcome;
    }
    //---> reject: reverse the changes
    for (std::vector<triangle*>::const_iterator it = v_triangles.begin(); it != v_triangles.end(); ++it) {
        (*it)->ReverseConstantMesh_Copy();
    }
    for (std::vector<links*>::const_iterator it = v_links.begin(); it != v_links.end(); ++it) {
        (*it)->ReverseConstantMesh_Copy();
        (*it)->GetNeighborLink1()->ReverseConstantMesh_Copy();
    }
    for (std::vector<links*>::iterator it = m_AffectedLinks.begin(); it != m_AffectedLinks.end(); ++it) {
        (*it)->Reverse_InteractionEnergy();
        (*it)->Reverse_VFInteractionEnergy();
    }
    if (p_vertex->GetVertexType() == 1) {
        p_vertex->GetPrecedingEdgeLink()->ReverseConstantMesh_Copy();
    }
    p_vertex->ReverseConstantMesh_Copy();
    p_vertex->Reverse_VFsBindingEnergy();
    for (std::vector<vertex*>::const_iterator it = v_neighbours.begin(); it != v_neighbours.end(); ++it) {
        (*it)->ReverseConstantMesh_Copy();
        (*it)->Reverse_VFsBindingEnergy();
    }
    return outcome;
}
// face angles are checked after the move, in CheckFacesAfterAVertexMove
bool VertexMoveKernel::VertexMoveIsFine(vertex* p_vertex, const Vec3D& dx) {

    double new_x = p_vertex->GetXPos() + dx(0);
    double new_y = p_vertex->GetYPos() + dx(1);
    double new_z = p_vertex->GetZPos() + dx(2);
    if (new_x >= (*m_pBox)(0)) {
        new_x = new_x - (*m_pBox)(0);
    } else if (new_x < 0) {
        new_x = new_x + (*m_pBox)(0);
    }
    if (new_y >= (*m_pBox)(1)) {
        new_y = new_y - (*m_pBox)(1);
    } else if (new_y < 0) {
        new_y = new_y + (*m_pBox)(1);
    }
    if (new_z >= (*m_pBox)(2)) {
        new_z = new_z - (*m_pBox)(2);
    } else if (new_z < 0) {
        new_z = new_z + (*m_pBox)(2);
    }
    const double mindist2 = *m_pMinLength2;
    const double maxdist2 = *m_pMaxLength2;
    const std::vector<vertex*>& v_neighbours = p_vertex->GetVNeighbourVertex();
    for (std::vector<vertex*>::const_iterator it = v_neighbours.begin(); it != v_neighbours.end(); ++it) {
        double dist2 = p_vertex->SquareDistanceOfAVertexFromAPoint(new_x, new_y, new_z, *it);
        if (dist2 < mindist2 || dist2 > maxdist2)
            return false;
    }
    //---> the target voxel may differ from the current one
    Voxel<vertex>* p_voxel = p_vertex->GetVoxel();
    int new_nx = int(new_x / p_voxel->GetXSideVoxel((*m_pBox)(0)));
    int new_ny = int(new_y / p_voxel->GetYSideVoxel((*m_pBox)(1)));
    int new_nz = int(new_z / p_voxel->GetZSideVoxel((*m_pBox)(2)));
    int i = Voxel<int>::Convert2LocalVoxelIndex(new_nx, p_voxel->GetXIndex(), p_voxel->GetXNoVoxel());
    int j = Voxel<int>::Convert2LocalVoxelIndex(new_ny, p_voxel->GetYIndex(), p_voxel->GetYNoVoxel());
    int k = Voxel<int>::Convert2LocalVoxelIndex(new_nz, p_voxel->GetZIndex(), p_voxel->GetZNoVoxel());
    if (i > 1 || i < -1 || j > 1 || j < -1 || k > 1 || k < -1) {
        std::cout << " ---> warning: the object might moved more than one voxel " << std::endl;
        return false;
    }
    Voxel<vertex>* new_pvox = p_voxel->GetANeighbourCell(i, j, k);
    for (int n = -1; n < 2; n++)
    for (int m = -1; m < 2; m++)
    for (int s = -1; s < 2; s++) {
        std::vector<vertex*> CV = new_pvox->GetANeighbourCell(n, m, s)->GetContentObjects();
        for (std::vector<vertex*>::iterator it = CV.begin(); it != CV.end(); ++it) {
            if (*it != p_vertex && p_vertex->SquareDistanceOfAVertexFromAPoint(new_x, new_y, new_z, *it) < mindist2)
                return false;
        }
    }
    return true;
}
bool VertexMoveKernel::CheckFacesAfterAVertexMove(vertex* p_vertex) {
    const std::vector<links*>& v_links = p_vertex->GetVLinkList();
    for (std::vector<links*>::const_iterator it = v_links.begin(); it != v_links.end(); ++it) {
        if (!(*it)->CheckFaceAngleWithMirrorFace(*m_pMinAngle) || !(*it)->CheckFaceAngleWithNextEdgeFace(*m_pMinAngle)) {
            return false;
        }
    }
    return true;
}
void VertexMoveKernel::GetEdgesWithInteractionChange(vertex* p_vertex, std::vector<links*>& edges) {
    edges.clear();
    m_CandidateLinks.clear();
    const std::vector<vertex*>& v_neighbours = p_vertex->GetVNeighbourVertex();
    for (std::vector<vertex*>::const_iterator it = v_neighbours.begin(); it != v_neighbours.end(); ++it) {
        if ((*it)->VertexOwnInclusion() || p_vertex->GetNumberOfVF() != 0) {
            const std::vector<links*>& n_links = (*it)->GetVLinkList();
            m_CandidateLinks.insert(m_CandidateLinks.end(), n_links.begin(), n_links.end());
            if ((*it)->m_VertexType == 1) {
                m_CandidateLinks.push_back((*it)->m_pPrecedingEdgeLink);
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
