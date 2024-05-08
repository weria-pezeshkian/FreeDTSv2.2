

#include <stdio.h>
#include "EvolveVerticesByMetropolisAlgorithm.h"
#include "State.h"

EvolveVerticesByMetropolisAlgorithm::EvolveVerticesByMetropolisAlgorithm(State *pState) :
                        m_pState(pState),
                        MESH(*(pState->GetMesh())),
                        AbstractSimulation(*(pState->GetSimulation())){

}
EvolveVerticesByMetropolisAlgorithm::~EvolveVerticesByMetropolisAlgorithm(){

}
void EvolveVerticesByMetropolisAlgorithm::Initialize(){
    
    m_NumberOfAttemptedMoves = 0;
    m_AcceptedMoves = 0;

    return;
}
bool EvolveVerticesByMetropolisAlgorithm::EvolveOneStep(int step){
 

    int no_surf_v = m_pSurfV.size();
    int no_edge_v = m_pEdgeV.size();
    int no_steps_edge = no_edge_v*m_NumberOfMovePerStep_Edge;
    int no_steps_surf = no_surf_v*m_NumberOfMovePerStep_Surf;

  for (int i = 0; i< no_steps_surf;i++) {
    
      int r_vid = m_pState->GetRandomNumberGenerator()->IntRNG(no_surf_v);
      vertex *pvertex = m_pSurfV[r_vid];
      if( m_FreezGroupName == pvertex->GetGroupName()){
          continue;
      }
      double dx=1-2*(m_pState->GetRandomNumberGenerator()->UniformRNG(1.0));
      double dy=1-2*(m_pState->GetRandomNumberGenerator()->UniformRNG(1.0));
      double dz=1-2*(m_pState->GetRandomNumberGenerator()->UniformRNG(1.0));
      
      if(!m_pState->GetBoundary()->MoveHappensWithinTheBoundary(dx,dy,dz, pvertex)){
          continue;
      }
      
      double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
      if(EvolveOneVertex(step, pvertex, m_DR*dx, m_DR*dy, m_DR*dz,thermal)){
          m_AcceptedMoves++;
      }
      m_NumberOfAttemptedMoves++;
    }
    
    for (int i = 0; i< no_steps_edge;i++) {
      
        int r_vid = m_pState->GetRandomNumberGenerator()->IntRNG(no_edge_v);
        vertex *pvertex = m_pEdgeV[r_vid];
        if( m_FreezGroupName == pvertex->GetGroupName()){
            continue;
        }
        double dx=1-2*(m_pState->GetRandomNumberGenerator()->UniformRNG(1.0));
        double dy=1-2*(m_pState->GetRandomNumberGenerator()->UniformRNG(1.0));
        double dz=1-2*(m_pState->GetRandomNumberGenerator()->UniformRNG(1.0));
        if(!m_pState->GetBoundary()->MoveHappensWithinTheBoundary(dx,dy,dz, pvertex)){
            continue;
        }
        
        double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);

        if(EvolveOneVertex(step, pvertex, m_DR*dx, m_DR*dy, m_DR*dz,thermal)){
            m_AcceptedMoves++;
        }
        m_NumberOfAttemptedMoves++;
    }
    
    return true;
}
bool EvolveVerticesByMetropolisAlgorithm::EvolveOneVertex(int step, vertex *pvertex, double dx, double dy, double dz,double temp){
    
    double old_energy = 0;
    double new_energy = 0;

//---> first checking if all the distances will be fine if we move the vertex
    if(!VertexMoveIsFine(pvertex,dx,dy,dz,m_MinLength2,m_MaxLength2))
        return 0;
    
//---> lets get some variables before moving
    // --- get the initial total area and volume if they are needed.
     double old_Tvolume = 0;
     double old_Tarea = 0;
     double old_Tcurvature = 0;
    if(m_pState->GetVolumeCoupling()->GetDerivedDefaultReadName()!="NoCoupling"){
        //GetRingVolumeOfVertex(pvertex);
    }
    //--- obtain vertices energy terms
    old_energy = pvertex->GetEnergy();
    const std::vector<vertex *>& vNeighbourV = pvertex->GetVNeighbourVertex();
    for (std::vector<vertex *>::const_iterator it = vNeighbourV.begin() ; it != vNeighbourV.end(); ++it){
        
        old_energy+=(*it)->GetEnergy();
    }
    //-- obtain link interaction energies
    
    //--- Move the vertex;
    //--- update triangles normal
    
    //--- check the angle
    
    //obtain new energy
    
    
    // see if it matches
    
    
    
    
    return true;
}
//---> this does not check the angle of the faces. Because this should be done after the move:
//waste of calculation if we do ith before the move. Unless, we store the values. That also not good because move could get rejected.
bool EvolveVerticesByMetropolisAlgorithm::VertexMoveIsFine(vertex* pvertex, double dx,double dy, double dz,  double mindist2, double maxdist2){
//--->  vertex new position if get accepted
    
        double new_x = pvertex->GetXPos() + dx;
        double new_y = pvertex->GetYPos() + dy;
        double new_z = pvertex->GetZPos() + dz;
    
//--->  let check the distances with the nighbours
    std::vector <vertex *> npvertex = pvertex->GetVNeighbourVertex();
    for (std::vector<vertex *>::iterator it = npvertex.begin() ; it != npvertex.end(); ++it){
        double dist2 = (*it)->SquareDistanceOfAVertexFromAPoint(new_x, new_y, new_z, *it);
            if(dist2<mindist2 || dist2>maxdist2)
            return false;
    }
      
//---> now check it within the voxel cells
//---> lets get the object voxel, note: the voxel could be different from the associated one as it is moving

        //-- obtain the object (vertex) new cell id, with respect to the current cell
        int i = int(new_x/pvertex->GetVoxel()->GetXSideVoxel())-pvertex->GetVoxel()->GetXIndex();
        int j = int(new_y/pvertex->GetVoxel()->GetYSideVoxel())-pvertex->GetVoxel()->GetYIndex();
        int k = int(new_z/pvertex->GetVoxel()->GetZSideVoxel())-pvertex->GetVoxel()->GetZIndex();
        //-- check if it has moved too far
        if(i > 1 || i < -1 || j > 1 || j < -1 || k > 1 || k < -1) {
            std::cout << " ---> warning: the object might moved more than one voxel " << std::endl;
            return false;
        }
        Voxel<vertex>* new_pvox = pvertex->GetVoxel()->GetANeighbourCell(i, j, k);
    
        for(int n=-1;n<2;n++)
        for(int m=-1;m<2;m++)
        for(int s=-1;s<2;s++){
            std::vector <vertex *> CV = new_pvox->GetANeighbourCell(n, m, s)->GetContentObjects();
            for (std::vector<vertex *>::iterator it = CV.begin() ; it != CV.end(); ++it){
                if(*it!=pvertex){
                    if((*it)->SquareDistanceOfAVertexFromAPoint(new_x, new_y, new_z, *it)<mindist2)
                        return false;
                }
                
            }
        } ///   for(int s=-1;s<2;s++){
        
 
    return true;
}
// finding the distance of the current vertex from the pv2; also considering the pbc conditions
bool EvolveVerticesByMetropolisAlgorithm::CheckFacesAfterAVertexMove(double &minangle, vertex* p_vertex) {
    std::vector<links*> linkList = p_vertex->GetVLinkList();
    for (std::vector<links*>::iterator it = linkList.begin(); it != linkList.end(); ++it) {
        links* link = *it;
        if (!link->CheckFaceAngleWithMirrorFace(minangle) || !link->CheckFaceAngleWithNextEdgeFace(minangle)) {
            return false;
        }
    }
    return true;
}


//========
// this function can be deleted any time; it is for test cases only
double  EvolveVerticesByMetropolisAlgorithm::SystemEnergy()
{
    /*
    MESH* m_pMESH = m_pState->m_pMesh;
    std::vector<vertex *> ActiveV = m_pMESH->m_pActiveV;
    std::vector<triangle *> pActiveT = m_pMESH->m_pActiveT;
    std::vector<links *> mLink = m_pMESH->m_pHL;
    std::vector<links *>  pEdgeL = m_pMESH->m_pEdgeL;
    std::vector<vertex *> EdgeV  = m_pMESH->m_pEdgeV;
    double en = 0;
    

    for (std::vector<triangle *>::iterator it = pActiveT.begin() ; it != pActiveT.end(); ++it)
        (*it)->UpdateNormal_Area(m_pBox);
    
    
    for (std::vector<links *>::iterator it = mLink.begin() ; it != mLink.end(); ++it)
    {
        (*it)->UpdateNormal();
        (*it)->UpdateShapeOperator(m_pBox);
    }
    
    
    for (std::vector<vertex *>::iterator it = ActiveV.begin() ; it != ActiveV.end(); ++it)
        (m_pState->CurvatureCalculator())->SurfVertexCurvature(*it);

    //====== edge links should be updated
    for (std::vector<links *>::iterator it = pEdgeL.begin() ; it != pEdgeL.end(); ++it)
            (*it)->UpdateEdgeVector(m_pBox);

    for (std::vector<vertex *>::iterator it = EdgeV.begin() ; it != EdgeV.end(); ++it)
            (m_pState->CurvatureCalculator())->EdgeVertexCurvature(*it);
    
    en=m_pEnergyCalculator->TotalEnergy(ActiveV,mLink);
    //en=en+ m_pEnergyCalculator->TotalEnergy(EdgeV,pEdgeL);
   
    
    return en;
     */
    return 0;
}
std::string EvolveVerticesByMetropolisAlgorithm::CurrentState(){
    
    std::string state = AbstractVertexPositionIntegrator::GetBaseDefaultReadName() +" = "+ GetDerivedDefaultReadName();
    return state;
}
