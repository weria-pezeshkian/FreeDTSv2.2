

#include <stdio.h>
#include "AlexanderMoveByMetropolisAlgorithm.h"
#include "State.h"

AlexanderMoveByMetropolisAlgorithm::AlexanderMoveByMetropolisAlgorithm(){

}
AlexanderMoveByMetropolisAlgorithm::~AlexanderMoveByMetropolisAlgorithm(){
    
}
bool AlexanderMoveByMetropolisAlgorithm::Initialize(State *pState){
    
    m_pState    = pState;
    m_pBox      = pState->GetMesh()->GetBox();
  //  m_pminAngle = &(m_pState->m_MinFaceAngle);
 //   m_pLmin2    = &(m_pState->m_MinVerticesDistanceSquare);
  //  m_pLmax2    = &(m_pState->m_MaxLinkLengthSquare);
  //  m_pBeta     = &(m_pState->m_Beta);

    
    return true;
}
bool AlexanderMoveByMetropolisAlgorithm::EvolveOneStep(int step){
 
    const std::vector<links *>& pEdges = m_pState->GetMesh()->GetRightL();

    int no_edges = pEdges.size();
    int no_steps = no_edges*m_NumberOfMovePerStep;
    
  for (int i = 0; i< no_steps;i++) {
    
      int r_lid = m_pState->GetRandomNumberGenerator()->IntRNG(no_edges);
      links *p_link = pEdges[r_lid];

      double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
      if(FlipOneEdge(step, p_link,thermal)){
          m_AcceptedMoves++;
      }
      m_NumberOfAttemptedMoves++;
    }

    
    return true;
}
bool AlexanderMoveByMetropolisAlgorithm::FlipOneEdge(int step, links *p_link, double temp){
    int moveresult = 0;
    double old_energy = 0;
    double new_energy = 0;

//---> first checking if all the distances will be fine if we move the vertex
    if(!EdgeCanBeFliped(p_link,*m_pLmin2,*m_pLmax2))
        return 0;
    
//---> lets get some variables before moving
    // --- get the initial total area and volume if they are needed.
    // this section could be centrilized. Somewhere the total volume and total area and dA should be one
   /* double old_volume = 0;
    double old_total_area = 0;
    double old_total_delta_area = 0;
    if((m_pState->GetVolumeCoupling())->GetState()==true  || (m_pState->GetTotalAreaCoupling())->GetState()==true)
    {
        if((m_pState->GetVolumeCoupling())->GetState()==true ){
            old_volume = m_pState->GetVolumeCoupling()->GetTotalVolume();
            old_total_area = m_pState->GetVolumeCoupling()->GetTotalArea();
        }
        else if((m_pState->GetTotalAreaCoupling())->GetState()==true){
            old_total_area = m_pState->GetTotalAreaCoupling()->GetTotalArea();
        }
    }
    if(m_pState->GetGlobalCurvature()->GetState()==true){

    }*/
    //----
    //--- obtain vertices energy terms

    //-- obtain link interaction energies
    
    
    
    return moveresult;
}
// this function can be deleted any time; it is for test cases only
double  AlexanderMoveByMetropolisAlgorithm::SystemEnergy()
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


// finding the distance of the current vertex from the pv2; also considering the pbc conditions
bool AlexanderMoveByMetropolisAlgorithm::CheckFacesAfterFlip(double &minangle, links* p_links) {
   
    /*std::vector<links*> linkList = p_vertex->GetVLinkList();
    for (std::vector<links*>::iterator it = linkList.begin(); it != linkList.end(); ++it) {
        links* link = *it;
        if (!link->CheckFaceAngleWithMirrorFace(minangle) || !link->CheckFaceAngleWithNextEdgeFace(minangle)) {
            return false;
        }
    }*/
    return true;
}
bool AlexanderMoveByMetropolisAlgorithm::EdgeCanBeFliped(links *pedge, double mindist2, double maxdist2){
    
    
    
    return true;
}
std::string AlexanderMoveByMetropolisAlgorithm::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    return state;
}
