

#include <stdio.h>
#include "AlexanderMoveByMetropolisAlgorithm.h"
#include "State.h"

AlexanderMoveByMetropolisAlgorithm::AlexanderMoveByMetropolisAlgorithm(){

}
AlexanderMoveByMetropolisAlgorithm::~AlexanderMoveByMetropolisAlgorithm(){
    
}
bool AlexanderMoveByMetropolisAlgorithm::Initialize(State *pState){
    
    m_pState    = pState;
    m_pBox      = (pState->m_pMesh)->m_pBox;
    m_pminAngle = &(m_pState->m_MinFaceAngle);
    m_pLmin2    = &(m_pState->m_MinVerticesDistanceSquare);
    m_pLmax2    = &(m_pState->m_MaxLinkLengthSquare);
    m_pBeta     = &(m_pState->m_Beta);
    *m_RateOfVMovePerStep = 1;

    
    return true;
}
bool AlexanderMoveByMetropolisAlgorithm::EvolveOneStep(int step){
 
  //  for (int i=0;i<;i++) {
    
    double dx;//=1-2*Random1.UniformRNG(1.0);            // Inside a cube with the side length of R
    double dy;//=1-2*Random1.UniformRNG(1.0);
    double dz;//=1-2*Random1.UniformRNG(1.0);
    double thermal;//=Random1.UniformRNG(1.0);
    vertex *pvertex;
    EvolveOneVertex(step, pvertex, dx, dy, dz,thermal);
   // }
    
    return true;
}
int AlexanderMoveByMetropolisAlgorithm::EvolveOneVertex(int step, vertex *pvertex, double dx, double dy, double dz,double temp){
    int moveresult = 0;
    double old_energy = 0;
    double new_energy = 0;

//---> first checking if all the distances will be fine if we move the vertex
    if(!VertexMoveIsFine(pvertex,dx,dy,dz,*m_pLmin2,*m_pLmax2))
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
    old_energy=pvertex->GetEnergy();
    std::vector <vertex *> vNeighbourV = pvertex->GetVNeighbourVertex();
    for (std::vector<vertex *>::iterator it = vNeighbourV.begin() ; it != vNeighbourV.end(); ++it){
        old_energy+=(*it)->GetEnergy();
    }
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
bool AlexanderMoveByMetropolisAlgorithm::CheckFacesAfterAVertexMove(double &minangle, vertex* p_vertex) {
    std::vector<links*> linkList = p_vertex->GetVLinkList();
    for (std::vector<links*>::iterator it = linkList.begin(); it != linkList.end(); ++it) {
        links* link = *it;
        if (!link->CheckFaceAngleWithMirrorFace(minangle) || !link->CheckFaceAngleWithNextEdgeFace(minangle)) {
            return false;
        }
    }
    return true;
}
