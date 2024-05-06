

#include <stdio.h>
#include "State.h"
#include "InclusionPoseUpdateByMetropolisAlgorithm.h"

InclusionPoseUpdateByMetropolisAlgorithm::InclusionPoseUpdateByMetropolisAlgorithm (){
    
}
InclusionPoseUpdateByMetropolisAlgorithm::~InclusionPoseUpdateByMetropolisAlgorithm() {
    
}
bool InclusionPoseUpdateByMetropolisAlgorithm::Initialize (State *pState)
{
    m_pState = pState;
   // m_Beta = m_pState->m_Beta;
    return true;
}
bool InclusionPoseUpdateByMetropolisAlgorithm::EvolveOneStep(int step){
    
    const std::vector<inclusion *>& pAllInclusion = m_pState->GetMesh()->GetInclusion();

    int no_incs = pAllInclusion.size();
    int no_steps_angle = no_incs*m_NumberOfMovePerStep_Angle;
    int no_steps_kawa = no_incs* m_NumberOfMovePerStep_Kawasaki;
    
  for (int i = 0; i< no_steps_kawa;i++) {
    
      int r_inc_id = m_pState->GetRandomNumberGenerator()->IntRNG(no_incs);
      inclusion *p_inc = pAllInclusion[r_inc_id];
      double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
      
      std::vector <links*> linklist=p_inc->Getvertex()->GetVLinkList();
      int m = m_pState->GetRandomNumberGenerator()->IntRNG(linklist.size());
      links *d_link = linklist[m];
      
      if(KawasakiMove(step, p_inc,d_link, thermal)){
          m_AcceptedMoves++;
      }
      m_NumberOfAttemptedMoves++;
    }
    for (int i = 0; i< no_steps_angle;i++) {
      
        int r_inc_id = m_pState->GetRandomNumberGenerator()->IntRNG(no_incs);
        inclusion *p_inc = pAllInclusion[r_inc_id];
        double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
        double dx=1-2*(m_pState->GetRandomNumberGenerator()->UniformRNG(1.0));
        double dy=1-2*(m_pState->GetRandomNumberGenerator()->UniformRNG(1.0));
        
        if(RotationMove(step, p_inc,thermal, m_DR*dx, m_DR*dy)){
            m_AcceptedMoves++;
        }
        m_NumberOfAttemptedMoves++;
      }
    
    return true;
}
bool InclusionPoseUpdateByMetropolisAlgorithm::KawasakiMove(int step, inclusion* p_inc, links * d_links, double thermal) {

    
    return true;
}
bool InclusionPoseUpdateByMetropolisAlgorithm::RotationMove(int step, inclusion *p_inc, double dx, double dy, double thermal) {

  
    return true;
}










