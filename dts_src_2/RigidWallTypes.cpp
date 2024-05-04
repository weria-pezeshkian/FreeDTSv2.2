
#include "State.h"
#include "RigidWallTypes.h"

TwoFlatParallelWall::TwoFlatParallelWall(State* pState, double thickness,  char direction){
    
    m_pState = pState;
    m_HalfThickness = thickness/2;
    m_MidPlane = 0;
    if(direction=='X'){
        m_Element = 0;
    }
    else if(direction=='Y'){
        m_Element = 1;
    }
    if(direction=='Z'){
        m_Element = 2;
    }
    else {
        *(m_pState->GetTimeSeriesLog()) << "---> such a direction for wall is wrong \n";
    }
}
TwoFlatParallelWall::~TwoFlatParallelWall(){
    
}
void TwoFlatParallelWall::Initialize() {

    m_MidPlane = 0;
    const std::vector<vertex *>& pAllVertices = m_pState->GetMesh()->GetSurfV();
    for (std::vector<vertex *>::const_iterator it = pAllVertices.begin(); it != pAllVertices.end(); ++it) {
        m_MidPlane+=((*it)->GetPos())(m_Element);
    }
    m_MidPlane = m_MidPlane/double(pAllVertices.size());

}
bool TwoFlatParallelWall::MoveHappensWithinTheBoundary(double dx, double dy, double dz, vertex* p_ver){
    
    Vec3D dr(dx,dy,dz);
    dr=dr+p_ver->GetPos();
    
    if(dr(m_Element)>m_MidPlane + m_HalfThickness || dr(m_Element)>m_MidPlane - m_HalfThickness)
        return false;

    return true;
}
