


#include "CouplingtoFixedGlobalCurvature.h"
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 This object affects the simulation results in every step. It is a global coupling, penalize the deviations of total mean curvature  from a target value
 It can also model the energy associated with the differnce in the inner and outer monolayer area
 The energy coupling is as below
 E=k/(2A)*(sum(2h-gc0)av)^2
 E=k/(2Ah^2)*(DA-DA0)^2
 */
CouplingtoFixedGlobalCurvature::CouplingtoFixedGlobalCurvature(VAHGlobalMeshProperties *VHA, State *pstate, double Gkappa, double GlobalC0) :
                AbstractGlobalCurvature(VHA, pstate),
                m_K(Gkappa/2.0),
                m_gC0(GlobalC0) {

}

CouplingtoFixedGlobalCurvature::~CouplingtoFixedGlobalCurvature()
{
    
}
// A function to initialize total area, total curvature and energy asscoiated with the coupling
void CouplingtoFixedGlobalCurvature::Initialize() {

    double dh=(m_TotalCurvature-m_gC0*m_TotalArea);
    m_Energy =m_K/(m_TotalArea)*dh*dh;
}
// when a vertex moves or a link flips or any other changes, we can see the changes in the total area and total mean curavture:
// this function can tell us how much this changes cost energy but it does not update the change since it can be rejected.
double CouplingtoFixedGlobalCurvature::CalculateEnergyChange(double DA, double DC)
{
    double de=0;
    double dh=(m_TotalCurvature+DC-m_gC0*(m_TotalArea+DA));
    double e= m_K/((m_TotalArea+DA))*dh*dh;  // new energy
    
    de=e-m_Energy;   // de = newenergy -oldenergy
    
    return de;
}
// this function update the changes if the move get accetped
void CouplingtoFixedGlobalCurvature::UpdateEnergyChange(double DA, double DC)
{
    double de=0;
    
    m_TotalArea+= DA;
    m_TotalCurvature+= DC;
    double dh=(m_TotalCurvature-m_gC0*m_TotalArea);
    m_Energy= m_K*dh*dh/(m_TotalArea);
 
}
std::string CouplingtoFixedGlobalCurvature::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    return state;
}
