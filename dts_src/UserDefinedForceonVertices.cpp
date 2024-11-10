#include "UserDefinedForceonVertices.h"


UserDefinedForceonVertices::UserDefinedForceonVertices(std::string input) {

    m_Inputs = input;
}
UserDefinedForceonVertices::~UserDefinedForceonVertices() {
    
}
double UserDefinedForceonVertices::Energy_of_Force(vertex *pv, Vec3D dx) {
    
    
    Vec3D Force = GetForce(pv);
    return -Vec3D::dot(dx , Force);
}
Vec3D UserDefinedForceonVertices::GetForce(vertex *v) // gives force in the local coordinate
{
    Vec3D f;
///====

//== write you code here

    std::cout<<" user is applying force \n";

    return f;
}
std::string UserDefinedForceonVertices::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    state = state + m_Inputs;
    return state;
}

