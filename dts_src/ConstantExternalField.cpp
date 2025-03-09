#include "ConstantExternalField.h"

// E = -k*(Field*Inc_direction)^2 --> we could make it k1*(Field*Inc_direction)-k2*(Field*Inc_direction)^2

ConstantExternalField::ConstantExternalField(double k, double x, double y, double z) {
    // Initialize the field strength
    m_FieldStrength = k;

    // Calculate the magnitude of the vector (x, y, z)
    double norm = sqrt(x*x + y*y + z*z);

    // Normalize the vector (x, y, z) to obtain the field direction
    Vec3D tem(x, y, z);
    tem = tem * (1.0 / norm); // Scale the vector by the reciprocal of its magnitude to normalize it

    // Assign the normalized vector as the field direction
    m_FieldDirection = tem;
    
    
    m_CouplingFunction = CouplingType_1;
}
ConstantExternalField::~ConstantExternalField()
{
    
}
double ConstantExternalField::GetCouplingEnergy(vertex *pvertex) {
    // Check if the vertex owns an inclusion
    if (!pvertex->VertexOwnInclusion())
        return 0; // If not, return zero coupling energy
    
    if(m_FieldStrength==0)
        return 0; // return zero coupling energy

    // Get the local direction of the inclusion
    Vec3D LD = pvertex->GetInclusion()->GetLDirection();

    // Transform the local direction to global direction using the transfer matrix
    Vec3D GD = pvertex->GetL2GTransferMatrix() * LD;


    // Calculate and return the coupling energy
    // Note: The negative sign indicates that the field tends to minimize the angle with the field direction
    return -m_FieldStrength * m_CouplingFunction(Vec3D::dot(GD , m_FieldDirection));
}
double ConstantExternalField::CouplingType_1(const double &Cos){
    return Cos*Cos;
}
double ConstantExternalField::CouplingType_2(const double &Cos){
    return Cos;
}
std::string ConstantExternalField::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    state += " "+Nfunction::D2S(m_FieldStrength) +" "+ Nfunction::D2S(m_FieldDirection(0))+" "+ Nfunction::D2S(m_FieldDirection(1))+" "+ Nfunction::D2S(m_FieldDirection(2));
    return state;
}
