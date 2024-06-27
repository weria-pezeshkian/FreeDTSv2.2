
#include <stdio.h>
#include "VectorField.h"
#include "vertex.h"

VectorField::VectorField()  {

}
VectorField::VectorField(InclusionType *inctype, double x, double y) : m_IncType(inctype) {
    
    m_LDirection(0) = x;
    m_LDirection(1) = y;
    m_LDirection(2) = 0;
    m_LDirection.normalize();
    m_MembraneBindingEnergy = 0;
}
VectorField::~VectorField() {
    
}
void VectorField::UpdateLocalDirection(const Vec3D & lo_dir) {
    m_LDirection = lo_dir;
    return;
}
void VectorField::UpdateGlobalDirection(const Vec3D & lg_dir) {
    m_GDirection = lg_dir;
    return;
}
void VectorField::UpdateMembraneBindingEnergy(const double &en){
    
    m_MembraneBindingEnergy = en;
    return;
}
void VectorField::Add2MembraneBindingEnergy(const double &en){
    
    m_MembraneBindingEnergy += en;
    return;
}
double VectorField::CalculateMembraneBindingEnergy(vertex *p_vertex){
    
    double en = 0;
    
if(p_vertex->GetVertexType() == 0) {
    double c1 = p_vertex->GetP1Curvature();
    double c2 = p_vertex->GetP2Curvature();

    double mean_times2 = c1 + c2;
    double gussian = c1 * c1;
    double area = p_vertex->GetArea();
    
    double k0 = m_IncType->ITk;
    double kg = m_IncType->ITkg;
    double k1 = m_IncType->ITk1;
    double k2 = m_IncType->ITk2;
    double c0 = m_IncType->ITc0;
    double cp10 = m_IncType->ITc1;
    double cn20 = m_IncType->ITc2;
    //--- kappa/2*(2H-c0)^2-kgK+
    double ev = k0 * (mean_times2-c0) * (mean_times2-c0) - kg * gussian;
    en += ev * area;
    
//--> if the k1 and k2 are zero, we do not need to calculate the rest
    if(k1 != 0 || k2 != 0){  // k1/2(cp-cp0)^2+k2/2(cn-cn0)^2
    
        double Cos = m_LDirection(0);
        double Sin = m_LDirection(1);
        double Cp = c1 * Cos * Cos + c2 * Sin * Sin;
        double Cn = c1 * Cos * Cos + c2 * Sin * Sin;
        double Delta_Cp = Cp - cp10;
        double Delta_Cn = Cn - cn20;
        en += (k1 * Delta_Cp * Delta_Cp + k2 * Delta_Cn * Delta_Cn) * area;
    }
}
else if( p_vertex->GetVertexType() == 1) {
    
    double geo_c = p_vertex->GetGeodesicCurvature();
    double norm_c = p_vertex->GetNormalCurvature();
    double length = p_vertex->GetLength();
    
        double lambda = m_IncType->ITelambda;
        double kg     = m_IncType->ITekg;
        double kn     = m_IncType->ITekn;
        double cn0    = m_IncType->ITecn;
        
        en += (lambda + kg*geo_c*geo_c)*length;
        
        if(kn != 0) {
            
            double Cos = m_LDirection(0);
            kn = kn*Cos*Cos;
            double Delta_norm_c = (norm_c-cn0);
            en += (kn*Delta_norm_c*Delta_norm_c)*length;
        }
    
}
else{
    
    std::cout<<" error-> this is unexpected 903 \n";
}
    return en;
}


