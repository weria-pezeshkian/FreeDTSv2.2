#include "Constant_NematicForce.h"


Constant_NematicForce::Constant_NematicForce()
{
    m_F0 =0;
    m_Fd= 1;   // in the direction of the protein
    m_Fn= 0;    // in the direction of surface normal
    m_En= 0;

}
Constant_NematicForce::~Constant_NematicForce()
{
    
}
void Constant_NematicForce::Initialize()
{
    double size = sqrt(m_Fd*m_Fd+m_Fn*m_Fn);
    m_Fd=m_Fd/size;
    m_Fn=m_Fn/size;
}
double Constant_NematicForce::Energy_of_Force(vertex *pv, Vec3D dx)
{
    m_En = 0;
    if(pv->VertexOwnInclusion()==true && m_F0!=0)
    {
        Tensor2  G2L = pv->GetG2LTransferMatrix();
        Vec3D ldx = G2L*dx;
        Vec3D incd = (pv->GetInclusion())->GetLDirection();
        Vec3D nv = pv->GetNormalVector();
        m_En = -m_Fd*m_F0*(incd.dot(ldx,incd));  // in the local space
        m_En+= -m_Fn*m_F0*(incd.dot(dx,nv));   // in the golbal space

    }
    
    return m_En;
}


