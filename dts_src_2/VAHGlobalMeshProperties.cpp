
#include "VAHGlobalMeshProperties.h"
#include "State.h"
VAHGlobalMeshProperties::VAHGlobalMeshProperties(State* pState){
    m_pState = pState;
}
VAHGlobalMeshProperties::~VAHGlobalMeshProperties() {
    
}
double VAHGlobalMeshProperties::VolumeofTrianglesAroundVertex(vertex *pv)
{
    double vol=0;
    std::vector <triangle *> pvT=pv->GetVTraingleList();
    for (std::vector<triangle *>::iterator it = pvT.begin() ; it != pvT.end(); ++it){
        vol+=SingleTriangleVolume(*it);
    }
    return vol;
}
//==========================================================
void VAHGlobalMeshProperties::Initialize(std::vector<triangle *> pTriangle){
    
    double V=0.0;
    double A=0.0;

    for (std::vector<triangle *>::iterator it = pTriangle.begin() ; it != pTriangle.end(); ++it)
    {
        V+=SingleTriangleVolume(*it);
        A+=(*it)->GetArea();
    }
    m_TotalVolume = V;
    m_TotalArea = A;
    
    return;
}
double VAHGlobalMeshProperties::SingleTriangleVolume(triangle *pt){
    double vol=0;
    vertex* pv= pt->GetV1();
    double area= pt->GetArea();
    Vec3D Norm=pt->GetNormalVector();
    Vec3D R (pv->GetVXPos(),pv->GetVYPos(),pv->GetVZPos());
    
    // If the system is PBC broken, then the volume will be wrong .....
    {
        vertex* pv2= pt->GetV2();
        Vec3D R2 (pv2->GetVXPos(),pv2->GetVYPos(),pv2->GetVZPos());
        R2=R2-R;
        vertex* pv3= pt->GetV3();
        Vec3D R3 (pv3->GetVXPos(),pv3->GetVYPos(),pv3->GetVZPos());
        R3=R3-R;
        
        if(R2.dot(R2,R2)>4 || R3.dot(R3,R3)>4)
        {
            Nfunction f;
            std::string sms="---> Error, the system crossed the PBC while using volume coupling; use a large box .. ";
            std::cout<<sms<<std::endl;
            exit(0);
        }
    }
    vol=area*(R.dot(Norm,R))/3.0;
    return vol;
}


