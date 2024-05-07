


#include "HarmonicPotentialBetweenTwoGroups.h"
#include "Nfunction.h"
#include "State.h"
HarmonicPotentialBetweenTwoGroups::HarmonicPotentialBetweenTwoGroups(State* pState, double K, double l0, std::string group1,std::string group2,double nx,double ny,double nz) : MESH(*(pState->GetMesh())) {
    m_K = K/2;
    m_Group1Name = group1;
    m_Group2Name = group2;
    m_Direction(0) = nx;
    m_Direction(1) = ny;
    m_Direction(2) = nz;
    m_L0 = l0;
    m_Energy = 0;
    
    // We should prepare the system here
    //1) Make groups and find there com
    // create a function to update group com if a vertex moves, it should check if the vertex exist in the group and if yes, update the com.
    // We should also reject the move or accept it, need more functions
}
HarmonicPotentialBetweenTwoGroups::~HarmonicPotentialBetweenTwoGroups(){
    
}
                                                                                                                                                                
void HarmonicPotentialBetweenTwoGroups::CalculateEnergy(int step)
{


    Vec3D *pBox = (m_pGroup1.at(0))->GetBox();
    Vec3D PX1 = m_Group1COG;
    Vec3D PX2 = m_Group2COG;
    Vec3D Dist;

    //std::cout<<" distance before "<<PX1(2)-PX2(2)<<"\n";
    for (int i=0;i<3;i++) {
        double dist=PX1(i)-PX2(i);
        Dist(i)=dist;
    }
    //std::cout<<" distance after "<<Dist(2)<<"\n";

    double dist=m_Direction(2)*Dist(2)*Dist(2)+m_Direction(0)*Dist(0)*Dist(0)+m_Direction(1)*Dist(1)*Dist(1);
    m_Dist=sqrt(dist);
    m_Dist = dist;
    m_Energy = m_K*(dist-m_L0)*(dist-m_L0);
}
void HarmonicPotentialBetweenTwoGroups::MovingVertex(vertex* v, Vec3D Dx)
{
    if(v->GetGroupName() == m_Group1Name)
    {
        m_Group1COG = m_Group1COG + Dx*(1/double(m_pGroup1.size()));
    }
    if(v->GetGroupName() == m_Group2Name)
    {
        m_Group2COG = m_Group2COG + Dx*(1/double(m_pGroup2.size()));
    }

}
void HarmonicPotentialBetweenTwoGroups::RejectMovingVertex(vertex* v, Vec3D Dx)
{

    if(v->GetGroupName() == m_Group1Name)
    {
        m_Group1COG = m_Group1COG - Dx*(1/double(m_pGroup1.size()));
    }
    if(v->GetGroupName() == m_Group2Name)
    {
        m_Group2COG = m_Group2COG - Dx*(1/double(m_pGroup2.size()));
    }
}
Vec3D HarmonicPotentialBetweenTwoGroups::COMVertexGroup(std::vector<vertex *> ver)
{
    Vec3D com;
    double x=0;
    double y=0;
    double z=0;
    for (std::vector<vertex *>::iterator it = ver.begin() ; it != ver.end(); ++it)
    {
        x+=(*it)->GetVXPos();
        y+=(*it)->GetVYPos();
        z+=(*it)->GetVZPos();
    }
    com(0)=x/ver.size();
    com(1)=y/ver.size();
    com(2)=z/ver.size();

    return com;
}


bool HarmonicPotentialBetweenTwoGroups::Initialize() {
    
    if (m_Groups.find(m_Group1Name) != m_Groups.end() && m_Groups.find(m_Group2Name) != m_Groups.end() ) {
        
        m_pGroup1 = m_Groups.at(m_Group1Name);
        m_pGroup2 = m_Groups.at(m_Group2Name);
        
    } else {
        std::cout<<"---error--> groups for "<<GetDefaultReadName()<<" does not exist.\n";
        return false;
    }
    m_Group1COG = COMVertexGroup(m_pGroup1);
    m_Group2COG = COMVertexGroup(m_pGroup2);
    Vec3D PX1 = m_Group1COG;
    Vec3D PX2 = m_Group2COG;
    Vec3D Dist;
    
    for (int i=0;i<3;i++)
    {
        double dist=PX1(i)-PX2(i);
        if(dist>(*m_pBox)(i)/2.0)
        {
            dist=(*m_pBox)(i)-dist;
        }
        else if(dist<0)
        {
            dist=-dist;
        }
        Dist(i)=dist;
    }
    m_L0=m_Direction(0)*Dist(0)*Dist(0)+m_Direction(1)*Dist(1)*Dist(1)+m_Direction(2)*Dist(2)*Dist(2);
    m_L0=sqrt(m_L0);
    
    return true;
}
