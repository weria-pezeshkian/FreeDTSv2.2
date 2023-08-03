#include "OpenEdgeEvolutionWithConstantVertex.h"
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian

 Edge treatment, this is a new development since June 2023;
 What it does:
 
 // 1. We need to give error if there is a edge; we cannot have osmotic pressure ...
 
 
 */
OpenEdgeEvolutionWithConstantVertex::OpenEdgeEvolutionWithConstantVertex()
{

}
OpenEdgeEvolutionWithConstantVertex::OpenEdgeEvolutionWithConstantVertex(bool state, MESH* pmesh,int rate, double lambda, double k1,double k2)
{
    m_pMESH = pmesh;
    m_State = state;
    m_Rate = rate;
    m_Lambda = lambda;
    m_K1 = k1;
    m_k2 = k2;
    m_WholeSize = 0;
    
}
OpenEdgeEvolutionWithConstantVertex::~OpenEdgeEvolutionWithConstantVertex()
{
    
}
void OpenEdgeEvolutionWithConstantVertex::Initialize()
{
  // creating some trinagle and links to later create from them
    int Nt = m_pMESH->m_pEdgeV.size();
    m_pBox = m_pMESH->m_pBox;
    // due to mirro links 2*(Nt-3)+Nt
    for (int i=0;i<3*Nt-6;i++)
    {
        int lid;
        links teml(lid);
        m_GhostL.push_back(teml);
    }
    for (int i=0;i<Nt-2;i++)
    {
        int tid;
        triangle temt(tid);
        m_GhostT.push_back(temt);
    }
    for (std::vector<links>::iterator it = m_GhostL.begin() ; it != m_GhostL.end(); ++it)
        m_pGhostL.push_back(&(*it));
    
    for (std::vector<triangle>::iterator it = m_GhostT.begin() ; it != m_GhostT.end(); ++it)
        m_pGhostT.push_back(&(*it));
    
    // calculate the total links size; not sure yet.
    std::vector<links*> edgelink = m_pMESH->m_pEdgeL;
    for (std::vector<links*>::iterator it = edgelink.begin() ; it != edgelink.end(); ++it)
    {
        (*it)->UpdateEdgeVector(m_pBox);
        m_WholeSize+=(*it)->GetEdgeSize();
    }
    
    // end
}
void OpenEdgeEvolutionWithConstantVertex::MC_Move(RNG* rng, double lmin, double lmax, double minangle)
{
    
    double createorkill = rng->UniformRNG(1.0);
    double thermal = rng->UniformRNG(1.0);
     
#if DEBUG_MODE == Enabled
    std::cout<<" mess_d: MC_Move Function  for the open edge: Starts \n";
#endif
    
if(createorkill<0.5)
{
        // an atempt to create a link                       //       v3
                                                            //      /   \
                                                          //     v1-----v2
    
    

    
    double oldedge = 0;
    
    if((m_pMESH->m_pEdgeL).size()>4)
    {
        int n = rng->IntRNG(m_pMESH->m_pEdgeV.size());
        vertex *v1 =m_pMESH->m_pEdgeV.at(n);
    
        if(Linkisvalid(v1,lmin, lmax, minangle)==true)
             links *l1 = CreateALink(v1);
    }
    else if((m_pMESH->m_pEdgeL).size()==4)
    {
            int n = rng->IntRNG(m_pMESH->m_pEdgeV.size());
            vertex *v1 =m_pMESH->m_pEdgeV.at(n);
        if(Linkisvalid(v1,lmin, lmax, minangle)==true)
       {
            links* l2 = v1->m_pEdgeLink;
            vertex *v3 = l2->GetV2();
            links* l1 = v3->m_pEdgeLink;
            vertex *v2 = l1->GetV2();  // we need v2, and the only way we can find it
            links *newL1 = CreateALink(v1);
            links *mnewL1 = CreateALink(v2);

           RemoveFromLinkList(newL1,m_pMESH->m_pEdgeL);
           RemoveFromVertexList(v1,m_pMESH->m_pEdgeV);
           RemoveFromVertexList(v2,m_pMESH->m_pEdgeV);
           RemoveFromLinkList(mnewL1,m_pMESH->m_pEdgeL);
           AddtoVertexList(v1,m_pMESH->m_pSurfV);
           AddtoVertexList(v2,m_pMESH->m_pSurfV);
           AddtoLinkList(newL1,m_pMESH->m_pMHL);
           AddtoLinkList(mnewL1,m_pMESH->m_pHL);

       }

    }
    else if((m_pMESH->m_pEdgeL).size()==0)
    {
        
    }
    else
    {
        std::cout<<"error---> it should not happen, this is not expected \n";
    }
}
else
{
        // an atempt to kill a link
    if(m_pMESH->m_pEdgeL.size()!=0)
    {
    
            int n = rng->IntRNG(m_pMESH->m_pEdgeL.size());
            links *plink =m_pMESH->m_pEdgeL.at(n);
            vertex *v1 = plink->GetV1();
            vertex *v2 = plink->GetV2();
    
        if(v1->GetVLinkList().size()>1 && v2->GetVLinkList().size()>1)
            KillALink(plink);
    }
}
  
#if DEBUG_MODE == Enabled
    std::cout<<" mess_d: MC_Move Function  for the open edge: Ends \n";
#endif

}
void OpenEdgeEvolutionWithConstantVertex::RemoveFromLinkList(links* z, std::vector<links*> &vect)
{
    vect.erase(std::remove(vect.begin(), vect.end(), z), vect.end());
}
void OpenEdgeEvolutionWithConstantVertex::RemoveFromVertexList(vertex* z, std::vector<vertex*> &vect)
{

    vect.erase(std::remove(vect.begin(), vect.end(), z), vect.end());
}
void OpenEdgeEvolutionWithConstantVertex::RemoveFromTriangleList(triangle* z, std::vector<triangle*> &vect)
{

    vect.erase(std::remove(vect.begin(), vect.end(), z), vect.end());
}
// this can be written as template
void OpenEdgeEvolutionWithConstantVertex::AddtoLinkList(links* z, std::vector<links*> &vect)
{
    vect.push_back(z);
}
void OpenEdgeEvolutionWithConstantVertex::AddtoVertexList(vertex* z, std::vector<vertex*> &vect)
{
    vect.push_back(z);
}
void OpenEdgeEvolutionWithConstantVertex::AddtoTriangleList(triangle* z, std::vector<triangle*> &vect)
{
    vect.push_back(z);
}
links* OpenEdgeEvolutionWithConstantVertex::CreateALink(vertex *v1)
{
    // an atempt to create a link                       //       v3
                                                        //      /   \
                                                      //     v1-----v2

    links* l2 = v1->m_pEdgeLink;
    vertex *v3 = l2->GetV2();
    links* l1 = v3->m_pEdgeLink;
    vertex *v2 = l1->GetV2();

    links *outlink;
    double diff_energy = 0;
    
    //if(exp(-diff_energy)>thermal )
    {
        // we added the triangle
        if(m_pGhostT.size()<1)
        {
            std::cout<<" error 882---> we do not have enough storage for the new trinagles \n";
            exit(0);
        }
        triangle *t = m_pGhostT.at(0);
        t->UpdateVertex(v1,v2,v3);
        RemoveFromTriangleList(t,m_pGhostT);
        AddtoTriangleList(t,m_pMESH->m_pActiveT);
        
        // add the trinagle to the vertcies list
        v1->AddtoTraingleList(t);
        v2->AddtoTraingleList(t);
        v3->AddtoTraingleList(t);
        
        // v3 is no longer an edge vertex
        RemoveFromVertexList(v3,m_pMESH->m_pEdgeV);
        AddtoVertexList(v3,m_pMESH->m_pSurfV);
        
        // we create three links, two mirror and one edge
        if(m_pGhostL.size()<3)
        {
            std::cout<<" error 883---> we do not have enough storage for the new links \n";
            exit(0);
        }
        RemoveFromLinkList(l1,m_pMESH->m_pEdgeL);
        RemoveFromLinkList(l2,m_pMESH->m_pEdgeL);
        
        links *l0 = m_pGhostL[0];
        links *ml1 = m_pGhostL[1];
        links *ml2 = m_pGhostL[2];
        RemoveFromLinkList(ml1,m_pGhostL);
        RemoveFromLinkList(ml2,m_pGhostL);
        RemoveFromLinkList(l0,m_pGhostL);
        AddtoLinkList(l0,m_pMESH->m_pEdgeL);
        AddtoLinkList(ml1,m_pMESH->m_pActiveL);
        AddtoLinkList(ml2,m_pMESH->m_pActiveL);

        v1->m_pEdgeLink = l0;
        l0->UpdateV(v1,v2,v3);
        l0->UpdateTriangle(t);

        
        ml1->UpdateV(v2,v3,v1);
        ml2->UpdateV(v3,v1,v2);
        ml1->UpdateTriangle(t);
        ml2->UpdateTriangle(t);
        
        ml1->UpdateMirrorLink(l1);
        ml2->UpdateMirrorLink(l2);
        l1->UpdateMirrorLink(ml1);
        l2->UpdateMirrorLink(ml2);
        ml1->UpdateNeighborLink1(ml2);
        ml1->UpdateNeighborLink2(l0);
        ml2->UpdateNeighborLink1(l0);
        ml2->UpdateNeighborLink2(ml1);
        l0->UpdateNeighborLink1(ml1);
        l0->UpdateNeighborLink2(ml2);
        
        // adding ml1 and ml2 and l1 and l2 to m_pMHL and m_pHL
        AddtoLinkList(l1,m_pMESH->m_pMHL);
        AddtoLinkList(l2,m_pMESH->m_pMHL);
        AddtoLinkList(ml1,m_pMESH->m_pHL);
        AddtoLinkList(ml2,m_pMESH->m_pHL);
        
        
        outlink = l0;

}
 
    return outlink;
}
void OpenEdgeEvolutionWithConstantVertex::KillALink(links *plink)
{
    double diff_energy = 0;
    // get the energy
    // ...
    // ...
    //if(exp(-diff_energy)>thermal )
    {

        triangle *t = plink->GetTriangle();
        vertex *v1 = plink->GetV1();
        vertex *v2 = plink->GetV2();
        vertex *v3 = plink->GetV3();
        links *l1 = plink->GetNeighborLink1();
        links *l2 = plink->GetNeighborLink2();
        links *ml1 = l1->GetMirrorLink();
        links *ml2 = l2->GetMirrorLink();
        // make the plink a ghost
        RemoveFromLinkList(plink,m_pMESH->m_pEdgeL);
        AddtoLinkList(plink,m_pGhostL);
        // adding the two mirror into the edge
        AddtoLinkList(ml1,m_pMESH->m_pEdgeL);
        AddtoLinkList(ml2,m_pMESH->m_pEdgeL);
        RemoveFromLinkList(ml1,m_pMESH->m_pHL);  // too expensive; I should find a better way
        RemoveFromLinkList(ml2,m_pMESH->m_pHL);  // too expensive
        RemoveFromLinkList(ml1,m_pMESH->m_pMHL);  // too expensive
        RemoveFromLinkList(ml2,m_pMESH->m_pMHL);  // too expensive
        
        // the two edge lose mirror
        ml1->UpdateMirrorFlag(false);
        ml2->UpdateMirrorFlag(false);
        v1->m_pEdgeLink = ml2;
        v3->m_pEdgeLink = ml1;
        
        // now two  links need to be removed as well
        RemoveFromLinkList(l1,m_pMESH->m_pActiveL);  // too expensive
        RemoveFromLinkList(l2,m_pMESH->m_pActiveL);  // too expensive
        RemoveFromLinkList(l1,m_pMESH->m_pHL);  // too expensive
        RemoveFromLinkList(l2,m_pMESH->m_pHL);  // too expensive
        RemoveFromLinkList(l1,m_pMESH->m_pMHL);  // too expensive
        RemoveFromLinkList(l2,m_pMESH->m_pMHL);  // too expensive
        AddtoLinkList(l1,m_pGhostL);
        AddtoLinkList(l2,m_pGhostL);
        // we need to also remove them from m_pHL and m_pMHL
        

       // removing the edge from vertex list
        v1->RemoveFromLinkList(plink);
        v3->RemoveFromLinkList(l2);
        v2->RemoveFromLinkList(l1);
        // v3 becomes an edge vertex
        RemoveFromVertexList(v3,m_pMESH->m_pSurfV);  // too expensive
        AddtoVertexList(v3,m_pMESH->m_pEdgeV);
    
        // the trinagle will be killed
        RemoveFromTriangleList(t,m_pMESH->m_pActiveT); // too expensive
        AddtoTriangleList(t,m_pGhostT);
        // remove the trinagle from all the vertcies list
        v1->RemoveFromTraingleList(t);
        v2->RemoveFromTraingleList(t);
        v3->RemoveFromTraingleList(t);
    }
}
bool OpenEdgeEvolutionWithConstantVertex::Linkisvalid(vertex *v1, double lmin, double lmax, double minangle)
{
    
    
    // check if the new link length is within the allowed range and also if the angle of the new trinagule is fine with respect to the two other trinagules
    bool isvalid = true;
    
    links* l2 = v1->m_pEdgeLink;
    vertex *v3 = l2->GetV2();
    links* l1 = v3->m_pEdgeLink;
    vertex *v2 = l1->GetV2();
    
    triangle t(0,v1,v2,v3);
    Vec3D *Box = v1->GetBox();
    t.UpdateNormal_Area(Box);
    Vec3D n1 = t.GetNormalVector();
    Vec3D n2 = (l2->GetTriangle())->GetNormalVector();
    Vec3D n3 = (l1->GetTriangle())->GetNormalVector();

    Vec3D P1(v1->GetVXPos(),v1->GetVYPos(),v1->GetVZPos());
    Vec3D P2(v2->GetVXPos(),v2->GetVYPos(),v2->GetVZPos());
    Vec3D DP = P2-P1;
    
    for (int i=0;i<3;i++)
    if(fabs(DP(i))>(*Box)(i)/2.0)
    {
        if(DP(i)<0)
            DP(i)=(*Box)(i)+DP(i);
        else if(DP(i)>0)
            DP(i)=DP(i)-(*Box)(i);
    }
    
    double dist2 = P1.dot(DP,DP);
    
    
    
    if(dist2<lmin || dist2>lmax)
        isvalid = false;
    if( n1.dot(n1,n2)<minangle)
        isvalid = false;
    
    if( n1.dot(n1,n3)<minangle)
        isvalid = false;

    return isvalid;
    
}
