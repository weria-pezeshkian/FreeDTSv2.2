#include <chrono>
#include "Three_Edge_Scission.h"
#include "State.h"

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian

 Edge treatment, this is a new development since June 2023;
 What it does:
 
 // 1. We need to give error if there is a edge; we cannot have osmotic pressure ...
 // 2. We can create a hole closer function; instead of invoking the making links twice 
 
 */
Three_Edge_Scission::Three_Edge_Scission()
{

}
Three_Edge_Scission::Three_Edge_Scission(int period, State *pState)
{

    m_pState = pState;
    m_Beta =   m_pState->m_Beta;
    m_pEnergyCalculator = pState->GetEnergyCalculator();
}
bool Three_Edge_Scission::MCMove(double * TotalEnergy, double temp, GenerateCNTCells *pGenCNT )
{
    
    // finding the pair list for cutting the trinagles
    std::vector<pair_pot_triangle> pair_list  = FindPotentialTriangles(m_pMESH);

    for (std::vector<pair_pot_triangle>::iterator it = pair_list.begin() ; it != pair_list.end(); ++it)
    {
        if(DoAScission(*it))
        {
            std::cout<<" cut was good \n";
        }
        else
        {
            std::cout<<" we faild to cut \n";
        }
    }

    return false;
}
Three_Edge_Scission::~Three_Edge_Scission()
{
    
}
void Three_Edge_Scission::initialize()
{
  // creating some trinagle and links to later create from them
    
    std::cout<<"  Three_Edge_Scission initalized \n";
    
    
    m_pMESH = m_pState->m_pMesh;
    m_pBox = m_pMESH->m_pBox;
    // ==== 
    int lid = 2*((m_pMESH->m_pMHL).size())+(m_pMESH->m_pEdgeL).size();
    int tid = (m_pMESH->m_pActiveT).size() ;

    
    for (int i=0;i<30;i++)
    {
        links teml(lid);
        m_GhostL.push_back(teml);
        lid++;
    }
    for (int i=0;i<30;i++)
    {
        triangle temt(tid);
        m_GhostT.push_back(temt);
        tid++;
    }
    for (std::vector<links>::iterator it = m_GhostL.begin() ; it != m_GhostL.end(); ++it)
        m_pGhostL.push_back(&(*it));
    
    for (std::vector<triangle>::iterator it = m_GhostT.begin() ; it != m_GhostT.end(); ++it)
        m_pGhostT.push_back(&(*it));
    


}
void Three_Edge_Scission::RemoveFromLinkList(links* z, std::vector<links*> &vect)
{
    vect.erase(std::remove(vect.begin(), vect.end(), z), vect.end());
}
void Three_Edge_Scission::RemoveFromVertexList(vertex* z, std::vector<vertex*> &vect)
{
    vect.erase(std::remove(vect.begin(), vect.end(), z), vect.end());
}
void Three_Edge_Scission::RemoveFromTriangleList(triangle* z, std::vector<triangle*> &vect)
{
    vect.erase(std::remove(vect.begin(), vect.end(), z), vect.end());
}
// this can be written as template
void Three_Edge_Scission::AddtoLinkList(links* z, std::vector<links*> &vect)
{
    vect.push_back(z);
}
void Three_Edge_Scission::AddtoVertexList(vertex* z, std::vector<vertex*> &vect)
{
    vect.push_back(z);
}
void Three_Edge_Scission::AddtoTriangleList(triangle* z, std::vector<triangle*> &vect)
{
    vect.push_back(z);
    return;
}

//== second method
std::vector<pair_pot_triangle> Three_Edge_Scission::FindPotentialTriangles(MESH* mesh){
    
    int id = 0;
    std::vector<pot_triangle> list;
    std::vector<links*>  all_link =  mesh->m_pActiveL;
    for (std::vector<links*>::iterator il1 = all_link.begin() ; il1 != all_link.end(); il1++){
        vertex *pv1 = (*il1)->GetV1();
        vertex *pv2 = (*il1)->GetV2();
        std::vector<links*>  v2_l = pv2->GetVLinkList();
        for (std::vector<links*>::iterator il2 = v2_l.begin() ; il2 != v2_l.end(); il2++){
            vertex *pv3 = (*il2)->GetV2();
            std::vector<links*>  v3_l = pv3->GetVLinkList();
                for (std::vector<links*>::iterator il3 = v3_l.begin() ; il3 != v3_l.end(); il3++){
                    if(pv1==(*il3)->GetV2() && (*il3)->GetV3()!=pv2 && ((*il3)->GetMirrorLink())->GetV3()!=pv2 ){
                        if( pv2->m_VertexType==0 && pv3->m_VertexType==0 && pv1->m_VertexType==0){
                            // pv1, pv2, pv3 are our vertices
                                                    
                            pot_triangle PotT;
                            PotT.id = id; id++;
                            PotT.cid = -1;
                            PotT.pv1 = pv1; PotT.pv2 = pv2; PotT.pv3 = pv3;
                            PotT.pl1= (*il1); PotT.pl2= (*il2); PotT.pl3= (*il3);
                            list.push_back(PotT);
                            //std::cout<<pv1->GetVID()<<"   "<<pv2->GetVID()<<"   "<<pv3->GetVID()<<"   \n";
                            all_link.erase(std::remove(all_link.begin(), all_link.end(), (*il2)), all_link.end());
                            all_link.erase(std::remove(all_link.begin(), all_link.end(), (*il3)), all_link.end());
                            all_link.erase(std::remove(all_link.begin(), all_link.end(), (*il2)->GetMirrorLink()), all_link.end());
                            all_link.erase(std::remove(all_link.begin(), all_link.end(), (*il3)->GetMirrorLink()), all_link.end());
                            all_link.erase(std::remove(all_link.begin(), all_link.end(), (*il1)->GetMirrorLink()), all_link.end());

                        }
                    }
            }
        }
    }
//===================
    std::vector <pair_pot_triangle> pair_list;
    int idpair=0;
for (int i=0;i<list.size();i++)
{
    bool isgood = false;
    for (int j=i+1;j<list.size();j++)
    {
        if(connected_2pot_triangles(list[i],list[j]))
        {
            pair_pot_triangle temp;
            temp.PT1 = list[i];
            temp.PT2 = list[j];
            temp.id = idpair;
            idpair++;
            pair_list.push_back(temp);
        }
    }
}
   //==== shoudl be removed
    for (int i=0;i<pair_list.size();i++)
    {
        pot_triangle p1 = (pair_list[i]).PT1;
        pot_triangle p2 = (pair_list[i]).PT2;

        (p1.pv1)->UpdateGroup(1+(p1.pv1)->GetGroup());
        (p1.pv2)->UpdateGroup(1+(p1.pv2)->GetGroup());
        (p1.pv3)->UpdateGroup(1+(p1.pv3)->GetGroup());
        (p2.pv1)->UpdateGroup(2+(p2.pv1)->GetGroup());
        (p2.pv2)->UpdateGroup(2+(p2.pv2)->GetGroup());
        (p2.pv3)->UpdateGroup(2+(p2.pv3)->GetGroup());
    }
    //===
    
    return pair_list;
}
bool Three_Edge_Scission::connected_2pot_triangles(pot_triangle potT1, pot_triangle potT2)
{
    if(potT2.cid != -1 || potT1.cid != -1)
        return false;
    if(potT2.pv1 == potT1.pv1 || potT2.pv1 == potT1.pv2 || potT2.pv1 == potT1.pv3)
        return false;
    if(potT2.pv2 == potT1.pv1 || potT2.pv2 == potT1.pv2 || potT2.pv2 == potT1.pv3)
        return false;
    if(potT2.pv3 == potT1.pv1 || potT2.pv3 == potT1.pv2 || potT2.pv3 == potT1.pv3)
        return false;
    
    
    int test1 = 0;
    int test2 = 0;
    int test3 = 0;
    std::vector <vertex *> nv1 = (potT1.pv1)->GetVNeighbourVertex();
    std::vector <vertex *> nv2 = (potT1.pv2)->GetVNeighbourVertex();
    std::vector <vertex *> nv3 = (potT1.pv3)->GetVNeighbourVertex();
    int connect[3][3] = {0};
    for (std::vector<vertex*>::iterator it = nv1.begin() ; it != nv1.end(); it++){
        
        if(potT2.pv1 == (*it) || potT2.pv1 == (*it) || potT2.pv1 == (*it))
            test1 = true;
        if(potT2.pv1 == (*it))
            connect[0][0]=1;
        if(potT2.pv1 == (*it))
            connect[1][0]=1;
        if(potT2.pv1 == (*it))
            connect[2][0]=1;
        
    }
    for (std::vector<vertex*>::iterator it = nv2.begin() ; it != nv2.end(); it++){
       
        if(potT2.pv1 == (*it) || potT2.pv1 == (*it) || potT2.pv1 == (*it))
            test2 = true;
        if(potT2.pv1 == (*it))
            connect[0][1]=1;
        if(potT2.pv1 == (*it))
            connect[1][1]=1;
        if(potT2.pv1 == (*it))
            connect[2][1]=1;
    }
    for (std::vector<vertex*>::iterator it = nv3.begin() ; it != nv3.end(); it++){
        
        if(potT2.pv1 == (*it) || potT2.pv1 == (*it) || potT2.pv1 == (*it))
            test3 = true;
        if(potT2.pv1 == (*it))
            connect[0][2]=1;
        if(potT2.pv1 == (*it))
            connect[1][2]=1;
        if(potT2.pv1 == (*it))
            connect[2][2]=1;
    }

    
    if(test1==true && test2==true && test3==true )
    {
        for (int i=0;i<3;i++){
            int row = 0;
            for (int j=0;j<3;j++){
                row+= connect[j][i];
            }
            if(row==0)
                return false;
        }
        
        potT2.cid = potT1.id;
        potT1.cid = potT2.id;
        return true;
    }
    
    return false;
}
bool Three_Edge_Scission::DoAScission(pair_pot_triangle pair)
{
    pot_triangle p1 = pair.PT1;
    pot_triangle p2 = pair.PT2;
    
    // find all the links and trinagles that is shared between these two
    // find the right orinatation of the two trinagles
    // create two new triangles
    
    return true;
}
bool Three_Edge_Scission::DoAFussion(pair_pot_triangle pair)
{
    pot_triangle p1 = pair.PT1;
    pot_triangle p2 = pair.PT2;
    
    
    return false;
}
