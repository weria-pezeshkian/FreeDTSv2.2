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
    std::vector<triangle> t  = FindThreeEdgedLoop(m_pMESH);
    std::cout<<"  Three_Edge_Scission move "<<t.size()<<"\n";


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
std::vector<triangle> Three_Edge_Scission::FindThreeEdgedLoop(MESH* mesh){
    
    std::vector<triangle> outlist;
    
    std::vector<vertex*>  allv =  mesh->m_pActiveV;
    for (std::vector<vertex*>::iterator it = allv.begin() ; it != allv.end(); ++it){
        std::vector<vertex*>  n1_v = (*it)->GetVNeighbourVertex();
        if((*it)->m_VertexType==0)
        for (std::vector<vertex*>::iterator it1 = n1_v.begin() ; it1 != n1_v.end(); ++it1){
            std::vector<vertex*>  n2_v = (*it1)->GetVNeighbourVertex();
            for (std::vector<vertex*>::iterator it2 = n2_v.begin() ; it2 != n2_v.end(); ++it2){
                std::vector<vertex*>  n3_v = (*it2)->GetVNeighbourVertex();
                for (std::vector<vertex*>::iterator it3 = n2_v.begin() ; it3 != n2_v.end(); ++it3){
                    std::vector<vertex*>  n4_v = (*it3)->GetVNeighbourVertex();
                    if((*it)->GetVID()==(*it3)->GetVID()){
                        // check if the triangle exist
                        bool t_exsit = false;
                        std::vector <links *> llist = (*it)->GetVLinkList();
                        for (std::vector<links*>::iterator itl1 = llist.begin() ; itl1 != llist.end(); ++itl1){
                            if((*itl1)->GetV3()==(*it2) || (*itl1)->GetV3()==(*it1)){
                                t_exsit = true;
                            }
                        }
                        if(t_exsit==false && (*it1)->m_VertexType==0 && (*it2)->m_VertexType==0)
                        {
                            // (*it), (*it1), (*it2) are our vertices
                            triangle T (0, (*it), (*it1), (*it2));
                            outlist.push_back(T);
                        }
                    }

                }
            }
        }
    }

    
    return outlist;
}
