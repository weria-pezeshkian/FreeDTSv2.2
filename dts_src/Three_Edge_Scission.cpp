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
    if(m_GhostT.size()<4 || m_pGhostL.size()<4){
        std::cout<<" --->note: the number of the links and trinagles in repository is not enough, restart the simulations \n";
        exit(0);
    }
    
   // auto start = std::chrono::steady_clock::now();

    // finding the pair list for cutting the trinagles
    std::vector<pair_pot_triangle> pair_list  = FindPotentialTriangles(m_pMESH);

    for (std::vector<pair_pot_triangle>::iterator it = pair_list.begin() ; it != pair_list.end(); ++it)
    {
        double enew = 0;
        double eold = ((it->PT1).pv1)->GetEnergy();
        eold+= ((it->PT1).pv2)->GetEnergy();
        eold+= ((it->PT1).pv3)->GetEnergy();
        eold+= ((it->PT2).pv1)->GetEnergy();
        eold+= ((it->PT2).pv2)->GetEnergy();
        eold+= ((it->PT2).pv3)->GetEnergy();

        if(DoAScission(*it))
        {
            std::cout<<" cut was good \n";
            (m_pState->CurvatureCalculator())->SurfVertexCurvature((it->PT1).pv1);
            (m_pState->CurvatureCalculator())->SurfVertexCurvature((it->PT1).pv2);
            (m_pState->CurvatureCalculator())->SurfVertexCurvature((it->PT1).pv3);
            (m_pState->CurvatureCalculator())->SurfVertexCurvature((it->PT2).pv1);
            (m_pState->CurvatureCalculator())->SurfVertexCurvature((it->PT2).pv2);
            (m_pState->CurvatureCalculator())->SurfVertexCurvature((it->PT2).pv3);
            
            enew+= (m_pState->GetEnergyCalculator())->SingleVertexEnergy((it->PT1).pv1);
            enew+= (m_pState->GetEnergyCalculator())->SingleVertexEnergy((it->PT1).pv2);
            enew+= (m_pState->GetEnergyCalculator())->SingleVertexEnergy((it->PT1).pv3);
            enew+= (m_pState->GetEnergyCalculator())->SingleVertexEnergy((it->PT2).pv1);
            enew+= (m_pState->GetEnergyCalculator())->SingleVertexEnergy((it->PT2).pv2);
            enew+= (m_pState->GetEnergyCalculator())->SingleVertexEnergy((it->PT2).pv3);
            
            double de = enew - eold;
            double *glo_energy=&(m_pState->m_TotEnergy);
            (*glo_energy)=(*glo_energy)+de;
        }
        else
        {
            std::cout<<" we faild to cut \n";
        }
    }
   // auto end = std::chrono::steady_clock::now();
   // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
   // std::cout << "Time taken: " << duration.count() << " milliseconds" << std::endl;
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
                            // this could be made better if we create our out delet function and delete them all at once 
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
bool Three_Edge_Scission::CorrectOrientation(pot_triangle p1,pot_triangle p2)
{
    links* ml1 = (p1.pl1)->GetMirrorLink();
    if((p1.pl1)->GetV3()==p2.pv1 || (p1.pl1)->GetV3()==p2.pv2 || (p1.pl1)->GetV3()==p2.pv3){
        // Orientation is correct, no change is needed
        return true;
    }
    else if(ml1->GetV3()==p2.pv1 || ml1->GetV3()==p2.pv2 || ml1->GetV3()==p2.pv3){

        // Orientation is not correct, we reverse it
        links* ml2 = (p1.pl2)->GetMirrorLink();
        links* ml3 = (p1.pl3)->GetMirrorLink();
        p1.pl1 = ml1;
        p1.pl2 = ml2;
        p1.pl3 = ml3;
        vertex *v1 = p1.pv1;
        vertex *v2 = p1.pv2;
        p1.pv1 = v2;
        p1.pv2 = v1;

        return true;
    }
    else {
        return false;
    }
    
    return true;
}
bool Three_Edge_Scission::DoAScission(pair_pot_triangle pair)
{
    if(m_pGhostT.size()<2){
        std::cout<<" ---> not enough reserved trinagles \n";
        return false;
    }
    pot_triangle p1 = pair.PT1;
    pot_triangle p2 = pair.PT2;
    
  if (CorrectOrientation(p1,p2) && CorrectOrientation(p2,p1)){

      //== Generate the first triangle
      CreateATriangleFromAPotentialTriangle(p1);

      //== Generate the second triangle
      CreateATriangleFromAPotentialTriangle(p2);

      //=== remove links and trinagles and vertices
      std::vector <links *> nv1_links = (p1.pv1)->GetVLinkList();
      for (std::vector<links*>::iterator it = nv1_links.begin() ; it != nv1_links.end(); it++){
      
          if((*it)->GetV2()==(p2.pv1) || (*it)->GetV2()==(p2.pv2) || (*it)->GetV2()==(p2.pv3))
          {
              RemoveFromLinkList((*it),m_pMESH->m_pActiveL);
              RemoveFromLinkList((*it),m_pMESH->m_pHL);
              RemoveFromLinkList((*it),m_pMESH->m_pMHL);
              RemoveFromLinkList((*it)->GetMirrorLink(),m_pMESH->m_pActiveL);
              RemoveFromLinkList((*it)->GetMirrorLink(),m_pMESH->m_pHL);
              RemoveFromLinkList((*it)->GetMirrorLink(),m_pMESH->m_pMHL);
              RemoveFromTriangleList((*it)->GetTriangle(),m_pMESH->m_pActiveT);
              //== specific to the vertex
              (p1.pv1)->RemoveFromLinkList(*it);
              (p1.pv1)->RemoveFromTraingleList((*it)->GetTriangle());
              (p1.pv1)->RemoveFromTraingleList(((*it)->GetMirrorLink())->GetTriangle());

              if((*it)->GetV2()==(p2.pv1))
              {
                  (p2.pv1)->RemoveFromLinkList((*it)->GetMirrorLink());
                  (p2.pv1)->RemoveFromTraingleList((*it)->GetTriangle());
                  (p2.pv1)->RemoveFromTraingleList(((*it)->GetMirrorLink())->GetTriangle());
              }
              else if((*it)->GetV2()==(p2.pv2))
              {
                  (p2.pv2)->RemoveFromLinkList((*it)->GetMirrorLink());
                  (p2.pv2)->RemoveFromTraingleList((*it)->GetTriangle());
                  (p2.pv2)->RemoveFromTraingleList(((*it)->GetMirrorLink())->GetTriangle());

              }
              else if((*it)->GetV2()==(p2.pv3))
              {
                  (p2.pv3)->RemoveFromLinkList((*it)->GetMirrorLink());
                  (p2.pv3)->RemoveFromTraingleList((*it)->GetTriangle());
                  (p2.pv3)->RemoveFromTraingleList(((*it)->GetMirrorLink())->GetTriangle());

              }
          }
      }
      std::vector <links *> nv2_links = (p1.pv2)->GetVLinkList();
      for (std::vector<links*>::iterator it = nv2_links.begin() ; it != nv2_links.end(); it++){
      
          if((*it)->GetV2()==(p2.pv1) || (*it)->GetV2()==(p2.pv2) || (*it)->GetV2()==(p2.pv3))
          {
              RemoveFromLinkList((*it),m_pMESH->m_pActiveL);
              RemoveFromLinkList((*it),m_pMESH->m_pHL);
              RemoveFromLinkList((*it),m_pMESH->m_pMHL);
              RemoveFromLinkList((*it)->GetMirrorLink(),m_pMESH->m_pActiveL);
              RemoveFromLinkList((*it)->GetMirrorLink(),m_pMESH->m_pHL);
              RemoveFromLinkList((*it)->GetMirrorLink(),m_pMESH->m_pMHL);
              RemoveFromTriangleList((*it)->GetTriangle(),m_pMESH->m_pActiveT);
              //== specific to the vertex
              (p1.pv2)->RemoveFromLinkList(*it);
              (p1.pv2)->RemoveFromTraingleList((*it)->GetTriangle());
              (p1.pv2)->RemoveFromTraingleList(((*it)->GetMirrorLink())->GetTriangle());

              if((*it)->GetV2()==(p2.pv1))
              {
                  (p2.pv1)->RemoveFromLinkList((*it)->GetMirrorLink());
                  (p2.pv1)->RemoveFromTraingleList((*it)->GetTriangle());
                  (p2.pv1)->RemoveFromTraingleList(((*it)->GetMirrorLink())->GetTriangle());
              }
              else if((*it)->GetV2()==(p2.pv2))
              {
                  (p2.pv2)->RemoveFromLinkList((*it)->GetMirrorLink());
                  (p2.pv2)->RemoveFromTraingleList((*it)->GetTriangle());
                  (p2.pv2)->RemoveFromTraingleList(((*it)->GetMirrorLink())->GetTriangle());

              }
              else if((*it)->GetV2()==(p2.pv3))
              {
                  (p2.pv3)->RemoveFromLinkList((*it)->GetMirrorLink());
                  (p2.pv3)->RemoveFromTraingleList((*it)->GetTriangle());
                  (p2.pv3)->RemoveFromTraingleList(((*it)->GetMirrorLink())->GetTriangle());

              }
          }
      }
      std::vector <links *> nv3_links = (p1.pv3)->GetVLinkList();
      for (std::vector<links*>::iterator it = nv3_links.begin() ; it != nv3_links.end(); it++){
      
          if((*it)->GetV2()==(p2.pv1) || (*it)->GetV2()==(p2.pv2) || (*it)->GetV2()==(p2.pv3))
          {
              RemoveFromLinkList((*it),m_pMESH->m_pActiveL);
              RemoveFromLinkList((*it),m_pMESH->m_pHL);
              RemoveFromLinkList((*it),m_pMESH->m_pMHL);
              RemoveFromLinkList((*it)->GetMirrorLink(),m_pMESH->m_pActiveL);
              RemoveFromLinkList((*it)->GetMirrorLink(),m_pMESH->m_pHL);
              RemoveFromLinkList((*it)->GetMirrorLink(),m_pMESH->m_pMHL);
              RemoveFromTriangleList((*it)->GetTriangle(),m_pMESH->m_pActiveT);
              
              //== specific to the vertex
              (p1.pv3)->RemoveFromLinkList(*it);
              (p1.pv3)->RemoveFromTraingleList((*it)->GetTriangle());
              (p1.pv3)->RemoveFromTraingleList(((*it)->GetMirrorLink())->GetTriangle());

              if((*it)->GetV2()==(p2.pv1))
              {
                  (p2.pv1)->RemoveFromLinkList((*it)->GetMirrorLink());
                  (p2.pv1)->RemoveFromTraingleList((*it)->GetTriangle());
                  (p2.pv1)->RemoveFromTraingleList(((*it)->GetMirrorLink())->GetTriangle());
              }
              else if((*it)->GetV2()==(p2.pv2))
              {
                  (p2.pv2)->RemoveFromLinkList((*it)->GetMirrorLink());
                  (p2.pv2)->RemoveFromTraingleList((*it)->GetTriangle());
                  (p2.pv2)->RemoveFromTraingleList(((*it)->GetMirrorLink())->GetTriangle());

              }
              else if((*it)->GetV2()==(p2.pv3))
              {
                  (p2.pv3)->RemoveFromLinkList((*it)->GetMirrorLink());
                  (p2.pv3)->RemoveFromTraingleList((*it)->GetTriangle());
                  (p2.pv3)->RemoveFromTraingleList(((*it)->GetMirrorLink())->GetTriangle());

              }
          }
      }
     /*
      std::cout<<" this should be deleted \n";
      p2 = pair.PT1;
      p1 = pair.PT2;
      
      std::vector <triangle *> linkss = (p2.pv2)->GetVTraingleList();
      for (std::vector<triangle*>::iterator it = linkss.begin() ; it != linkss.end(); it++){
      
          if((*it)->GetV1()==(p1.pv1) || (*it)->GetV1()==(p1.pv2) || (*it)->GetV1()==(p1.pv3))
          {
              std::cout<<" 1this should not happen \n";
          }
          if((*it)->GetV2()==(p1.pv1) || (*it)->GetV2()==(p1.pv2) || (*it)->GetV2()==(p1.pv3))
          {
              std::cout<<" 2this should not happen \n";
          }
          if((*it)->GetV3()==(p1.pv1) || (*it)->GetV3()==(p1.pv2) || (*it)->GetV3()==(p1.pv3))
          {
              std::cout<<" 3this should not happen \n";
          }
          
      }*/
      
      //================
      (p1.pv1)->RemoveFromNeighbourVertex(p2.pv1);
      (p1.pv1)->RemoveFromNeighbourVertex(p2.pv2);
      (p1.pv1)->RemoveFromNeighbourVertex(p2.pv3);
      (p1.pv2)->RemoveFromNeighbourVertex(p2.pv1);
      (p1.pv2)->RemoveFromNeighbourVertex(p2.pv2);
      (p1.pv2)->RemoveFromNeighbourVertex(p2.pv3);
      (p1.pv3)->RemoveFromNeighbourVertex(p2.pv1);
      (p1.pv3)->RemoveFromNeighbourVertex(p2.pv2);
      (p1.pv3)->RemoveFromNeighbourVertex(p2.pv3);
      
      (p2.pv1)->RemoveFromNeighbourVertex(p1.pv1);
      (p2.pv1)->RemoveFromNeighbourVertex(p1.pv2);
      (p2.pv1)->RemoveFromNeighbourVertex(p1.pv3);
      (p2.pv2)->RemoveFromNeighbourVertex(p1.pv1);
      (p2.pv2)->RemoveFromNeighbourVertex(p1.pv2);
      (p2.pv2)->RemoveFromNeighbourVertex(p1.pv3);
      (p2.pv3)->RemoveFromNeighbourVertex(p1.pv1);
      (p2.pv3)->RemoveFromNeighbourVertex(p1.pv2);
      (p2.pv3)->RemoveFromNeighbourVertex(p1.pv3);

    //=== update links shape Operator
      (p1.pl1)->UpdateShapeOperator(m_pBox);
      (p1.pl2)->UpdateShapeOperator(m_pBox);
      (p1.pl3)->UpdateShapeOperator(m_pBox);
      (p2.pl1)->UpdateShapeOperator(m_pBox);
      (p2.pl2)->UpdateShapeOperator(m_pBox);
      (p2.pl3)->UpdateShapeOperator(m_pBox);
      
    }
    else
    {
        std::cout<<" not expected \n";
    }


    

    // find all the links and trinagles that is shared between these two
    // find the right orinatation of the two trinagles
    // create two new triangles
    
    return true;
}
bool Three_Edge_Scission::CreateATriangleFromAPotentialTriangle(pot_triangle p1)
{
    triangle *gt1 = m_pGhostT[m_pGhostT.size()-1];
    triangle TmT((m_pMESH->m_pActiveT).size(),p1.pv1,p1.pv2,p1.pv3);
    *gt1 = TmT;
    (p1.pl1)->UpdateTriangle(gt1);
    (p1.pl2)->UpdateTriangle(gt1);
    (p1.pl3)->UpdateTriangle(gt1);
    p1.pl1->UpdateNeighborLink1(p1.pl2);
    p1.pl1->UpdateNeighborLink2(p1.pl3);
    p1.pl2->UpdateNeighborLink1(p1.pl3);
    p1.pl2->UpdateNeighborLink2(p1.pl1);
    p1.pl3->UpdateNeighborLink1(p1.pl1);
    p1.pl3->UpdateNeighborLink2(p1.pl2);
    p1.pl1->UpdateV3(p1.pv3);
    p1.pl2->UpdateV3(p1.pv1);
    p1.pl3->UpdateV3(p1.pv2);
    (m_pMESH->m_pActiveT).push_back(gt1);
    m_pGhostT.pop_back();
    p1.pv1->AddtoTraingleList(gt1);
    p1.pv2->AddtoTraingleList(gt1);
    p1.pv3->AddtoTraingleList(gt1);
    gt1->UpdateNormal_Area(m_pBox);   //trinagule normal and area should be obtained

    return true;
}
bool Three_Edge_Scission::DoAFussion(pair_pot_triangle pair)
{
    pot_triangle p1 = pair.PT1;
    pot_triangle p2 = pair.PT2;
    
    
    return false;
}
void Three_Edge_Scission::RemoveFromLinkList(links* z, std::vector<links*> &vect)
{
    vect.erase(std::remove(vect.begin(), vect.end(), z), vect.end());
}
void Three_Edge_Scission::RemoveFromTriangleList(triangle* z, std::vector<triangle*> &vect)
{
    vect.erase(std::remove(vect.begin(), vect.end(), z), vect.end());
}
void Three_Edge_Scission::AddtoLinkList(links* z, std::vector<links*> &vect)
{
    vect.push_back(z);
}
void Three_Edge_Scission::AddtoTriangleList(triangle* z, std::vector<triangle*> &vect)
{
    vect.push_back(z);
}
