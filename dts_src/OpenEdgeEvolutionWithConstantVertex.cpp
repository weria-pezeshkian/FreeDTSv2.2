#include "OpenEdgeEvolutionWithConstantVertex.h"
#include "State.h"
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian

 Edge treatment, this is a new development since June 2023;
 
 What it does:
 
 // 1. We need to give error if there is a edge; we cannot have osmotic pressure ...
 // 2. We can create a hole closer function; instead of invoking the making links twice 
 
 */

OpenEdgeEvolutionWithConstantVertex::OpenEdgeEvolutionWithConstantVertex(int period, double rate, State *pState ) :
                                                    m_pState(pState),
                                                    m_pEdgeL(pState->GetMesh()->GetEdgeL()),
                                                    m_pGhostL(pState->GetMesh()->GetGhostL()),
                                                    m_pGhostT(pState->GetMesh()->GetGhostT()),
                                                    m_pActiveT(pState->GetMesh()->GetActiveT()),
                                                    m_pRightL(pState->GetMesh()->GetRightL()),
                                                    m_pLeftL(pState->GetMesh()->GetLeftL()),
                                                    m_pActiveL(pState->GetMesh()->GetActiveL()),
                                                    m_pSurfV(pState->GetMesh()->GetSurfV()),
                                                    m_pEdgeV(pState->GetMesh()->GetEdgeV()),
                                                    m_Period(period),
                                                    m_NumberOfMovePerStep(rate),
                                                    m_Beta(pState->GetSimulation()->GetBeta()),
                                                    m_DBeta(pState->GetSimulation()->GetDBeta()),
                                                    m_MinLength2(pState->GetSimulation()->GetMinL2()),
                                                    m_MaxLength2(pState->GetSimulation()->GetMaxL2()),
                                                    m_MinAngle(pState->GetSimulation()->GetMinAngle())
{
    m_EdgeSize = m_pEdgeL.size();
}
OpenEdgeEvolutionWithConstantVertex::~OpenEdgeEvolutionWithConstantVertex(){
    
}
void OpenEdgeEvolutionWithConstantVertex::Initialize(){
    
    m_pBox = m_pState->GetMesh()->GetBox();
    m_pMesh = m_pState->GetMesh();
    m_NumberOfAttemptedMoves = 0;
    m_AcceptedMoves = 0;
    
    m_EdgeSize = m_pEdgeL.size();

    return;
}
bool OpenEdgeEvolutionWithConstantVertex::Move(int step) {

    if(step % m_Period != 0 || m_pEdgeL.size() == 0 )
        return false;
   
    int N = m_pEdgeL.size();
    N = int(m_NumberOfMovePerStep*double(N));
    for (int i = 0; i< N ;i++) {

        if(m_pEdgeL.size() == 0 )
            break;
        if(MCAttemptedToAddALink()){
            m_AcceptedMoves++;
        }
        if(MCAttemptedToRemoveALink()){
            m_AcceptedMoves++;
        }
        m_NumberOfAttemptedMoves++;
        m_NumberOfAttemptedMoves++;
    }
    m_EdgeSize = m_pEdgeL.size();

    return true;
}
bool OpenEdgeEvolutionWithConstantVertex::MCAttemptedToRemoveALink(){
    if( m_pEdgeL.size() == 0 )
        return false;
    
    int n = m_pState->GetRandomNumberGenerator()->IntRNG(m_pEdgeL.size());
    links *plink = m_pEdgeL[n];
    
    vertex *v1 = plink->GetV1();
    vertex *v2 = plink->GetV2();
    vertex *v3 = plink->GetV3();

    if(v1->GetVLinkList().size() < 2 || v2->GetVLinkList().size() < 2 || v3->GetVertexType() == 1)
        return false;
        
    double eold = 0;
    double enew = 0;
    
        // calculate the old energy
        eold = v1->GetEnergy();
        eold += v2->GetEnergy();
        eold += v3->GetEnergy();

        //
        {
        std::vector <links *> nvl1 = v1->GetVLinkList();
        for (std::vector<links *>::iterator it = nvl1.begin() ; it != nvl1.end(); ++it)
            eold += 2*(*it)->GetIntEnergy();
        
        std::vector <links *> nvl2 = v2->GetVLinkList();
        for (std::vector<links *>::iterator it = nvl2.begin() ; it != nvl2.end(); ++it)
            eold += 2*(*it)->GetIntEnergy();
        
        std::vector <links *> nvl3 = v3->GetVLinkList();
        for (std::vector<links *>::iterator it = nvl3.begin() ; it != nvl3.end(); ++it)
            eold += 2*(*it)->GetIntEnergy();
            
            //--- this link does not exist in the nb of any vertices
            eold += 2*(v1->m_pPrecedingEdgeLink)->GetIntEnergy();

            // because these two links were counted two times
            eold -= 2*(plink->GetNeighborLink1())->GetIntEnergy();
            eold -= 2*(plink->GetNeighborLink2())->GetIntEnergy();
        }
        //== this is a new piece to use global inclusion direction as a constant variable during the move: note global direction should be projected on the new local plane
        if(v1->VertexOwnInclusion()){
            if(!v1->GetInclusion()->UpdateGlobalDirectionFromLocal())
                return false;
        }
        if(v2->VertexOwnInclusion()){
            if(!v2->GetInclusion()->UpdateGlobalDirectionFromLocal())
                return false;
        }
        if(v3->VertexOwnInclusion()){
            if(!v3->GetInclusion()->UpdateGlobalDirectionFromLocal())
                return false;
        }

        // we kill a link and update the geomotry
        KillALink(plink);

        // updating the local inc direction from the global one
        if (v1->VertexOwnInclusion() || v2->VertexOwnInclusion() || v3->VertexOwnInclusion()) {
            Vec3D LD1, LD2, LD3;
            if(v1->VertexOwnInclusion()){
                LD1 = (v1->GetInclusion())->GetGDirection();
                LD1 = (v1->GetG2LTransferMatrix())*LD1;
                if(LD1.isbad()){
                    CreateALink(v1);
                    return false;
                }
                LD1(2) = 0;
                LD1.normalize();
            }
            if(v2->VertexOwnInclusion()){
                LD2 = (v2->GetInclusion())->GetGDirection();
                LD2 = (v2->GetG2LTransferMatrix())*LD2;
                if(LD2.isbad()){
                    CreateALink(v1);
                    return false;
                }
                LD2(2) = 0;
                LD2.normalize();
            }
            if(v3->VertexOwnInclusion()){
                LD3 = (v3->GetInclusion())->GetGDirection();
                LD3 = (v3->GetG2LTransferMatrix())*LD3;
                if(LD3.isbad()){
                    CreateALink(v1);
                    return false;
                }
                LD3(2) = 0;
                LD3.normalize();
            }
//-- this should happen at the end, otherwise this early links get bad number
            if(v1->VertexOwnInclusion())
                v1->GetInclusion()->UpdateLocalDirection(LD1);
            if(v2->VertexOwnInclusion())
                v2->GetInclusion()->UpdateLocalDirection(LD2);
            if(v3->VertexOwnInclusion())
                v3->GetInclusion()->UpdateLocalDirection(LD3);
        }
 
        // new energy
               enew = m_pState->GetEnergyCalculator()->SingleVertexEnergy(v1);
               enew += m_pState->GetEnergyCalculator()->SingleVertexEnergy(v2);
               enew += m_pState->GetEnergyCalculator()->SingleVertexEnergy(v3);
        
        
        {
        std::vector <links *> nvl1 = v1->GetVLinkList();
        for (std::vector<links *>::iterator it = nvl1.begin() ; it != nvl1.end(); ++it)
            enew += m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(*it);
            
        std::vector <links *> nvl2 = v2->GetVLinkList();
        for (std::vector<links *>::iterator it = nvl2.begin() ; it != nvl2.end(); ++it)
            enew += m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(*it);
            
        std::vector <links *> nvl3 = v3->GetVLinkList();
        for (std::vector<links *>::iterator it = nvl3.begin() ; it != nvl3.end(); ++it)
            enew += m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(*it);
        }
        enew += m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(v1->GetPrecedingEdgeLink());
        
        double diff_energy = enew - eold;
        double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
    
        //if(double(NS)/double(NL+1)*(exp(-m_Beta*DE)>thermal ))
    if(exp( -m_Beta * diff_energy + m_DBeta) > thermal) {

            m_pState->GetEnergyCalculator()->AddToTotalEnergy(diff_energy);
        }
        else{
            // reject the move
                CreateALink(v1);
                if(v1->VertexOwnInclusion())
                {
                    Tensor2  G2L = v1->GetG2LTransferMatrix();
                    Vec3D LD = G2L*((v1->GetInclusion())->GetGDirection());
                    if(fabs(LD(2))>0.00001){
                        std::cout<<" something is wrong here, this should not happen. vector should be on the plane \n";
                    }
                    (v1->GetInclusion())->UpdateLocalDirection(LD);
                }
                if(v2->VertexOwnInclusion())
                {
                    Tensor2  G2L = v2->GetG2LTransferMatrix();
                    Vec3D LD = G2L*((v2->GetInclusion())->GetGDirection());
                    if(fabs(LD(2))>0.00001){
                        std::cout<<" something is wrong here, this should not happen. vector should be on the plane \n";
                    }
                    (v2->GetInclusion())->UpdateLocalDirection(LD);
                }
                if(v3->VertexOwnInclusion())
                {
                    Tensor2  G2L = v3->GetG2LTransferMatrix();
                    Vec3D LD = G2L*((v3->GetInclusion())->GetGDirection());
                    if(fabs(LD(2))>0.00001){
                        std::cout<<" something is wrong here, this should not happen. vector should be on the plane \n";
                    }
                    (v3->GetInclusion())->UpdateLocalDirection(LD);
                }
            
            {
            double e = m_pState->GetEnergyCalculator()->SingleVertexEnergy(v1);
            e = m_pState->GetEnergyCalculator()->SingleVertexEnergy(v2);
            e = m_pState->GetEnergyCalculator()->SingleVertexEnergy(v3);
            
            std::vector <links *> nvl1 = v1->GetVLinkList();
            for (std::vector<links *>::iterator it = nvl1.begin() ; it != nvl1.end(); ++it)
                e=m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(*it);
                
            std::vector <links *> nvl2 = v2->GetVLinkList();
            for (std::vector<links *>::iterator it = nvl2.begin() ; it != nvl2.end(); ++it)
                e=m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(*it);
                
            std::vector <links *> nvl3 = v3->GetVLinkList();
            for (std::vector<links *>::iterator it = nvl3.begin() ; it != nvl3.end(); ++it)
                e=m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(*it);
            
            e=m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(v1->m_pPrecedingEdgeLink);
            }
           // std::cout<<"killed rejected \n";

        }
        // int energy;
        //double eint = TwoInclusionsInteractionEnergy(l0);
    
    return true; // MCAttemptedToRemoveALink(){
} // MCAttemptedToRemoveALink(){

bool OpenEdgeEvolutionWithConstantVertex::MCAttemptedToAddALink(){
    
    // an atempt to create a link                       //       v3
                                                        //      /   \
                                                      //     v1-----v2
    
    
    if(m_pEdgeL.size() < 3 && m_pEdgeL.size() != 0){
        std::cout<<" error--> 123 should not happen \n";
        return false;
    }

    // select an edge vertex
    int n = m_pState->GetRandomNumberGenerator()->IntRNG(m_pEdgeV.size());
    vertex *v1 = m_pEdgeV[n];
    
    if( !Linkisvalid(v1) ){
        return false;
    }
    
    double eold = 0;
    double enew = 0;
    
    links* l2 = v1->m_pEdgeLink;
    vertex *v3 = l2->GetV2();
    links* l1 = v3->m_pEdgeLink;
    vertex *v2 = l1->GetV2();
    
//----   for constant global direction type of moves
    if(v1->VertexOwnInclusion()){
        if(!(v1->GetInclusion())->UpdateGlobalDirectionFromLocal())
            return false;
    }
    if(v2->VertexOwnInclusion()){
        if(!(v2->GetInclusion())->UpdateGlobalDirectionFromLocal())
            return false;
    }
    if(v3->VertexOwnInclusion()){
        if(!(v3->GetInclusion())->UpdateGlobalDirectionFromLocal())
            return false;
    }
    
    eold += v1->GetEnergy();
    eold += v2->GetEnergy();
    eold += v3->GetEnergy();

    //=== inclusion interaction energy; it is multipled by two because only half is assinged to each edge
    std::vector <links *> nvl1 = v1->GetVLinkList();
    for (std::vector<links *>::iterator it = nvl1.begin() ; it != nvl1.end(); ++it)
        eold += 2*(*it)->GetIntEnergy();

    std::vector <links *> nvl2 = v2->GetVLinkList();
    for (std::vector<links *>::iterator it = nvl2.begin() ; it != nvl2.end(); ++it)
        eold += 2*(*it)->GetIntEnergy();

    std::vector <links *> nvl3 = v3->GetVLinkList();
    for (std::vector<links *>::iterator it = nvl3.begin() ; it != nvl3.end(); ++it)
        eold += 2*(*it)->GetIntEnergy();

    eold += 2*(v1->m_pPrecedingEdgeLink)->GetIntEnergy();
    // create a link (this also updates the gemotry)
     links *newlink = CreateALink(v1);
  //  return true;
    //----   for constant global direction type of moves
    if (v1->VertexOwnInclusion() || v2->VertexOwnInclusion() || v3->VertexOwnInclusion()) {
        Vec3D LD1, LD2, LD3;
        if(v1->VertexOwnInclusion()){
            LD1 = (v1->GetInclusion())->GetGDirection();
            LD1 = (v1->GetG2LTransferMatrix())*LD1;
            if(LD1.isbad()){
                KillALink(newlink);
                return false;
            }
            LD1(2) = 0;
            LD1.normalize();
        }
        if(v2->VertexOwnInclusion()){
            LD2 = (v2->GetInclusion())->GetGDirection();
            LD2 = (v2->GetG2LTransferMatrix())*LD2;
            if(LD2.isbad()){
                //== reject the move
                KillALink(newlink);
                return false;
            }
            LD2(2) = 0;
            LD2.normalize();
        }
        if(v3->VertexOwnInclusion()){
            LD3 = (v3->GetInclusion())->GetGDirection();
            LD3 = (v3->GetG2LTransferMatrix())*LD3;
            if(LD3.isbad()){
                //== reject the move
                KillALink(newlink);
                return false;
            }
            LD3(2) = 0;
            LD3.normalize();
        }
//-- this should happen at the end, otherwise this early links get bad number
        if(v1->VertexOwnInclusion())
        (v1->GetInclusion())->UpdateLocalDirection(LD1);
        if(v2->VertexOwnInclusion())
        (v2->GetInclusion())->UpdateLocalDirection(LD2);
        if(v3->VertexOwnInclusion())
        (v3->GetInclusion())->UpdateLocalDirection(LD3);
    }
    
    enew = m_pState->GetEnergyCalculator()->SingleVertexEnergy(v1);
    enew += m_pState->GetEnergyCalculator()->SingleVertexEnergy(v2);
    enew += m_pState->GetEnergyCalculator()->SingleVertexEnergy(v3);

    //=== inclusion interaction energy
    for (std::vector<links *>::iterator it = nvl1.begin() ; it != nvl1.end(); ++it)
        enew += m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(*it);
    for (std::vector<links *>::iterator it = nvl2.begin() ; it != nvl2.end(); ++it)
        enew += m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(*it);
    for (std::vector<links *>::iterator it = nvl3.begin() ; it != nvl3.end(); ++it)
        enew += m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(*it);

    enew += m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(v1->m_pPrecedingEdgeLink);
        // the new created link
    enew += m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(newlink);

    double diff_energy = enew - eold;
    double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);

    if(exp(-m_Beta  *diff_energy + m_DBeta ) > thermal ) {
        
        m_pState->GetEnergyCalculator()->AddToTotalEnergy(diff_energy);
    }
    else
    {
        KillALink(newlink);
            if(v1->VertexOwnInclusion()){
                Tensor2  G2L = v1->GetG2LTransferMatrix();
                Vec3D LD = G2L*((v1->GetInclusion())->GetGDirection());
                if(fabs(LD(2))>0.00001){
                    std::cout<<" something is wrong here, this should not happen. vector should be on the plane \n";
                }
                (v1->GetInclusion())->UpdateLocalDirection(LD);
            }
            if(v2->VertexOwnInclusion()){
                Tensor2  G2L = v2->GetG2LTransferMatrix();
                Vec3D LD = G2L*((v2->GetInclusion())->GetGDirection());
                if(fabs(LD(2))>0.00001){
                    std::cout<<" something is wrong here, this should not happen. vector should be on the plane \n";
                }
                (v2->GetInclusion())->UpdateLocalDirection(LD);
            }
            if(v3->VertexOwnInclusion()){
                Tensor2  G2L = v3->GetG2LTransferMatrix();
                Vec3D LD = G2L*((v3->GetInclusion())->GetGDirection());
                if(fabs(LD(2))>0.00001){
                    std::cout<<" something is wrong here, this should not happen. vector should be on the plane \n";
                }
                (v3->GetInclusion())->UpdateLocalDirection(LD);
            }
        double e = m_pState->GetEnergyCalculator()->SingleVertexEnergy(v1);
        e = m_pState->GetEnergyCalculator()->SingleVertexEnergy(v2);
        e = m_pState->GetEnergyCalculator()->SingleVertexEnergy(v3);
        {
        for (std::vector<links *>::iterator it = nvl1.begin() ; it != nvl1.end(); ++it)
            e = m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(*it);
        for (std::vector<links *>::iterator it = nvl2.begin() ; it != nvl2.end(); ++it)
            e = m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(*it);
        for (std::vector<links *>::iterator it = nvl3.begin() ; it != nvl3.end(); ++it)
            e = m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(*it);
        }
        e = m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(v1->m_pPrecedingEdgeLink);

    }
    
    return true;
}
// an atempt to create a link                       //        v3                      v3
// between v1 and v2                                // ml3  /   \  ml2---> ml3,l3  // T \\  l2, ml2   ml3(v1,v3)  ml2(v3,v2)
                                                  //     v1      v2               v1---->v2
links* OpenEdgeEvolutionWithConstantVertex::CreateALink(vertex *v1) {           //     l1
// l1 will be created and
    
    if(m_pGhostT.size()<1 && m_pGhostL.size()<3){
        std::cout<<" error 882---> we do not have enough storage for the new trinagle and links \n";
        exit(0);
    }
    
    links* ml3 = v1->GetEdgeLink();
    vertex *v3 = ml3->GetV2();
    links* ml2 = v3->GetEdgeLink();
    vertex *v2 = ml2->GetV2();
    links *l1 = m_pGhostL[0];
    links *l2 = m_pGhostL[1];
    links *l3 = m_pGhostL[2];
    triangle *tre = m_pGhostT[0];  // one trinagule is created

//---- create the triangle
    tre->UpdateVertex(v1,v2,v3);
    m_pGhostT.erase(m_pGhostT.begin());
    AddtoTriangleList(tre, m_pActiveT);
    l1->UpdateTriangle(tre);
    l2->UpdateTriangle(tre);
    l3->UpdateTriangle(tre);
    v1->AddtoTraingleList(tre);
    v2->AddtoTraingleList(tre);
    v3->AddtoTraingleList(tre);

//----- update the links
      //-- update l1 and l2 that have became a surface edge
    RemoveFromLinkList(ml2, m_pEdgeL);
    RemoveFromLinkList(ml3, m_pEdgeL);
    AddtoLinkList(ml2, m_pRightL);
    AddtoLinkList(ml3, m_pRightL);
    ml2->m_LinkType = 0;
    ml3->m_LinkType = 0;
    l2->m_LinkType = 0;
    l3->m_LinkType = 0;
    l1->m_LinkType = 1;
    ml2->UpdateMirrorFlag(true);
    ml3->UpdateMirrorFlag(true);
    l2->UpdateMirrorFlag(true);
    l3->UpdateMirrorFlag(true);
    l1->UpdateMirrorFlag(false);
    ml2->UpdateMirrorLink(l2);
    ml3->UpdateMirrorLink(l3);
    l2->UpdateMirrorLink(ml2);
    l3->UpdateMirrorLink(ml3);
    AddtoLinkList(l2, m_pLeftL);   // since we just added ml1 to m_pRightL
    AddtoLinkList(l3, m_pLeftL);
    m_pGhostL.erase(m_pGhostL.begin());   // three times for l1,l2,l1
    m_pGhostL.erase(m_pGhostL.begin());   // three times for l1,l2,l1
    m_pGhostL.erase(m_pGhostL.begin());    // three times for l1,l2,l1
    AddtoLinkList(l1, m_pEdgeL);
    AddtoLinkList(l1, m_pActiveL);
    AddtoLinkList(l2, m_pActiveL);
    AddtoLinkList(l3, m_pActiveL);

    l1->UpdateNeighborLink1(l2);
    l1->UpdateNeighborLink2(l3);
    l2->UpdateNeighborLink1(l3);
    l2->UpdateNeighborLink2(l1);
    l3->UpdateNeighborLink1(l1);
    l3->UpdateNeighborLink2(l2);
    l1->UpdateV(v1,v2,v3);
    l2->UpdateV(v2,v3,v1);
    l3->UpdateV(v3,v1,v2);
//----update vertices
    //----update v3 that has become a surf
    RemoveFromVertexList(v3, m_pEdgeV);   // only this vertex becames a surf
    AddtoVertexList(v3, m_pSurfV);
    v3->m_VertexType = 0;                          // meaning it is not an edge vertex anymore
    //----update other vertices
    v1->m_pEdgeLink = l1;
    v2->m_pPrecedingEdgeLink = l1;
    v1->AddtoNeighbourVertex(v2);
    v2->AddtoNeighbourVertex(v1);
    v1->AddtoLinkList(l1);
    v2->AddtoLinkList(l2);
    v3->AddtoLinkList(l3);
    
 //----- now we should update the geomtry of the affected v, l, t
    tre->UpdateNormal_Area(m_pBox);   //trinagule normal and area should be found
   // l2->UpdateNormal();   // normal of the links with mirror should be updated
   // l3->UpdateNormal();   // normal of the links with mirror should be updated

    // their mirror will be updated by the function within
    l1->UpdateEdgeVector(m_pBox);   // l1 is an edge link, we need only the edge vector and length
    l2->UpdateShapeOperator(m_pBox);   // this link is now a surface link and should be found a shape operator
    l3->UpdateShapeOperator(m_pBox);   // this link is now a surface link and should be found a shape operator
    
    (m_pState->GetCurvatureCalculator())->UpdateEdgeVertexCurvature(v1);  // v1 is still an edge vertex
    (m_pState->GetCurvatureCalculator())->UpdateEdgeVertexCurvature(v2);  // // v2 is still an edge vertex
    (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(v3);  // v3 is now a surface vertex

    return l1;
}

//         v3                    v3
//  l3  // T1 \\ l2  -->  ml3  /   \  ml2
//      v1-l1- v2             v1    v2
bool OpenEdgeEvolutionWithConstantVertex::KillALink(links *l1)
{
        triangle *tri = l1->GetTriangle();
        vertex *v1 = l1->GetV1();
        vertex *v2 = l1->GetV2();
        vertex *v3 = l1->GetV3();
        links *l2 = l1->GetNeighborLink1();
        links *l3 = l2->GetNeighborLink1();
        links *ml2 = l2->GetMirrorLink();
        links *ml3 = l3->GetMirrorLink();
    
//--- remove the triangle
    RemoveFromTriangleList(tri, m_pActiveT);
    AddtoTriangleList(tri,m_pGhostT);
          //--- remove the trinagle from all the vertcies list
    v1->RemoveFromTraingleList(tri);
    v2->RemoveFromTraingleList(tri);
    v3->RemoveFromTraingleList(tri);
    
//---- remove the links l1, l2 and l3 ; note the mirors remain allive but will become an edge link
    RemoveFromLinkList(l1, m_pActiveL);  // too expensive
    RemoveFromLinkList(l2, m_pActiveL);  // too expensive
    RemoveFromLinkList(l3, m_pActiveL);  // too expensive
    RemoveFromLinkList(l2, m_pRightL);  // too expensive
    RemoveFromLinkList(l3, m_pRightL);  // too expensive
    RemoveFromLinkList(l2, m_pLeftL);  // too expensive
    RemoveFromLinkList(l3, m_pLeftL);  // too expensive
    RemoveFromLinkList(l1, m_pEdgeL);    // this link does not exist in m_pMHL and m_pHL
    AddtoLinkList(l1,m_pGhostL);
    AddtoLinkList(l2,m_pGhostL);
    AddtoLinkList(l3,m_pGhostL);
    v1->RemoveFromLinkList(l1);
    v2->RemoveFromLinkList(l2);
    v3->RemoveFromLinkList(l3);
    
//--- convert the links into edge links ml2 and ml3
    RemoveFromLinkList(ml2, m_pRightL);  // too expensive; I should find a better way
    RemoveFromLinkList(ml3, m_pRightL);  // too expensive
    RemoveFromLinkList(ml2, m_pLeftL);  // too expensive
    RemoveFromLinkList(ml3, m_pLeftL);  // too expensive
    ml2->UpdateMirrorFlag(false);
    ml3->UpdateMirrorFlag(false);
    ml2->m_LinkType = 1;
    ml3->m_LinkType = 1;
    AddtoLinkList(ml2, m_pEdgeL);    // adding the two mirror into the edge
    AddtoLinkList(ml3, m_pEdgeL);   // adding the two mirror into the edge

//--- convert the three vertices into edge vertex,
        //--- only v2 needs an update: v1 and v3 are already an edge vertex.
    RemoveFromVertexList(v3, m_pSurfV);  // too expensive
    AddtoVertexList(v3, m_pEdgeV);
    v3->m_VertexType = 1; //
       //-- v1 and v2 are no longer connected
    v1->RemoveFromNeighbourVertex(v2);
    v2->RemoveFromNeighbourVertex(v1);
       // now we updated the vertex edge links and preceding link; note: v3 edge links will not be changed and also v1  preceding link will not be changed
        v1->m_pEdgeLink = ml3;
        v3->m_pEdgeLink = ml2;
        v3->m_pPrecedingEdgeLink = ml3;
        v2->m_pPrecedingEdgeLink = ml2 ;
//--- now we should update the geomtry of the affected v, l
        ml2->UpdateEdgeVector(m_pBox);   // edge vector should be updated
        ml3->UpdateEdgeVector(m_pBox);   // edge vector should be updated
  
        (m_pState->GetCurvatureCalculator())->UpdateEdgeVertexCurvature(v1);  // v1 is still an edge vertex
        (m_pState->GetCurvatureCalculator())->UpdateEdgeVertexCurvature(v2);  // // v2 is still an edge vertex
        (m_pState->GetCurvatureCalculator())->UpdateEdgeVertexCurvature(v3);  // v3 is now an edge vertex
    
    return true;
    
}
double  OpenEdgeEvolutionWithConstantVertex::SystemEnergy() {
    double en = 0;
    /*
    std::vector<vertex *> ActiveV =  m_pSurfV;
    std::vector<triangle *> pActiveT =  m_pActiveT;
    std::vector<links *> mLink =  m_pHL;
    std::vector<links *>  pEdgeL =  m_pEdgeL;
    std::vector<vertex *> EdgeV  =  m_pEdgeV;
    

    for (std::vector<triangle *>::iterator it = pActiveT.begin() ; it != pActiveT.end(); ++it)
        (*it)->UpdateNormal_Area(m_pBox);
    
    
    for (std::vector<links *>::iterator it = ( m_pHL).begin() ; it != ( m_pHL).end(); ++it){
        (*it)->UpdateNormal();
        (*it)->UpdateShapeOperator(m_pBox);
    }
    
    for (std::vector<vertex *>::iterator it = ( m_pSurfV).begin() ; it != ( m_pSurfV).end(); ++it)
        (m_pState->GetCurvatureCalculator())->UpdateSurfVertexCurvature(*it);

    //====== edge links should be updated
    for (std::vector<links *>::iterator it = pEdgeL.begin() ; it != pEdgeL.end(); ++it)
            (*it)->UpdateEdgeVector(m_pBox);

    for (std::vector<vertex *>::iterator it = EdgeV.begin() ; it != EdgeV.end(); ++it)
            (m_pState->GetCurvatureCalculator())->UpdateEdgeVertexCurvature(*it);
    
    en = m_pState->GetEnergyCalculator()->CalculateAllLocalEnergy();
    */
    return en;
}
bool OpenEdgeEvolutionWithConstantVertex::Linkisvalid(vertex *v1) {
    // Check if the new link length is within the allowed range and if the angle of the new triangle
    // is acceptable with respect to the two other triangles.
    //           va -- v3----vb
    //            \n2 /n1\n3/
    //             v1--- v2
    //                DP
    
    // check if the new link length is within the allowed range and also if the angle of the new trinagule is fine with respect to the two other trinagules
    links* l2   = v1->GetEdgeLink();
    vertex *v3  = l2->GetV2();
    links* l1   = v3->GetEdgeLink();
    vertex *v2  = l1->GetV2();


    Vec3D DP = v2->GetPos() - v1->GetPos();
    
    for (int i=0;i<3;i++)
    if(fabs(DP(i)) > (*m_pBox)(i)/2.0)
    {
        if(DP(i) < 0)
            DP(i) = (*m_pBox)(i)+DP(i);
        else if(DP(i) > 0)
            DP(i) = DP(i)-(*m_pBox)(i);
    }
    double dist2 = DP.dot(DP,DP);
   // std::cout<<dist2<<"\n";
    if(dist2 < m_MinLength2 || dist2 > m_MaxLength2){
        return false;
    }
    
    triangle t(0,v1,v2,v3);
    t.UpdateNormal_Area(m_pBox);
    Vec3D n1 = t.GetNormalVector();
    Vec3D n2 = l2->GetTriangle()->GetNormalVector();
    Vec3D n3 = l1->GetTriangle()->GetNormalVector();
    
    if( n1.dot(n1,n2) < m_MinAngle){
        return false;
    }
    if( n1.dot(n1,n3) < m_MinAngle){
        return false;
    }
    
    return true;
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
std::string OpenEdgeEvolutionWithConstantVertex::CurrentState(){
    
    std::string state = AbstractOpenEdgeEvolution::GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    state = state + +" "+ Nfunction::D2S(m_Period) +" "+ Nfunction::D2S(m_NumberOfMovePerStep);

    return state;
}

