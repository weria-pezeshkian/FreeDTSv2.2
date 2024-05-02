


#include "CurvatureByShapeOperatorType1.h"
#include "Tensor2.h"
CurvatureByShapeOperatorType1::CurvatureByShapeOperatorType1(){
    
}
CurvatureByShapeOperatorType1::~CurvatureByShapeOperatorType1(){

}
bool CurvatureByShapeOperatorType1::VertexCurvature(vertex *pvertex){
 
    if(pvertex->GetVertexType()==0){
        return SurfVertexCurvature(pvertex);
    }
    else if(pvertex->GetVertexType()==1){
        return EdgeVertexCurvature(pvertex);
    }
    else{
        std::cout<<" vertex type "<<pvertex->GetVertexType()<<" is unrecognized \n";
        return false;
    }
    
    return true;
}
bool CurvatureByShapeOperatorType1::SurfVertexCurvature(vertex * pvertex){
    
    std::vector<triangle *> Ntr=pvertex->GetVTraingleList();
    double Area=0.0;
    Vec3D Normal;
    for (std::vector<triangle *>::iterator it = Ntr.begin() ; it != Ntr.end(); ++it)
    {
        const Vec3D& Nv = (*it)->GetAreaVector();
        Normal=Normal+Nv;
        Area+=(*it)->GetArea();
    }
    Area=Area/3.0;
    if(Area<=0){
        std::string sms=" ---> error: vertex has a negetive or zero area \n";
        std::cout<<sms<<"\n";
        return false;
    }
    double normalsize=Normal.norm();
    Normal=Normal*(1.0/normalsize);
    pvertex->UpdateNormal_Area(Normal,Area);
        ///=======
    //=== Shape Operator
    //========
    Tensor2  SV;
    Tensor2 IT('I');
    
    Tensor2 P=IT-IT.makeTen(Normal);

    Tensor2 Pt=P.Transpose(P);
    std::vector<links *> NLinks=pvertex->GetVLinkList();

    // what if the edge vertex has one trinagle?
    for (std::vector<links *>::iterator it = NLinks.begin() ; it != NLinks.end(); ++it)
    {
       if((*it)->GetMirrorFlag())
       {
           Vec3D ve=(*it)->GetNormal();
           double we=ve.dot(Normal,ve);
           Vec3D Be=(*it)->GetBe();
           double he=(*it)->GetHe();
           Vec3D Se = P*Be;

           // ff should be 1 but just for sake of numerical errors
           double ff=Se.norm();
           if(ff==0)
           {
               std::cout<<"-----> Error: projection is zero error"<<"\n";
               exit(0);
           }
           else
           {
               Se=Se*(1.0/ff);
           }
        
           Tensor2 Q=P.makeTen(Se);
        
           SV=SV+(Q)*(we*he);
       }//if((*it)->GetMirrorFlag())
    }

    ///=============
    //==== Find Curvature and local frame
    //=============
    
    Tensor2 Hous = Householder(Normal);
    Tensor2 LSV;// Local SV
    LSV=(Hous.Transpose(Hous))*(SV*Hous);    // LSV is the curvature matrix in the local frame, the it is a 2*2 minor matrix since all the 3 component are zero.

    double b=LSV(0,0)+LSV(1,1);
    double c=LSV(0,0)*LSV(1,1)-LSV(1,0)*LSV(0,1);


    double delta=b*b-4*c;
    double c1,c2;
    if(delta>0.0){
        delta=sqrt(delta);
        c1=b+delta;
        c1=0.5*c1;      // c1 always will be larger then c2
        c2=b-delta;
        c2=0.5*c2;
        
    }
    else if (fabs(delta)<0.0001){
        
        c1=0.5*b;
        c2=c1;

    }
    else{
        c1=100;
        c2=100;
        std::cout<<"WARNING: faild to find curvature on vertex "<<pvertex->GetVID()<<"  because delta is "<<delta<<"  c1 and c2 are set to 100 \n";
        std::cout<<" if you face this too much, you should stop the job and .... \n";
    }

    Tensor2 EigenvMat('O');
    
    double p=LSV(0,0);
    double q=LSV(0,1);


    double size=sqrt(q*q+(c1-p)*(c1-p));                   // The Eigenvectors can be calculated using this equation LSV*R=c1*R
    EigenvMat(0,0)=q/size;                                  // only one of them needs to be calculated, one is normal vector and the other is perpendicular to first one
    EigenvMat(1,0)=(c1-p)/size;
    EigenvMat(0,1)=-EigenvMat(1,0);
    EigenvMat(1,1)=EigenvMat(0,0);
    EigenvMat(2,2)=1;
        
    Tensor2 TransferMatLG=Hous*EigenvMat;   /// This matrix transfers vectors from local coordinate to global coordinate
    Tensor2 TransferMatGL=TransferMatLG.Transpose(TransferMatLG);   /// This matrix transfers vectors from Global coordinate to local coordinate
    pvertex->UpdateL2GTransferMatrix(TransferMatLG);
    pvertex->UpdateG2LTransferMatrix(TransferMatGL);

    c1=c1/Area;
    c2=c2/Area;
    pvertex->UpdateP1Curvature(c1);
    pvertex->UpdateP2Curvature(c2);

    return true;
}
bool CurvatureByShapeOperatorType1::EdgeVertexCurvature(vertex * pvertex)
{
    // first we obtain the vertex area and normal. Area is not important here as the vertex is an edge vertex
    std::vector<triangle *> Ntr=pvertex->GetVTraingleList();
    
    double Area=0.0;
    Vec3D Normal;
    for (std::vector<triangle *>::iterator it = Ntr.begin() ; it != Ntr.end(); ++it){
        const Vec3D& Nv = (*it)->GetAreaVector();
        Normal=Normal+Nv;
        Area+=(*it)->GetArea();
    }
    Area=Area/3.0;
    // just in case; this should never happen unless the inputs are wrong
    if(Area<=0){
        std::cout<<Ntr.size()<<"\n";
        std::string sms=" error----> bad area for vertex \n";
        std::cout<<sms<<"\n";
        return false;
    }
    double normalsize=Normal.norm();
    Normal=Normal*(1.0/normalsize);
    pvertex->UpdateNormal_Area(Normal,Area);
    

    
    // the shape of the system                          //         v
                                                        //     l1 / \ l2
                                                 
    links* link1 = pvertex->m_pPrecedingEdgeLink;
    links* link2 = pvertex->m_pEdgeLink;

    double l1 = link1->m_EdgeSize;
    double l2 = link2->m_EdgeSize;
    double vlenght = 0.5*(l1+l2);
    pvertex->m_VLength = vlenght;
    Vec3D L1=link1->m_EdgeVector;
    Vec3D L2=link2->m_EdgeVector;
    L1 = L1*(1/l1);
    L2 = L2*(1/l2);
    Vec3D Norm = (L1-L2)*(1.0/vlenght);    // dT/ds = -Norm ; the size is curvature
    Vec3D Tv = Norm*Normal;  // T at each vertex
    Tv = Tv*(1/Tv.norm());
    Vec3D P = Normal*Tv;
    double cn = Norm.dot(Norm,Normal);
    double cg = Norm.dot(Norm,P);
    pvertex->m_Geodesic_Curvature = cg;
    pvertex->m_Normal_Curvature =cn;
    Tensor2 TransferMatGL(Tv,P,Normal);     // P1,P2,N is for other surfaces
    Tensor2 TransferMatLG=TransferMatGL.Transpose(TransferMatGL);
    pvertex->UpdateL2GTransferMatrix(TransferMatLG);
    pvertex->UpdateG2LTransferMatrix(TransferMatGL);
    
    return true;
    
}
Tensor2 CurvatureByShapeOperatorType1::Householder(Vec3D N){
    
    Tensor2 Hous;
    Vec3D Zk;
    Zk(2)=1.0;
    Zk=Zk+N;
    Zk=Zk*(1.0/Zk.norm());
    
    Tensor2 I('I');
    Tensor2 W=Hous.makeTen(Zk);
    Hous=(I-W*2)*(-1);
    
    return Hous;
}
Vec3D CurvatureByShapeOperatorType1::Calculate_Vertex_Normal(vertex *pvertex){
    // first we obtain the vertex area and normal.
    std::vector<triangle *> Ntr=pvertex->GetVTraingleList();
    
    double Area=0.0;
    Vec3D Normal;
    for (std::vector<triangle *>::iterator it = Ntr.begin() ; it != Ntr.end(); ++it)
    {
        const Vec3D& Nv = (*it)->GetAreaVector();
        Normal=Normal+Nv;
        Area+=(*it)->GetArea();
    }
    Area=Area/3.0;
    // just in case; this should never happen unless the inputs are wrong
    if(Area<=0)
    {
        std::string sms=" error----> vertex has a zero or negetive area \n";
        std::cout<<sms<<"\n";
        exit(0);
    }
    double normalsize=Normal.norm();
    Normal=Normal*(1.0/normalsize);
    pvertex->UpdateNormal_Area(Normal,Area);
    

    return Normal;
}




