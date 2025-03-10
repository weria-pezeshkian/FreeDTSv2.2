

#include <time.h>
#include "Energy.h"
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 Energy of a single vertex
 Energy of a link (when connected vertices has inclusions)
 Energy of the whole system
 Energy change of a vertex move
 */
Energy::Energy()
{
    
}
Energy::Energy(Inclusion_Interaction_Map * pint)
{
   m_pInt=pint;
   m_Kappa = pint->m_BendingRigidity;
   m_KappaG = pint->m_GaussianRigidity;
   m_mem_c0 = pint->m_Spontaneous_Curvature;
   m_Membrane_model_parameters = pint->m_Membrane_model_parameters;   // this is useless for now
   m_NO_Membrane_model_parameters = m_Membrane_model_parameters.size();
    
   //== edge parameter;
   m_Lambda = pint->m_Lambda;
   m_KnEdge = pint->m_KnEdge;
   m_KgEdge = pint->m_KgEdge;
   m_KvaEdge = pint->m_KvaEdge;
   m_av0Edge = pint->m_av0Edge;
  //== vertex area energy
   m_Kva = pint->m_Kva;
   m_av0 = pint->m_av0;
    
    m_FieldDirection = pint->m_FieldDirection;
    m_FieldStrength = pint->m_FieldStrength;
}
Energy::~Energy()
{
    
}
double Energy::SingleVertexEnergy(vertex *pv)
{
    double Energy=0.0;
    
if(pv->m_VertexType==0)
{
    std::vector<double> Curve=pv->GetCurvature();
    double mean=(Curve.at(0)+Curve.at(1));
    double gussian=(Curve.at(0))*(Curve.at(1));
    double area=pv->GetArea();


        if(pv->VertexOwnInclusion()==true)
        {
        inclusion *inc=pv->GetInclusion();
        InclusionType *inctype = inc->GetInclusionType();
        double k0 = inctype->ITk;
        double kg = inctype->ITkg;
        double k1 = inctype->ITk1;
        double k2 = inctype->ITk2;
        double c0 = inctype->ITc0;
        double c1 = inctype->ITc1;
        double c2 = inctype->ITc2;

            double H=(mean-c0);
            k0=k0/2.0;
            Energy+=(k0*H*H-kg*gussian)*area;         /// this means that inclsuion overwrite the vertex bending rigidity
            if(k1!=0 || k2!=0)
            {
            Vec3D LD=inc->GetLDirection();
            double Cos=LD(0);
            double Sin=LD(1);
            double C1=Curve.at(0)*Cos*Cos+Curve.at(1)*Sin*Sin;
            double C2=Curve.at(1)*Cos*Cos+Curve.at(0)*Sin*Sin;
            double H1=(C1-c1);
            double H2=(C2-c2);
            Energy+=(k1*H1*H1+k2*H2*H2)*area/2;

            }
            if(m_FieldStrength!=0)
            {
                Vec3D LD=inc->GetLDirection();
                Vec3D GD = (pv->GetL2GTransferMatrix())*LD;
                double Cangle = GD.dot(GD,m_FieldDirection);
                Energy+=-m_FieldStrength*Cangle*Cangle;
            }
        }
        else
        {
            Energy=(m_Kappa*(mean-m_mem_c0)*(mean-m_mem_c0)-m_KappaG*gussian)*area;
        }

    
    // energy for area
    if(m_Kva!=0)
    Energy+=m_Kva*(area-m_av0)*(area-m_av0);
}
else
{
    Energy = SingleEdgeVertexEnergy(pv);

}

    pv->UpdateEnergy(Energy);


#if TEST_MODE == Enabled
    if(isnan(Energy)==1)
    std::cout<<Energy<<" energy is a problem \n";
#endif
    return Energy;
}
double Energy::SingleEdgeVertexEnergy(vertex *pv)
{
    
    double Energy = 0;
    
    double gc = pv->m_Geodesic_Curvature;
    double nc = pv->m_Normal_Curvature;
    double length = pv->m_VLength;
    double area=pv->GetArea();

    
    //== this should be changed
    //m_KnEdge; m_KgEdge;
    //====

    if(pv->VertexOwnInclusion()==true)
    {
        inclusion *inc=pv->GetInclusion();
        InclusionType *inctype = inc->GetInclusionType();
        
        double lambda = inctype->ITelambda;
        double kg = inctype->ITekg;
        double kn = inctype->ITekn;
        double cn0 = inctype->ITecn;
            Energy=(lambda+ +m_KgEdge*gc*gc+m_KnEdge*nc*nc)*length;         /// this means that inclsuion overwrite the vertex bending rigidity
        if(kn!=0)  // this is not a general model, we need to modify it later
        {
            Vec3D LD=inc->GetLDirection();
            double Cos=LD(0);
            kn=kn*Cos*Cos;
            double H1=(nc-cn0);
            Energy+=(kn*H1*H1)*length;
        }
        if(m_FieldStrength!=0)
        {
            Vec3D LD=inc->GetLDirection();
            Vec3D GD = (pv->GetL2GTransferMatrix())*LD;
            double Cangle = GD.dot(GD,m_FieldDirection);
            Energy+=-m_FieldStrength*Cangle*Cangle;
        }
                
    }
    else
    {
        Energy = m_Lambda+m_KgEdge*gc*gc+m_KnEdge*nc*nc;
        Energy=Energy*length;
    }

    // energy for area
    if(m_KvaEdge!=0)
    Energy+=m_KvaEdge*(length-m_av0Edge)*(length-m_av0Edge);
    return Energy;
}

double Energy::TotalEnergy(std::vector<vertex *> pVeretx, std::vector<links *> pLinks)
{
    double E=0.0;
    
    for (std::vector<vertex *>::iterator it = pVeretx.begin() ; it != pVeretx.end(); ++it)
    {
        E+=SingleVertexEnergy(*it);
    }
    for (std::vector<links *>::iterator it = pLinks.begin() ; it != pLinks.end(); ++it)
    {
        E+=TwoInclusionsInteractionEnergy(*it);
    }
    
    return E;
}
double Energy::Energy_OneLinkFlip(links * plinks)
{
 double E=0.0;
    
    links *mirorl=plinks->GetMirrorLink();
    vertex *v1=plinks->GetV1();
    vertex *v2=plinks->GetV2();
    vertex *v3=plinks->GetV3();
    vertex *v4=mirorl->GetV3();

    E+=SingleVertexEnergy(v1);
    E+=SingleVertexEnergy(v2);
    E+=SingleVertexEnergy(v3);
    E+=SingleVertexEnergy(v4);

  return E;
    

}
double Energy::Energy_OneVertexMove(vertex * pVeretx){
// moving a vertex leads to changes in the energy of the vertex and its beighbouring vertices. The function give the total energy of such system. 
    double E=0.0;
    std::vector<vertex *> NpVer=pVeretx->GetVNeighbourVertex(); /// Get The vertexs on the ring
    E+=SingleVertexEnergy(pVeretx);
    for (std::vector<vertex *>::iterator it = NpVer.begin() ; it != NpVer.end(); ++it){
        E+=SingleVertexEnergy(*it);
    }

    return E;
}
double Energy::TwoInclusionsInteractionEnergy(links * lp)
{

    vertex * v1 = lp->GetV1();
    vertex * v2 = lp->GetV2();
    double e=0;
    
    bool has1=v1->VertexOwnInclusion();
    bool has2=v2->VertexOwnInclusion();
    
    if(has1==true && has2==true)
    {
        int id1=(v1->GetInclusion())->GetInclusionTypeID();
        int id2=(v2->GetInclusion())->GetInclusionTypeID();
        PairInt pair_ab = m_pInt->GetPairInt(id1,id2);
        std::vector <double> ff = pair_ab.Varibale;
        int FunctionType  = pair_ab.FunctionType;
        
            double theta = 0;
         if(FunctionType == 0)
         {
                e=0;
         }
          else if(FunctionType == 1)
          {
            if(ff.at(2)!=0)
            theta = Geo_Theta(v1,v2);
            e=InteractionFunction(ff.at(0), ff.at(1),ff.at(2),theta);
          }
          else if(FunctionType == 2)
          {
              m_Angle3D = 0;
              m_Angle2D = 0;
              e = F2(v1,v2,ff);
          }
         else if(FunctionType == 10)
         {
             m_Angle3D = 0;
             m_Angle2D = 0;
             e = F10(v1,v2,ff);
         }
        else if(FunctionType == 11)
        {
            m_Angle3D = 0;
            m_Angle2D = 0;
            e = F11(v1,v2,ff);
        }
         else
         {
             std::cout<<"---> Error: Unregognized function typeid -->"<<FunctionType<<std::endl;
             exit(0);
         }

    }
    
    if(lp->GetMirrorFlag()==true)      //// this check is not needed for now, for later developments
    {
        lp->UpdateIntEnergy(e/2.0);
        (lp->GetMirrorLink())->UpdateIntEnergy(e/2.0);
    }
    else
    {
        // divided by 2 here is also fine as everywhere will be muliplied by two again
        lp->UpdateIntEnergy(e/2.0);
    }
    

    
    return e;
}
double Energy::InteractionFunction(double N, double A, double B, double theta)
{
    double e = 0;


        e=1+cos(double(N)*theta);
        e=-A+B*e;

    return e;
}
double Energy::Geo_Theta(vertex *v1, vertex *v2)
{

    double theta;
        
        Vec3D X1(v1->GetVXPos(),v1->GetVYPos(),v1->GetVZPos());
        Vec3D X2(v2->GetVXPos(),v2->GetVYPos(),v2->GetVZPos());
        Vec3D geodesic_dir=(X2-X1);
        Vec3D *pBox=v1->GetBox();
        
        for (int i=0;i<3;i++)
        {
            if(fabs(geodesic_dir(i))>(*pBox)(i)/2)
            {
                if(geodesic_dir(i)<0)
                    geodesic_dir(i)=geodesic_dir(i)+(*pBox)(i);
                else if(geodesic_dir(i)>0)
                    geodesic_dir(i)=geodesic_dir(i)-(*pBox)(i);
            }
        }

        
        Vec3D y1=(v1->GetG2LTransferMatrix())*geodesic_dir;
        y1(2)=0;
        y1=y1*(1/(y1.norm()));
        Vec3D y2=(v2->GetG2LTransferMatrix())*geodesic_dir;
        y2(2)=0;
        y2=y2*(1/(y2.norm()));
        Vec3D n(0,0,1);
        Vec3D d1 = (v1->GetInclusion())->GetLDirection();
        Vec3D d2 = (v2->GetInclusion())->GetLDirection();
        double cos1 = y1.dot(y1,d1);
        double sin1 = n.dot(n*y1,d1);
        double cos2 = y1.dot(y2,d2);
        double sin2 = n.dot(n*y2,d2);
        double S_an = sin1*sin2+cos1*cos2;
         theta=acos(S_an);
        
    
    
    return theta;
}
double Energy::F2(vertex *v1, vertex *v2, std::vector<double> var)
{
    /// F = -e0+e1*cosN(phi-phi0)+e2*(l/l0-1)cos(beta-beta0)
    double E = 0;
    m_Angle2D = Geo_Theta(v1,v2);
    Vec3D N1 = v1->GetNormalVector();
    Vec3D N2 = v2->GetNormalVector();
    double beta = acos(N1.dot(N1,N2));
    
    
    double e0 = var.at(0);
    double e1 = var.at(1);
    double N = var.at(2);
    double theta0 = var.at(3);
    double e2 = var.at(4);
    double beta0 = var.at(5);
    //std::cout<< e0 <<"  "<< e1 <<"  "<< e2 <<"  "<< N <<"  "<<alpha<<"  "<< theta0 <<" \n ";

    theta0 = theta0/180.0*3.14;
    beta0 = beta0/180.0*3.14;
    
    //====== obtain the orinatation of the beta
    Vec3D P1 (v1->GetVXPos(),v1->GetVYPos(),v1->GetVZPos()) ;
    Vec3D P2 (v2->GetVXPos(),v2->GetVYPos(),v2->GetVZPos()) ;
    double l = (P2-P1).norm();
    Vec3D l1 = N2-N1+(P2-P1)*(1.0/l);
    if(l1.norm()<1)
        beta=-beta;
    
    E = -e0-e1*cos(N*(m_Angle2D-theta0))+e2*(beta-beta0)*(beta-beta0);
    
    return E;
}
double Energy::F10(vertex *v1, vertex *v2, std::vector<double> var)
{
    /// F = -e0+e1*cosN(phi-phi0)+e2*exp(-alpha(theta-theta0))
    double E = 0;
    
    double e0 = var.at(0);
    double e1 = var.at(1);
    double N = var.at(2);
    double phi0 = var.at(3);
    double e2 = var.at(4);
    double alpha = var.at(5);
    double theta0 = var.at(6);
    
    
    //std::cout<< e0 <<"  "<< e1 <<"  "<< e2 <<"  "<< N <<"  "<<alpha<<"  "<< theta0 <<" \n ";

    Vec3D X1(v1->GetVXPos(),v1->GetVYPos(),v1->GetVZPos());
    Vec3D X2(v2->GetVXPos(),v2->GetVYPos(),v2->GetVZPos());
    Vec3D geodesic_dir=(X2-X1);
    Vec3D *pBox=v1->GetBox();
    
    for (int i=0;i<3;i++)
    {
        if(fabs(geodesic_dir(i))>(*pBox)(i)/2)
        {
            if(geodesic_dir(i)<0)
                geodesic_dir(i)=geodesic_dir(i)+(*pBox)(i);
            else if(geodesic_dir(i)>0)
                geodesic_dir(i)=geodesic_dir(i)-(*pBox)(i);
        }
    }
    
    
    Vec3D y1=(v1->GetG2LTransferMatrix())*geodesic_dir;
    y1(2)=0;
    y1=y1*(1/(y1.norm()));
    Vec3D y2=(v2->GetG2LTransferMatrix())*geodesic_dir;
    y2(2)=0;
    y2=y2*(1/(y2.norm()));
    Vec3D n(0,0,1);
    Vec3D d1 = (v1->GetInclusion())->GetLDirection();
    Vec3D d2 = (v2->GetInclusion())->GetLDirection();
    double cos1 = y1.dot(y1,d1);
    double sin1 = n.dot(n*y1,d1);
    double cos2 = y1.dot(y2,d2);
    double sin2 = n.dot(n*y2,d2);
    double S_an = sin1*sin2+cos1*cos2;
    m_Angle2D=acos(S_an);
    
    
    Vec3D gd1 = (v2->GetL2GTransferMatrix())*d1;
    Vec3D gd2 = (v2->GetL2GTransferMatrix())*d2;
    
    m_Angle2D  = acos(n.dot(gd1,gd2));
    theta0 = theta0/180.0*3.14;
    
    E = e0+e1*cos(N*(m_Angle2D-phi0))+e2*exp(-alpha*(m_Angle2D-theta0)*(m_Angle2D-theta0));
    E = -E;
    
    return E;
}
double Energy::F11(vertex *v1, vertex *v2, std::vector<double> var)
{
    /// F = e0+e1*cosN(phi-phi0)+e2*(3*(m1.r)(m2.r2)-m1.m2)
    
    /// m1 = Dir1+DeltaD
    double E = 0;
    
    double e0 = var.at(0);
    double e1 = var.at(1);
    double N = var.at(2);
    double phi0 = var.at(3);
    double e2 = var.at(4);
    double Q0 = var.at(5);
    
    
    /*
      /
     /Q0
     -------->
     */
    
    
    //std::cout<< e0 <<"  "<< e1 <<"  "<< e2 <<"  "<< N <<"  "<<alpha<<"  "<< theta0 <<" \n ";
    
    Vec3D X1(v1->GetVXPos(),v1->GetVYPos(),v1->GetVZPos());
    Vec3D X2(v2->GetVXPos(),v2->GetVYPos(),v2->GetVZPos());
    Vec3D geodesic_dir=(X2-X1);
    Vec3D *pBox=v1->GetBox();
    
    for (int i=0;i<3;i++)
    {
        if(fabs(geodesic_dir(i))>(*pBox)(i)/2)
        {
            if(geodesic_dir(i)<0)
                geodesic_dir(i)=geodesic_dir(i)+(*pBox)(i);
            else if(geodesic_dir(i)>0)
                geodesic_dir(i)=geodesic_dir(i)-(*pBox)(i);
        }
    }
    double GEOLength = geodesic_dir.norm();
    Vec3D GEOUnit = geodesic_dir*(1.0/GEOLength);   /// geo direction
    
    Vec3D y1=(v1->GetG2LTransferMatrix())*geodesic_dir;
    y1(2)=0;
    y1=y1*(1/(y1.norm()));
    Vec3D y2=(v2->GetG2LTransferMatrix())*geodesic_dir;
    y2(2)=0;
    y2=y2*(1/(y2.norm()));
    Vec3D n(0,0,1);
    Vec3D d1 = (v1->GetInclusion())->GetLDirection();
    Vec3D d2 = (v2->GetInclusion())->GetLDirection();
    double cos1 = y1.dot(y1,d1);
    double sin1 = n.dot(n*y1,d1);
    double cos2 = y1.dot(y2,d2);
    double sin2 = n.dot(n*y2,d2);
    double S_an = sin1*sin2+cos1*cos2;
    m_Angle2D=acos(S_an);
    
    Vec3D q0(0,0,tan(Q0));
    Vec3D gd1 = d1+q0;
    Vec3D gd2 = d2+q0;
    gd1 = gd1*(1/(gd1.norm()));
    gd2 = gd2*(1/(gd2.norm()));
    gd1 = (v1->GetL2GTransferMatrix())*gd1;
    gd2 = (v2->GetL2GTransferMatrix())*gd2;
    
    double Q1=gd1.dot(gd1,GEOUnit);
    double Q2=gd2.dot(gd2,GEOUnit);
    

    


    
    E = e0+e1*cos(N*(m_Angle2D-phi0))+e2*(3*Q1*Q2-n.dot(gd1,gd2));
    E = -E;
    
    return E;
}
