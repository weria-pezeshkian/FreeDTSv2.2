#if !defined(AFX_VAHGlobalMeshProperties_H_224B21B8_C13C_2248_BF23_124095086233__INCLUDED_)
#define AFX_VAHGlobalMeshProperties_H_224B21B8_C13C_2248_BF23_124095086233__INCLUDED_
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 this class is created in version 1.2 to couple the system to a force
for taraget VAHGlobalMeshProperties.
*/
#include "SimDef.h"


class State;
class vertex;
class triangle;
class VAHGlobalMeshProperties  {
public:
    VAHGlobalMeshProperties(State* pstate);
    virtual ~VAHGlobalMeshProperties();


       inline double GetTotalVolume()             const             {return m_TotalVolume;}
       inline double GetTotalArea()               const             {return m_TotalArea;}
       inline double GetTotalMeanCurvature()      const             {return m_TotalCurvature;}

public:
    

    
    
    void   Initialize();   //
    ///
protected:
    double GetRingVolumeOfVertex(vertex * pVeretx);   /// this does not mean anything outside of this code
    void CalculateRingVolumeOfVertex(vertex * pVeretx, double &vol, double &area);   //
    void CalculateRingCurvatureOfVertex(vertex * pVeretx, double &curve, double &area);   //

    ///
private:
    double CalculateSingleTriangleVolume(triangle * ptriangle);   // and this one
    double CalculateSingleTriangleArea(triangle * ptriangle);   //
    double CalculateSingleVertexArea(triangle * vertex);   //

    //=====
    
protected:
    State* m_pState;
    double m_TotalVolume;
    double m_TotalArea;
    int m_TotalCurvature;    //  Delta A = h*m_TotalCurvature = h* Sum [2H_vA_v]
    


};


#endif
