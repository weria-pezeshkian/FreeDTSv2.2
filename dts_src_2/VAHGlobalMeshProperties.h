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
    ~VAHGlobalMeshProperties();


       inline double GetTotalVolume()             const     {return m_TotalVolume;}
       inline double GetTotalArea()               const          {return m_TotalArea;}
       inline double GetTotalMeanCurvature()      const            {return m_TotalH;}

public:
    

    
    
    void Initialize(std::vector<triangle *> pTriangle);   ///
    double VolumeofTrianglesAroundVertex(vertex * pVeretx);   /// this does not mean anything outside of this code
    double SingleTriangleVolume(triangle * ptriangle);   /// and this one

    //=====
    
    
    
private:
    double m_TotalVolume;
    double m_TotalArea;
    int m_TotalH;
    
    State* m_pState;


};


#endif
