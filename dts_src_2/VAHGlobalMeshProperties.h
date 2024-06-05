#if !defined(AFX_VAHGlobalMeshProperties_H_224B21B8_C13C_2248_BF23_124095086233__INCLUDED_)
#define AFX_VAHGlobalMeshProperties_H_224B21B8_C13C_2248_BF23_124095086233__INCLUDED_
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 this class is created in version 1.2 to centerlized global variables like voume, area and total curvature.
*/
#include "SimDef.h"

class vertex;
class triangle;
class links;
class State;
class VAHGlobalMeshProperties  {
public:
    VAHGlobalMeshProperties();
    virtual ~VAHGlobalMeshProperties();


       inline double GetTotalVolume()                 const             {return m_TotalVolume;}
       inline double GetTotalArea()                   const             {return m_TotalArea;}
       inline double GetTotalMeanCurvature()          const             {return m_TotalCurvature;}
       inline bool   GetCalculateVAH()                const             {return m_CalculatedGlobalVariable;}

public:
    void CalculateAVertexRingContributionToGlobalVariables(vertex *p_vertex, double &vol, double &area, double& curvature);
    void CalculateALinkTrianglesContributionToGlobalVariables(links *p_link, double &vol, double &area, double& curvature);
    void Initialize(State* pState);

protected:
    void UpdateCalculatedGlobalVariable(){
        m_CalculatedGlobalVariable = true;
        return;
    }
    double CalculateSingleTriangleVolume(triangle *pTriangle);

    
protected:
    double m_TotalVolume;
    double m_TotalArea;
    int m_TotalCurvature;    //  Delta A = h*m_TotalCurvature = h* Sum [2H_vA_v]
    bool m_CalculatedGlobalVariable;

private:
    State *m_pState;
};


#endif
