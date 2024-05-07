#if !defined(AFX_CouplingTotalAreaToHarmonicPotential_H_334B21B8_C13C_2248_BF23_124095086233__INCLUDED_)
#define AFX_CouplingTotalAreaToHarmonicPotential_H_334B21B8_C13C_2248_BF23_124095086233__INCLUDED_
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 this class try to force the membrane area to a targeted size
*/
#include "SimDef.h"
#include "vertex.h"
#include "triangle.h"
#include "links.h"
#include "AbstractTotalAreaCoupling.h"
class CouplingTotalAreaToHarmonicPotential : public  AbstractTotalAreaCoupling {
public:
    CouplingTotalAreaToHarmonicPotential(VAHGlobalMeshProperties *VHA, State *pstate, double gamma,  double K0);
    ~CouplingTotalAreaToHarmonicPotential();



    


public:
    void Initialize(std::vector<triangle *> &pTriangle);   ///
    double CalculateEnergyChange(double oldarea,  double newarea);
    void UpdateArea(double oldarea, double newarea);

    inline  std::string GetDerivedDefaultReadName()  {return "HarmonicPotential";}
    inline static std::string GetDefaultReadName() {return "HarmonicPotential";}
    
private:
    double AreaofTrianglesAroundVertex(vertex * pVeretx);   ///
    double m_TotalArea;
    int m_NoEQStep;
    double m_K0;
    double m_A0;
    double m_NT;  // number of the triangles in the system
    double m_Gamma;
    double m_6SQPI;   /// 1/6pi^1/2

};


#endif
