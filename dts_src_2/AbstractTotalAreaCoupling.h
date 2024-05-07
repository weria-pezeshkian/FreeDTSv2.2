#if !defined(AFX_AbstractTotalAreaCoupling_H)
#define AFX_AbstractTotalAreaCoupling_H
#include <iostream>
#include "VAHGlobalMeshProperties.h"
// Define a base class with a virtual function
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a base class for global Curvature and how energy should change
========================================================
*/
class  AbstractTotalAreaCoupling : public VAHGlobalMeshProperties {
public:
    AbstractTotalAreaCoupling(VAHGlobalMeshProperties *VHA) : VAHGlobalMeshProperties(*VHA) {
        
    }
    virtual ~ AbstractTotalAreaCoupling(){
        
    }
    virtual  void Initialize(std::vector<triangle *> &pTriangle) = 0;
    virtual  void UpdateArea(double oldarea, double newarea) = 0;
    virtual  double CalculateEnergyChange(double oldarea,  double newarea) = 0;
    virtual inline std::string GetDerivedDefaultReadName() = 0;
    
    inline static std::string GetBaseDefaultReadName() {return "TotalAreaCoupling";}

};
//---- a class for no box change
class NoTotalAreaCoupling : public  AbstractTotalAreaCoupling {
public:
    NoTotalAreaCoupling(VAHGlobalMeshProperties *VHA) : AbstractTotalAreaCoupling(VHA) {
        
    }
    ~NoTotalAreaCoupling(){
        
    }
    virtual inline std::string GetDerivedDefaultReadName()  {return "NoTotalAreaCoupling";}
    
    void Initialize(std::vector<triangle *> &pTriangle){
    return;}
    void UpdateArea(double oldarea, double newarea){
        return;
    }
    double CalculateEnergyChange(double oldarea,  double newarea){
        return 0;
    }
    
};

#endif
