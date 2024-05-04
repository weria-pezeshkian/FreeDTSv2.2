#if !defined(AFX_AbstractGlobalCurvature_H)
#define AFX_AbstractGlobalCurvature_H
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
class  AbstractGlobalCurvature  : public VAHGlobalMeshProperties {
public:
    AbstractGlobalCurvature(VAHGlobalMeshProperties *VHA, State *pstate) : VAHGlobalMeshProperties(pstate) {
        
    }
    virtual ~ AbstractGlobalCurvature(){
        
    }
    virtual  double GetEnergy() = 0;
    virtual  bool GetState() = 0;
    virtual  void Initialize(std::vector<vertex *> &Ver) = 0;
    virtual  void UpdateEnergyChange(double delta_area, double delta_curvature) = 0;
    virtual  double CalculateEnergyChange(double delta_area, double delta_curvature) = 0;
    
    virtual inline std::string GetDerivedDefaultReadName() = 0;
    inline static std::string GetBaseDefaultReadName() {return "GlobalCurvature";}
    

};
//---- a class for no box change
class NoGlobalCurvature : public  AbstractGlobalCurvature {
public:
    NoGlobalCurvature(VAHGlobalMeshProperties *VHA, State *pstate) : AbstractGlobalCurvature(VHA, pstate) {
        
    }
    ~NoGlobalCurvature(){
        
    }
    
    inline bool GetState()                           {return false;} // if the coupling is active
    inline double GetEnergy()                        {return 0;} // Only part of energy asscoiated with this class
    virtual inline std::string GetDerivedDefaultReadName()  {return "NoGlobalCurvature";}

    
    void Initialize(std::vector<vertex *> &Ver){
     
        return;
    }
    void UpdateEnergyChange(double delta_area, double delta_curvature){
        return;
    }
    double CalculateEnergyChange(double delta_area, double delta_curvature){
        return 0;
    }
    
};

#endif
