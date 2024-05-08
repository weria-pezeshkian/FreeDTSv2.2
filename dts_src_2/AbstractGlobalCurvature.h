#ifndef ABSTRACT_GLOBAL_CURVATURE_H
#define ABSTRACT_GLOBAL_CURVATURE_H

#include <iostream>
#include "VAHGlobalMeshProperties.h"

/*
 * AbstractGlobalCurvature: A base class for global curvature and energy changes.
 * Developed in 2024 by Weria Pezeshkian (weria.pezeshkian@gmail.com)
 * Copyright (c) Weria Pezeshkian
 */

class  AbstractGlobalCurvature  : public VAHGlobalMeshProperties {
public:
    AbstractGlobalCurvature(VAHGlobalMeshProperties *VHA, State *pstate) : VAHGlobalMeshProperties(pstate), m_Energy(0) {
        
    }
    virtual ~ AbstractGlobalCurvature(){
        
    }
    inline double GetEnergy()  {return m_Energy;}
    
    virtual  void Initialize() = 0;
    virtual  void UpdateEnergyChange(double delta_area, double delta_curvature) = 0;
    virtual  double CalculateEnergyChange(double delta_area, double delta_curvature) = 0;
    virtual std::string CurrentState() = 0;
    virtual inline std::string GetDerivedDefaultReadName() = 0;
    inline static std::string GetBaseDefaultReadName() {return "GlobalCurvature";}

protected:
    double m_Energy;

};
//---- a class for no box change
class NoGlobalCurvature : public  AbstractGlobalCurvature {
public:
    NoGlobalCurvature(VAHGlobalMeshProperties *VHA, State *pstate) : AbstractGlobalCurvature(VHA, pstate) {}
    ~NoGlobalCurvature(){ }
    inline std::string GetDerivedDefaultReadName()  {return "NoGlobalCurvature";}
    void Initialize(){return;}
    void UpdateEnergyChange(double delta_area, double delta_curvature){return;}
    double CalculateEnergyChange(double delta_area, double delta_curvature){return 0;}
    std::string CurrentState(){
        
        std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
        return state;
    }
};

#endif
