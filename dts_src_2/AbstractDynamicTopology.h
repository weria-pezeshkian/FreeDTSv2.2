#if !defined(AFX_DynamicTopology_H)
#define AFX_DynamicTopology_H
#include <iostream>
#include "RNG.h"

// Define a base class with a virtual function
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a base class for changing the box
========================================================
*/
class AbstractDynamicTopology {
public:
    AbstractDynamicTopology();
    virtual ~AbstractDynamicTopology();

    virtual bool MCMove(int step, double* totalenergy, RNG *rng, Voxelization<vertex>* p_Allvoxel) = 0;
    virtual void initialize() = 0;
    
    virtual inline std::string GetDerivedDefaultReadName() {return "";}
    inline static std::string GetBaseDefaultReadName() {return "DynamicTopology";}

};
//---- a class for no box change
class ConstantTopology : public AbstractDynamicTopology {
public:
    ConstantTopology(){
        
    }
    ~ConstantTopology(){
        
    }

    inline std::string GetDerivedDefaultReadName() {return "ConstantTopology";}
    void initialize(){
        return;
    }
    bool MCMove(int step, double * TotalEnergy, RNG *rng, Voxelization<vertex>* p_Allvoxel){
        return false;
    }
};

#endif
