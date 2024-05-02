#if !defined(AFX_DynamicBox_H)
#define AFX_DynamicBox_H
#include <iostream>
#include "Voxelization.h"
// Define a base class with a virtual function
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a base class for changing the box
========================================================
*/
class AbstractDynamicBox {
public:
    AbstractDynamicBox(){
        
    }
    virtual ~AbstractDynamicBox(){
        
    }

    virtual bool MCMoveBoxChange(double dx, double* totalenergy, double temp, int step, Voxelization<vertex>* p_Allvoxel) = 0;
    virtual void initialize() = 0;
    virtual  bool GetCNTCondition() = 0;
    virtual int GetTau() = 0;
    virtual inline  std::string GetDerivedDefaultReadName()  {return "";}

    
    inline static std::string GetBaseDefaultReadName()  {return "Dynamic_Box";}
    
    inline double GetDR()                               {return m_DR;}
    void UpdateDR(double dr){
        m_DR = dr;
        return;
    }


private:
    double m_DR;
    
};
//---- a class for no box change
class NoBoxChange : public AbstractDynamicBox {
public:
    NoBoxChange(){
        
    }
    ~NoBoxChange(){
        
    }

    inline  std::string GetDerivedDefaultReadName()  {return "ConstantBox";}

    bool GetCNTCondition(){
        return false;
    }
    int GetTau(){
        return 0;
    }
    void initialize(){
        return;
    }
    bool MCMoveBoxChange(double dx, double * TotalEnergy, double temp, int step, Voxelization<vertex>* p_Allvoxel){
        return false;
    }

};

#endif
