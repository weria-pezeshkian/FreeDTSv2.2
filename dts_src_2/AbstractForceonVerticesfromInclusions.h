#if !defined(AFX_AbstractForceonVerticesfromInclusions_H)
#define AFX_AbstractForceonVerticesfromInclusions_H
#include <iostream>

// Define a Abstract class with a virtual function and some main function
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a Abstract class for Force on Vertices by Inclusions.
========================================================
*/
class State;
class  AbstractForceonVerticesfromInclusions {
public:
    AbstractForceonVerticesfromInclusions(){
        
    }
    virtual ~ AbstractForceonVerticesfromInclusions(){
        
    }
    virtual double Energy_of_Force(vertex *p, Vec3D dx) = 0;
    virtual inline std::string GetDerivedDefaultReadName() = 0;
    
    inline static std::string GetBaseDefaultReadName() {return "InclusionInducedForceOnVertex";}
    
private:
    
};
//---- a class for no box change
class NoForce : public AbstractForceonVerticesfromInclusions {
public:
    NoForce(){
        
    }
    ~NoForce(){
        
    }
    
    virtual inline std::string GetDerivedDefaultReadName()  {return "NoForce";}

    
    double Energy_of_Force(vertex *p, Vec3D dx){
        return 0;
    }


};
#endif
