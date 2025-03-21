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
    virtual Vec3D Inclusion_Force(vertex *p) = 0;
    virtual inline std::string GetDerivedDefaultReadName() = 0;
    virtual std::string CurrentState() = 0;
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
    
    inline std::string GetDerivedDefaultReadName()  {return "No";}
    inline static std::string GetDefaultReadName()  {return "No";}

    
    double Energy_of_Force(vertex *p, Vec3D dx){
        return 0;
    }
    Vec3D Inclusion_Force(vertex *p){
        Vec3D f(0,0,0);
        return f;
    }
    
    std::string CurrentState(){
        
        std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
        return state;
    }

};
#endif
