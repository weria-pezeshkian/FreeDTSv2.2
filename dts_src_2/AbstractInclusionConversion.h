#if !defined(AFX_AbstractInclusionConversion_H)
#define AFX_AbstractInclusionConversion_H
#include <iostream>

// Define a Abstract class with a virtual function and some main function
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a Abstract class for curvature calculations.
========================================================
*/
class State;
class  AbstractInclusionConversion {
public:
    AbstractInclusionConversion(){
        
    }
    virtual ~ AbstractInclusionConversion(){
        
    }
    virtual inline bool GetState()= 0;
    virtual void Initialize(State *pstate) = 0;
    virtual inline std::string GetDerivedDefaultReadName() = 0;
    
    inline static std::string GetBaseDefaultReadName() {return "InclusionConversion";}

    
private:
    
};
//---- a class for no box change
class NoInclusionConversion : public AbstractInclusionConversion {
public:
    NoInclusionConversion(){
        
    }
    ~NoInclusionConversion(){
        
    }
    virtual inline std::string GetDerivedDefaultReadName() {return "NoInclusionConversion";}
    inline bool GetState(){
        return false;
    }

    
    void Initialize(State *pstate){
        return;
    }

};
#endif
