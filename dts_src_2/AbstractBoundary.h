#if !defined(AFX_Boundary_H)
#define AFX_Boundary_H
#include <iostream>

// Define a base class with a virtual function
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a base class for rigid wall boundry
========================================================
*/
class  AbstractBoundary {
public:
    AbstractBoundary(){
        
    }
    virtual ~ AbstractBoundary(){
        
    }
    virtual void Initialize(std::vector <vertex *> &Apv);
    virtual bool MoveHappensWithinTheBoundary(int step, double x, double y, double z, vertex* v) = 0;

    virtual inline std::string GetDerivedDefaultReadName() {return "";}
    inline static std::string GetBaseDefaultReadName() {return "Boundary";}
    

};
//---- a class for no box change
class PBCBoundary : public  AbstractBoundary {
public:
    PBCBoundary(){
        
    }
    ~PBCBoundary(){
        
    }

    inline std::string GetDerivedDefaultReadName() {return "PBC";}
    void Initialize(std::vector <vertex *> &Apv){
     
        return;
    }
    bool MoveHappensWithinTheBoundary(int step, double x, double y, double z, vertex* v){
        return true;
    }
    
};

#endif
