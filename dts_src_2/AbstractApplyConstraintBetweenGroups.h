#if !defined(AFX_AbstractApplyConstraintBetweenGroups_H)
#define AFX_AbstractApplyConstraintBetweenGroups_H
#include <iostream>

// Define a Abstract class with a virtual function and some main function
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a Abstract class for Applying Constraint Between Groups
========================================================
*/
class State;
class  AbstractApplyConstraintBetweenGroups {
public:
    AbstractApplyConstraintBetweenGroups(){
        
    }
    virtual ~ AbstractApplyConstraintBetweenGroups(){
        
    }
    virtual bool Initialize() = 0;
    virtual inline double GetEnergy()   = 0;
    virtual inline double GetDistance() = 0;
    
    virtual inline std::string GetDerivedDefaultReadName() {return "";}
    inline static std::string GetBaseDefaultReadName() {return "ConstraintBetweenGroups";}

    
private:
    
};
//---- a class for no constraint
class NoConstraint : public AbstractApplyConstraintBetweenGroups {
public:
    NoConstraint(){
        
    }
    ~NoConstraint(){
        
    }
    inline std::string GetDerivedDefaultReadName() {return "NoConstraint";}
    bool Initialize(){return true;}
    inline double GetEnergy()   {return 0;}
    inline double GetDistance()  {return 0;}

};
#endif
