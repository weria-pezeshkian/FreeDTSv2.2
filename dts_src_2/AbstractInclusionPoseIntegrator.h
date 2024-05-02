#if !defined(AFX_AbstractInclusionPoseIntegrator_H)
#define AFX_AbstractInclusionPoseIntegrator_H
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
class  AbstractInclusionPoseIntegrator {
public:
    AbstractInclusionPoseIntegrator(){
        m_UpdateRate_Angle = 1;
        m_UpdateRate_Kawasaki = 1;
    }
    virtual ~ AbstractInclusionPoseIntegrator(){
        
    }
    virtual bool Initialize(State *pState) = 0;
    virtual bool EvolveOneStep(int step) = 0;
    virtual inline std::string GetDerivedDefaultReadName() = 0;
    inline static std::string GetBaseDefaultReadName() {return "InclusionPoseIntegrator";}
    
    inline double GetRate_Kawasaki()                       {return m_UpdateRate_Kawasaki;}
    inline double GetRate_Angle()                          {return m_UpdateRate_Angle;}

    
    void SetMoveRate(double rate_angle, double rate_kawa){
        m_UpdateRate_Angle = rate_angle;
        m_UpdateRate_Kawasaki = rate_kawa;

        return;
    }
    
private:
    double m_UpdateRate_Angle;
    double m_UpdateRate_Kawasaki;

};

#endif
