#if !defined(AFX_AbstractInclusionPoseIntegrator_H)
#define AFX_AbstractInclusionPoseIntegrator_H
#include <iostream>
#include "SimDef.h"

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
        m_NumberOfMovePerStep_Angle = 1;
        m_NumberOfMovePerStep_Kawasaki = 1;
        m_DR = 0.1;
    }
    virtual ~ AbstractInclusionPoseIntegrator(){
        
    }
    virtual bool Initialize(State *pState) = 0;
    virtual bool EvolveOneStep(int step) = 0;
    virtual std::string CurrentState() = 0;
    virtual inline std::string GetDerivedDefaultReadName() = 0;
    inline static std::string GetBaseDefaultReadName() {return "InclusionPoseIntegrator";}
    


    
    void SetMoveRate(double rate_angle, double rate_kawa){
        m_NumberOfMovePerStep_Angle = rate_angle;
        m_NumberOfMovePerStep_Kawasaki = rate_kawa;

        return;
    }
    
protected:
    double m_NumberOfAttemptedMoves;
    double m_AcceptedMoves;
    double m_NumberOfMovePerStep_Angle;   // how many updates should be made per step
    double m_NumberOfMovePerStep_Kawasaki;   // how many updates should be made per step
    double m_DR;

};

#endif
