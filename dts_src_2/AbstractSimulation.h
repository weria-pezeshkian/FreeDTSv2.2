#if !defined(AFX_AbstractSimulation_H)
#define AFX_AbstractSimulation_H
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
class  AbstractSimulation {
public:
    AbstractSimulation() : m_CenteringFrequently(0), m_Initial_Step(1), m_Final_Step(10), m_Beta(1.0), m_DBeta (0.0),
                           m_MinLength2(1),  m_MaxLength2(3), m_MinAngle(-5){
    }
    
    /*AbstractSimulation(AbstractSimulation *pSim) : AbstractSimulation(pSim), m_CenteringFrequently(0), m_Initial_Step(1), m_Final_Step(10), m_Beta(1.0), m_DBeta (0.0) {

    }*/
    virtual ~ AbstractSimulation(){
        
    }
    virtual void Initialize() = 0;
    virtual bool do_Simulation() = 0;
    
    inline int GetInitialStep()                         {return m_Initial_Step;}
    
    virtual inline std::string GetDerivedDefaultReadName() = 0;
    inline static std::string GetBaseDefaultReadName() {return "Integrator_Type";}
    
protected:
    inline int GetBoxCentering()                        {return m_CenteringFrequently;}
    inline int GetFinalStep()                           {return m_Final_Step;}
    inline int GetBeta()                                {return m_Beta;}


public:
    void SetCentering(int centering){
        m_CenteringFrequently = centering;
        return;
    }
    void UpdateFinalStep(int final_step){
        m_Final_Step = final_step;
        return;
    }
    void UpdateInitialStep(int ini_step){
        m_Initial_Step = ini_step;
        return;
    }
    void SetBeta(double beta, double dbeta){
        m_Beta = beta;
        m_DBeta = dbeta;
        return;
    }
private:
    int m_CenteringFrequently; // how often centering the system in the box
    int m_Initial_Step;
    int m_Final_Step;

protected:
    double m_Beta;
    double m_DBeta;
    double m_MinLength2;
    double m_MaxLength2;
    double m_MinAngle;
    
};
#endif
