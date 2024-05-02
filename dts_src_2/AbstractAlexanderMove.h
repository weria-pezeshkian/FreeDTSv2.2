#if !defined(AFX_AbstractAlexanderMove_H)
#define AFX_AbstractAlexanderMove_H
#include <iostream>

// Define a base class with a virtual function and some main function
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a base class for curvature calculations.
========================================================
*/
class State;
class  AbstractAlexanderMove {
public:
    AbstractAlexanderMove(){
        m_Rate = 1;
    }
    virtual ~ AbstractAlexanderMove(){
        
    }
    virtual bool Initialize(State *pState) = 0;
    virtual bool EvolveOneStep(int step) = 0;
    virtual inline std::string GetDerivedDefaultReadName() {return "";}
    
    inline static std::string GetBaseDefaultReadName() {return "AlexanderMove";}
    inline double GetRate()                       {return m_Rate;}

    void SetMoveRate(double rate){
        m_Rate = rate;
        return;
    }
    
private:
    double m_Rate;

};

#endif
