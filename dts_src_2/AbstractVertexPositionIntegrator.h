#if !defined(AFX_AbstractVertexPositionIntegrator_H)
#define AFX_AbstractVertexPositionIntegrator_H

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
class  AbstractVertexPositionIntegrator {
public:
    AbstractVertexPositionIntegrator(){
        m_FreezGroupName = "";
        m_DR = 0.05;
        m_UpdateRate_Surf = 1;
        m_UpdateRate_Edge = 1;
    }
    virtual ~ AbstractVertexPositionIntegrator(){
        
    }
    virtual bool Initialize(State *pState) = 0;
    virtual bool EvolveOneStep(int step) = 0;

    inline std::string GetFreezGroupName()              {return m_FreezGroupName;}
    inline double GetDR()                               {return m_DR;}
    inline double GetRate_Surf()                       {return m_UpdateRate_Surf;}
    inline double GetRate_Edge()                       {return m_UpdateRate_Edge;}
    virtual inline std::string GetDerivedDefaultReadName() = 0;
    
    inline static std::string GetBaseDefaultReadName() {return "VertexPositionIntegrator";}
    
    //---- update
    void UpdateFreezGroupName(std::string name_freezgroup){
        m_FreezGroupName = name_freezgroup;
        return;
    }
    void UpdateDR(double dr){
        m_DR = dr;
        return;
    }
    void SetMoveRate(double update_rateSur, double update_rateEdge ){
        m_UpdateRate_Surf = update_rateSur;
        m_UpdateRate_Edge = update_rateEdge;
        return;
    }
    
private:
    std::string m_FreezGroupName;
    double m_DR;
    double m_UpdateRate_Surf;
    double m_UpdateRate_Edge;

};

#endif
