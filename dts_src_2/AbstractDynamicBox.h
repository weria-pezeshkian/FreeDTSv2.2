#if !defined(AFX_DynamicBox_H)
#define AFX_DynamicBox_H
#include <iostream>
// Define a base class with a virtual function
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a base class for changing the box
========================================================
*/
class AbstractDynamicBox {
public:
    AbstractDynamicBox(){
        m_DR = 0.01;
    }
    virtual ~AbstractDynamicBox(){
        
    }

    virtual bool ChangeBoxSize(int step) = 0;
    virtual void Initialize() = 0;
    virtual std::string CurrentState() = 0;

    
    virtual inline  std::string GetDerivedDefaultReadName()=0;
    inline static std::string GetBaseDefaultReadName()  {return "Dynamic_Box";}
    
    
    
    inline double GetDR()                               {return m_DR;}
    void UpdateDR(double dr){
        m_DR = dr;
        return;
    }

protected:
    double m_NumberOfAttemptedMoves;
    double m_AcceptedMoves;
    double m_DR;
    
};
//---- a class for no box change
class NoBoxChange : public AbstractDynamicBox {
public:
    NoBoxChange(){
        
    }
    ~NoBoxChange(){
        
    }

    inline  std::string GetDerivedDefaultReadName()  {return "ConstantBox";}

    void Initialize(){
        return;
    }
    bool ChangeBoxSize(int step){
        return false;
    }
    std::string CurrentState(){
        
        std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
        return state;
    }
};

#endif
