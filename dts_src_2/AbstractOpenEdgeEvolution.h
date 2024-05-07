#if !defined(AFX_AbstractOpenEdgeEvolution_H)
#define AFX_AbstractOpenEdgeEvolution_H
#include <iostream>
#include "RNG.h"

// Define a base class with a virtual function for different open edge treatment algorthems
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a base class for changing the box
========================================================
*/
class AbstractOpenEdgeEvolution {
public:
    AbstractOpenEdgeEvolution(){
        
    }
    virtual ~AbstractOpenEdgeEvolution(){
        
    }
    virtual bool Move(int step) = 0;
    virtual void Initialize() = 0;
    
    virtual inline std::string GetDerivedDefaultReadName() = 0;
    inline static std::string GetBaseDefaultReadName() {return "OpenEdgeEvolution";}
    
protected:
    double m_NumberOfAttemptedMoves;
    double m_AcceptedMoves;

};
//---- a class for no edge change
class NoEvolution : public AbstractOpenEdgeEvolution {
public:
    NoEvolution(){
        
    }
    ~NoEvolution(){
        
    }

    void Initialize(){
        return;
    }
    bool Move(int step){
        return false;
    }
    inline std::string GetDerivedDefaultReadName()  {return "NoEvolution";}


};

#endif
