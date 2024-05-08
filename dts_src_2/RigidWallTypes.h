#if !defined(AFX_FreeDTSRigidWallType_H)
#define AFX_FreeDTSRigidWallType_H
#include <iostream>
#include "AbstractBoundary.h"
// Define a base class with a virtual function
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a base class for rigid wall boundry
========================================================
*/
//---- a class for no box change
class State;
class TwoFlatParallelWall : public  AbstractBoundary {
public:
    TwoFlatParallelWall(State* pState, double thickness, char direction);
    ~TwoFlatParallelWall();
    inline std::string GetDerivedDefaultReadName() {return "TwoFlatParallelWall";}
    inline static std::string GetDefaultReadName() {return "TwoFlatParallelWall";}

    void Initialize();
    bool MoveHappensWithinTheBoundary(double x, double y, double z, vertex* v);
    std::string CurrentState();

private:
    State* m_pState;
    double m_HalfThickness;
    double m_MidPlane;
    int m_Element;  // 0 X, 1, Y, 2 Z

    
};

#endif
