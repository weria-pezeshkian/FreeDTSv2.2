#if !defined(AFX_OpenEdgeEvolution_H)
#define AFX_OpenEdgeEvolution_H
#include <iostream>
#include "GenerateCNTCells.h"
#include "RNG.h"

// Define a base class with a virtual function
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a base class for changing the box
========================================================
*/
class OpenEdgeEvolution {
public:
    OpenEdgeEvolution();
    virtual ~OpenEdgeEvolution();

    virtual void MC_Move(RNG* rng, double lmin, double lmax, double maxangle) = 0;
    virtual void Initialize() = 0;
    int m_Rate;


};
//---- a class for no box change
class NoEvolution : public OpenEdgeEvolution {
public:
    NoEvolution();
    ~NoEvolution();

    void Initialize();
    void MC_Move(RNG* rng, double lmin, double lmax, double maxangle);


private:

    int m_Rate;
};

#endif
