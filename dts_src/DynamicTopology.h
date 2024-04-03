#if !defined(AFX_DynamicTopology_H)
#define AFX_DynamicTopology_H
#include <iostream>
#include "GenerateCNTCells.h"

// Define a base class with a virtual function
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a base class for changing the box
========================================================
*/
class DynamicTopology {
public:
    DynamicTopology();
    virtual ~DynamicTopology();

    virtual bool MCMove(double* totalenergy, double temp, GenerateCNTCells* pGenCNT) = 0;
    virtual void initialize() = 0;
    

};
//---- a class for no box change
class ConstantTopology : public DynamicTopology {
public:
    ConstantTopology();
    ~ConstantTopology();

    
    void initialize();
    bool MCMove(double * TotalEnergy, double temp, GenerateCNTCells *pGenCNT );
};

#endif
