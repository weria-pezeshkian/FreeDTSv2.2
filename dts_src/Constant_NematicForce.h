#if !defined(AFX_Constant_NematicForce_H_334B21B8_INCLUDED_)
#define AFX_Constant_NematicForce_H_334B21B8_INCLUDED_
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 2024
*/
#include "SimDef.h"
#include "vertex.h"

class Constant_NematicForce
{
public:
    Constant_NematicForce();
    ~Constant_NematicForce();



    


public:
    double m_F0;
    double m_Fd;   // in the direction of the protein
    double m_Fn;    // in the direction of surface normal
    double m_En;
public:
    double Energy_of_Force(vertex *p, Vec3D dx);
    void  Initialize();
    //=====
    
    
    
private:


};


#endif
