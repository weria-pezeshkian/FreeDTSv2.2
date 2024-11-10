#if !defined(AFX_UserDefinedForceonVertices_H_334B21B8_INCLUDED_)
#define AFX_UserDefinedForceonVertices_H_334B21B8_INCLUDED_
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 2024
*/
#include "SimDef.h"
#include "vertex.h"
#include "links.h"
#include "AbstractForceonVertices.h"
class UserDefinedForceonVertices : public AbstractForceonVertices{
public:
    UserDefinedForceonVertices(std::string inputs);
    ~UserDefinedForceonVertices();
    double Energy_of_Force(vertex *p, Vec3D dx);
    std::string CurrentState();
    inline  std::string GetDerivedDefaultReadName()  {return "User";}
    inline static std::string GetDefaultReadName() {return "User";}
    
private:
    Vec3D GetForce(vertex *pv);
    std::string m_Inputs;

};


#endif
