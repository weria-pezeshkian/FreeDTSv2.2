#if !defined(AFX_ConstantExternalField_H_334B21B8_INCLUDED_)
#define AFX_ConstantExternalField_H_334B21B8_INCLUDED_
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 2024
*/
#include "vertex.h"
#include "AbstractExternalFieldOnVectorFields.h"

class ConstantExternalField : public AbstractExternalFieldOnVectorFields {

public:
    ConstantExternalField(double k, double x, double y, double z);
    ~ConstantExternalField();

    double GetCouplingEnergy(vertex *pvertex);
private:
    Vec3D m_FieldDirection;
    double m_FieldStrength;

};


#endif
