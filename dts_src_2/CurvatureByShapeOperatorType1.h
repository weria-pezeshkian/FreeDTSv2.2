#if !defined(AFX_Curvature_H_524B21B8_C13C_2248_BF23_124095086233__INCLUDED_)
#define AFX_Curvature_H_524B21B8_C13C_2248_BF23_124095086233__INCLUDED_

#include "SimDef.h"
#include "vertex.h"
#include "triangle.h"
#include "links.h"
#include "AbstractCurvature.h"

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 An class to obtain curvature of a single vertex
 This class only give the correct answer if the area of each triangle has been calculated correclty.
 */
class CurvatureByShapeOperatorType1 : public AbstractCurvature{
public:
    
    CurvatureByShapeOperatorType1();
	 ~CurvatureByShapeOperatorType1();

public:
    
    bool SurfVertexCurvature(vertex *p);
    bool EdgeVertexCurvature(vertex *p);
    bool VertexCurvature(vertex *p);

private:
    Tensor2 Householder(Vec3D N);
    Vec3D Calculate_Vertex_Normal(vertex *p);
};


#endif
