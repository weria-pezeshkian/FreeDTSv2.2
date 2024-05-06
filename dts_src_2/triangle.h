#if !defined(AFX_triangle_H_6Q4B21B8_C13C_5648_BF23_124095086233__INCLUDED_)
#define AFX_triangle_H_6Q4B21B8_C13C_5648_BF23_124095086233__INCLUDED_

#include "SimDef.h"
#include "Vec3D.h"
/*
 * File: triangle.h
 * Author: Weria Pezeshkian (weria.pezeshkian@gmail.com)
 * Description: Definition of the Triangle class, representing a geometric triangle in 3D space.
 *              This class encapsulates functionality related to triangles, such as storing vertex pointers,
 *              calculating area and normal vectors, and updating triangle properties.
 */

class vertex;   // Forward declaration

class triangle {

public:
    triangle(int id);
	triangle(int id, vertex *v1, vertex *v2, vertex *v3);
	 ~triangle();

        inline const int GetTriID()          const  {return m_ID;}
        inline vertex *GetV1()                  	{return m_V1;}
        inline vertex *GetV2()                  	{return m_V2;}
        inline vertex *GetV3()                  	{return m_V3;}
        inline double GetArea()                  	{return m_Area;}
        inline double GetVolume()                   {return m_Volume;}
    	inline Vec3D GetAreaVector()                {return m_AreaVector;}
        inline Vec3D GetNormalVector()              {return m_Normal;}
    
//--- only needed for visualization
        inline bool GetRepresentation()             {return m_Representation;}

public:
void UpdateRepresentation(bool); 	/// this is for visulaization output and does not effect the simulation
void UpdateNormal_Area(Vec3D *Box);     // this might be send to other classes
void UpdateNormal_Area(Vec3D& norm, double& area);     // A function to uopdate normal and area

void UpdateVertex(vertex *v1,vertex *v2,vertex *v3); // If a link flips the triangle changes its vertices, this function do the job
void UpdateVolume(double vol);
void UpdateID(int id); // this should not be used for active trinagles

private:
    vertex *m_V1;
    vertex *m_V2;
    vertex *m_V3;
  int m_ID;
  Vec3D m_Normal;
  Vec3D m_AreaVector;
  double m_Area;
double m_Volume;

    bool m_Representation;


};


#endif
