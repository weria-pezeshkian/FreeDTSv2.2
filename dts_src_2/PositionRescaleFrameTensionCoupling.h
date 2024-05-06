#if !defined(AFX_PositionRescaleFrameTensionCoupling_H_DD4B21B8_C13C_5648_BF23_344095086239__INCLUDED_)
#define AFX_PositionRescaleFrameTensionCoupling_H_DD4B21B8_C13C_5648_BF23_344095086239__INCLUDED_


#include "SimDef.h"
#include "links.h"
#include "Vec3D.h"
#include "Energy.h"
#include "HarmonicPotentialBetweenTwoGroups.h"
#include "triangle.h"
#include "MESH.h"
#include "AbstractDynamicBox.h"

class State;
class PositionRescaleFrameTensionCoupling : public AbstractDynamicBox {
public:
    
    PositionRescaleFrameTensionCoupling(int period, double , State *pState);
	~PositionRescaleFrameTensionCoupling();

    void Initialize();
    bool ChangeBoxSize(int step);

    inline  std::string GetDerivedDefaultReadName()  {return "IsotropicFrameTension";}
    inline  static std::string GetDefaultReadName()  {return "IsotropicFrameTension";}

private:
    
    bool VertexMoveIsFine(double dx,double dy, double mindist2, double maxdist2);



    
private:
    double m_SigmaP;
    int m_Period;
    Vec3D *m_pBox;
    State *m_pState;
    double *m_pLmin2;
    double *m_pLmax2;
    double *m_pminAngle;

};


#endif
