#if !defined(AFX_HarmonicPotentialBetweenTwoGroups_H_334B21B8_D13C_2248_QF23_124095086255__INCLUDED_)
#define AFX_HarmonicPotentialBetweenTwoGroups_H_334B21B8_D13C_2248_QF23_124095086255__INCLUDED_
#include "SimDef.h"
#include "MESH.h"
#include "AbstractApplyConstraintBetweenGroups.h"
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 This object is to couple the system to a harmonic potential between two groups, and also change the ...
 */
class State;
class HarmonicPotentialBetweenTwoGroups : public AbstractApplyConstraintBetweenGroups, public MESH {
public:
    HarmonicPotentialBetweenTwoGroups(State* pState, double K, double l0, std::string group1,std::string group2,double nx,double ny,double nz);
    ~HarmonicPotentialBetweenTwoGroups();

    inline double GetEnergy()                           {return m_Energy;}
    inline double GetDistance()                         {return m_Dist;}

public:
    
    bool Initialize();
    void CalculateEnergy(int step);
    void MovingVertex(vertex* v, Vec3D Dx);
    void RejectMovingVertex(vertex* v, Vec3D Dx);
    
    inline  std::string GetDerivedDefaultReadName() {return "HarmonicPotentialBetweenTwoGroups";}
    inline static std::string GetDefaultReadName() {return "HarmonicPotentialBetweenTwoGroups";}
    std::string CurrentState();
    //=====
private:  
    Vec3D COMVertexGroup(std::vector<vertex *>);
    
    
private:
    double m_K;
    std::string m_Group1Name;
    std::string m_Group2Name;
    std::vector<vertex *> m_pGroup1;
    std::vector<vertex *> m_pGroup2;
    Vec3D m_Direction;
    double m_L0;
    Vec3D m_Group1COG;
    Vec3D m_Group2COG;
    double m_Energy;


    double m_Dist;

    





};


#endif
