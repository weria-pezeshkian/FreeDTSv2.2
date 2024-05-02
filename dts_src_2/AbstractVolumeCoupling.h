#if !defined(AFX_AbstractVolumeCoupling_H)
#define AFX_AbstractVolumeCoupling_H
#include <iostream>

// Define a base class with a virtual function
/*
=======================================================
 developed 2024 by Weria
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This class is a base class for changing the box
========================================================
*/
class AbstractVolumeCoupling {
public:
    AbstractVolumeCoupling();
    virtual ~AbstractVolumeCoupling();

    virtual inline bool GetState()= 0;
    virtual inline double GetTotalVolume() = 0;
    virtual inline double GetTotalArea() = 0;
    virtual void Initialize(std::vector<triangle *> pTriangle) = 0;
    virtual double SingleTriangleVolume(triangle * ptriangle) = 0;
    virtual double GetEnergyChange(int s, double oa, double ov, double na, double nv) = 0;
    virtual void UpdateArea_Volume(double oa, double ov, double na, double nv) = 0;
    virtual inline std::string GetDerivedDefaultReadName() = 0;
    
    inline static std::string GetBaseDefaultReadName() {return "VolumeCoupling";}

};
//---- a class for no box change
class NoCoupling : public AbstractVolumeCoupling {
public:
    NoCoupling(){
        
    }
    ~NoCoupling(){
        
    }
    inline bool GetState()                   {return false;}
    inline double GetTotalVolume()                  {return 0;}
    inline double GetTotalArea()                  {return 0;}
    virtual inline std::string GetDerivedDefaultReadName() = {return "NoVolumeCoupling";}

    
    void    Initialize(std::vector<triangle *> pTriangle)       {return;}
    double  VolumeofTrianglesAroundVertex(vertex * pVeretx)     {return 0;}
    double  GetEnergyChange(int s, double oa, double ov, double na, double nv) {return 0;}
    double  Energy(double volume, double area, double a)        {return 0;}
    void    UpdateArea_Volume(double oa, double ov, double na, double nv) {return ;}
    double  SingleTriangleVolume(triangle * ptriangle)          {return 0;}



};

#endif
