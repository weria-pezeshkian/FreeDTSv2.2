#ifndef INCLUSION_H
#define INCLUSION_H

#include "SimDef.h"
#include "Vec3D.h"
#include "InclusionType.h"

/*
 * @brief inclusion class.
 *
 * This class models an inclusion object.
 *
 * Author: Weria Pezeshkian (weria.pezeshkian@gmail.com)
 * Copyright (c) Weria Pezeshkian
 */


class vertex;
class inclusion {

public:
    
    inclusion(int id, InclusionType *inctype);
	 ~inclusion();

    inline const int GetID()                               const  {return m_ID;}
    inline vertex* Getvertex()                                    {return m_pvertex;}
    inline Vec3D GetLDirection()                                  {return m_LDirection;}
    inline Vec3D GetGDirection()                                  {return m_GDirection;}
    inline InclusionType *GetInclusionType()                      {return m_IncType;}
    InclusionType *m_IncType;


  void Updatevertex(vertex * );
  void UpdateLocalDirection(const Vec3D & lo_dir);
  void UpdateGlobalDirection(const Vec3D & lg_dir);
  bool UpdateGlobalDirectionFromLocal();
  bool UpdateLocalDirectionFromGlobal();


private:
    int m_ID;
    Vec3D m_LDirection;      /// its direction in the local frame
    Vec3D m_GDirection;       /// its direction in the global frame
    vertex *m_pvertex;

};


#endif
