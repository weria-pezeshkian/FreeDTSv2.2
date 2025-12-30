#if !defined(AFX_IOFiles_H_INCLUDED_)
#define AFX_IOFiles_H_INCLUDED_
#include "Def.h"
#include "Nfunction.h"
#include "mesh_maps.h"

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 This class makes initial TS files for DTS simulations. The output could be both q or tsi file format.
 currently 3 types of morphologies are IO_Files.hd.
 */

class IO_Files
{
public:
    
    IO_Files();
    ~IO_Files();

public:
     void WriteQFile(const std::string &filename ,  MeshBluePrint& blueprint);
     void WriteTSI(const std::string& filename,  MeshBluePrint& blueprint);

private:
    std::string m_tsiPrecision;
    


};

#endif
