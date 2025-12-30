

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 */
#include <fstream>
#include "IO_Files.h"
#include "Vec3D.h"

IO_Files::IO_Files()
{
    Nfunction f;
    m_tsiPrecision = f.Int_to_String(18)+"."+f.Int_to_String(10);
}
IO_Files::~IO_Files()
{
    
}

void IO_Files::WriteQFile(const std::string &filename ,  MeshBluePrint &blueprint)
{
    int pres=10;
    std::ofstream output;
    output.open(filename.c_str());
    
    
    output<<std::fixed;
    output<<std::setprecision( pres )<<(blueprint.simbox)(0)<<"   "<<(blueprint.simbox)(1)<<"   "<<(blueprint.simbox)(2)<<"   \n";
    output<<(blueprint.bvertex).size()<<"\n";
    
    
    for (std::vector<Vertex_Map>::iterator it = (blueprint.bvertex).begin() ; it != (blueprint.bvertex).end(); ++it)
        output<<std::setprecision( pres )<<it->id<<"  "<<it->x<<"  "<<it->y<<"  "<<it->z<<"  "<<it->domain<<std::endl;

    output<< (blueprint.btriangle).size()<<"\n";
    for (std::vector<Triangle_Map>::iterator it = (blueprint.btriangle).begin() ; it != (blueprint.btriangle).end(); ++it)
        output<<it->id<<"  "<<it->v1<<"   "<<it->v2<<"  "<<it->v3<<"  0  "<<std::endl;
    
    output.close();
}
void IO_Files::WriteTSI(const std::string &filename ,  MeshBluePrint &blueprint)
{

    FILE * output;
    output = fopen(filename.c_str(), "w");
    std::string format = "%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf\n";
    const char* version="version 1.1";
    fprintf(output,"%s\n",version);
    //------
    const char* box="box";
    fprintf(output,"%s%18.10lf%18.10lf%18.10lf\n",box,(blueprint.simbox)(0),(blueprint.simbox)(1),(blueprint.simbox)(2));
    
    const char* ver="vertex";
    int size=(blueprint.bvertex).size();
    fprintf(output,"%s%20d\n",ver,size);
    format = "%5d%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf\n";
    for (std::vector<Vertex_Map>::iterator it = (blueprint.bvertex).begin() ; it != (blueprint.bvertex).end(); ++it)
        fprintf(output,format.c_str(),it->id,it->x,it->y,it->z);
  
    const char* tri="triangle";
    size = (blueprint.btriangle).size();
    fprintf(output,"%s%20d\n",tri,size);
   for (std::vector<Triangle_Map>::iterator it = (blueprint.btriangle).begin() ; it != (blueprint.btriangle).end(); ++it)
        fprintf(output,"%10d%10d%10d%10d\n",it->id,it->v1,it->v2,it->v3);
    
 
    const char* inc="inclusion";
    size = (blueprint.binclusion).size();
    fprintf(output,"%s%20d\n",inc,size);
     format = "%10d%10d%10d%"+m_tsiPrecision+"lf%"+m_tsiPrecision+"lf\n";
     for (std::vector<Inclusion_Map>::iterator it = (blueprint.binclusion).begin() ; it != (blueprint.binclusion).end(); ++it)
        fprintf(output,format.c_str(),it->id,it->tid,it->vid,it->x,it->y);

    fclose(output);
}

