

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 */
#include <fstream>
#include "Restart.h"
#include "State.h"
#include "CreateMashBluePrint.h"

Restart::Restart(){
    m_Period = 1000;
    m_RestartFileName = "";
}
Restart::Restart(State *pState)
{
    m_pState = pState;
    m_TEMFileName = "-1.res";
    m_Period = 1000;
    m_RestartFileName = ".res";
}
Restart::~Restart()
{
    
}
void Restart::UpdatePeriod(int period){
    
    m_Period = period;
    return;
}
void Restart::SetRestartFileName(){
    
    m_RestartFileName = m_pState->GetRunTag() +"."+ RestartExt;
    m_TEMFileName = m_pState->GetRunTag() +"-1."+ RestartExt;
    return;
}
bool Restart::UpdateRestartState(int step,  double r_vertex, double r_box){
        
    if(m_Period!=0 && step%m_Period == 0 ){
        WriteRestart(m_RestartFileName, step, m_pState->GetMesh(),r_vertex, r_box);
    }
    return true;
}
//=== Writing a restart file: this file is just the active [State] Object with an updated initial step
void Restart::WriteRestart(std::string &filename, int step, MESH * pmesh, double r, double rb){
    
    MeshBluePrint blueprint = pmesh->Convert_Mesh_2_BluePrint(pmesh);

    //=== we write the restart into tem restart file
    std::fstream Rfile;
    Rfile.open(m_TEMFileName.c_str(),std::ios::out | std::ios::binary);
    (Rfile).write((char *) &step, sizeof(int));    
    (Rfile).write((char *) &r, sizeof(double));
    (Rfile).write((char *) &rb, sizeof(double));
    Vec3D box = blueprint.simbox;
    (Rfile).write((char *) &(box(0)), sizeof(double));
    (Rfile).write((char *) &(box(1)), sizeof(double));
    (Rfile).write((char *) &(box(2)), sizeof(double));

    int size = (blueprint.bvertex).size();
    (Rfile).write((char *) &size, sizeof(int));
    for (std::vector<Vertex_Map>::iterator it = (blueprint.bvertex).begin() ; it != (blueprint.bvertex).end(); ++it)
        (Rfile).write((char *) &(*it), sizeof(Vertex_Map));
    size = (blueprint.btriangle).size();
    (Rfile).write((char *) &size, sizeof(int));
    for (std::vector<Triangle_Map>::iterator it = (blueprint.btriangle).begin() ; it != (blueprint.btriangle).end(); ++it)
        (Rfile).write((char *) &(*it), sizeof(Triangle_Map));
    size = (blueprint.binclusion).size();
    (Rfile).write((char *) &size, sizeof(int));
    for (std::vector<Inclusion_Map>::iterator it = (blueprint.binclusion).begin() ; it != (blueprint.binclusion).end(); ++it)
        (Rfile).write((char *) &(*it), sizeof(Inclusion_Map));

    Rfile.close();
    Rfile.flush();

    Nfunction::CopyBinaryFile(m_TEMFileName,m_RestartFileName, FreeDTS_BUFFERSIZE);

    remove(m_TEMFileName.c_str());

    return;
}
MeshBluePrint Restart::ReadFromRestart(const std::string& inputFilename, int& step, bool& readOk, double& rVertex, double& rBox) {
    readOk = false;
    std::string restartFilename = m_TEMFileName; // Default to temporary file name

    // Check if the temporary restart file exists
    std::ifstream tempFile(m_TEMFileName);
    if (!tempFile.good()) {
        // Use the provided input filename if the temporary file doesn't exist
        restartFilename = inputFilename;

        // Append the restart file extension if needed
        std::string ext = inputFilename.substr(inputFilename.find_last_of(".") + 1);
        if (ext != RestartExt) {
            restartFilename = restartFilename + "." + RestartExt;
        }
    }

    // Read from the determined restart file
    return ReadRestart(restartFilename, step, readOk, rVertex, rBox);
}
//=== Read a restart file and load to the  active [State] Object
MeshBluePrint Restart::ReadRestart(std::string filename , int &step, bool &readok, double &r_vertex, double &r_box) {

    MeshBluePrint blueprint;

    // We should make a check and see if the restart file is writtn correctly or not
    //=== just read the restart file and copy it into the current active [State] object
    std::fstream Rfile;
    Rfile.open(filename.c_str(), std::ios::in |std::ios::binary);
if(Rfile.is_open())
{
    Rfile.read((char *) &step, sizeof(int));

    double lx,ly,lz;
    Rfile.read((char *) &r_vertex, sizeof(double));
    Rfile.read((char *) &r_box, sizeof(double));
    Rfile.read((char *) &lx, sizeof(double));
    Rfile.read((char *) &ly, sizeof(double));
    Rfile.read((char *) &lz, sizeof(double));
    Vec3D box(lx,ly,lz);
    blueprint.simbox= box;
    std::vector<Vertex_Map> bvertex;       // a vector of all vertices (only the blueprint not the object) in the mesh
    std::vector<Triangle_Map> btriangle;   // a vector of all triangles (only the blueprint not the object) in the mesh
    std::vector<Inclusion_Map> binclusion; // a vector of all inclusions (only the blueprint not the object) in the mesh
    int size;
    Rfile.read((char *) &size, sizeof(int));
    for (int i=0;i<size;i++){
        Vertex_Map vmap;
        (Rfile).read((char *) &vmap, sizeof(Vertex_Map));
        bvertex.push_back(vmap);
    }
    Rfile.read((char *) &size, sizeof(int));
    for (int i=0;i<size;i++){
        Triangle_Map tmap;
        (Rfile).read((char *) &tmap, sizeof(Triangle_Map));
        btriangle.push_back(tmap);
    }
    Rfile.read((char *) &size, sizeof(int));
    for (int i=0;i<size;i++){
        Inclusion_Map incmap;
        (Rfile).read((char *) &incmap, sizeof(Inclusion_Map));
        binclusion.push_back(incmap);
    }
    blueprint.bvertex = bvertex;
    blueprint.btriangle = btriangle;
    blueprint.binclusion = binclusion;

    Rfile.close();

}
    readok = true;
    return blueprint;
}

