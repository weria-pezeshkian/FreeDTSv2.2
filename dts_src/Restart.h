#ifndef RESTART_H_INCLUDED
#define RESTART_H_INCLUDED

#include "SimDef.h"
#include "CreateMashBluePrint.h"
/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
This for reading and writing the restart file (a binary file)
 */
class State;
struct MESH;

class Restart {
public:
    Restart();
    Restart(State* pstate);
    ~Restart();
    
    bool UpdateRestartState(int step, double r_vertex, double r_box);  // Writing the restart file
    MeshBluePrint ReadFromRestart(const std::string& filename, int& step, bool& readok, double& r_vertex, double& r_box); // Reading a restart file
    void SetRestartFileName();
    void UpdatePeriod(int period);

private:
    void WriteRestart(std::string &filename, int step, MESH* pmesh, double r, double rb);
    MeshBluePrint ReadRestart(std::string filename, int& step, bool& readok, double& r_vertex, double& r_box); // Reading a restart file

    State* m_pState;
    std::string m_TEMFileName;
    std::string m_RestartFileName;
    int m_Period;
};

#endif // RESTART_H_INCLUDED
