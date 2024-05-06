

/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 */
#include <iostream>
#include "TimeSeriesLogInformation.h"
#include "State.h"
#include "Nfunction.h"
#include "SimDef.h"

TimeSeriesLogInformation::TimeSeriesLogInformation(){

}
TimeSeriesLogInformation::TimeSeriesLogInformation(State *pState){
    m_pState = pState;
}
TimeSeriesLogInformation::~TimeSeriesLogInformation(){
    if (m_TimeSeriesFile.is_open()) {
        m_TimeSeriesFile.flush();
        m_TimeSeriesFile.close();
     }
}
bool TimeSeriesLogInformation::WriteAStringIntoFile(std::string &str){
    
    m_TimeSeriesFile<<str;
    m_TimeSeriesFile<<std::endl;
    
    return true;
}
bool TimeSeriesLogInformation::OpenFile(bool clearfile){
    
    std::string filename = m_pState->GetRunTag();
    filename = filename + ".log";
    
    //-- if it is a restart simulation, then do not clear the file and check for the last line in the file
    if(!clearfile){
        m_TimeSeriesFile.open(filename.c_str(),std::fstream::app);
    }
    else{
        m_TimeSeriesFile.open(filename.c_str());
    }
    
    if (!m_TimeSeriesFile.is_open()) {
        std::cerr << "---> error: Unable to open file: " << filename << std::endl;
        return false;
    }
    
    WriteStateInfo();

    return true;
}
bool TimeSeriesLogInformation::FlushLogFile(){ // the energy file should be flushed first

    m_TimeSeriesFile.flush(); //

    return true;
}
void TimeSeriesLogInformation::WriteStateInfo(){
    
    std::vector<std::string> argument = m_pState->GetCommandLineArgument();
    for (std::vector<std::string>::iterator it = argument.begin() ; it != argument.end(); ++it){
        m_TimeSeriesFile<<(*it)<<"   ";
    }
    m_TimeSeriesFile<<std::endl;
    m_TimeSeriesFile<<";--------- this part can be used as an input.dts file ---------------------------------  "<<std::endl;
    /* m_TimeSeriesFile<<"Integrator = "<<m_Integrator<<std::endl;
     m_TimeSeriesFile<<"MC_Moves = "<<(m_MCMove.VertexMove)<<"  "<<(m_MCMove.LinkFlip)<<"  "<<m_MCMove.EdgeVertexMove<<"  "<<(m_MCMove.InclusionMove_Angle)<<"  "<<(m_MCMove.InclusionMove_Kawasaki)<<std::endl;
     m_TimeSeriesFile<<"Initial_Step = "<<m_Initial_Step<<std::endl;
     m_TimeSeriesFile<<"Final_Step = "<<m_Final_Step<<std::endl;
     m_TimeSeriesFile<<"MinfaceAngle = "<<m_MinFaceAngle<<std::endl;
     m_TimeSeriesFile<<"OutPutEnergy_periodic = "<<m_OutPutEnergy_periodic<<std::endl;
     m_TimeSeriesFile<<"Restart_periodic = "<<m_RESTART.restartPeriod<<std::endl;
     m_TimeSeriesFile<<"TopologyFile = "<<m_TopologyFile<<std::endl;
     m_TimeSeriesFile<<"Seed = "<<m_Seed<<std::endl;
     m_TimeSeriesFile<<"Kappa = "<<m_inc_ForceField.m_BendingRigidity<<"  "<<m_inc_ForceField.m_GaussianRigidity<<std::endl;
     m_TimeSeriesFile<<"Display_periodic = "<<m_Display_periodic<<std::endl;
     m_TimeSeriesFile<<"CNTCELL = "<<m_CNTCELL(0)<<" "<<m_CNTCELL(1)<<" "<<m_CNTCELL(2)<<" "<<std::endl;
     m_TimeSeriesFile<<"GeneralOutputFilename = "<<m_GeneralOutputFilename<<std::endl;
     m_TimeSeriesFile<<"Min_Max_LinkLenghtsSquare = "<<m_MinVerticesDistanceSquare<<"  "<<m_MaxLinkLengthSquare<<std::endl;
     m_TimeSeriesFile<<"OutPutTRJ_TSI = "<<(m_TRJTSI.tsiPeriod)<<"  "<<(m_TRJTSI.tsiPrecision)<<"  "<<(m_TRJTSI.tsiFolder_name)<<"  "<<std::endl;
     m_TimeSeriesFile<<"OutPutTRJ_BTS = "<<(m_TRJBTS.btsPeriod)<<"  "<<(m_TRJBTS.btsPrecision)<<"  "<<(m_TRJBTS.btsFile_name)<<"  "<<std::endl;
    std::string state = "off";
    */
    m_TimeSeriesFile<<";------------------------------------------  "<<std::endl;

}
