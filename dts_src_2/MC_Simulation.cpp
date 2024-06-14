

#ifdef _OPENMP
# include <omp.h>
#endif
#include <thread>
#include <iostream>
#include <string.h>
#include "MC_Simulation.h"
#include "State.h"
#include "SimDef.h"
/*
 List of skipped function due to lack of clarity based on current state of the code
 They need to be finished before calling this a new version.
 
 1) void TimeSeriesLogInformation::WriteStateInfo(){
 
 */






/*
 Weria Pezeshkian (weria.pezeshkian@gmail.com)
 Copyright (c) Weria Pezeshkian
 MC simulation class, runs mc simulation if it is defined in the input file.
 */
MC_Simulation::MC_Simulation(State *pState) : m_pState(pState) {

}
MC_Simulation::~MC_Simulation(){
    
}
void MC_Simulation::Initialize(){
    
    return;
}
bool MC_Simulation::do_Simulation(){
#if DEBUG_MODE == Enabled
    std::cout<<" do_Simulation function is starting  \n";
#endif
    
//---> Voxelize the mesh for the first time; this should be done before any calculation
    m_pState->GetVoxelization()->Voxelize(m_pState->GetMesh()->GetActiveV());
    
#if DEBUG_MODE == Enabled
    std::cout<<" system has been voxelaized  \n";
#endif

//----> checking if the mesh is good, within the bond of the simulation type. For here, it should be within
        //CheckMesh();
 
//--- before simualtion lets have a frame of the initial system
        m_pState->GetVisualization()->WriteAFrame(0);
   // time_t startTime;
   // time(&startTime);
#if DEBUG_MODE == Enabled
    std::cout<<" We have reached simulation run loop!  \n";
#endif
for(int step = GetInitialStep(); step <= GetFinalStep(); step++){
        
//---> centering the simulation box
    if(GetBoxCentering()!=0 && step%GetBoxCentering()==0){
        m_pState->GetMesh()->CenterMesh();
        m_pState->GetVoxelization()->ReassignMembersToVoxels(m_pState->GetMesh()->GetActiveV());
    } // [ if(GetBoxCentering()!=0 && step%GetBoxCentering()==0)]

//---> Run standard Integrators
        //--- run the vertex position update
        m_pState->GetVertexPositionUpdate()->EvolveOneStep(step); // we may need the final step as well to check if the update of move size should be done
        //--- run the link flip update
        m_pState->GetAlexanderMove()->EvolveOneStep(step);
        //--- run the inclusion update
        m_pState->GetInclusionPoseUpdate()->EvolveOneStep(step);

//----> Run the supplementary integrators
       //--- update the box side
         m_pState->GetDynamicBox()->ChangeBoxSize(step); // we may need the final step as well to check if the update of move size should be done
       //--- update the mesh topology
       //m_pState->GetDynamicTopology()->MCMove();
        //--- update edge of mesh open edge
         m_pState->GetOpenEdgeEvolution()->Move(step);
        //---- convert inclusions
        m_pState->GetInclusionConversion()->Exchange(step);
    
//----> write files
    //--- write visulaization frame
    m_pState->GetVisualization()->WriteAFrame(step);
    //--- write non-binary trejectory e.g., tsi, tsg
    m_pState->GetNonbinaryTrajectory()->WriteAFrame(step,m_pState->GetRunTag());
    //--- write binary trejectory e.g., bts
    m_pState->GetBinaryTrajectory()->WriteAFrame(step);
    //--- write into time seri file, e.g., energy, volume ...
    m_pState->GetTimeSeriesDataOutput()->WriteTimeSeriesDataOutput(step);
    //--- write check point for the state
    m_pState->GetRestart()->UpdateRestartState(step, m_pState->GetVertexPositionUpdate()->GetDR(), m_pState->GetDynamicBox()->GetDR());

//----> print info about the simulation, e.g., rate,
   // time_t currentTime;
   // time(&currentTime);
    if (step%50 == 0) {
        std::cout<<"Step = "<<step<<"/"<<GetFinalStep()<<std::flush;
        std::cout << std::fixed << std::setprecision(3);
        std::cout<<" Rates: "<<std::flush;
        std::cout<<" vertex move = "<<m_pState->GetVertexPositionUpdate()->GetAcceptanceRate(true)<<std::flush;
        std::cout<<"; alexander move = "<<m_pState->GetAlexanderMove()->GetAcceptanceRate(true)<<std::flush;
        std::cout<<"; inclusion move = "<<m_pState->GetInclusionPoseUpdate()->GetAcceptanceRate(true)<<std::flush;
        if(m_pState->GetDynamicBox()->GetDerivedDefaultReadName() != "No")
        std::cout<<"; Box Move = "<<m_pState->GetDynamicBox()->GetAcceptanceRate(true)<<std::flush;
        std::cout << '\r';
        std::cout << "\033[K";
    }

} // for(int step=GetInitialStep(); step<GetFinalStep(); step++)
   

    /*m_pState->GetCurvatureCalculator()->print();
    m_pState->GetCurvatureCalculator()->print();*/
    m_pState->GetCurvatureCalculator()->Initialize();
    double Final_energy = m_pState->GetEnergyCalculator()->CalculateAllLocalEnergy();
    double energy_leak = Final_energy - m_pState->GetEnergyCalculator()->GetEnergy();
    std::cout << std::fixed << std::setprecision(4);
    if(fabs(energy_leak) > 0.0001){
        
        std::cout<<"---> possible source of code error: energy leak... "<<energy_leak<<" with real energy of "<<Final_energy<<"  and stored energy of "<<m_pState->GetEnergyCalculator()->GetEnergy()<<"\n";
    }
        
    return true;
}
std::string MC_Simulation::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ GetDerivedDefaultReadName();
    state = state + "\n Min_Max_Lenghts = "+Nfunction::D2S(m_MinLength2)+" "+Nfunction::D2S(m_MaxLength2);
    state = state + "\n MinfaceAngle = "+Nfunction::D2S(m_MinAngle);
    state = state + "\n Temprature = "+Nfunction::D2S(m_Beta)+" "+Nfunction::D2S(m_DBeta);
    state = state + "\n Box_Centering_F = "+Nfunction::D2S(m_CenteringFrequently);
    state = state + "\n Set_Steps = "+Nfunction::D2S(m_Initial_Step)+" "+Nfunction::D2S(m_Final_Step);
    
    return state;
}



