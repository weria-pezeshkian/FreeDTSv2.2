

#ifdef _OPENMP
# include <omp.h>
#endif
#include <thread>
#include <iostream>
#include <string.h>
#include "MC_Simulation.h"
#include "State.h"
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
MC_Simulation::MC_Simulation(State *pState){
    m_pState = pState;
}
MC_Simulation::~MC_Simulation(){
    
}
void MC_Simulation::Initialize(){
    
    return;
}
bool MC_Simulation::do_Simulation(){
    
//--- before doing anything lets have a frame of the initial system
       m_pState->GetVisualization()->WriteAFrame(-1);
    
//---> Voxelize the mesh for the first time; this should be done before any calculation
        m_pState->GetVoxelization()->Voxelize(m_pState->GetMesh()->GetActiveV());
    
//----> checking if the mesh is good, within the bond of the simulation type. For here, it should be within
        //CheckMesh();
    
//---> Update curvature and energy of the system
    //-- update system curvature
    // m_pState->GetCurvatureCalculator()->UpdateSystemCurvature()
    //-- update system energy
    // m_pState->GetEnergyCalculator()->UpdateSystemEnergy()
    

   // open input files (except bts as it has opened before)
   // m_pState->
//--- before simualtion lets have a frame of the initial system
        m_pState->GetVisualization()->WriteAFrame(0);
    time_t startTime;
    time(&startTime);
for(int step = GetInitialStep(); step < GetFinalStep(); step++){
        
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
       //m_pState->GetDynamicBox()->MCMoveBoxChange(); // we may need the final step as well to check if the update of move size should be done
       //--- update the mesh topology
       //m_pState->GetDynamicTopology()->MCMove();
        //--- update edge of mesh open edge
        // m_pState->GetOpenEdgeEvolution()->MC_Move();
        //---- convert inclusions
        //m_pState->GetInclusionConversion()->ActiveExchange();
    
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
    time_t currentTime;
    time(&currentTime);
    if (currentTime - startTime >= 30) {
        std::cout<<"Step = "<<step<<"/"<<GetFinalStep()<<std::flush;
        std::cout << std::fixed << std::setprecision(2);
        std::cout<<"Rates: "<<std::flush;
        std::cout<<"vertex move = "<<m_pState->GetVertexPositionUpdate()->GetAcceptanceRate(true)<<std::flush;
        std::cout << '\r';
        std::cout << "\033[K";
    }

} // for(int step=GetInitialStep(); step<GetFinalStep(); step++)
   
    
//---> check energy to see if there is any energy leak in the calculations
       //if(!energygood(totalenergy)){
            //std::cout<<"  energy is ... "<<std::endl;
         //  return false;
      //  }
        
    return true;
}



