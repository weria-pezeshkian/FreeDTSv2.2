/*
Author: Weria Pezeshkian (weria.pezeshkian@gmail.com && weria.pezeshkian@nbi.ku.dk)
Copyright (c) Weria Pezeshkian

Description:
    This class checks the name of the executable, although it's not a critical task for the current version.
*/

#ifdef _OPENMP
    #include <omp.h>
#endif

#include <vector>
#include <string>
#include "SimDef.h"
#include "Job.h"
#include "State.h"
#include "RNG.h"

/*
Description:
    This class handles task distribution and allows for execution.

Parameters:
    argument (std::vector<std::string>&): Vector containing input arguments.

*/
Job::Job(const std::vector<std::string> &argument) {
    // Extract executable name from the argument list
    std::string ex_name = Nfunction::SubstringFromRight(argument[0], '/');
    
    // Check if the executable name matches the expected name
    if (ex_name != EXE_NAME) { // EXE_NAME is defined in the SimDef.h file
        std::cout << "--> unrecognized executable name ---> " << ex_name << " :( " << " it should be " << EXE_NAME << std::endl;
        exit(0);
    }
#ifndef _OPENMP
    // Perform a normal simulation on a single CPU if OpenMP is not enabled
#if DEBUG_MODE == Enabled
    std::cout<<" Job started  \n";
#endif
    State T_state(argument);
#if DEBUG_MODE == Enabled
    std::cout<<" stated created  \n";
#endif
    T_state.Initialize();
#if DEBUG_MODE == Enabled
    std::cout<<" state ini  \n";
#endif
    T_state.GetSimulation()->do_Simulation();
#if DEBUG_MODE == Enabled
    std::cout<<" sim is done  \n";
#endif
#else
//---> constract an State object
    State T_state(argument);
//---> get parallel tempering data in the input.dts file
    Parallel_Tempering parallel_tempering_data = T_state.m_Parallel_Tempering;
    
//---> here is one openmp is on but still want to perform one single simulation
    if (!parallel_tempering_data.State) {
        T_state.Initialize();
        T_state.GetSimulation()->do_Simulation();
    }
//**********************************************//
    else { // run parallel tempering simulations
//  Parallel_Tempering  = on PT_steps  PT_minbeta    PT_maxbeta
   
        //ParallelReplicaSimulation Simulation(.....);
        
        double PT_minbeta = IPTdata.PT_minbeta;          // from inputfile
        double PT_maxbeta = IPTdata.PT_maxbeta;            // from inputfile
        int PT_steps = IPTdata.PT_steps;                    // from inputfile
        std::string gfile = T_state.GetRunTag();              // general output file name
    // std::string tsifile = (T_state.m_TRJTSI).tsiFolder_name;
    //    RNG Random1(T_state.m_Seed);
        int itime = T_state.GetSimulation()->GetFinalStep();
        int etime = T_state.GetSimulation()->GetFinalStep();
        int num_threads = T_state.GetThreads_Number();
        int exchange_step_length = (etime-itime)/PT_steps;
        omp_set_num_threads(num_threads);         // setting the number of thread
        std::vector<double> threads_energy;      // a container to store energy of each thread just (i) is the thread id
        std::vector<int>  tempid_thread_id;   // a map for temp_id to thread_id, temp = temperatures.at(temp_id)
        std::vector<double> betas;        // the 1/temprature of each temp_id
        for (int i=0;i<num_threads;i++){
            
                threads_energy.push_back(0);
                betas.push_back(PT_minbeta + double(i)*(PT_maxbeta - PT_minbeta)/double(num_threads-1));
                tempid_thread_id.push_back(i);

        } // for (int i=0;i<num_threads;i++)
#pragma omp parallel if(Parallel_Tempering)
{
     State ReplicaState(argument);
     int Thread_ID = omp_get_thread_num();
     int Thread_num = omp_get_num_threads();
    if(Thread_num!=num_threads) {
        std::cout<<" error ---> requested thread number is availabe. \n";
        exit(0);
    } //     if(Thread_num!=num_threads) {
    
    if(num_threads>1){
        ReplicaState.GetSimulation()->GetBeta(PT_minbeta + double(Thread_ID)*(PT_maxbeta - PT_minbeta)/double(Thread_num-1));
    }
    else {
        std::cout<<" error--> Parallel_Tempering cannot be done with a single thread. "<<std::endl;
        exit(0);
    }
    
    if(Thread_ID == Thread_num-1){
      //  ReplicaState.m_Targeted_State = true;
    }
    else {
      //  ReplicaState.m_Targeted_State = false;
    }

     //============================
    //========: Runing the simulations
    for (int pt_step = 0; pt_step<exchange_step_length;pt_step++) {
        ReplicaState.GetSimulation()->UpdateInitialStep (itime + pt_step*PT_steps);
        ReplicaState.GetSimulation()->UpdateFinalStep (itime + (pt_step+1)*PT_steps);

        ReplicaState.GetSimulation()->do_Simulation();

#pragma omp critical //(filling) not sure if it is needed
        
        threads_energy[Thread_ID] = ReplicaState.GetEnergyCalculator()->GetEnergy();   // get the latest energy of the systems
        
#pragma omp barrier
//change temprature
#pragma omp single
// set the tempratur eof each state
        for (int c=0;c<betas.size()-1;c++)
        {
            double b1 = betas[c];
            double b2 = betas[c+1];
            int t1 = tempid_thread_id[c];
            int t2 = tempid_thread_id[c+1];
            double e1 = threads_energy[t1];
            double e2 = threads_energy[t2];
            double cra = (e2-e1)*(b2-b1);  // should be checked
            double ran = Random1.UniformRNG(1.0); //should be checked
            if(exp(cra)>ran)
            {
                // P = min(1,exp([E_i-Ej]*(1/Ti-1/Tj))) = min(1,exp([E_i-E_j]*(1/Ti-1/Tj)))
                tempid_thread_id[c] = t2;
                tempid_thread_id[c+1] = t1;
            }
        }
#pragma omp critical //(filling) not sure if it is needed
        std::vector<int>::iterator it= find(tempid_thread_id.begin(), tempid_thread_id.end(), Thread_ID);
        int temid_of_thread = std::distance(tempid_thread_id.begin(), it);
        if(temid_of_thread<betas.size()){
            
            ReplicaState.GetSimulation()->SetBeta(betas[temid_of_thread]);
            if(temid_of_thread == num_threads-1){
             //   ReplicaState.m_Targeted_State = true;
            }
            else {
             //   ReplicaState.m_Targeted_State = false;
            }
        }
        else {
             std::cout<<"error-->3627 should not happen\n";
            exit(0);
         }
                        
                    }
#pragma omp barrier
        // S.m_GeneralOutputFilename = gfile+ +"_temp_" + f.Int_to_String(S.m_Beta);
        // (S.m_TRJTSI).tsiFolder_name =  tsifile+"_temp_" +f.Int_to_String(S.m_Beta);
                }
                
            }
} // end pragma omp parallel
}  // end of Parallel_Tempering

#endif
#if DEBUG_MODE == Enabled
std::cout<<" End of Job class \n";
#endif

}
Job::~Job() {
    
}




