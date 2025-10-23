

#include <stdio.h>
#include "EquilibriumInclusionExchangeByChemicalPotential.h"
#include "Nfunction.h"
#include "State.h"
 /*
  * @file EquilibriumInclusionExchangeByChemicalPotential.cpp
  * @brief Implementation of the EquilibriumInclusionExchangeByChemicalPotential class methods.
  *
  * This file contains the implementation of the methods declared in the EquilibriumInclusionExchangeByChemicalPotential class.
  * The class is responsible for exchanging inclusion types between two states based on active algorthem.
  * The exchange is not done based on the energetic of the states. It is active, but energy get updated after the exchange
  * It initializes the inclusion exchange process, performs exchanges at defined intervals, and manages energy calculations.
  *
  * @author Weria Pezeshkian (weria.pezeshkian@gmail.com)
  */

EquilibriumInclusionExchangeByChemicalPotential::EquilibriumInclusionExchangeByChemicalPotential(State *pstate, int period, double rate, double chemicalpotential, std::string t_name1, std::string t_name2) :
                        m_Period(period),
                        m_Rate(rate),
                        m_Mu(chemicalpotential),
                        m_TypeName_1(t_name1),
                        m_TypeName_2(t_name2),
                        m_Beta(pstate->GetSimulation()->GetBeta()),
                        m_DBeta(pstate->GetSimulation()->GetDBeta()),
                        m_N2(0),
                        m_N1(0)
{
    m_NumberOfMoves = 0;
}
EquilibriumInclusionExchangeByChemicalPotential::~EquilibriumInclusionExchangeByChemicalPotential() {
    
}
void EquilibriumInclusionExchangeByChemicalPotential::Initialize(State *pstate) {
    
    const std::vector<inclusion *>& pAllInclusion = m_pState->GetMesh()->GetInclusion();
    for (std::vector<inclusion *>::const_iterator it = pAllInclusion.begin() ; it != pAllInclusion.end(); ++it){
        if((*it)->m_IncType->ITName == m_TypeName_1){
            m_pSubInc.push_back(*it);
            m_pIncType1 = (*it)->m_IncType;
            m_N1++;
        }
        else if((*it)->m_IncType->ITName == m_TypeName_2){
            m_pSubInc.push_back(*it);
            m_pIncType2 = (*it)->m_IncType;
            m_N2++;
        }
    } // end for (std::vector<inclusion *>::const_iterator it = pAllInclusion.begin() ; it != pAllInclusion.end(); ++it)
    
    m_N = m_N1 + m_N2;
    
    
    m_NumberOfMoves = int (m_Rate * double(m_N));
}
bool EquilibriumInclusionExchangeByChemicalPotential::Exchange(int step){
    
    
    if( step%m_Period !=0 )
        return false;
    
    // in the old version it was as: double eta = 2*double(m_N-m_DeltaN0)/double(m_N+m_DeltaN0)-1;
    

    
    for (int i = 0; i< m_NumberOfMoves;i++) {
      
        int r_inc_id = m_pState->GetRandomNumberGenerator()->IntRNG(m_N);
        inclusion *p_inc = m_pSubInc[r_inc_id];
        double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);

        if(TryForOneInclusion(p_inc,thermal)){


        }
      }
    
    return true;
}
bool EquilibriumInclusionExchangeByChemicalPotential::TryForOneInclusion(inclusion * p_inc, double thermal){
    
    
    double old_energy = 0.0;
    double new_energy = 0.0;

    vertex * pver = p_inc->Getvertex();      // Current vertex of the inclusion
    double en = pver->GetEnergy();

    InclusionType* inctype = p_inc->m_IncType;

    if(p_inc->m_IncType == m_pIncType1){
        p_inc->m_IncType  = m_pIncType2;
    }
    else{
        p_inc->m_IncType  = m_pIncType1;
    }
 
    /*std::vector<links*> Affected_links = GetEdgesWithInteractionChange(d_links);
    for (std::vector<links *>::iterator it = Affected_links.begin() ; it != Affected_links.end(); ++it){
        (*it)->Copy_InteractionEnergy();
        old_energy += 2*(*it)->GetIntEnergy();
    }*/
    
    
    
    //--- now, lets update energy
  /*  new_energy += m_pState->GetEnergyCalculator()->SingleVertexEnergy(tver); //
    for (std::vector<links *>::iterator it = n_edges.begin() ; it != n_edges.end(); ++it){
            new_energy += m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(*it);
    }
    */
    
    //pver->UpdateInclusion((tver->GetInclusion()));


    double diff_energy = new_energy - old_energy;
    double U = m_Beta * diff_energy - m_DBeta;
    //---> accept or reject the move
    if(U <= 0 || exp(-U) > thermal ) {
        // move is accepted
       // (m_pState->GetEnergyCalculator())->AddToTotalEnergy(diff_energy);
        return true;
    }
    else{
        
        if(p_inc->m_IncType == m_pIncType1){
            p_inc->m_IncType  = m_pIncType2;
        }
        else{
            p_inc->m_IncType  = m_pIncType1;
        }
        
    }


/*for (std::vector<inclusion *>::iterator it = m_pSubInc.begin() ; it != m_pSubInc.end(); ++it){
        
        
        
    }
    
        double old_energy = 0;
        double new_energy = 0;
        
        bool exchange_happend = false;
        int delta_N = m_N1 - m_N2;
        double phi = double(delta_N-m_Delta_N0);

        vertex* tver = (*it)->Getvertex();
        InclusionType* inctype = (*it)->m_IncType;
        double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);
        
        if((*it)->m_IncType == m_pIncType1){ // transition from 1->2
            double Prob = m_Epsilon / static_cast<double>(m_N) / (1.0 + exp(-m_Gama * phi));

            if( Prob > thermal ) {
                (*it)->m_IncType = m_pIncType2;
                m_N1--;
                m_N2++;
                exchange_happend = true;
            }
        } // if((*it)->m_IncType->ITName == m_TypeName_1)
        else { // transition from 2->1
            double Prob = m_Epsilon / static_cast<double>(m_N) / (eta + exp(m_Gama * phi));

            if(Prob > thermal) {
                (*it)->m_IncType  = m_pIncType1;
                m_N1++;
                m_N2--;
                exchange_happend = true;
            }
        }
        if(exchange_happend){ // note, while the exchange has happend, energies are not updated. so old energy can be obtained.
            
            old_energy += tver->GetEnergy();
            std::vector<links *> n_edges = tver->GetVLinkList();     // links that their interaction energies may change
            
            for (std::vector<links *>::iterator it = n_edges.begin() ; it != n_edges.end(); ++it){
                old_energy += 2*((*it)->GetIntEnergy());
            }
            //--- now, lets update energy
            new_energy += m_pState->GetEnergyCalculator()->SingleVertexEnergy(tver); //
            for (std::vector<links *>::iterator it = n_edges.begin() ; it != n_edges.end(); ++it){
                    new_energy += m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(*it);
            }
            //--- update the total elastic energy of the system
            double T_en = m_pState->GetEnergyCalculator()->GetEnergy();
            T_en += new_energy - old_energy;
            m_pState->GetEnergyCalculator()->UpdateTotalEnergy(T_en);
            m_ActiveEnergy += (new_energy - old_energy);
            
            m_NumberOfAttemptedMoves++;
            m_AcceptedMoves++;
            
        } // if(exchange_happend)
        else { // otherwise
            m_NumberOfAttemptedMoves++;
        }

        
    } //for (m_pSubInc.begin())
    */
    
    return true;
}
std::string EquilibriumInclusionExchangeByChemicalPotential::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    state = state + " "+ Nfunction::D2S(m_Period);
    state = state + " "+ Nfunction::D2S(m_Rate);
    state = state + " "+ Nfunction::D2S(m_Mu);
    state = state + " "+ m_TypeName_1;
    state = state + " "+ m_TypeName_2;

    return state;
}
