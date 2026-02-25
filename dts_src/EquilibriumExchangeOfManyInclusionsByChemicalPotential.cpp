

#include <stdio.h>
#include "EquilibriumExchangeOfManyInclusionsByChemicalPotential.h"
#include "FactoryInclusionConversionMethod.h"
#include "Nfunction.h"
#include "State.h"
/*
 * @file EquilibriumExchangeOfManyInclusionsByChemicalPotential.cpp
 * @brief Implements the EquilibriumExchangeOfManyInclusionsByChemicalPotential class methods.
 *
 * This file provides the definitions of the methods declared in the
 * EquilibriumExchangeOfManyInclusionsByChemicalPotential class. The class handles
 * the exchange of inclusion types between two states based on their chemical
 * potentials. Exchanges are performed according to the energetics of the states,
 * ensuring an equilibrium process.
 *
 * Key responsibilities include initializing the inclusion exchange process,
 * performing exchanges at specified intervals, and managing associated
 * energy calculations.
 *
 * @author Weria Pezeshkian (weria.pezeshkian@gmail.com)
 */

EquilibriumExchangeOfManyInclusionsByChemicalPotential::EquilibriumExchangeOfManyInclusionsByChemicalPotential(std::string info) :                   
                        m_Input_Data(info),
                        m_N2(0),
                        m_N1(0)
{
    
    std::vector<std::string> data = Nfunction::Split(m_Input_Data);
    
    // Prevent out-of-bounds access
    if (data.size() < 3) {
        std::cout << "---> error: insufficient input data for "
                  << this->GetDefaultReadName() << ".\n";
        exit(-1);
    }

    // Use correct parsing function (double instead of int)
    //m_EP = Nfunction::String_to_Double(data[0]);
   // m_R0 = Nfunction::String_to_Double(data[1]);
    
   // m_NumberOfMoves = 0;
    
    /*m_Rate(rate),
    m_Mu(chemicalpotential),
    m_TypeName_1(t_name1),
    m_TypeName_2(t_name2),*/

}
EquilibriumExchangeOfManyInclusionsByChemicalPotential::~EquilibriumExchangeOfManyInclusionsByChemicalPotential() {
    
}
void EquilibriumExchangeOfManyInclusionsByChemicalPotential::Initialize() {
    
    
    m_Beta = &(m_pState->GetSimulation()->GetBeta());
    m_DBeta = &(m_pState->GetSimulation()->GetDBeta());
    
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
bool EquilibriumExchangeOfManyInclusionsByChemicalPotential::Exchange(int step){
    
    if( step%m_Period !=0 )
        return false;
        
    for (int i = 0; i< m_NumberOfMoves;i++) {
      
        int r_inc_id = m_pState->GetRandomNumberGenerator()->IntRNG(m_N);
        inclusion *p_inc = m_pSubInc[r_inc_id];
        double thermal = m_pState->GetRandomNumberGenerator()->UniformRNG(1.0);

        if(TryForOneInclusion(p_inc,thermal)){


        }
      }
    
    return true;
}
bool EquilibriumExchangeOfManyInclusionsByChemicalPotential::TryForOneInclusion(inclusion * p_inc, double thermal){
    
    
    double old_energy = 0.0;
    double new_energy = 0.0;
    int n1 = 0;
    int n2 = 0;
    InclusionType* Type1;  // temp to keep type1
    InclusionType* Type2;
    double chem = m_Mu ;
    vertex * pver = p_inc->Getvertex();      // Current vertex of the inclusion
    double en = pver->GetEnergy();
    old_energy = en;
    std::vector<links *> Affected_links = pver->GetVLinkList();// links that their interaction energies may change

    for (std::vector<links *>::iterator it = Affected_links.begin() ; it != Affected_links.end(); ++it){
        (*it)->Copy_InteractionEnergy();
        old_energy += 2*(*it)->GetIntEnergy();
    }
    
    /// make the move
    InclusionType* inctype = p_inc->m_IncType;
    if(p_inc->m_IncType == m_pIncType1){
        p_inc->m_IncType  = m_pIncType2;
        n1 = m_N1;
        n2 = m_N2;
        Type1 = m_pIncType1;
        Type2 = m_pIncType2;
        chem = -m_Mu;   // because m_Mu is chemical potentail of m_N1
    }
    else{
        p_inc->m_IncType  = m_pIncType1;
        n1 = m_N2;
        n2 = m_N1;
        Type1 = m_pIncType2;
        Type2 = m_pIncType1;
    }
    
    //--- now, lets update energy
    
    
    new_energy = m_pState->GetEnergyCalculator()->SingleVertexEnergy(pver); //
    for (std::vector<links *>::iterator it = Affected_links.begin() ; it != Affected_links.end(); ++it){
            new_energy += m_pState->GetEnergyCalculator()->TwoInclusionsInteractionEnergy(*it);
    }
    
    double diff_energy = new_energy - old_energy - chem ;
    double pre_factor = n1/(m_N - n1 + 1);
    double U = *m_Beta * diff_energy - *m_DBeta;
    
    //---> accept or reject the move
    if(pre_factor * exp(-U) > thermal ) {
        
        // move is accepted
        (m_pState->GetEnergyCalculator())->AddToTotalEnergy(diff_energy);
        if(p_inc->m_IncType == m_pIncType2){
            m_N1--;
            m_N2++;
        }
        else{
            m_N1++;
            m_N2--;
        }
        return true;
    }
    else {
            p_inc->m_IncType  = Type1;
            pver->UpdateEnergy(en);
        for (std::vector<links *>::iterator it = Affected_links.begin() ; it != Affected_links.end(); ++it){
            (*it)->Reverse_InteractionEnergy();
        }
        return false;

    }

    
    return true;
}
std::string EquilibriumExchangeOfManyInclusionsByChemicalPotential::CurrentState(){
    
    std::string state = GetBaseDefaultReadName() +" = "+ this->GetDerivedDefaultReadName();
    state = state + " "+ m_Input_Data;

    return state;
}
static class EquilibriumExchangeOfManyInclusionsByChemicalPotentialRegister {

    // Static create function for the factory
    static AbstractInclusionConversion* Create(std::istream& input)
    {
        std::string info;
        getline(input,info);

        // Construct and return a new object
        return new EquilibriumExchangeOfManyInclusionsByChemicalPotential(info);
    }

public:
    // Constructor automatically registers the class with the factory
    EquilibriumExchangeOfManyInclusionsByChemicalPotentialRegister() {
        FactoryInclusionConversionMethod::Instance().Register(
            "EquilibriumExchangeOfManyInclusions", // string identifier for the factory
            Create                      // the creator function
        );
    }

} EquilibriumExchangeOfManyInclusionsByChemicalPotentialRegisterObject; // static instance triggers registration
