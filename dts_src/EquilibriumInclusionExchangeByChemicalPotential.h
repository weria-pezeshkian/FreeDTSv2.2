#ifndef EQUILIBRIUM_INCLUSION_EXCHANGE_BY_CHEMICAL_POTENTIAL_H
#define EQUILIBRIUM_INCLUSION_EXCHANGE_BY_CHEMICAL_POTENTIAL_H

#include <string>
#include <vector>

#include "SimDef.h"
#include "AbstractInclusionConversion.h"
#include "inclusion.h"
class State;

/**
 * @class EquilibriumInclusionExchangeByChemicalPotential
 * @brief Performs equilibrium inclusion conversion based on chemical potential.
 *
 * This class implements an exchange algorithm that adjusts inclusion types
 * according to a target chemical potential between two states.
 *
 * @author
 * Weria Pezeshkian
 */

class EquilibriumInclusionExchangeByChemicalPotential : public AbstractInclusionConversion
{
public:
    /**
     * @brief Constructor
     * @param period   Exchange cycle period (simulation steps)
     * @param rate       rate, how many
     * @param mu          Chemical potential difference
     * @param type1    Name of the first inclusion type
     * @param type2    Name of the second inclusion type
     */
    EquilibriumInclusionExchangeByChemicalPotential(State *pstate, int period,
                                                    double rate,
                                                    double mu,
                                                    std::string type1,
                                                    std::string type2);

    ~EquilibriumInclusionExchangeByChemicalPotential() override;

    void Initialize(State* pState) override;
    bool Exchange(int step) override;

    std::string CurrentState() override;

    static inline std::string GetDefaultReadName() { return "EquilibriumExchange"; }
    inline std::string GetDerivedDefaultReadName() override {
        return GetDefaultReadName();
    }

    
    bool TryForOneInclusion(inclusion * pinc, double thermal);

private:
    int m_Period;                   ///< Exchange cycle period
    double m_Mu;                    ///< Chemical potential
    double m_Rate;                    ///< Chemical potential

    int m_N1;                       ///< Number of type1 inclusions
    int m_N2;                       ///< Number of type2 inclusions
    int m_N;                        ///< Total number of inclusions

    int m_NumberOfMoves;
    std::string m_TypeName_1;        ///< Type name 1
    std::string m_TypeName_2;        ///< Type name 2

    std::vector<inclusion*> m_pSubInc; ///< Inclusions available for exchange

    InclusionType* m_pIncType1;     ///< Pointer to type 1 definition
    InclusionType* m_pIncType2;     ///< Pointer to type 2 definition
    double &m_Beta;
    double &m_DBeta;
    State* m_pState;                ///< Pointer to simulation state
};

#endif // EQUILIBRIUM_INCLUSION_EXCHANGE_BY_CHEMICAL_POTENTIAL_H
