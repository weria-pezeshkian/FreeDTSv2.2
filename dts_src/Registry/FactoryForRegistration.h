#if !defined(AFX_GenericFactoryTemplate_H)
#define AFX_GenericFactoryTemplate_H

#include "GenericFactoryTemplate.h"

using FactoryAlexanderMove = Factory<AbstractAlexanderMove>;

using FactoryOpenEdgeEvolutionMethod = Factory<AbstractOpenEdgeEvolution>;

using FactoryEnergyFunctions  = Factory<AbstractEnergy>;

using FactorySimulationScheme  = Factory<AbstractSimulation>;


#endif
