#if !defined(AFX_GenericFactoryTemplate_H)
#define AFX_GenericFactoryTemplate_H

#include "GenericFactoryTemplate.h"

class AbstractAlexanderMove;
class AbstractInclusionPoseIntegrator;
class AbstractOpenEdgeEvolution;
class AbstractEnergy;
class AbstractSimulation;
class AbstractVisualizationFile;
class AbstractApplyConstraintBetweenGroups;
class AbstractInclusionPoseIntegrator;
class AbstractCurvature;

using FactoryAlexanderMove = Factory<AbstractAlexanderMove>;
using FactoryOpenEdgeEvolutionMethod = Factory<AbstractOpenEdgeEvolution>;
using FactoryEnergyFunctions  = Factory<AbstractEnergy>;
using FactorySimulationScheme  = Factory<AbstractSimulation>;
using FactoryVisualizationFile  = Factory<AbstractVisualizationFile>;
using FactoryConstraintBetweenGroups  = Factory<AbstractApplyConstraintBetweenGroups>;
using FactoryInclusionPoseIntegrator  = Factory<AbstractInclusionPoseIntegrator>;
using FactoryCurvatureMethod  = Factory<AbstractCurvature>;

#endif
