#if !defined(AFX_GenericFactoryTemplate_H)
#define AFX_GenericFactoryTemplate_H

#include "GenericFactoryTemplate.h"  // will be only included 

class AbstractAlexanderMove;
class AbstractInclusionPoseIntegrator;
class AbstractOpenEdgeEvolution;
class AbstractEnergy;
class AbstractSimulation;
class AbstractVisualizationFile;
class AbstractApplyConstraintBetweenGroups;
class AbstractInclusionPoseIntegrator;
class AbstractCurvature;
class AbstractExternalFieldOnInclusions;
class AbstractVectorFieldsRotationMove;
class AbstractInclusionConversion;
class AbstractNonbondedInteractionBetweenVertices;

using FactoryAlexanderMove = Factory<AbstractAlexanderMove>;
using FactoryOpenEdgeEvolutionMethod = Factory<AbstractOpenEdgeEvolution>;
using FactoryEnergyFunctions  = Factory<AbstractEnergy>;
using FactorySimulationScheme  = Factory<AbstractSimulation>;
using FactoryVisualizationFile  = Factory<AbstractVisualizationFile>;
using FactoryConstraintBetweenGroups  = Factory<AbstractApplyConstraintBetweenGroups>;
using FactoryInclusionPoseIntegrator  = Factory<AbstractInclusionPoseIntegrator>;
using FactoryCurvatureMethod  = Factory<AbstractCurvature>;
using FactoryExternalFieldOnInclusions  = Factory<AbstractExternalFieldOnInclusions>;
using FactoryVectorFieldsRotationMove  = Factory<AbstractVectorFieldsRotationMove>;
using FactoryInclusionConversionMethod = Factory<AbstractInclusionConversion>;
using FactoryNonbondedInteractionBetweenVertices = Factory<AbstractNonbondedInteractionBetweenVertices>;



#endif
