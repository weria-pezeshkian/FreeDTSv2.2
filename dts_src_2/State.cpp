
#include <fstream>
#include <time.h>
#include "MESH.h"
#include "CreateMashBluePrint.h"
#include "State.h"
State::State(){
    
}
State::State(std::vector<std::string> argument) :
      //============ Initialization of all files
      m_pTimeSeriesDataOutput(new TimeSeriesDataOutput(this)),  // Initialize TimeSeriesDataOutput
      m_pTimeSeriesLogInformation(new TimeSeriesLogInformation(this)),  // Initialize TimeSeriesLogInformation
      m_pRestart(new Restart(this)),  // Initialize Restart
      m_pNonbinaryTrajectory(new Traj_tsi(this)),  // Initialize Traj_tsi
      m_pVisualizationFile(new WritevtuFiles(this)),  // Initialize WritevtuFiles
      m_pBinaryTrajectory(new NoFile),  // Initialize BinaryTrajectory

      //---- Initialize constraint components
      m_pVAHCalculator(new VAHGlobalMeshProperties(this)),
      m_pVolumeCoupling(new NoCoupling(m_pVAHCalculator,this)),
      m_pCoupleGlobalCurvature(new NoGlobalCurvature(m_pVAHCalculator,this)),
      m_pTotalAreaCoupling(new NoTotalAreaCoupling(m_pVAHCalculator, this)),
      m_pForceonVerticesfromInclusions(new NoForce),  // Initialize ForceonVerticesfromInclusions
      m_pExternalFieldOnVectorFields(new NoExternalField),  // Initialize ExternalFieldOnVectorFields
      m_pApplyConstraintBetweenGroups(new NoConstraint),  // Initialize ApplyConstraintBetweenGroups
      m_pBoundary(new PBCBoundary),  // Initialize Boundary

      //---- Initialize integrators
      m_pVertexPositionIntegrator(new EvolveVerticesByMetropolisAlgorithm),  // Initialize VertexPositionIntegrator
      m_pAlexanderMove(new AlexanderMoveByMetropolisAlgorithm),  // Initialize AlexanderMove
      m_pInclusionPoseIntegrator(new InclusionPoseUpdateByMetropolisAlgorithm),  // Initialize InclusionPoseIntegrator
     
      //---- Initialize supplementary integrators
      m_pDynamicBox(new NoBoxChange),  // Initialize DynamicBox
      m_pDynamicTopology(new ConstantTopology),  // Initialize DynamicTopology
      m_pOpenEdgeEvolution(new NoEvolution),  // Initialize OpenEdgeEvolution
      m_pInclusionConversion(new NoInclusionConversion),  // Initialize InclusionConversion

      //--- Initialize accessory objects
      m_pCurvatureCalculations(new CurvatureByShapeOperatorType1),  // Initialize CurvatureCalculations
      m_RandomNumberGenerator(new RNG(1234)),  // Initialize RandomNumberGenerator
      m_pEnergyCalculator(new Energy(this)),  // Initialize EnergyCalculator

      m_pSimulation(new MC_Simulation(this)),  // Initialize Simulation
      m_pParallel_Replica(new ParallelReplicaData(false)),

     //--- Local state member variables
      m_Argument(argument),  // Set Argument
      m_TopologyFile("topology.top"),  // Set TopologyFile
      m_InputFileName("Input.dts"),  // Set InputFileName
      m_IndexFileName("Index.inx"),  // Set IndexFileName
      m_GeneralOutputFilename("dts"),  // Set GeneralOutputFilename
      m_Total_no_Threads(1),  // Set Total_no_Threads
      m_Targeted_State(true),  // Set Targeted_State
      m_RestartFileName(""),  // Set RestartFileName
      m_MinFaceAngle(-0.5),  // Set MinFaceAngle
      m_MinVerticesDistanceSquare(1.0),  // Set MinVerticesDistanceSquare
      m_MaxLinkLengthSquare(3.0)  // Set MaxLinkLengthSquare

{ // start of the "class State" constructor
//-- Find input files name (input.tsi, top.top||top.tsi, index.inx, restart.res)
        if (!ExploreArguments(argument)) {
            exit(0);
     }
        ReadInputFile(m_InputFileName);  // Reading input file, i.e., file.dts
        ExploreArguments(argument);      // Update last state
}
State::~State()
{
    delete m_pTimeSeriesDataOutput;
    delete m_pTimeSeriesLogInformation;
    delete m_pRestart;
    delete m_pNonbinaryTrajectory;
    delete m_pVisualizationFile;
    delete m_pBinaryTrajectory;
    delete m_pVAHCalculator;
    delete m_pVolumeCoupling;
    delete m_pCoupleGlobalCurvature;
    delete m_pTotalAreaCoupling;
    delete m_pForceonVerticesfromInclusions;
    delete m_pExternalFieldOnVectorFields;
    delete m_pApplyConstraintBetweenGroups;
    delete m_pBoundary;
    delete m_pVertexPositionIntegrator;
    delete m_pAlexanderMove;
    delete m_pInclusionPoseIntegrator;
    delete m_pDynamicBox;
    delete m_pDynamicTopology;
    delete m_pOpenEdgeEvolution;
    delete m_pInclusionConversion;
    delete m_RandomNumberGenerator;
    delete m_pEnergyCalculator;
    delete m_pSimulation;
    delete m_pParallel_Replica;
}
bool State::ExploreArguments(std::vector<std::string> &argument){
    /*
     * ExploreArguments Function
     * Description:
     * This function parses the command-line arguments passed to the program and updates the State object accordingly.
     * The function iterates through the vector of arguments, processing each flag and its corresponding value.
     *
     * Arguments:
     * - argument: A vector of strings containing command-line arguments and their values.
     *
     * Flags and Corresponding Actions:
     * - "-h" or "--help": Displays a help message and returns false.
     * - "-in": Sets the input file name in the State object.
     * - "-top": Sets the topology file name in the State object.
     * - "-defout": Sets the general output filename in the State object.
     * - "-b": Updates the initial step value in the simulation object.
     * - "-e": Updates the final step value in the simulation object.
     * - "-seed": Updates the seed value in the simulation object.
     * - "-restart": Sets the restart filename in the State object.
     * - "-nt": Updates the total number of threads in the State object.
     * - "-ndx": Sets the index filename in the State object.
     *
     * Return Value:
     * - True if all arguments are successfully processed, false otherwise.
     *
     * Additional Notes:
     * - The function checks for unknown command-line arguments and displays an error message.
     * - It provides a usage message if an unknown flag is encountered.
     * - String_to_Int function is used to convert string arguments to integer values.
     */
    const std::string HELP_FLAG        = "-h";
    const std::string LONG_HELP_FLAG   = "--help";
    const std::string INPUT_FLAG       = "-in";
    const std::string TOPOLOGY_FLAG    = "-top";
    const std::string INI_STEP_FLAG    = "-b";
    const std::string FINAL_STEP_FLAG  = "-e";
    const std::string SEED_FLAG        = "-seed";
    const std::string RESTART_FLAG     = "-restart";
    const std::string INDEX_FLAG       = "-ndx";
    const std::string DEFOUT_FLAG      = "-defout";
    const std::string THREAD_FLAG      = "-nt";

    for (size_t i=1;i<argument.size();i=i+2)
    {
        const std::string& flag = argument[i];
        
        if (flag == HELP_FLAG || flag == LONG_HELP_FLAG) {
            HelpMessage(); // write a help message
            return false;
        }
        else if(flag == INPUT_FLAG) {
            
            m_InputFileName = argument[i+1];
        }
        else if(flag == TOPOLOGY_FLAG) {
            
            m_TopologyFile = argument[i+1];
        }
        else if(flag == DEFOUT_FLAG)
        {
            m_GeneralOutputFilename = argument[i+1];
        }
        else if(flag == INI_STEP_FLAG){
            
            m_pSimulation->UpdateInitialStep(Nfunction::String_to_Int(argument[i+1]));
        }
        else if(flag == FINAL_STEP_FLAG){
            
            m_pSimulation->UpdateFinalStep(Nfunction::String_to_Int(argument[i+1]));
        }
        else if(flag == SEED_FLAG){
            
            m_RandomNumberGenerator = new RNG(Nfunction::String_to_Int(argument[i+1]));
        }
        else if(flag == RESTART_FLAG){
            
            m_RestartFileName = argument[i+1];
        }
        else if(flag == THREAD_FLAG){
            
            m_Total_no_Threads = Nfunction::String_to_Int(argument[i+1]);
        }
        else if(flag == INDEX_FLAG){   // get the index file and also check if the file exist
            m_IndexFileName = argument[i+1];
        }
        else
        {
            std::cout << "--> error: unknown commandline argument! "<<flag<<std::endl;
            std::cout << "    For more information and tips run "<< EXE_NAME <<" -h"<<std::endl;
            return false;
        }
    }

    return true;
}
bool State::ReadInputFile(std::string file)
{
    // enforcing correct file extension:  check simdef file for value of InExt
    std::string ext = file.substr(file.find_last_of(".") + 1);
    if(ext!=InExt)
    file = file + "." + InExt;
    if (Nfunction::FileExist(file)!=true)
    {
        std::cout<<"----> Error: the input file with the name "<<file<< " does not exist "<<std::endl;
        return false;
    }

    std::ifstream input;
    input.open(file.c_str());
    std::string firstword,rest,str;

    while (true)
    {
        input>>firstword;
        if(input.eof())
            break;
        
        
        if(firstword == "Initial_Step")
        {
            int ini;
            input>>str>>ini;
            m_pSimulation->UpdateInitialStep(ini);
            getline(input,rest);
        }
        else if(firstword == "Final_Step")
        {
            int fi;
            m_pSimulation->UpdateFinalStep(fi);
            getline(input,rest);
        }
        else if(firstword == "Box_Centering_F")
        {
            double rate;
            input>>str>>rate;
            m_pSimulation->SetCentering(rate);
            getline(input,rest);
        }
        else if(firstword == "Simulation")
        {
            std::string type;
            input >> str >> type;
            if(type=="MC"){
                m_pSimulation = new MC_Simulation(this);
            }
            getline(input,rest);
        }
        else if(firstword == "MC_Moves")
        {
            getline(input, str);
            std::vector<std::string> data = Nfunction::split(str);

            // Check if the correct number of moves is provided
            const int expectedMoves = 5;
            if (data.size() != expectedMoves + 1) { // +1 to account for the first element is just = sign
                std::cout << "---> error: the defined MC moves does not cover all the moves. It should be " << expectedMoves << " numbers" << std::endl;
                exit(0);
            }
            m_pVertexPositionIntegrator->SetMoveRate(Nfunction::String_to_Double(data[1]), Nfunction::String_to_Double(data[2]));
            m_pAlexanderMove->SetMoveRate(Nfunction::String_to_Double(data[3]));
            m_pInclusionPoseIntegrator->SetMoveRate(Nfunction::String_to_Double(data[4]),Nfunction::String_to_Double(data[5]));

        }
        else if(firstword == "FreezingAGroup")
        {
            input>>str>>str;
            m_pVertexPositionIntegrator->UpdateFreezGroupName(str);
            getline(input,rest);
        }
        else if(firstword=="InclusionConversion")
        {
            std::string type;
            input >> str >> type;
            if(type == "ActiveTwoStateInclusion"){
                std::string inctype1, inctype2;
                double ep1,ep2,persentage,gama;
                input>>inctype1>>inctype2>>ep1>>ep2>>persentage>>gama;
                m_pInclusionConversion = new ActiveTwoStateInclusion(ep1, ep2, persentage, gama, inctype1,inctype2);

            }
            else{
                std::cout<<" unknown InclusionConversion method "<<std::endl;
            }
                
            getline(input,rest);
        }
        else if(firstword == AbstractBoundary::GetBaseDefaultReadName()) // " Boundary "
        {
            std::string type;
            input >> str >> type;
            if(type == TwoFlatParallelWall::GetDefaultReadName()){   // " "TwoFlatParallelWall" "
                double thickness;
                char direction;
                input>>str>>thickness>>direction;
                m_pBoundary = new TwoFlatParallelWall(this, thickness,direction);
            }
            else {
                std::cout<<" unknown Boundary type: "<<type<<std::endl;
            }
            getline(input,rest);
        }
        else if(firstword == "ForceFromInclusions")
        {
            std::string type;
            input >> str >> type;
            if(type == "Constant_NematicForce"){
                double f0; //(d,p,n)
                input>>str>>f0;
                m_pForceonVerticesfromInclusions = new Constant_NematicForce(f0);
            }
            else {
                std::cout<<" unknown ForceFromInclusions method "<<std::endl;
            }
            getline(input,rest);
        }
        else if(firstword == m_pDynamicBox->GetBaseDefaultReadName())
        {
            
            std::string type;
            int period = 0;
            double force = 0;

            input >> str >> type >> period >> force;

            if (type == "SingleSideConstantForce") {
                m_pDynamicBox = new DynamicBoxSide(period, force, this);
            }
            else if (type == "PositionRescale_FrameTension") {
                m_pDynamicBox = new PositionRescaleFrameTensionCoupling(period, force,this);
            }

            // Consume remaining input line
            getline(input, rest);
        }
        else if(firstword == "Dynamic_Topology")
        {
            
            std::string type;
            int period = 0;

            input >> str >> type >> period;

            if (type == "Scission_By3Edges") {
                
                m_pDynamicTopology = new Three_Edge_Scission(period, this);
            }
            // Consume remaining input line
            getline(input, rest);
        }
        else if(firstword == "OpenEdgeEvolution")
        {
            // OpenEdgeEvolution =  EvolutionWithConstantVertex PT_steps
            std::string type;
            int period = 0;
            input >> str >> type >> period;
            getline(input,rest);
            if (type == "EvolutionWithConstantVertex") {
                
                m_pOpenEdgeEvolution = new OpenEdgeEvolutionWithConstantVertex(period, this);
            }
            else if (type == "EvolutionWithConstantVertex_B") {
                
                std::cout<<"---> warning: this algorithm for edge treatment has errors. It is up to you to use it \n";
                m_pOpenEdgeEvolution = new OpenEdgeEvolutionWithConstantVertex_B(period, this);
            }
            else {
                
                std::cout<<"---> error: unknown Open Edge Evolution type: "<<type<<"\n";
            }
        }
        else if(firstword == "GlobalCurvature")
        {
            std::string type;
            double k,gc0;
            input>>str>>type>>k>>gc0;
            if(type == "CoupledToHarmonicPotential") {
                m_pCoupleGlobalCurvature = new CouplingtoFixedGlobalCurvature(k,gc0);
            }
            else {
                std::cout<<"---> error: this algorithm for global curvature is unknown \n";
                std::cout<<"----------- no global curvature coupling with be applied \n";
            }
            getline(input,rest);
            
        }
        else if(firstword == "Volume_Constraint")
        {
            
            std::string type;
            input >> str >> type;
            if(type == "SecondOrderCoupling"){
                int eqsteps = 0;
                double dp,k,vt;
                input>>eqsteps>>dp>>k>>vt;
                m_pVolumeCoupling = new CmdVolumeCouplingSecondOrder(eqsteps, dp,  k, vt);
            }
            else if(type == "Osmotic_Pressure"){
                int eqsteps = 0;
                double gamma,p0;
                input>>eqsteps>>gamma>>p0;
                m_pVolumeCoupling = new Apply_Osmotic_Pressure(eqsteps,gamma,p0);
            }
            getline(input,rest);
        }
        else if(firstword == "TotalAreaCoupling"){
            
            std::string type;
            int eq_steps;
            double gamma , k0;
            input>>str>>type;

            if(type=="CouplingTotalAreaToHarmonicPotential"){
                input>>eq_steps>>gamma>>k0;
                m_pTotalAreaCoupling = new CouplingTotalAreaToHarmonicPotential(eq_steps,gamma,k0);
            }
            getline(input,rest);
        }
        else if(firstword == "ConstraintBetweenGroups")
        {
            std::string type;
            input >> str >> type;
            if(type == "HarmonicPotentialBetweenTwoGroups"){
                // ConstraintBetweenGroups  = HarmonicPotentialBetweenTwoGroups 10 0.1 2000 Group1 Group2 0 1 1
                double k,dx,rate,nx,ny,nz;
                std::string g1,g2;
                input>>k>>dx>>rate>>g1>>g2>>nx>>ny>>nz;
                m_pApplyConstraintBetweenGroups = new HarmonicPotentialBetweenTwoGroups(k, dx,rate,g1,g2,nx,ny,nz);
            }
            else{
                std::cout<<" unknown type of ConstraintBetweenGroups "<<std::endl;
            }
            getline(input,rest);

        }
        else if(firstword == "TotalAreaCoupling"){
            std::string type;
            input>>str>>type;
            if(type == "Spherical"){
                //m_pBoundary
            }
        }
        else if(firstword == "MinfaceAngle")
        {
            input>>str>>m_MinFaceAngle;
            getline(input,rest);
        }
        else if(firstword == "OutPutEnergy_periodic")
        {
            int period;
            input>>str>>period;
            m_pTimeSeriesDataOutput->UpdatePeriod(period);
            getline(input,rest);
        }
        else if(firstword == "Restart_periodic")
        {
            int period;
            input>>str>>period;
            m_pRestart->UpdatePeriod(period);
            getline(input,rest);
        }
        else if(firstword == "TopologyFile")
        {
            input>>str>>m_TopologyFile;
            getline(input,rest);
        }
        else if(firstword == "Seed")
        {
            int seed;
            input>>str>>seed;
            m_RandomNumberGenerator = new RNG(seed);
            getline(input,rest);
        }
        else if(firstword == "Kappa")
        {
            double k,kg,c0;
            input>>str>>k>>kg>>c0;
            m_pEnergyCalculator->SetSurfRigidity(k,kg,c0);
            getline(input,rest);
        }
        else if(firstword == "Edge_Parameters")
        {
            double lambda,kg,kn;
            input>>str>>lambda>>kg>>kn;
            m_pEnergyCalculator->SetEdgeRigidity(lambda,kg,kn);
            //=== send it to force field class
            getline(input,rest);
        }
        else if(firstword == "ConstantField") {
            double k,x,y,z;
            input>>str>>k>>x>>y>>z;
            m_pExternalFieldOnVectorFields = new ConstantExternalField(k,x,y,z);
            getline(input,rest);

        }
        else if(firstword == "VertexArea")
        {
            double a,b,c,d;
            input>>str>>a>>b>>c>>d;
            double A01 = (1+2*b)*sqrt(3)/2.0;   // selecting b between 0-1
            double A02 = sqrt(1+2*d);   // selecting d between 0-1

            if(b<0 || b>1 || d<0 || d>1){
                std::cout<<"---> error in constant vertex area; gamma is bad; it should be a double number between 0-1 \n";
                exit(0);
            }
            m_pEnergyCalculator->SetSizeCoupling (a,A01,c,A02);
        }
        else if(firstword == "Display_periodic")
        {
            int period;
            input>>str>>period;
            m_pVisualizationFile->SetPeriod(period);
            getline(input,rest);
        }
        else if(firstword == "Voxel_Size"){
            double lx,ly,lz;
            input>>str>>lx>>ly>>lz;
            m_pVoxelization->UpdateVoxelSize(lx, ly,  lz);
            getline(input,rest);
        }
        else if(firstword == "GeneralOutputFilename")
        {
            input>>str>>m_GeneralOutputFilename;
            getline(input,rest);
        }
        else if(firstword == "Min_Max_LinkLenghtsSquare")
        {
            input>>str>>m_MinVerticesDistanceSquare>>m_MaxLinkLengthSquare;
            getline(input,rest);
        }
        else if(firstword == Traj_tsi::GetDefaultReadName()){ // "OutPutTRJ_TSI"
            
            int period;
            std::string tsiPrecision, tsiFolder_name;
            input>>str>>period>>tsiPrecision>>tsiFolder_name;
            m_pNonbinaryTrajectory  = new Traj_tsi(this, period, tsiFolder_name, tsiPrecision );
        }
        else if(firstword == "Temprature"){
            
            double beta,delta_beta;
            input>>str>>beta>>delta_beta;
            m_pSimulation->SetBeta(beta,delta_beta);

        }
        else if(firstword == "Parallel_Tempering")
        {
            // Parallel_Tempering  = on PT_steps  PT_minbeta    PT_maxbeta
            std::string state;
            input>>str>>state>>(m_pParallel_Replica->PT_steps)>>(m_pParallel_Replica->PT_minbeta)>>(m_pParallel_Replica->PT_maxbeta);
            getline(input,rest);
            if(state=="on"|| state=="yes"|| state=="On"|| state=="ON"|| state=="Yes"|| state=="YES")
                m_pParallel_Replica->State = true;
            else
                m_pParallel_Replica->State = false;
        }
        else if(firstword == BTSFile::GetDefaultReadName() ){ // "OutPutTRJ_BTS"
            int periodic, precision;
            std::string filename;
            input>>str>>periodic>>precision>>filename;
            m_pBinaryTrajectory = new BTSFile(this,periodic,precision,filename);
            getline(input,rest);

        }
        else if(firstword == "INCLUSION")
        {
            break;
        }
        else
        {
            if(firstword.at(0)!=';')
            {
                std::cout<<"Error: bad keyword in the input file *** "<<firstword<<" ***\n";
                return false;
            }
            getline(input,rest);
        }
    }
    
    // read inclusion
    ReadInclusionType(input);
    
    input.close();
    
    return true;
}
void State::HelpMessage(){
    
    std::cout<<"----------------------copyright: Weria Pezeshkian------------------"<<"\n";
    std::cout<<"---------------------- weria.pezeshkian@gmail.com ------------------"<<"\n";
    std::cout<<"------------------------------ FreeDTS -----------------------------"<<"\n";
    std::cout<<" =================================================================  \n";
    std::cout<<"---------------------version "<<SoftWareVersion<<" ------------------"<<"\n";
    std::cout<<" =================================================================  \n";
    // Print example command
    std::cout << "Example: DTS -f Input.dts -top topol.q -restart res.res" << std::endl;
    std::cout << "=============================================================================" << std::endl;
    
    // Print options
    std::cout <<"     "<<std::setw(20) << std::left << "Option" << std::setw(20) << "Type" << std::setw(20) << "Default" << "Description" << std::endl;
    std::cout <<"----------------------------------------------------------------------------------" << std::endl;
    std::cout <<"     "<< std::setw(20) << std::left << "-in" << std::setw(20) << "string" << std::setw(20) << "Input.dts" << "Input file name" << std::endl;
    std::cout <<"     "<< std::setw(20) << std::left << "-top" << std::setw(20) << "string" << std::setw(20) << "topology.top" << "Topology file name" << std::endl;
    std::cout <<"     "<< std::setw(20) << std::left << "-b" << std::setw(20) << "int" << std::setw(20) << "1" << "Initial time step" << std::endl;
    std::cout <<"     "<< std::setw(20) << std::left << "-e" << std::setw(20) << "int" << std::setw(20) << "10" << "Final time step" << std::endl;
    std::cout <<"     "<< std::setw(20) << std::left << "-seed" << std::setw(20) << "int" << std::setw(20) << "36723" << "Random number seed" << std::endl;
    std::cout <<"     "<< std::setw(20) << std::left << "-defout" << std::setw(20) << "string" << std::setw(20) << "dts" << "Output file prefix" << std::endl;
    std::cout <<"     "<< std::setw(20) << std::left << "-ndx" << std::setw(20) << "string" << std::setw(20) << "Index.inx" << "Index file name" << std::endl;
    std::cout <<"     "<< std::setw(20) << std::left << "-nt" << std::setw(20) << "int" << std::setw(20) << "1" << "Total number of threads" << std::endl;
    std::cout <<"     "<< std::setw(20) << std::left << "-restart" << std::setw(20) << "string" << std::setw(20) << "NO" << "Restart file name" << std::endl;
    std::cout <<"====================================================================================" << std::endl;
    
    return;
}
bool State::Initialize(){
//-----> Get the mesh
        // Create a MeshBluePrint object
        CreateMashBluePrint Create_BluePrint;
        MeshBluePrint mesh_blueprint;
//---->  Check if the restart file name is provided
        bool restartReadSuccess = false;
        if (!m_RestartFileName.empty()) {
            int step;
            double r_vertex;
            double r_box;

            // Attempt to open the restart file
            std::cout << "---> Note: attempting to open the restart file: " << m_RestartFileName << std::endl;

            // Read the restart file and update the State object to that state
            mesh_blueprint = m_pRestart->ReadFromRestart(m_RestartFileName, step, restartReadSuccess, r_vertex, r_box);

            // Check if the restart file was successfully read
            if (restartReadSuccess) {
                std::cout << "---> Note: Restart file was successfully read" << std::endl;
                m_pVertexPositionIntegrator->UpdateDR(r_vertex);
                m_pDynamicBox->UpdateDR(r_box);
                m_pSimulation->UpdateInitialStep(step+1);
                //----- open log file
                m_pTimeSeriesLogInformation->OpenFile(false);
                //---- open the binar trajectory
                m_pBinaryTrajectory->OpenFile(false, 'w');
            }
            else{
                // If failed to read restart file, generate MeshBluePrint from input topology file
                std::cout << "---> Note: Failed to read restart file, using topology file" << std::endl;
            }
        }
        if(!restartReadSuccess){ // this for also a situation where m_RestartFileName is empty
            mesh_blueprint = Create_BluePrint.MashBluePrintFromInput_Top(m_InputFileName, m_TopologyFile);
            
            //----- open log file
            m_pTimeSeriesLogInformation->OpenFile(false);
            
            // Open binary trajectory file for writing
            m_pBinaryTrajectory->OpenFile(true, 'w');

            // Open folder for non-binary trajectory
            m_pNonbinaryTrajectory->OpenFolder();
            
            // open folder for Visualization
            if(m_pVisualizationFile->GetDerivedDefaultReadName()=="VTUFileFormat"){
               m_pVisualizationFile->OpenFolder();
            }

        }

        // Generate mesh from the mesh blueprint
        m_Mesh.GenerateMesh(mesh_blueprint);
    
        // Set the pointer to the mesh object
        m_pMesh = &m_Mesh;

        // Update group from index file
        m_pMesh->UpdateGroupFromIndexFile(m_IndexFileName);
    
//============ activate of all inputs and data structures
        m_pRestart->SetRestartFileName();
    
//----> calaculate curvature for all vertices
        m_pCurvatureCalculations->Initialize(this);

    //----> set some easy access for Vertex Position Integrator
        m_pVertexPositionIntegrator->Initialize(this);

//----> boundry of the simulations
        m_pBoundary->Initialize();
    
//----> energy class
    //---> to get interaction energies
     m_pEnergyCalculator->Initialize(m_InputFileName);
    //---> to update each vertex and edge energy. Up to now May 2024, edge energy is not zero when both vertices has inclusions
     m_pEnergyCalculator->CalculateAllLocalEnergy();
    
    /*
    //---- Initialize constraint components
     m_pVolumeCoupling
    m_pCoupleGlobalCurvature
    m_pTotalAreaCoupling               
    m_pForceonVerticesfromInclusions
    m_pApplyConstraintBetweenGroups
    
    // Initialize integrators
    m_pAlexanderMove                    = new AlexanderMoveByMetropolisAlgorithm;
    m_pInclusionPoseIntegrator          = new InclusionPoseUpdateByMetropolisAlgorithm;
   
    // Initialize supplementary integrators
    m_pDynamicBox                       = new NoBoxChange;
    m_pDynamicTopology                  = new ConstantTopology;
    m_pOpenEdgeEvolution                = new NoEvolution;
    m_pInclusionConversion              = new NoInclusionConversion;
    
    */
    
        m_pSimulation->Initialize();
    
    return true;
}
bool State::ReadInclusionType(std::ifstream& input) {
    std::string firstword, rest, str1, str2, TypeNames;
    int N, TypeID, NoType;
    double Kappa, KappaG, KappaP, KappaL, C0, C0P, C0N;

    // Store inclusion types in a vector
    std::vector<InclusionType> all_InclusionType;

    // Add a default inclusion type
    InclusionType emptyIncType;
    all_InclusionType.push_back(emptyIncType);

    // Read the header line
    input >> str1 >> NoType >> str2;
    getline(input, rest);
    getline(input, rest); // Discard the header line

    // Check if the header line indicates inclusion type definition
    if (str1 == "Define" || str1 == "define" || str1 == "DEFINE") {
        for (int i = 0; i < NoType; i++) {
            input >> N >> TypeNames >> Kappa >> KappaG >> KappaP >> KappaL >> C0 >> C0P >> C0N;
            std::string edgedata;
            getline(input, edgedata);
            std::vector<std::string> edge_data = Nfunction::split(edgedata);

            // Parse edge data if available
            double lam = 0, ekg = 0, ekn = 0, ecn = 0;
            if (edge_data.size() >= 4) {
                lam = Nfunction::String_to_Double(edge_data[0]);
                ekg = Nfunction::String_to_Double(edge_data[1]);
                ekn = Nfunction::String_to_Double(edge_data[2]);
                ecn = Nfunction::String_to_Double(edge_data[3]);
            }

            // Create inclusion type and add to vector
            InclusionType incType(TypeNames, i + 1, N, Kappa/2, KappaG, KappaP/2, KappaL/2, C0, C0P, C0N, lam, ekg, ekn, ecn);
            all_InclusionType.push_back(incType);
        }
    }

    // Set inclusion types in the mesh object
    m_Mesh.m_InclusionType = all_InclusionType;

    // Set pointers to inclusion types
    m_Mesh.m_pInclusionType.clear();
    for (size_t i = 0; i < m_Mesh.m_InclusionType.size(); ++i) {
        m_Mesh.m_pInclusionType.push_back(&m_Mesh.m_InclusionType[i]);
    }

    return true;
}

