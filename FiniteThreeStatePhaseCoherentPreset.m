function qkdInput = FiniteThreeStatePhaseCoherentPreset(lift)
% FiniteThreeStatePhaseCoherent A preset for three-state protocol where
% Alice sends phase coherent states and we use flag state-squasher 
% Bob's side. 

qkdInput = QKDSolverInput();

%% Parameters

lossdB = linspace(0,20,21);
lossEta = 10.^(-lossdB/10);
qkdInput.addScanParameter("eta", num2cell(lossEta));

qkdInput.addFixedParameter("protocolTime",3600); % protocol run time in seconds
qkdInput.addFixedParameter("laserRepetitionRate",3*10^9); % Assume we have a  3 GHz source

qkdInput.addFixedParameter("pzAlice",0.8);


qkdInput.addFixedParameter("basisChoiceBeamSplittingRatio",0.2); %to X

qkdInput.addFixedParameter("f",1.16); %effiecienct of error-correction.

qkdInput.addFixedParameter("photonNumberCutOff",1);

qkdInput.addFixedParameter("darkCountRate", 0);

% finite size parameters.

qkdInput.addFixedParameter("alphabet", 2); % encoding alphabet size; for qubits, this is 2

qkdInput.addFixedParameter("lift",lift); % decides whether we use postselection lift or not.

% x in g_{n,x}; vac, qubit, 8 flags.
switch lift
    case "Postselection-no-improvement"
        liftDimension = 2^2*11^4; % No block-diagonal structure, qubit, no d_B^2 improvement
    case "Postselection-generic"
        liftDimension = 22^2; % No block-diagonal structure; qubit
    case "Postselection-blockdiagonal"
        liftDimension = 4^2+9*2^2; % Block-diagonal qubit. See Section V.A. for explanation.
    case "Postselection-decoy-generic"
        liftDimension = 33^2*900; % No block-diagonal structure; decoy with photon-number cut-off of 8.
    case "Postselection-decoy"
        liftDimension = (6^2+9*3^2)*30; % Block-diagonal structure; decoy with photon-number cut-off of 8.
    otherwise
        liftDimension = 1; % Placeholder value so that the code runs for non-PS cases.
end

qkdInput.addFixedParameter("liftDimension",liftDimension);



log2targetepsSec = -12 / log10(2);
log2targetepsCor = -12 / log10(2);  % heuristic choice. 
qkdInput.addFixedParameter("log2targetepsSec",log2targetepsSec);
qkdInput.addFixedParameter("log2targetepsCor",log2targetepsCor);

epsBiasFactor = 1/2;
qkdInput.addFixedParameter("epsBiasFactor",epsBiasFactor); %adjusts bias between (epsPA), and (epsAT)

qkdInput.addFixedParameter("ptest", 0.05); % ptest*N signals used for testing
 
descriptionModule = QKDDescriptionModule(@ThreeStateDescriptionFunc);
qkdInput.setDescriptionModule(descriptionModule);

% channel model. 
channelModule = QKDChannelModule(@ThreeStateFlagStateChannelFunc);
qkdInput.setChannelModule(channelModule);

% Key rate module performs squashing and decoy analysis
keyRateOptions = struct();
keyRateOptions.decoyTolerance = 1e-14;
keyRateOptions.decoySolver = "mosek";
keyRateOptions.decoyForceSep = true;
keyMod = QKDKeyRateModule(@FiniteThreeStateKeyRateFunc, keyRateOptions);
qkdInput.setKeyRateModule(keyMod);

optimizerMod = QKDOptimizerModule(@coordinateDescentFunc,struct("verboseLevel",0),struct("verboseLevel",0));
qkdInput.setOptimizerModule(optimizerMod);

% math solver options
mathSolverOptions = struct();
mathSolverOptions.initMethod = 1;
mathSolverOptions.maxIter = 50;
mathSolverOptions.blockDiagonal = true;
mathSolverMod = QKDMathSolverModule(@FW2StepSolver,mathSolverOptions,mathSolverOptions);
qkdInput.setMathSolverModule(mathSolverMod);

% global options
qkdInput.setGlobalOptions(struct("errorHandling",3,"verboseLevel",1,"cvxSolver","mosek","cvxPrecision","default"));
