function [keyRate, modParser, debugInfo] = FiniteThreeStateKeyRateFunc(params,options,mathSolverFunc,mathSolverOptions,debugInfo)
% FiniteThreeStateKeyRateFuncRenyi A finite size key rate function for the
%  Three state BB84 protocol.
%
% Input parameters:
% * dimA: dimension of Alice's system.
% * dimB: dimension of Bob's system.
% * f: error correction effiency. A 1 means we are correcting at the
%   Shannon limit. A more practical value would be around 1.16.
% * observablesJoint: The joint observables from Alice and Bob's
%   measurments. These are organized in a 4x6 table. The
%   observables must be hermitian and each must be the size dimA*dimB by
%   dimA*dimB. The observables assume the spaces are ordered A \otimes B.
%   They also should be positive semi-definite and should sum to identity,
%   but this is hard to check because of machine precision issues.
% * protocolTime: Time in seconds for which protocol runs.
% * laserRepetitionRate: Frequency of laser being used in Hz.
% * lift: Type of Coherent lift being used.
% * distance: Distance over which QKD is being performed.
% * eps: (finite size) a struct containing four epsilons used in finite key
%   rate analysis: PE (parameter estimation), bar (smoothing min entropy),
%   EC (error correction), PA (privacy amplification)
% * ptest: the fraction of signals sent that are used for testing.
%   Typically should be fairly low. This cuts into key rate, but also
%   entirely destroys key rate if it is too low.
% * expectationsJoint : Joint expectations of Bob's outcome
%   Alice sent a given signal.
% * ObservablesJoint : POVMs corresponding to Alice sending a given signal
%   and Bob obtaining a given outcome,
% * log2targetepsSec : target values of log2(targetepsSecurity)
% * log2targetepsSCor: target values of log2(targetepsCorrectness)
% * epsBiasFactor : bias between epsPA and epsAT  for the IID calculations
% * liftdimension : the dimension of the postselection list. Corresponds to
%   x in g_{n,x}
% * protocolTime : total protocol run time
% * laserRepetitionRate : Repetition rate of Alice's  laser.
% * distance : distance between Alice and Bob (in metres)
%
% Outputs:
% * keyrate: Key rate of the QKD protocol.
% Options:
% * verboseLevel: (global option) See makeGlobalOptionsParser for details.
% * errorHandling: (global option) See makeGlobalOptionsParser for details.
% DebugInfo:
% * krausSum: sum_i K^\dagger_i*K_i. For a completely positive trace
%   non-increasing map, this sum should be <=I. 
%
% See also QKDKeyRateModule, PM46DescriptionFunc, makeGlobalOptionsParser
arguments
    params (1,1) struct
    options (1,1) struct
    mathSolverFunc (1,1) function_handle
    mathSolverOptions (1,1) struct
    debugInfo (1,1) DebugInfo
end

optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.parse(options);
options = optionsParser.Results;


modParser = moduleParser(mfilename);

modParser.addRequiredParam("observablesJoint",@(x) allCells(x,@ishermitian));

modParser.addRequiredParam("expectationsJoint",@mustBeProbDist);
modParser.addAdditionalConstraint(@isEqualSize,["observablesJoint","expectationsJoint"]);

modParser.addRequiredParam("krausOps", @isCPTNIKrausOps);
modParser.addRequiredParam("keyProj", @(x) mustBeAKeyProj(x));

modParser.addRequiredParam("dimA",@mustBeInteger);
modParser.addRequiredParam("dimB", @mustBeInteger);
modParser.addAdditionalConstraint(@mustBePositive,"dimA")
modParser.addAdditionalConstraint(@mustBePositive,"dimB")
modParser.addAdditionalConstraint(@observablesAndDimensionsMustBeTheSame,["observablesJoint","dimA","dimB"])

modParser.addRequiredParam("announcementsA")
modParser.addRequiredParam("announcementsB")
modParser.addRequiredParam("keyMap",@(x)mustBeA(x,"KeyMapElement"))

modParser.addRequiredParam("f", @(x) mustBeGreaterThanOrEqual(x,1));
modParser.addRequiredParam("rhoA",@isDensityOperator);
modParser.addRequiredParam("alphabet", @(x) mustBeInteger(x));

%% finite key analysis parameters
modParser.addRequiredParam("ptest", @(x) mustBeInRange(x, 0, 1));

modParser.addRequiredParam("log2targetepsSec");
modParser.addRequiredParam("log2targetepsCor");
modParser.addRequiredParam("epsBiasFactor", @(x) mustBeInRange(x, 0, 1));
modParser.addRequiredParam("liftDimension");

modParser.addRequiredParam("protocolTime", @(x) mustBeGreaterThan(x, 0));
modParser.addRequiredParam("laserRepetitionRate", @(x) mustBeGreaterThan(x, 0));
modParser.addRequiredParam("lift");
modParser.addRequiredParam("distance");

modParser.addOptionalParam("blockDimsA", nan);
modParser.addOptionalParam("blockDimsB", nan);
modParser.addAdditionalConstraint(@(x,y) blockDimsMustMatch(x,y),["blockDimsA","dimA"]);
modParser.addAdditionalConstraint(@(x,y) blockDimsMustMatch(x,y),["blockDimsB","dimB"]);
modParser.addAdditionalConstraint(@(blockDimsA,blockDimsB) ~xor(all(isnan(blockDimsA),"all"),all(isnan(blockDimsB),"all")),["blockDimsA","blockDimsB"]);

modParser.parse(params);

params = modParser.Results;

%% extra parameters

%% parse parameters
modParser.parse(params);

params = modParser.Results;


debugMathSolver = debugInfo.addLeaves("mathSolver");
mathSolverInput = struct();


%%% we start computations %%%


[deltaLeak, gains] = errorCorrectionCost(params.announcementsA,params.announcementsB, params.expectationsJoint ,params.keyMap,params.f); 
debugInfo.storeInfo("deltaLeak",deltaLeak);
debugInfo.storeInfo("gains",gains);





%% finite size calculations.


time = params.protocolTime;
laserRepetitionRate = params.laserRepetitionRate;
distance = params.distance;

log2targetepsSec = params.log2targetepsSec;
log2targetepsCor = params.log2targetepsCor;

liftDimension = params.liftDimension;

switch params.lift
    case 'IID'

        N=laserRepetitionRate*time; % Only limited by laser repetition rate

        logCost = 0;
        params.logCost = logCost;

        logepsilon.EC = log2targetepsCor; % failure probability for error-correction
        
        logepsilon.AT = log2targetepsSec +log2(1-params.epsBiasFactor); % failure probability for acceptance test 
        logepsilon.PA = log2targetepsSec+ log2(params.epsBiasFactor); % failure probability for privacy amplification
        keyCost = 0;

    case 'Sequential-IID'

        if distance == 0
            N = laserRepetitionRate*time;
        else
            repRate=min(laserRepetitionRate*time,1.5*10^8/distance); % Eq. (41) from https://arxiv.org/pdf/2301.11340.pdf
            N = repRate*time;
        end

        logCost = 0;
        params.logCost = logCost;
            
        logepsilon.EC = log2targetepsCor; % failure probability for error-correction
        
        logepsilon.AT = log2targetepsSec +log2(1-params.epsBiasFactor); % failure probability for acceptance test 
        logepsilon.PA = log2targetepsSec+ log2(params.epsBiasFactor); % failure probability for privacy amplification
        keyCost = 0;

    otherwise % For postselection technique
        
        N=laserRepetitionRate*time; % Only limited by laser repetition rate
        
        logCostOld = (liftDimension-1)*log2(N+1);
        logCostNew = (liftDimension-1)*(log2(exp(1))+log2(N)+(liftDimension-1)/(N*log(2)) - log2(liftDimension-1));
        if (logCostNew > logCostOld)
            warning("The new log cost is HIGHER \n");
        end
        logCost = min(logCostNew,logCostOld);
        params.logCost = logCost;

        logepsilon.EC = log2targetepsCor; %failure probability for error-correction

        %Assuming \sqrt{8 \epsSec} = \epsilontilde / 2 = targetEpsSec /
        %2g_{n,x} as stated in the paper.
                            
        logepsilon.PA = 2*(log2targetepsSec - logCost) + log2(params.epsBiasFactor)-5;
        logepsilon.AT =2*(log2targetepsSec - logCost) +log2(1-params.epsBiasFactor)-5;
        keyCost = 4*logCost - 2*log2targetepsSec; %Assuming \sqrt{8 \epsSec} = \epsilontilde / 2 = targetEpsSec / 2g_{n,x}
end

params.keyCost = keyCost;
params.logepsilon = logepsilon;

debugInfo.storeInfo("numberOfSignals",N)
params.N = N;

m = N*params.ptest; % received testing signals

% we construct muball for entrywise.
muBall = sqrt( ( log2(2*numel(params.expectationsJoint))-logepsilon.AT )/(2*m)  );
params.muBall = muBall;

debugInfo.storeInfo("muBall",muBall);



%Now, we add constraints. Recall that the POVMs correspond to Bob and Alice


Constraints = arrayfun(@(index)InequalityConstraint(...
    params.observablesJoint{index},params.expectationsJoint(index) - muBall ,...
    params.expectationsJoint(index)+muBall), 1:numel(params.observablesJoint));

mathSolverInput.inequalityConstraints = [Constraints];




%% Translate for math solver
mathSolverInput.krausOps = params.krausOps;
mathSolverInput.keyProj = params.keyProj;
mathSolverInput.rhoA = params.rhoA;
if mathSolverOptions.blockDiagonal
    mathSolverInput.blockDimsA = params.blockDimsA;
    mathSolverInput.blockDimsB = params.blockDimsB;
end


[relEnt,~] = mathSolverFunc(mathSolverInput,mathSolverOptions, debugMathSolver);


[keyRate, keyLength] = finiteKeyRate(relEnt, deltaLeak, params);

if isfield(debugMathSolver.info,"relEntStep2Linearization")
    relEntStep2Linearization = debugMathSolver.info.relEntStep2Linearization; 
    
    keyRateStep2Linearization = finiteKeyRate( relEntStep2Linearization, deltaLeak, params);
    if options.verboseLevel>=2
        fprintf("Key rate using step 2 linearization intial value: %e\n",max(keyRateStep2Linearization,0))
    end
    debugInfo.storeInfo("keyRateRelEntStep2Linearization",keyRateStep2Linearization)
    
end




if options.verboseLevel>=1
    fprintf("Key rate: %e\n",keyRate);
end



debugInfo.storeInfo("logCost", logCost);
debugInfo.storeInfo("keyCost", keyCost);

debugInfo.storeInfo("keyRate", keyRate);
debugInfo.storeInfo("keyLength", keyLength);
debugInfo.storeInfo("relEnt",relEnt);
debugInfo.storeInfo("relEntStep2", debugMathSolver.info.relEntStep2Linearization);


end








%%%%%%%%%%%  FINITE FADING CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [keyRate, keyLength, debugInfo] = finiteKeyRate(relEnt, deltaLeak, params)
%computes the finite size keyrate.

d = params.alphabet; %size of key register.
ptest = params.ptest;
pgen = 1 - ptest;
logepsilon = params.logepsilon;
N = params.N;
n=pgen*N;

kappa = sqrt(-logepsilon.PA)/log2(d+1);
alpha =  1+kappa/sqrt(n);




keyCost = params.keyCost / N;




CorrectionTerm = ( kappa * sqrt(n)/N )*(log2(d+1))^2; % n_k(alpha-1) log^2(d_Z+1)

privacyAmplification = (alpha/kappa)*(-2-logepsilon.PA+2/alpha)* sqrt(n) / N; % alpha / (\alpha-1) * (log(1/4epsPA) + 2/alpha)
ECLeakage = n/N*deltaLeak + ceil(-logepsilon.EC)/N; 

keyRate = pgen*max(relEnt,0) - CorrectionTerm - ECLeakage - privacyAmplification - keyCost; 
keyLength = keyRate*N;

end



















function observablesAndDimensionsMustBeTheSame(observables,dimA,dimB)
if ~allCells(observables,@(x) size(x,1) == dimA*dimB)
    throwAsCaller(MException("BasicKeyRateFunc:ObservablesAndDimensionsMustBeTheSame","The Observables must have the same dimensions as Alice and Bob multiplied together."));
end
end

function mustBeSizedLikeAnnouncements(jointExpectations,announcementsA,announcementsB)
if ~isequal(size(jointExpectations),[numel(announcementsA),numel(announcementsB)])
    throwAsCaller(MException("BasicKeyRateFunc:jointKeyDoesNotHaveSameSizeAsAnnouncements",...
        "The joint key distribution must have size numel(announcementsA) by numel(announcementsB)."))
end
end

function mustBeProbDistCell(input)
mustBeProbDist([input{:}])
end


function eachRowMustBeAProbDist(expectationsConditional)

% get the dimensions of the conditional expectations. Then based on that
% pick a strategy to handle it
dimExpCon = size(expectationsConditional);

errorID ="BasicBB84_WCPKeyRateFunc:InvalidRowsAreNotProbDists";
errorTXT = "A row in the conditional distribution is not a valid probability distribution.";

if numel(dimExpCon) == 2 % Matlab's minimum number of dimensions is 2.
    % The array is 2d and the slicing is easy
    for index = 1:dimExpCon(1)
        if~isProbDist(expectationsConditional(index,:))
           throwAsCaller(MException(errorID,errorTXT));
        end
    end
else
    % We have some tricky slicing to do for 3 plus dimensions.
    % We need to index the first dimension and the combination of
    % dimensions 3 and up. The second dimension will just use :.
    maskedDims = [dimExpCon(1),prod(dimExpCon(3:end))];

    for index = 1:prod(maskedDims)
        vecIndex = ind2subPlus(maskedDims,index);
        if ~isProbDist(expectationsConditional(vecIndex(1),:,vecIndex(2)))
            throwAsCaller(MException(errorID,errorTXT));
        end
    end
end
end



