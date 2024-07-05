function [newParams,modParser] = ThreeStateDescriptionFunc(params, options, debugInfo)


% ThreeStateDescriptionFunc A description function for the threestate
% protocol with the Flag state squashed. 

% Input parameters:
% * pz: The probability that Alice sends signals in the Z-basis.
% * photonNumberCutOff: Bob's photon number cutoff
% * basisChoiceBeamSplittingRatio : beam-splitting ratio towards the monitoring line (X-basis)
% 
% Output parameters:
% * observablesJoint: The joint observables for Alice and Bob's measurement
%   of the signals.
% * dimA: dimension of Alice's system.
% * dimB: dimension of Bob's system.
% * probSignalsA: Probability of Alice selecting a signal to send to Bob.
%   In this protocol, it is half the probability of the basis choice.
% * rhoA: Alice's reduced density matrix for prepare-and-measure based
%   protocols.
% * POVMA: Alice's set of POVM operators which she measures her state with
%   in the source replacement scheme.
% * POVMB: Bob's set of POVM operators which he measures his state with.
% * announcementsA: Alice's announcements for each of her POVM operators.
%   Can be integers or strings.
% * announcementsB: Bob's announcements for each of his POVM operators.
%   Can be integers or strings.
% * keyMap: An array of KeyMap objects that contain pairs of accepted
%   announcements and an array dictating the mapping of Alice's measurement
%   outcome to key bits (May be written with Strings).
% * krausOps: A cell array of matrices. The Kraus operators that form the G
%   map on Alice and Bob's joint system. These should form a completely
%   postive trace non-increasing linear map. Each Kraus operator must be
%   the same size.
% * keyProj:  A cell array of projection operators that extract the key
%   from G(\rho). These projection operators should sum to identity. This
%   map is often called Z.
% Options:
% * none
% DebugInfo:
% * krausSum: sum_i K^\dagger_i*K_i which should be <= identity for
%   a CPTNI map.
%
% Reviewed by Devashish Tupkary 2023/09/18
% See also QKDDescriptionModule, makeGlobalOptionsParser
arguments
    params (1,1) struct
    options (1,1) struct
    debugInfo (1,1) DebugInfo
end

%% options parser
%Parsing technical options for the module
optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.parse(options);
options = optionsParser.Results;

%% module parser
%Parsing parameters for the module
modParser = moduleParser(mfilename);
modParser.addRequiredParam("pzAlice",@(x) mustBeInRange(x,0,1));

modParser.addRequiredParam("basisChoiceBeamSplittingRatio",@(x) mustBeInRange(x,0,1));


modParser.addRequiredParam("photonNumberCutOff" , @mustBeInteger);
modParser.addAdditionalConstraint(@mustBePositive,"photonNumberCutOff");



modParser.parse(params)
params = modParser.Results;

%% simple setup
%Setup parameters for later use
newParams = struct();

%z-basis choice
pzAlice = params.pzAlice;

% Alice sends one of three possible pure states. 
dimA = 3;

N = params.photonNumberCutOff;
nondiagdimB = (N+1)*(N+2)/2;  % <=BobCutOff-photon dimension
diagdimB = 8;   % flag-state dimension (number of coarse-grained 'click' events) : 5 single clicks, no click, cross click, and rest.
dimB = nondiagdimB + diagdimB; % Bob's total dimension

newParams.dimA = dimA;
newParams.dimB = dimB;

%% generate rhoA


signalStates = {zket(2,1),zket(2,2),1/sqrt(2)*(zket(2,1)+zket(2,2))};
probSignalsA = [pzAlice/2,pzAlice/2,(1-pzAlice)]; %probability of Alice sending each signal.
newParams.probSignalsA = probSignalsA;

rhoA = zeros(dimA,dimA);

for i = 1:numel(signalStates)
    for j = 1:numel(signalStates)
        overlap = signalStates{i}'*signalStates{j};
        rhoA = rhoA + zket(dimA,i)*(zket(dimA,j)')*overlap*sqrt(probSignalsA(i)*probSignalsA(j));
    end
end



newParams.rhoA = rhoA;

%% joint obserables
POVMsA = {diag(zket(dimA,1)),diag(zket(dimA,2)),diag(zket(dimA,3))};



% This generates the POVMs for the 3-state protocol with 3 detectors (which
% contains one additional detector than the protocol we consider.) These
% POVMs are squashed to Bob's smaller space but without flags. 
BPOVMs = COWPOVM(N, 1, params.basisChoiceBeamSplittingRatio); 

% BPOVM is a cell array POVM, where the second index is equal to
% the total number of detectors clicking + 1. The first index tells you
% the specific click pattern when that many detectors click.

% This maps the clicks in the '+' detector to the no-click event.
% Thus, it constructing the POVMs for the 3-state protocol with a
% missing '+' detector from the POVMs of the 3-state protocol that contains
% the '+' detector. 
BPOVM = CoarseGrainedPOVMs(BPOVMs); 

% This coarse-grains all multi-click events to cross click and no cross
% click. We also add the flags to the POVMs.
t = params.basisChoiceBeamSplittingRatio;
cutOffFactor = 1 - t^(N+1)-(1-t/4)^(N+1) + (3*t/4)^(N+1);
[POVMsB, ~] = MultiCoarseGrainedPOVMs(BPOVM, cutOffFactor);


newParams.POVMA = POVMsA; %we don't use this directly, but good to have.
newParams.POVMB = POVMsB;



% each POVM element is assigned an announcement made by their respective
% party
newParams.announcementsA = ["Z","Z","X"];
% we throw away vacuum, multi clicks, and all clicks in the '-' detector.
newParams.announcementsB = ["throw","throw","throw","throw","keep","keep","throw","throw"]; 

% For each pair of announcements Alice and Bob keep, we assign Alice's POVM
% elements a corresponding key dit value. We only produce a key when Alice
% uses the Z basis. Therefore the third value does not matter, since that corresponds to Alice's X basis
% outcome. 
newParams.keyMap = [KeyMapElement("Z","keep",[1,2,1])]; 

% add all joint observables
observablesJoint = cell(numel(POVMsA),numel(POVMsB));
for indexA = 1:numel(POVMsA)
    for indexB = 1:numel(POVMsB)
        observablesJoint{indexA,indexB} = kron(POVMsA{indexA},POVMsB{indexB});
    end
end

newParams.observablesJoint = observablesJoint;
debugInfo.storeInfo("observablesJoint", observablesJoint);


%% Kraus Ops (for G map)
% A: Alice's system, B: Bob's System, C: announcement register, R:
% key register.
% The Kraus operators are matrices from ABC \rightarrow RBC. We use results from 
% https://arxiv.org/abs/1905.10896 (Appendix A) to shrink the Kraus operators from outputing RABC
% to just RBC. 

%Bob's kraus operator for "keep".
KraussBob = sqrt(POVMsB{5}+POVMsB{6}); 
KraussAlice = zket(2,1)*zket(3,1)' + zket(2,2)*zket(3,2)';
krausOps = {kron(KraussAlice, KraussBob)};

krausSum = 0;
for index = 1:numel(krausOps)
    krausSum = krausOps{index}'*krausOps{index};
end
debugInfo.storeInfo("krausSum", krausSum);

%% Pinching map 
proj0 = kron(diag([1,0]),eye(dimB));
proj1 = kron(diag([0,1]),eye(dimB));
keyProj = {proj0,proj1};

%% set key map, kraus ops, and block diagonal structure in new parameters
newParams.krausOps = krausOps;
newParams.keyProj = keyProj;
newParams.blockDimsA = [dimA];
newParams.blockDimsB = [1,nondiagdimB-1,ones(1,diagdimB)]; % Vacuum, non-vacuum (this can be further broken up if BobCutoff > 1), flags



end

