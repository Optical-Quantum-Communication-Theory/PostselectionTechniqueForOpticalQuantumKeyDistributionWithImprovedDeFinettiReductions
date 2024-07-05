function [newParams, modParser]= ThreeStateFlagStateChannelFunc(params,options,debugInfo)

% ThreeStateFlagStateChannelFunc : A channel func for the three state
% protocol that computes the joint expectation values of Alice and Bob's
% measurement outcomes.

% Input parameters:
% * pzAlice: The probability that Alice sends signals in the Z-basis.
% * eta (1): the transmissivity of the quantum channel; Must be between 0 
%   and 1 inclusive.
% * photonNumberCutOff: Bob's photon number cutoff
% * basisChoiceBeamSplittingRatio : beam-splitting ratio towards the monitoring line (X-basis)
% Output parameters:
% * ADD
% Options:
% * None.
% DebugInfo:
% * None.
%
% 
% See also QKDChannelModule,  makeGlobalOptionsParser
arguments
    params (1,1) struct
    options (1,1) struct
    debugInfo (1,1) DebugInfo
end

%% options parser
optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.parse(options);
options = optionsParser.Results;

%% module parser
modParser = moduleParser(mfilename);

%Decoy intensities

%Z-basis choice
modParser.addRequiredParam("pzAlice", @(x) mustBeInRange(x, 0, 1));
modParser.addAdditionalConstraint(@isscalar,"pzAlice");

%Channel loss
modParser.addOptionalParam("eta", 1, @(x) mustBeInRange(x, 0, 1));
modParser.addAdditionalConstraint(@isscalar,"eta");

modParser.addRequiredParam("basisChoiceBeamSplittingRatio",@(x) mustBeInRange(x,0,1));

modParser.addRequiredParam("photonNumberCutOff",@mustBeInteger);


modParser.addRequiredParam("observablesJoint",@(x) allCells(x,@ishermitian));
modParser.addAdditionalConstraint(@mustBePositive,"photonNumberCutOff");

modParser.parse(params);

params = modParser.Results;


pzAlice = params.pzAlice;
t = params.basisChoiceBeamSplittingRatio;
eta = params.eta;
observablesJoint = params.observablesJoint;
distance = -10*log10(eta)/0.16*10^3; % Assuming ULL fibre with 0.16 dB/km loss. Unit is metres


%% Calculating the weight outside the subspace

%compute simulated rho for honest behaviour.

probSignalsA = [pzAlice/2,pzAlice/2,(1-pzAlice)];
ketPlus = (zket(2,1)+zket(2,2))/sqrt(2);

rhoAprime{1} = eta*blkdiag(0,diag(zket(2,1)))+(1-eta)*blkdiag(1,zeros(2)); %vacuum + qubit
rhoAprime{2} = eta*blkdiag(0,diag(zket(2,2)))+(1-eta)*blkdiag(1,zeros(2));
rhoAprime{3} = eta*blkdiag(0,ketPlus*ketPlus')+(1-eta)*blkdiag(1,zeros(2));

rhoAB = zeros(33);
for index = 1:numel(probSignalsA)
    rhoAB = rhoAB + probSignalsA(index)*kron(diag(zket(numel(probSignalsA),index)),...
        blkdiag(rhoAprime{index},zeros(8)));
end



expectationsJoint = zeros(size(params.observablesJoint));

for index = 1:numel(params.observablesJoint)
    expectationsJoint(index) = real(trace(params.observablesJoint{index}*rhoAB));
end

newParams.expectationsJoint = expectationsJoint;


newParams.distance = distance;

end