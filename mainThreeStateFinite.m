%pick the preset file to use 

clear;

% we plot 6 kinds of key rates..
liftType = {"IID","Sequential-IID","Postselection-generic",...
      "Postselection-blockdiagonal","Postselection-decoy","Postselection-no-improvement"};


for i=1:numel(liftType)
    qkdInput = FiniteThreeStatePhaseCoherentPreset(liftType{i});
    %run the QKDSolver with this input
    results = MainIteration(qkdInput);
    %save the results and preset to a file.
    filename = "FiniteThreeState"+liftType{i}+".mat";
    save(filename,"results","qkdInput");
end



%% plot the result
PSPlots(liftType);
