%% Bob's coarse-grained POVM operators for the 3-State protocol with multi-click coarse graining
% We coarse grain the POVMs to save space such that all multi-click POVMs
% are coarse grained into one POVM
%
% Input:
% * POVM   : POVMs as outputted by the CoarseGrainedPOVMs function
% 
% Output:
% * CPOVM : Coarse grained POVMs
%
% * CPOVMnoflag : Coarse grained POVMS without the flags
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright: Shlok Nahar                     Last Updated: 27th August 2020


function [CPOVM, CPOVMnoflag] = MultiCoarseGrainedPOVMs(POVM, cutOffFactor)
    
    CPOVM = cell(1, 8); % The first 6 are the vacuum and single-click events.
                        % The 8th is the cross-click POVM event.
                        % The 7th one is all remaining click patterns.
        
    CPOVM{1,1} = POVM{1,1};         % No-click event
    n = length(POVM{1,1});
    sum = POVM{1,1}; %sum over all single slicks
    for i = 1:1:5
        CPOVM{1,i+1} = POVM{i,2};   % Single-click events
        CPOVM{1,i+1}(abs(CPOVM{1,i+1})<eps) = 0;
        sum = sum + POVM{i,2};
    end


    % All clicks in the data line detectors while ignoring the middle click
    % in the minus detector..
    data = POVM{4,2} + POVM{5,2} + POVM{2,2} + ...  % single clicks
           POVM{6,3} + POVM{7,3} + POVM{10,3} + ... % double clicks
           POVM{9,4};                             % triple clicks
    
    % all clicks in the minus detector. 
    minus = POVM{1,2} + POVM{3,2} + ... % single clicks without {2,2} to avoid double counting
            POVM{1,3} + POVM{2,3} + POVM{5,3} + ... % double clicks
            POVM{1,4};                            % triple clicks

    CPOVM{1,8} = eye(n) - data - minus - POVM{1,1}; % cross click (while ignoring middle click in the minus detector)
    sum = sum + CPOVM{1,8};
    CPOVM{1,8}(abs(CPOVM{1,8})<eps) = 0;

    CPOVM{1,7} = eye(n) - sum;      % Multi-click event
    CPOVM{1,7}(abs(CPOVM{1,7})<eps) = 0;

    CPOVMnoflag = CPOVM;
    
%% Adding the flag-state portion of the POVMs

    totalflag = 8;
    currentflag = 1;
%     CPOVM{1,1} = blkdiag(CPOVM{1,1},zeros(totalflag));
    for i = 1:1:totalflag-1
        flagstate = zket(totalflag, currentflag)*zket(totalflag, currentflag)';
        CPOVM{i} = (CPOVM{i}+CPOVM{i}')/2;
%         CPOVM{i+1} = blkdiag(CPOVM{i+1},flagstate);
        CPOVM{i} = blkdiag(CPOVM{i},(1-cutOffFactor)*flagstate);
        currentflag = currentflag + 1;
    end
    % Different flag space for cross-click
    i = i+1;
    flagstate = zket(totalflag, currentflag)*zket(totalflag, currentflag)';
    CPOVM{i} = (CPOVM{i}+CPOVM{i}')/2;
%         CPOVM{i+1} = blkdiag(CPOVM{i+1},flagstate);
    CPOVM{i} = blkdiag(CPOVM{i},cutOffFactor*eye(totalflag)+(1-cutOffFactor)*flagstate);
end
