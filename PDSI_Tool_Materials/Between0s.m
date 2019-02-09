%% © 2012 Vanderbilt University %%

function [PV,PX1,PX2,PX3,PPe,X,BT] = Between0s(k,Z,X3,PX1,PX2,PX3,PPe,BT,X)

% This function is called when non-zero, non-one hundred PPe values occur
% between values of PPe = 0. When this happens, a possible abatement
% discontinues without ending the wet spell or drought. X should be
% assigned to PX3 for all months between, and including, the two instances
% of PPe = 0 (cf. Alley, 1984; Journal of Climate and Applied Meteorology, 
% Vol. 23, No. 7). To do this, backtrack up to the first instance of 
% PPe = 0 while setting X to PX3. 

% Since the possible abatement has ended, the drought or wet spell
% continues. Set PV, PX1, PX2, and PPe to 0. Calculate PX3 and set X = PX3.
% Set BT=3 in preparation for backtracking.
PV = 0;                         
PX1(k) = 0;                     
PX2(k) = 0;                     
PPe(k) = 0;                     
PX3(k) = 0.897 * X3 + (Z(k) / 3);
X(k) = PX3(k);
BT(k) = 3;

% In order to set all values of X between the two instances of PPe = 0, the
% first instance of PPe = 0 must be found. This "for" loop counts back 
% through previous PPe values to find the first instance where PPe = 0.
for count1 = k : -1 : 1
    if PPe(count1) == 0
        r = count1;
        break
    else
        continue
    end
end

% Backtrack from the current month where PPe = 0 to the last month where
% PPe = 0.
for count = k : -1 : r
    % Set X = PX3
    if BT(count) == 3
        X(count) = PX3(count);
        % If the end of the backtracking loop hasn't been reached, set the
        % BT value for the preceeding month to 3 to continue the
        % backtracking.
        if count ~= r
            BT(count-1) = 3;
        else
        end
        continue
    end
end

end