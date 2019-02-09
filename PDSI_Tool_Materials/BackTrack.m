%% © 2012 Vanderbilt University %%

function [X,BT] = BackTrack(k,PPe,PX1,PX2,PX3,X,BT)

% This function backtracks through previous PX1 and PX2 values.
% Backtracking occurs in two instances: (1) After the probability reaches 
% 100 and (2) When the probability is zero. In either case, the
% backtracking function works by backtracking through PX1 and PX2 until
% reaching a month where PPe = 0. Either PX1 or PX2 is assigned to X as the
% backtracking progresses.

% Backtracking occurs from either PPe(k) = 100 or PPe(k) = 0 to the first 
% instance in the previous record where PPe = 0. This "for" loop counts 
% back through previous PPe values to find the first instance 
% where PPe = 0.
% r is a variable used to mark the place of the last PPe = 0 before 
% PPe = 100.
for c = k : -1 : 1
    if PPe(c) == 0
        r = c;
        break
    else
        continue
    end
end

% Backtrack from either PPe = 100 or PPe = 0 to the last instance of 
% non-zero, non-one hundred probability.
for count = k : -1 : r + 1
    % When PPe(k) = 100 and |PX3| > 1 set X(k) = PX3(k).
    %                                                                       
    % Set the BT value of the previous month to either 1 or 2 based on the
    % sign of PX3(k). If PX3(k) is negative, a BT = 2 begins backtracking 
    % up X2 and vice versa.
    if PPe(count) == 100 && abs(PX3(count)) > 1
        X(count) = PX3(count);
        if PX3(count) < 0
            BT(count-1) = 2;
        else
            BT(count-1) = 1;
        end
        continue
    end
    
    % Everything below deals with months where PPe is not equal to 100. 
    % Based on the assigned BT value, start in either PX1 or PX2. If
    % that value is not 0, assign X and set the BT value for the preceeding
    % month to 1 if X = PX1 or 2 if X = PX2. If BT = 1 and PX1 = 0, assign 
    % X to PX2 and set the BT value for the preceeding month to 2 and vice
    % versa. Continue this process of backtracking up either PX1 or PX2
    % and switching when either PX1 or PX2 equals 0 or until the end of the
    % loop is reached.
    if BT(count) == 2
        if PX2(count) == 0
            X(count) = PX1(count);
            BT(count) = 1;
            BT(count-1) = 1;
            continue
        else
            X(count) = PX2(count);
            BT(count-1) = 2;
            continue
        end
    elseif BT(count) == 1
        if PX1(count) == 0
            X(count) = PX2(count);
            BT(count) = 2;
            BT(count-1) = 2;
            continue
        else
            X(count) = PX1(count);
            BT(count-1) = 1;
            continue
        end
    end
    
end

end