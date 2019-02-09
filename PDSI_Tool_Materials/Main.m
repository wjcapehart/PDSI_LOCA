%% © 2012 Vanderbilt University %%

function [PX1,PX2,PX3,X,BT] = Main(Z,k,PV,PPe,X1,X2,PX1,PX2,PX3,X,BT)

% This function calculates PX1 and PX2 and calls the backtracking loop. If 
% the absolute value of PX1 or PX2 goes over 1, that value becomes the new 
% PX3. 

% Calculate the current PX1 and PX2.
PX1(k) = max(0,0.897 * X1 + (Z(k) / 3));
PX2(k) = min(0,0.897 * X2 + (Z(k) / 3));

if PX1(k) >= 1 && PX3(k) == 0   
    % When PX1 >= 1 the wet spell becomes established. X is assigned as PX1
    % and PX3 = PX1. PX1 is set to zero after PX3 is set to PX1. BT is set
    % to 1 and the backtrack function is called to begin backtracking up
    % PX1.
    X(k) = PX1(k);              
    PX3(k) = PX1(k);            
    PX1(k) = 0;
    BT(k) = 1;
    [X,BT] = BackTrack(k,PPe,PX1,PX2,PX3,X,BT);                                                             
    return                                                                  
end

if PX2(k) <= -1 && PX3(k) == 0  
    % When PX2 <= -1 the drought becomes established. X is assigned as PX2
    % and PX3 = PX2. PX2 is set to zero after PX3 is set to PX2. BT is set
    % to 2 and the backtrack function is called to begin backtracking up
    % PX2.
    X(k) = PX2(k);              
    PX3(k) = PX2(k);            
    PX2(k) = 0;
    BT(k) = 2;
    [X,BT] = BackTrack(k,PPe,PX1,PX2,PX3,X,BT);                                                             
    return                                                                 
end

if PX3(k) == 0
    % When PX3 is zero and both |PX1| and |PX2| are less than 1, there is
    % no established drought or wet spell. X is set to whatever PX1 or PX2
    % value is not equal to zero. BT is set to either 1 or 2 depending on
    % which PX1 or PX2 value equals zero. The backtrack function is called
    % to begin backtracking up either PX1 or PX2 depending on the BT value.
    if PX1(k) == 0              
        X(k) = PX2(k);          
        BT(k) = 2;
        [X,BT] = BackTrack(k,PPe,PX1,PX2,PX3,X,BT);                                                                 
        return                                                              
    elseif PX2(k) == 0          
        X(k) = PX1(k);          
        BT(k) = 1;
        [X,BT] = BackTrack(k,PPe,PX1,PX2,PX3,X,BT);                                                          
        return                                                              
    end     
end

% There is no determined value to assign to X when PX3 ~= 0, 
% 0 <= PX1 < 1, and -1 < PX2 <= 0 so set X = PX3.   
X(k) = PX3(k);

end