%% © 2012 Vanderbilt University %%

function [table] = PDSI_Central(Z_all,count_loc,beg_row_col,end_row_col, ...
                   lat_col,yr_col,PET)

% NOTE:
% beg_row_col is a vector of the row numbers where the temperature data for
% each of the different locations begins.
% count_loc is the number of different locations in the input data.
% lat_col is the record of latitudes in degrees for all locations,
% arranged in a column vector, where a latitude is associated with each
% observation (i.e., month) of the record.

PET_col = reshape((PET'),length(PET)*12,1);
table = [];

for j = 1:count_loc
    
    Z = Z_all(beg_row_col(j):end_row_col(j),1);
    
    %% INITIALIZE PDSI AND PHDI CALCULATIONS
    
    V = 0;  % V is the sum of the Uw (Ud) values for the current and
            % previous months of an established dry (wet) spell and is
            % used in calculating the Pe value for a month.
    Pe = 0; % Pe is the probability that the current wet or dry spell
            % has ended in a month.
    X1 = 0; % X1 is the severity index value for an incipient wet spell
            % for a month.
    X2 = 0; % X2 is the severity index value for an incipient dry spell
            % for a month.
    X3 = 0; % X3 is the severity index value of the current established
            % wet or dry spell for a month.
    Q = 0;
    
    l(j) = (end_row_col(j) - beg_row_col(j) + 1);
    
    % BACTRACKING VARIABLES
    BT = zeros(l(j),1); % BT is the backtracking variable, and is 
                        % pre-allocated with zeros. Its value (1, 2, or 3)
                        % indicates which intermediate index (X1, X2, or
                        % X3) to backtrack up, selecting the associated 
                        % term (X1, X2, or X3) for the PDSI. NOTE: BT may
                        % be operationally left equal to 0, as it cannot
                        % be known in real time when an existing drought
                        % or wet spell may or may not be over.
    
    S = 10^4; % S is the rounding number, and it is used to round the final
              % calculations of the output variables to four decimal 
              % places.
    
    %% CALCULATE PDSI AND PHDI
    for k = 1:l(j)
        PX1(k) = 0; 
        PX2(k) = 0; 
        PX3(k) = 0; 
        PPe(k) = 0; 
        X(k) = 0; 
        Ze(k) = 0; % Ze is the soil moisture anomaly (Z) value that will 
                   % end the current established dry or wet spell in that 
                   % month and is used in calculating the Q value and 
                   % subsequently the Pe value for a month.
        Uw(k) = 0; % Uw is the effective wetness required in a month to end 
                   % the current established dry spell (drought).
        Ud(k) = 0; % Ud is the effective dryness required in a month to end 
                   % the current wet spell.
        
        if Pe == 100 || Pe == 0 % No abatement underway
            if abs(X3) <= 0.5 % Drought or wet spell ends
                PV = 0; % PV is the preliminary V value and is used in  
                        % operational calculations.
                PPe(k) = 0; % PPe is the preliminary Pe value and is used
                            % in operational calculations.
                PX3(k) = 0; % PX3 is the preliminary X3 value and is used
                            % in operational calculations.                
                [PX1,PX2,PX3,X,BT] = Main(Z,k,PV,PPe,X1,X2,PX1,PX2,PX3, ...
                                     X,BT);                                                          
            elseif X3 > 0.5 % Wet spell underway
                if Z(k) >= 0.15 % Wet spell intensifies
                    [PV,PX1,PX2,PX3,PPe,X,BT] = Between0s(k,Z,X3,PX1, ...
                                                PX2,PX3,PPe,BT,X);                                                  
                else % Wet spell starts to abate, and it may end.
                    [Ud,Q,PV,PPe,PX1,PX2,PX3,X,BT] = Function_Ud(k,Z,V, ...
                                                     Pe,PPe,PX1,PX2, ...
                                                     PX3,X1,X2,X3,X,BT);
                                                                                               
                end
            elseif X3 < -0.5 % Drought underway
                if Z(k) <= -0.15 % Drought intensifies 
                    [PV,PX1,PX2,PX3,PPe,X,BT] = Between0s(k,Z,X3,PX1, ...
                                                PX2,PX3,PPe,BT,X);                                                 
                else % Drought starts to abate, and it may end.
                    [Uw,Q,PV,PPe,PX1,PX2,PX3,X,BT] = Function_Uw(k,Z,V, ...
                                                     Pe,PPe,PX1,PX2, ...
                                                     PX3,X1,X2,X3,X,BT);
                                                                                          
                end
            end
        else % Abatement underway
            if X3 > 0 % Wet spell underway
                [Ud,Q,PV,PPe,PX1,PX2,PX3,X,BT] = Function_Ud(k,Z,V,Pe, ...
                                                 PPe,PX1,PX2,PX3,X1,X2, ...
                                                 X3,X,BT);                                                     
            else % Drought underway
                [Uw,Q,PV,PPe,PX1,PX2,PX3,X,BT] = Function_Uw(k,Z,V,Pe, ...
                                                 PPe,PX1,PX2,PX3,X1,X2, ...
                                                 X3,X,BT);
                                                                                            
            end
        end
        
        %% Assign V, Pe, X1, X2, and X3 for next month (k + 1)
        V = PV;
        Pe = PPe(k);
        X1 = PX1(k);
        X2 = PX2(k);
        X3 = PX3(k);
        
        %% ASSIGN X FOR CASES WHERE PX3 AND BT EQUAL ZERO
        % NOTE: This is a conflicting case that arises where X cannot be
        % assigned as X1, X2, or X3 in real time. Here 0 < PX1 < 1, 
        % -1 < PX2 < 0, and PX3 = 0, and it is not obvious which
        % intermediate index should be assigned to X. Therefore,
        % backtracking is used here, where BT is set equal to the next
        % month's BT value and X is assigned to the intermediate index
        % associated with that BT value.
        if k > 1
            if PX3(k-1) == 0 && BT(k-1) == 0
                for c = k -1 : -1 : 1
                    if BT(c) ~= 0 % Backtracking continues in a 
                                  % backstepping procedure up through the
                                  % first month where BT is not equal to 
                                  % zero.
                        r = c + 1; % r is the row number up through which
                                   % backtracking continues.
                        break
                    else
                        continue
                    end
                end
                for count0 = k - 1 : -1 : r
                    BT(count0) = BT(count0 + 1); % Assign BT to next 
                                                 % month's BT value.
                    if BT(count0) == 2
                        if PX2(count0) == 0 % If BT = 2, X = PX2 unless 
                                            % PX2 = 0, then X = PX1.
                            X(count0) = PX1(count0);
                            BT(count0) = 1;
                            continue 
                        else 
                            X(count0) = PX2(count0);
                            BT(count0) = 2;
                            continue
                        end 
                    elseif BT(count0) == 1
                        if PX1(count0) == 0 % If BT = 1, X = PX1 unless
                                            % PX1 = 0, then X = PX2.
                            X(count0) = PX2(count0); 
                            BT(count0) = 2;
                            continue 
                        else 
                            X(count0) = PX1(count0);
                            BT(count0) = 1;
                            continue
                        end 
                    else
                    end
                end
            end
        end
        
        % In instances where there is no established spell for the last 
        % monthly observation, X is initially assigned to 0. The code below 
        % sets X in the last month to greater of |PX1| or |PX2|. This 
        % prevents the PHDI from being inappropriately set to 0. 
        if k == l(j);
            if PX3(k) == 0 && X(k) == 0
                if abs(PX1(k)) > abs(PX2(k))
                    X(k) = PX1(k);
                else
                    X(k) = PX2(k);
                end
            end
        end
                
        
        % ROUND X1, X2, X3, Pe, V, PV, Q, X, PX1, PX2, PX3, PPe, Ud, Uw,
        % AND Ze TO FOUR DECIMAL PLACES.
        X1 = round(X1*S)/S;
        X2 = round(X2*S)/S;
        X3 = round(X3*S)/S;
        Pe = round(Pe*S)/S;
        V = round(V*S)/S;
        PV = round(PV*S)/S;
        Q = round(Q*S)/S;
        X(k) = round(X(k)*S)/S;
        PX1(k) = round(PX1(k)*S)/S;
        PX2(k) = round(PX2(k)*S)/S;
        PX3(k) = round(PX3(k)*S)/S;
        PPe(k) = round(PPe(k)*S)/S;
        Ud(k) = round(Ud(k)*S)/S;
        Uw(k) = round(Uw(k)*S)/S;
        Ze(k) = round(Ze(k)*S)/S;
        
    end
    
    %% ASSIGN PDSI VALUES
    % NOTE: 
    % In Palmer's effort to create a meteorological drough index (PDSI),
    % Palmer expressed the beginning and ending of dry (or wet) periods in
    % terms of the probability that the spell has started or ended (Pe). A
    % drought (wet spell) is definitely over when the probability reaches
    % or exceeds 100%, but the drought (wet spell) is considered to have
    % ended the first month when the probability becomes greater than 0%
    % and then continues to remain greater than 0% until it reaches 100% 
    % (cf. Palmer, 1965; US Weather Bureau Research Paper 45).
    PDSI = X';
    
    %% ASSIGN PHDI VALUES
    % NOTE:
    % There is a lag between the time that the drought-inducing
    % meteorological conditions end and the environment recovers from a
    % drought. Palmer made this distinction by computing a meteorological
    % drought index (described above) and a hydrological drought index. The
    % X3 term changes more slowly than the values of the incipient (X1 and
    % X2) terms. The X3 term is the index for the long-term hydrologic
    % moisture condition and is the PHDI.
    for s = 1:length(PX3)
        if PX3(s) == 0 % For calculation and program advancement purposes, 
                       % the PX3 term is sometimes set equal to 0. In such 
                       % instances, the PHDI is set equal to X (the PDSI), 
                       % which accurately reflects the X3 value.
            PHDI(s) = X(s);
        else
            PHDI(s) = PX3(s);
        end
    end
    
    lat_loc = lat_col(beg_row_col(j):end_row_col(j),1);
    yr_loc = yr_col(beg_row_col(j):end_row_col(j),1);
    PET_loc = PET_col(beg_row_col(j):end_row_col(j),1);
    Z_loc = Z_all(beg_row_col(j):end_row_col(j),1);
    
    % Create output table of variables for location j.
    table_loc = [lat_loc, yr_loc, PET_loc, Z_loc, PPe', PX1', PX2', ...
                  PX3', X', PDSI, PHDI'];
    % Vertically concatenate output tables for all locations.
    table = [table; table_loc];
    
    % Reset catalogued variables.
    X = [];
    PDSI = [];
    PHDI = [];
    PX1 = [];
    PX2 = [];
    PX3 = [];
    PPe = [];
    Ud = [];
    Uw = [];
    Ze = [];
    Z = [];
    
end

end

