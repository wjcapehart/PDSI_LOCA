%% © 2012 Vanderbilt University %%

function [Uw,Q,PV,PPe,PX1,PX2,PX3,X,BT] = Function_Uw(k,Z,V,Pe,PPe,PX1, ...
                                          PX2,PX3,X1,X2,X3,X,BT)

Uw(k) = Z(k) + 0.15; % In the case of an established drought, Palmer (1965) 
                     % notes that a value of Z = -0.15 will maintain an
                     % index of -0.50 from month to month. An established
                     % drought or wet spell is considered definitely over
                     % when the index reaches the "near normal" category
                     % which lies between -0.50 and +0.50. Therefore, any
                     % value of Z >= -0.15 will tend to end a drought.
PV = Uw(k) + max(V,0);
if PV <= 0 % During a drought, PV <= 0 implies PPe = 0 (i.e., the 
           % probability that the drought has ended returns to zero).                                                             
   Q = 0; 
   [PV,PX1,PX2,PX3,PPe,X,BT] = Between0s(k,Z,X3,PX1,PX2,PX3,PPe,BT,X);                                                                 
else
    Ze(k) = -2.691 * X3 - 1.5;
    if Pe == 100 
        Q = Ze(k); % Q is the total moisture anomaly required to end the 
                   % current drought.
    else                                                                    
        Q = Ze(k) + V;
    end
    PPe(k) = (PV / Q) * 100;
    if PPe(k) >= 100
        PPe(k) = 100;
        PX3(k) = 0;
    else
        PX3(k) = 0.897 * X3 + (Z(k) / 3);
    end
    [PX1,PX2,PX3,X,BT] = Main(Z,k,PV,PPe,X1,X2,PX1,PX2,PX3,X,BT);                                                                 
end

end