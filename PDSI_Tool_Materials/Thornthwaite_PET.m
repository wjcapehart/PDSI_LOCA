%% © 2012 Vanderbilt University %%

function [PET] = Thornthwaite_PET(T_F,lat_list,count_loc,beg_row_col, ...
          lat_col)

% The Thornthwaite_PET function calculates the potential evapotranspiration 
% using Thornthwaite's method (cf. (1) Thornthwaite, Wilm, et al., 1944; 
% Transactions of the American Geophysical Union, Vol. 25, pp. 683-693; 
% (2) Thornthwaite, 1948; Geographical Review, Vol. 38, No. 1, January
% 1948; (3) Thornthwaite and Mather, 1955; Publications in Climatology,
% Vol. 8, No. 1; and (4) Thornthwaite and Mather, 1957; Publications in 
% Climatology, Vol. 10, No. 3)

% NOTE:
% T_F is the temperature in degrees F. This data should be input as a
% column vector.
% lat_list is the list of latitudes in the input data. There should be one
% latitude for each location, and this informatio should be input as a
% column vector.
% count_loc is the number of different locations in the input data.
% beg_row_col is a vector of the row numbers where the temperature data for
% each of the different locations begins.
% lat_col is the record of latitudes in degrees for all locations, arranged
% in a column vector, where a latitude is associated with each observation
% (i.e., month) of the record.

%% ASSIGN VARIABLES

% days_mo is the number of days in each of the 12 months.
days_mo = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

% J is the Julian day of the middle of month i. Pre-allocate J for speed.
J = zeros(1,12);

for i = 1:12
    if i == 1
        J(i) = days_mo(i)/2;
    else
        J(i) = days_mo(i)/2 + sum(days_mo(1,1:i-1));
    end
end

%% CONVERT TEMPERATURE FROM DEGREES FAHRENHEIT (F) TO DEGREES CELSIUS (C)

% Haith and Shoemaker (cf. Haith and Shoemaker, 1987; Journal of the 
% American Water Resources Association; Vol. 23, No. 3, June 1987) set PET
% equal to zero on days for which the mean monthly temperature is less than
% or equal to zero. T_C is the temperature record in degrees Celsius where
% those temperatures that were less than zero are now zero, such that the
% subsequent calculation of the potential evapotranspiration (PET) renders
% a PET equal to zero.
T_C = (T_F - 32).*(5/9);
T_C(T_C < 0) = 0;

%% CALCULATE THE POTENTIAL EVAPOTRANSPIRATION FOR EACH LOCATION j

PET = [];

for j = 1:count_loc
    
    % lat_rad is the latitudes in radians for location j.
    lat_rad = lat_list(j)*2*pi/360;
    
    % TLAT is the negative tangent of the latitiude (where the latitude is 
    % in radians).
    TLAT = -tan(lat_rad);
    
    g = beg_row_col(j); % g is the row number of the temperature and
                        % precipitation data records (i.e., column vectors)
                        % where the data for location j begins.
    if j == count_loc
        d = length(lat_col); % d is the row number of the temperature and
                             % precipitation data records (i.e., column 
                             % vectors) where the data for location j ends.
    else
        d = beg_row_col(j + 1) - 1;
    end
    
    y = (d - g + 1)/12; % y is the number of years in the data record for 
                        % location j.
                  
    T_loc = T_C(g:d,1); % T_loc is the temperature record (in degrees C) 
                        % for location j, arranged as a row vector.
                        
    T = reshape((T_loc'),12,y)'; % T is the temperature (for location j) 
                                 % converted from degrees F to degrees C 
                                 % and reshaped so rows represent years 
                                 % and columns represent months.
    
    %% CALCULATE VARIABLES USED IN THE PET CALCULATION
    
    % These equations are different than those presented in 
    % the NCDC and CPC PDSI fortran codes. Both the NCDC and CPC use 
    % transformed equations. This code follows Allen et al. (1994a), Allen 
    % et al. (1994b), and Duffle and Beckman (1991).
        % (1) Allen et al., 1994a; "An Update for the Calculation of 
        % Reference Evapotranspiration," ICID Bulletin of the International 
        % Commission on Irrigation and Drainage; 
        % (2) Allen et al., 1994b; "An Update for the Definition of 
        % Reference Evapotranspiration," ICID Bulletin of the International 
        % Commission on Irrigation and Drainage; and
        % (3) Duffie and Beckman, 1974; "Solar  engineering thermal
        % processes," Wiley and Sons, NY. 
        
        
    % Calculate PHI, the solar declination angle on Julian day (J).
    PHI = .4093.*sin(((2.*(pi())./365).*J)-1.405);
    % Calculate w (omega), the sunset hour angle on Julian day (J).
    w = acos(TLAT.*PHI);
    % Calculate D, the number of daylight hours (day length) on Julian  
    % day (J).
    D = 24.*w./pi;
    % Calculate DL, the day length normalizer.
    DL = D./12;
    % Calculate DM, the month length normalizer.
    DM = days_mo./30;
    % Calculate h, the heat index for for month i of year k.
    h = ((T./5)).^1.514;

    % Calculate H, the yearly heat index for all years of the temperature 
    % record.
    H =  sum(h,2);
    
    % Calculate a, the PET exponent for year k.
    a = ((6.75*10^-7).*H.^3 - (7.71*10^-5).*H.^2 + (1.792*10^-2).*H + ...
        0.49239);

    %% CALCULATE PET FOR EACH MONTH OF EACH YEAR FOR LOCATION j 
    
    % These equations are different than those presented in 
    % the NCDC and CPC PDSI fortran codes. Both the NCDC and CPC use 
    % transformed equations (see Hobbins et al., 2008; Geophysical 
    % Research Letters, Vol. 23, L12403, doi:10.1029/2008GL033840; and  
    % Dai, 2011; Journal of Geophysical Research, Vol. 116, D12115, 
    % doi:10.1029/2010JD015541). 
    
    % These equations do not match the CPC and NCDC transformed equations,  
    % but instead follow Thornthwaite's untransformed equations.
    	% (1) Thornthwaite's initial proposal of the method: Wilm, et al., 
        % 1944; Transactions of the American Geophysical Union, Vol. 25, 
        % pp. 683-693; 
        % (2) Comprehensive methodology, PET calculations by temperature: 
        % Thornthwaite, 1948; Geographical Review, Vol. 38, No. 1, January 
        % 1948; 
        % (3) Modifications and instructions: Thornthwaite and Mather, 
        % 1955; Publications in Climatology, Vol. 8, No. 1; 
        % (4) Modifications and instructions: Thornthwaite and Mather, 
        % 1957; Publications in Climatology, Vol. 10, No. 3; and
        % (5) Detailed Methodology: Thornthwaite and Havens, April 1958; 
        % Monthly Weather Review. 
    % Thornthwaite's method for calculating the PET is used here (as 
    % opposed to the Penman-Monteith or Hamon methods). Nevertheless, PDSI
    % values should be insensitive to alternative methods (van der Schrier,
    % 2011). 
        % (1) G. van der Schrier et al., 2011; Journal of Geophysical 
        % Research, Vol. 116, D03106, doi:10.1029/2010JD015001
    % The PET is claculated here according to three temperature ranges: T
    % < 0C, 0C <= T < 26.5C, and T >= 26.5C (Thornthwaite, 1948; Willmott 
    % et al., 1985; Haith and Shoemaker, 1987).
        % (1) Willmott et al., 1985; Journal of Climatology, Vol 5, pp
        % 589-606; and
        % (2) Haith and Shoemaker, 1987; Journal of the 
        % American Water Resources Association; Vol. 23, No. 3, June 1987.
    
    % PET_mm is the potential evapotranspiration (in millimeters) for 
    % location j. Pre-allocate for speed.
    PET_mm = zeros(size(T));
    
    for k = 1:y % k is the counter for each year of the data on record for
                % each of the different locations.
        for i = 1:12 % i is the counter for each of the 12 months in a 
                     % year.
            % Calculate PET for Temperatures < 32 degrees F.
            if T(k,i) <= 0
                PET_mm(k,i) = 0;
            % Calculate PET for Temperatures <= 32 degrees F and < 80 
            % degress F.
            elseif T(k,i) > 0 && T(k,i) < 26.5
                PET_mm(k,i) = 16*DL(i)*DM(i)*(((10*T(k,i))/H(k,1))^a(k,1));
            % Calculate PET for Temperatures >= 80 degrees F.
            elseif T(k,i) >= 26.5
                PET_mm(k,i) = (-415.85 + 32.24*T(k,i) - 0.43* ...
                              (T(k,i)^2))*DL(i)*DM(i);
            end
        end
    end
    
    % Convert PET from millimeters to inches.
    PET_in = (PET_mm./25.4);
   % Catalogue PET for all locations.
    PET = [PET; PET_in];

                      
    %% NCDC PET Calculation %%
    %{
    %
    % If you would like to run the code using the NCDC equation for PDSI,
    % please comment out lines 81-192 and uncomment lines 195-260.
    %
    %
    load NCDC_soilcnst.txt; % Loads the text file that contains the NCDC
                            % constants used to calculate PET. 
    B = NCDC_soilcnst(:,3); % B is the NCDC equivalent of a, Thornthwaite's
                            % PET exponent; instead of a yearly exponent, 
                            % only one exponent is used for all years.
    h = NCDC_soilcnst(:,4); % Instead of a yearly heat index, only one heat
                            % index is used for all years.
    TLAT = NCDC_soilcnst(:,5); % TLAT is the negative tangent of the 
                               % latitude (where the latitude is in 
                               % radians).
    % For efficiency, move lines 201-208 above line 59.
    
    PHI = [-0.3865982; -0.2316132; -0.0378180; 0.1715539; 0.3458803; ...
           0.4308320; 0.3916645; 0.2452467; 0.0535511; -0.15583436; ...
           -0.3340551; -0.4310691];
    
    T_loc = T_F(g:d,1); % T_loc is the temperature record (in degrees C) 
                        % for location j, arranged as a row vector.
                        
    T = reshape((T_loc'),12,y)'; % T is the temperature (for location j) 
                                 % converted from degrees F to degrees C 
                                 % and reshaped so rows represent years 
                                 % and columns represent months.
    % PET_in is the potential evapotranspiration (in inches) for location
    % j. Pre-allocate for speed.                             
    PET_in = zeros(size(T));
    
    DUM = PHI.*TLAT(j);
    DK = atan(sqrt(1 - DUM.*DUM)./DUM);
    for i = 1:12
        if DK(i) < 0
            DK(i) = 3.141593 + DK(i);
        end
        DK(i) = (DK(i) + 0.0157)/1.57;
    end
    for k = 1:y % k is the counter for each year of the data on record for
                % each of the different locations.
        for i = 1:12 % i is the counter for each of the 12 months in a
                     % year.
            % Calculate PET for Temperatures <= 32 degrees F.
            if T(k,i) <= 32
                PET_in(k,i) = 0.0;
            % Calculate PET for Temperatures >= 80 degrees F.
            elseif T(k,i) >= 80
                PET_in(k,i) = (sin(T(k,i)/57.3 - 0.166) - 0.76)*DK(i);
                PET_in(k,i) = PET_in(k,i)*days_mo(i);
            % Calculate PET for Temperatures <= 32 degrees F and < 80 
            % degress F.
            else
                DUM = log(T(k,i) - 32);
                PET_in(k,i) = (exp(-3.863233 + B(j)*1.715598 - B(j)* ...
                               log(h(j)) + B(j)*DUM))*DK(i);
                PET_in(k,i) = PET_in(k,i)*days_mo(i);
            end
        end
    end
    
    % Catalogue PET for all locations.
    PET = [PET; PET_in];
    %}     
end

end