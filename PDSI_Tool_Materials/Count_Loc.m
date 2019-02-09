%% © 2012 Vanderbilt University %%

function  [beg_row_col,end_row_col,beg_row_mat,end_row_mat,count_loc, ...
           lat_list,lat_mat] = Count_Loc(lat_col)

% NOTE: lat_col is the record of latitudes in degrees for all locations, 
% arranged in a column vector, where a latitude is associated with each
% observation (i.e., month) of the record for all locations.

% lat_mat_rec is the record of latitudes in degrees for all locations, 
% arranged in a matrix, where a latitude is associated with each MONTH AND 
% YEAR of the record for all locations.
lat_mat_rec = (reshape(lat_col,12,length(lat_col)/12))';

% lat_mat is the record of latitudes in degrees for all locations arranged
% in a column vector, where a latitude is associated with each YEAR of the
% record for all locations.
lat_mat = lat_mat_rec(:,1);

% count1 is a counter that tracks the rows in the total data record (i.e., 
% column vector) where the temperature and precipitation data for each of 
% the different locations begins.
count1 = 1;

% beg_row_col is a vector of the row numbers where the temperature data for
% each of the different locations begins.
beg_row_col = count1;

for m = 2:length(lat_col)
    if lat_col(m) == lat_col(m-1)
        count1(m) = count1(m-1) + 1;
    else
        count1(m) = count1(m-1) + 1;
        beg_row_col = [beg_row_col; count1(m)];
    end
end

for i = 1:length(beg_row_col)
    if i == length(beg_row_col)
        end_row_col(i) = length(lat_col); % end_row_col is a vector of the 
                                          % row numbers where the 
                                          % temperature data for each of 
                                          % the different locations ends.
    else
        end_row_col(i) = beg_row_col(i + 1) - 1;
    end
end

% count2 is a counter that tracks the rows in the reshaped data record  
% (i.e., matrix) where the temperature and precipitation data for each of 
% the different locations begins.
count2 = 1;

% beg_row_mat is a vector of the row numbers where the temperature data for
% each of the different locations begins when the data is listed in a  
% matrix such that years represent rows and columns represent months.
beg_row_mat = count2;

for m = 2:length(lat_mat)
    if lat_mat(m) == lat_mat(m-1)
        count2(m) = count2(m-1) + 1;
    else
        count2(m) = count2(m-1) + 1;
        beg_row_mat = [beg_row_mat; count2(m)];
    end
end

for i = 1:length(beg_row_mat)
    if i == length(beg_row_mat)
        end_row_mat(i) = length(lat_mat); % end_row_col is a vector of the
                                          % row numbers where the
                                          % temperature data for each of
                                          % the different locations ends
                                          % when the data is listed in a 
                                          % matrix such that years 
                                          % represent rows and columns
                                          % represent months.
    else
        end_row_mat(i) = beg_row_mat(i + 1) - 1;
    end
end

% count_loc is the number of different locations in the total data record.
count_loc = length(beg_row_col);

for j = 1:count_loc
    m = beg_row_col(j);
    % lat_list is a vector of the latitudes in degrees for the different
    % locations in the total data record.
    lat_list(j) = lat_col(m);
end

end
