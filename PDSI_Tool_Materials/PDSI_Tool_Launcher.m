%% © 2012 Vanderbilt University %%

function varargout = PDSI_Tool_Launcher(varargin)
% PDSI_Tool_Launcher MATLAB code for PDSI_Tool_Launcher.fig
%      PDSI_Tool_Launcher, by itself, creates a new PDSI_Tool_Launcher or raises the existing
%      singleton*.
%
%      H = PDSI_Tool_Launcher returns the handle to a new PDSI_Tool_Launcher or the handle to
%      the existing singleton*.
%
%      PDSI_Tool_Launcher('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PDSI_Tool_Launcher.M with the given input 
%      arguments.
%
%      PDSI_Tool_Launcher('Property','Value',...) creates a new PDSI_Tool_Launcher or raises 
%      the existing singleton*.  Starting from the left, property value 
%      pairs are applied to the GUI before PDSI_Tool_Launcher_OpeningFcn gets called. 
%      An unrecognized property name or invalid value makes property 
%      application stop. All inputs are passed to PDSI_Tool_Launcher_OpeningFcn via 
%      varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PDSI_Tool_Launcher

% Last Modified by GUIDE v2.5 15-Apr-2013 14:05:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PDSI_Tool_Launcher_OpeningFcn, ...
                   'gui_OutputFcn',  @PDSI_Tool_Launcher_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before PDSI_Tool_Launcher is made visible.
function PDSI_Tool_Launcher_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PDSI_Tool_Launcher (see VARARGIN)

% Choose default command line output for PDSI_Tool_Launcher
handles.output = hObject;

set(handles.pet_buttongroup,'SelectionChangeFcn', ...
    @pet_buttongroup_SelectionChangeFcn);
set(handles.calyrs_buttongroup, 'SelectionChangeFcn', ...
    @calyrs_buttongroup_SelectionChangeFcn);
set(handles.radiobutton_thorn,'Value',1)
set(handles.radiobutton_noaa,'Value',1)
set(handles.help_text,'Visible','off')

W1 = {'Welcome to the MATLAB Palmer Drought Indices Calculator.'};
W2 = {'Please load your data, select your desired PET calculation method'};
W3 = {'and calibration period, and run the program.'};
W4 = {'Click the "Help" button for directions and answers to any questions.'};
dash = {'--------'};
space = {''};
set(handles.welcome_text,'String',[space;W1;dash;W2;W3;dash;W4;]);
H4 = {'------- Temperature and Precipitation Inputs -------'};
H5 = {'- The US Historical Climate Network (USHCN) has historical temperature and precipitation data available for US stations at http://cdiac.ornl.gov/epubs/ndp/ushcn/access.html.'};
H6 = {'- Data must be consecutive (i.e., no gaps) and combined into one .txt file. Do not include column headers.'};
H7 = {'- Data should be chronologically organized into four columns: Column 1 is the latitude in degrees, Column 2 is the year, Column 3 is temperature, and Column 4 is precipitation.'};
H8 = {'- Note that each temperature and precipitation observation must have a latitude and a year associated with it.'};
H9 = {'- If there are data for multiple locations, make sure that the data are arranged in the same order.'}; 
H10 = {'- Temperature and precipitation data should be input into the program in degrees Fahrenheit and inches, respectively.'};
H11 = {'-------- Available Water Capacity (AWC) Input --------'};
H12 = {'- Different locations have different field capacities. An AWC value for each location must be loaded into the program in inches.'};
H13 = {'- AWC data should be organized in a column in one .txt file. The AWCs should be organized in the same location order as the temperature and precipitation data.'};
H14 = {'- Note that only one AWC value is needed per location (i.e. the number of AWC values must equal the number of stations).'};
H15 = {'------- Calibration Period -------'};
H16 = {'- The calibration period is used to calculate the "Climatologically Appropriate for Existing Conditions" (CAFEC) precipitation for a location.'};
H17 = {'- The CAFEC is used to calculate weighting factors used in the Z-Index calculation (see Palmer, 1965).'};
H18 = {'- NOAA uses the period January 1931 to December 1990 as its calibration period (see Karl, 1986), and this option is provided.'};
H19 = {'- In the absence of long data records or to use a more comprehensive timespan, the option is also provided to use the full record as the calibration period.'};
H20 = {'------- Operation -------'};
H21 = {'- To operate the GUI, load your precipitation, temperature, and AWC data using the "Load Data" buttons, select your desired PET calculation method and calibration period, and click "Run".'};
H22 = {'------- Output Details -------'};
H23 = {'- Units for the output of the PET calculations are in inches. Main drought index outputs are PET, Z-Index, PDSI, and PHDI (see Palmer, 1965).'};
H24 = {'- Results are output into a text file (Palmer.txt) that is saved in the working folder of the current directory. The text file should be opened in Notepad or imported into Excel.'};
H25 = {''};
H26 = {'* For additional help please refer to README.pdf.'};
H27 = {'TO EXIT THIS HELP MENU, PRESS THE "HELP" BUTTON IN THE LOWER RIGHT CORNER'};
set(handles.help_text,'String',[H26;H25;H4;H25;H5;H6;H7;H8;H9;H10;H25;H11; ...
    H25;H12;H13;H14;H25;H15;H25;H16;H17;H18;H19;H25;H20;H25;H21;H25; ...
    H22;H25;H23;H24;H25;H27]);

global RB_Cal RB_PET 

% Set button values to defaults in case user doesn't press them.
RB_PET = 2;
RB_Cal = 1;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PDSI_Tool_Launcher wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PDSI_Tool_Launcher_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_data.
function pushbutton_data_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global DATA lat_col yr_col T_F P

D = uiimport('-file');
fprintf('%s\n%s %s.\n','Temperature and precipitation data successfully loaded.');
E = fieldnames(D);
DATA = getfield(D,char(E));

% Extract latitude, temperature, and precipitation data from large data
% file.
lat_col = DATA(1:end,1);
yr_col = DATA(1:end,2);
T_F = DATA(1:end,3);
P = DATA(1:end,4);

% --- Executes on button press in pushbutton_awc.
function pushbutton_awc_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_awc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global AWC

A = uiimport('-file');
fprintf('%s\n%s %s.\n','AWC data successfully loaded.');
W = fieldnames(A);
AWC_DATA = getfield(A,char(W));

% Extract AWC data from larger data file.
AWC = AWC_DATA(1:end,1);

% --- Executes on button press in togglebutton_help.
function togglebutton_help_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_help

button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
	set(handles.help_text,'Visible','on')
elseif button_state == get(hObject,'Min')
	set(handles.help_text,'Visible','off')
end

function pet_buttongroup_SelectionChangeFcn(hObject, eventdata)

global RB_PET

%retrieve GUI data, i.e. the handles structure
handles = guidata(hObject); 
 
switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'radiobutton_hamon'
      %execute this code when hamon_radiobutton is selected
      RB_PET = 1;
    case 'radiobutton_thorn'
      %execute this code when thorn_radiobutton is selected
      RB_PET = 2;
    otherwise
       % Code for when there is no match.
 
end
%updates the handles structure
guidata(hObject, handles);

function calyrs_buttongroup_SelectionChangeFcn(hObject, eventdata)

global RB_Cal

%retrieve GUI data, i.e. the handles structure
handles = guidata(hObject); 
 
switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'radiobutton_noaa'
      %execute this code when hamon_radiobutton is selected
      RB_Cal = 1;
    case 'radiobutton_full'
      %execute this code when thorn_radiobutton is selected
      RB_Cal = 2;
    otherwise
       % Code for when there is no match.
 
end
%updates the handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_run.
function pushbutton_run_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global RB_Cal RB_PET AWC lat_col yr_col T_F P 

[beg_row_col,end_row_col,beg_row_mat,end_row_mat,count_loc, ...
           lat_list,lat_mat] = Count_Loc(lat_col);

if RB_Cal == 1
    % Use the calibration period that the NCDC and CPC use - January 1931
    % through December 1990 (cf. Karl, 1986; Journal of Climate and Applied
    % Meteorology, Vol. 25, No. 1, January 1986).
    
    % NOTE:
    % THE ORDERING OF LOCATIONS MUST BE CONSISTENT ACROSS INPUTS!
    % yr_col is the record of years for all locations, arranged in a column
    % vector, where a year is associated with each observation (i.e., 
    % month) of the record for all locations.
    % beg_row_col is a vector of the row numbers where the temperature data
    % for each of the different locations begins.
    % beg_row_mat is a vector of the row numbers where the temperature data 
    % for each of the different locations begins when the data is listed in 
    % a matrix such that years represent rows and columns represent months.
    % count_loc is the number of different locations in the total data 
    % record.
    
    % yr_mat_rec is the record of years for all locations, arranged in a
    % matrix, where a year is associated with each MONTH of the record for 
    % all locations.
    yr_mat_rec = (reshape(yr_col,12,length(yr_col)/12))';
    
    % yr is the record of years for all locations, arranged in a column 
    % vector, where each year of the record is listed only once (i.e., one 
    % year for each year of observations) for all locations.
    yr = yr_mat_rec(:,1);
    
    % count1931 is a counter that tracks the rows in the total data record
    % (i.e. column vector) that holds the January 1931 temperature and
    % precipitation data for each of the different locations.
    count1931 = 1;
    
    % count1990 is a counter that tracks the rows in the total data record
    % (i.e., column vector) that holds the December 1990 temperature and
    % precipitation data for each of the different locations.
    count1990 = 1;
    
    % beg_yr_col is a vector of the row numbers that hold the January 1931
    % temperature and precipitation data for each of the different 
    % locations.
    beg_yr_col = [];
    
    % end_yr_col is a vector of the row numbers that the December 1990
    % temperature and precipitation data for each of the different 
    % locations.
    end_yr_col = [];
    
    for j = 1:count_loc
        for m = 1:(end_row_col(j) - beg_row_col(j) + 1)
            if m == 1
                v(m) = beg_row_col(j);
                if yr_col(v(m)) == 1931
                    count1931(m) = count1931(m);
                    count1990(m) = count1990(m);
                    beg_yr_col = [beg_yr_col; count1931(m)];
                elseif yr_col(v(m)) == 1990
                    count1931(m) = count1931(m);
                    count1990(m) = count1990(m);
                    end_yr_col = [end_yr_col; count1990(m)];
                else
                    count1931(m) = count1931(m);
                    count1990(m) = count1990(m);
                end
                continue
            elseif m == (end_row_col(j) - beg_row_col(j) + 1)
                v(m) = end_row_col(j);
            else
                v(m) = v(m-1) + 1;
            end
            if yr_col(v(m)) == 1931
                count1931(m) = count1931(m - 1) + 1;
                count1990(m) = count1990(m - 1) + 1;
                beg_yr_col = [beg_yr_col; count1931(m)];
            elseif yr_col(v(m)) == 1990
                count1931(m) = count1931(m - 1) + 1;
                count1990(m) = count1990(m - 1) + 1;
                end_yr_col = [end_yr_col; count1990(m)];
            else
                count1931(m) = count1931(m - 1) + 1;
                count1990(m) = count1990(m - 1) + 1;
            end
        end
        beg_calyr_col(j) = min(beg_yr_col);
        beg_yr_col = [];
        beg_calyr_mat(j) = (beg_calyr_col(j) - 1)/12 + 1;
        end_calyr_col(j) = max(end_yr_col);
        end_yr_col = [];
        end_calyr_mat(j) = end_calyr_col(j)/12;
    end
    
else
    % Use the entire period of record as the calibration period (cf. Karl, 
    % 1986; Journal of Climate and Applied Meteorology, Vol. 25, No. 1, 
    % January 1986).
    
    for j = 1:count_loc
        beg_calyr_col(j) = 1; 
        end_calyr_col(j) = (end_row_col(j) - beg_row_col(j) + 1);
        beg_calyr_mat(j) = (beg_calyr_col(j) - 1)/12 + 1;
        end_calyr_mat(j) = end_calyr_col(j)/12;
    end
        
end

if RB_PET == 1
    PET = Hamon_PET(T_F,lat_list,count_loc,beg_row_mat,lat_mat);
else
    PET = Thornthwaite_PET(T_F,lat_list,count_loc,beg_row_col, ...
          lat_col);
end

[ET,PR,R,RO,PRO,L,PL] = WaterBalance(AWC,PET,P,beg_row_col, ...
                        end_row_col,count_loc);

[Z_all] = Z_Index(P,PET,ET,PR,R,RO,PRO,L,PL,beg_row_mat,end_row_mat, ...
          count_loc,beg_calyr_mat,end_calyr_mat,beg_row_col,end_row_col);

[table] = PDSI_Central(Z_all,count_loc,beg_row_col,end_row_col, ...
          lat_col,yr_col,PET);

% Open the text file to which the table of values is to be written with 
% write permission.
fid = fopen('Palmer.txt','w');
fprintf(fid, '%10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\r\n','Latitude','Year','PET','Z-Index','PPe','PX1','PX2','PX3','X','PDSI','PHDI');
fprintf(fid, '%10.4f %10.0f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\r\n',(table'));
fclose(fid);

% If the tool successfully executes the Palmer calculations, then print 
% text, including the directory and folder to which the output text file 
% was written, to the command window.
currentFolder = pwd;
screen_text_l1 = 'Success! The output file - "Palmer.txt" -';
screen_text_l2 = 'is located in';
fprintf('%s\n%s %s.\n',screen_text_l1,screen_text_l2,currentFolder)
