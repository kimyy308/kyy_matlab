close all; clc; clear all;


locinfo.lonw=210;
locinfo.lone=270;
locinfo.lats=-20;
locinfo.latn=10;



datadir = '/Users/kimyy/Desktop/bgc_predictability/Linear_Inverse_model/data_nc';

% read LE
file_le = [datadir, filesep, 'LIM_input_LE_', ...
    'x1_', num2str(locinfo.lonw), ...
    '_x2_', num2str(locinfo.lone), ...
    '_y1_', num2str(locinfo.lats), ...
    '_y2_', num2str(locinfo.latn), '.nc'];

% Get info about the NetCDF file
LE_info = ncinfo(file_le);
% Create an empty struct to hold the variables
allVars = struct();
% Loop over each variable in the file
for i = 1:length(LE_info.Variables)
    varName = LE_info.Variables(i).Name;
    % Read the variable data
    LE_data = ncread(file_le, varName);
    % Store in struct with the variable name as a field
    LE_alldata.(varName) = LE_data;
end

refDate = datetime(1955, 1, 17);
timeDates = refDate + days(LE_alldata.time); % wrong
LE_alldata.year = repelem(1955:2020, 12)';
LE_alldata.month = repmat(1:12, 1, length(1955:2020))';


% Now allVars contains every variable as a field
disp(allVars);


% read SST obs
file_obsSST = [datadir, filesep, 'LIM_input_obs_TEMP_', ...
    'x1_', num2str(locinfo.lonw), ...
    '_x2_', num2str(locinfo.lone), ...
    '_y1_', num2str(locinfo.lats), ...
    '_y2_', num2str(locinfo.latn), '.nc'];

obsSST_info = ncinfo(file_obsSST);
allVars = struct();
for i = 1:length(obsSST_info.Variables)
    varName = obsSST_info.Variables(i).Name;
    obsSST_data = ncread(file_obsSST, varName);
    obsSST_alldata.(varName) = obsSST_data;
end

refDate = datetime(1870, 1, 1);
timeDates = refDate + days(obsSST_alldata.time);
obsSST_alldata.year = year(timeDates);
obsSST_alldata.month = month(timeDates);

disp(allVars);

% read SSH obs
file_obsSSH = [datadir, filesep, 'LIM_input_obs_SSH_', ...
    'x1_', num2str(locinfo.lonw), ...
    '_x2_', num2str(locinfo.lone), ...
    '_y1_', num2str(locinfo.lats), ...
    '_y2_', num2str(locinfo.latn), '.nc'];

obsSSH_info = ncinfo(file_obsSSH);
allVars = struct();
for i = 1:length(obsSSH_info.Variables)
    varName = obsSSH_info.Variables(i).Name;
    obsSSH_data = ncread(file_obsSSH, varName);
    obsSSH_alldata.(varName) = obsSSH_data;
end
obsSSH_alldata.time=obsSSH_alldata.time+365*5;

refDate = datetime(1988, 1, 17);
timeDates = refDate + days(obsSSH_alldata.time);
obsSSH_alldata.year = year(timeDates);
obsSSH_alldata.month = month(timeDates);

disp(allVars);


% read photoC_TOT_zint_100m obs
file_obsphotoC_TOT_zint_100m = [datadir, filesep, 'LIM_input_obs_photoC_TOT_zint_100m_', ...
    'x1_', num2str(locinfo.lonw), ...
    '_x2_', num2str(locinfo.lone), ...
    '_y1_', num2str(locinfo.lats), ...
    '_y2_', num2str(locinfo.latn), '.nc'];

obsphotoC_TOT_zint_100m_info = ncinfo(file_obsphotoC_TOT_zint_100m);
allVars = struct();
for i = 1:length(obsphotoC_TOT_zint_100m_info.Variables)
    varName = obsphotoC_TOT_zint_100m_info.Variables(i).Name;
    obsphotoC_TOT_zint_100m_data = ncread(file_obsphotoC_TOT_zint_100m, varName);
    obsphotoC_TOT_zint_100m_alldata.(varName) = obsphotoC_TOT_zint_100m_data;
end
obsphotoC_TOT_zint_100m_alldata.time=obsphotoC_TOT_zint_100m_alldata.time+365*5;

refDate = datetime(1993, 1, 17);
timeDates = refDate + days(obsphotoC_TOT_zint_100m_alldata.time);
obsphotoC_TOT_zint_100m_alldata.year = year(timeDates);
obsphotoC_TOT_zint_100m_alldata.month = month(timeDates);

disp(allVars);

save([datadir, filesep, 'LIM_input_mat_', 'x1_', num2str(locinfo.lonw), ...
    '_x2_', num2str(locinfo.lone), ...
    '_y1_', num2str(locinfo.lats), ...
    '_y2_', num2str(locinfo.latn), '.mat'], ...
    'LE_alldata', 'LE_info', 'obsSST_alldata', 'obsSST_info', 'obsSSH_alldata', 'obsSSH_info', 'obsphotoC_TOT_zint_100m_alldata', 'obsphotoC_TOT_zint_100m_info');
