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

% stime=1955;
stime=1850;
% etime=2020;
etime=2024;
refDate = datetime(stime, 1, 17);
timeDates = refDate + days(LE_alldata.time); % wrong
LE_alldata.year = repelem(stime:etime, 12)';
LE_alldata.month = repmat(1:12, 1, length(stime:etime))';


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

refDate = datetime(1945, 1, 17);
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

refDate = datetime(1895, 1, 17);
timeDates = refDate + days(obsphotoC_TOT_zint_100m_alldata.time);
obsphotoC_TOT_zint_100m_alldata.year = year(timeDates);
obsphotoC_TOT_zint_100m_alldata.month = month(timeDates);


% HBLT obs
file_obsHBLT = [datadir, filesep, 'LIM_input_obs_HBLT_', ...
    'x1_', num2str(locinfo.lonw), ...
    '_x2_', num2str(locinfo.lone), ...
    '_y1_', num2str(locinfo.lats), ...
    '_y2_', num2str(locinfo.latn), '.nc'];

obsHBLT_info = ncinfo(file_obsHBLT);
allVars = struct();
for i = 1:length(obsHBLT_info.Variables)
    varName = obsHBLT_info.Variables(i).Name;
    obsHBLT_data = ncread(file_obsHBLT, varName);
    obsHBLT_alldata.(varName) = obsHBLT_data;
end
obsHBLT_alldata.time=obsHBLT_alldata.time_counter;
obsHBLT_alldata.time=obsHBLT_alldata.time+365*5;

refDate = datetime(1958, 1, 17);
timeDates = refDate + seconds(obsHBLT_alldata.time);
obsHBLT_alldata.year = year(timeDates);
obsHBLT_alldata.month = month(timeDates);


% U010 obs
file_obsU010 = [datadir, filesep, 'LIM_input_obs_U010_', ...
    'x1_', num2str(locinfo.lonw), ...
    '_x2_', num2str(locinfo.lone), ...
    '_y1_', num2str(locinfo.lats), ...
    '_y2_', num2str(locinfo.latn), '.nc'];

obsU010_info = ncinfo(file_obsU010);
allVars = struct();
for i = 1:length(obsU010_info.Variables)
    varName = obsU010_info.Variables(i).Name;
    obsU010_data = ncread(file_obsU010, varName);
    obsU010_alldata.(varName) = obsU010_data;
end
obsU010_alldata.time=obsU010_alldata.time+365*5;
% refDate = datetime(1958, 1, 17);
% timeDates = refDate + hours(obsU010_alldata.time);
% obsU010_alldata.year = year(timeDates);
% obsU010_alldata.month = month(timeDates);

abc=repmat(1950:2024, [12 1]);
obsU010_alldata.year=abc(:);
obsU010_alldata.month=repmat(1:12, [1 900/12])';


% V010 obs
file_obsV010 = [datadir, filesep, 'LIM_input_obs_V010_', ...
    'x1_', num2str(locinfo.lonw), ...
    '_x2_', num2str(locinfo.lone), ...
    '_y1_', num2str(locinfo.lats), ...
    '_y2_', num2str(locinfo.latn), '.nc'];

obsV010_info = ncinfo(file_obsV010);
allVars = struct();
for i = 1:length(obsV010_info.Variables)
    varName = obsV010_info.Variables(i).Name;
    obsV010_data = ncread(file_obsV010, varName);
    obsV010_alldata.(varName) = obsV010_data;
end
obsV010_alldata.time=obsV010_alldata.time+365*5;
% refDate = datetime(1958, 1, 17);
% timeDates = refDate + hours(obsV010_alldata.time);
% obsV010_alldata.year = year(timeDates);
% obsV010_alldata.month = month(timeDates);

abc=repmat(1960:2024, [12 1]);
obsV010_alldata.year=abc(:);
obsV010_alldata.month=repmat(1:12, [1 780/12])';


disp(allVars);

save([datadir, filesep, 'LIM_input_mat_', 'x1_', num2str(locinfo.lonw), ...
    '_x2_', num2str(locinfo.lone), ...
    '_y1_', num2str(locinfo.lats), ...
    '_y2_', num2str(locinfo.latn), '.mat'], ...
    'LE_alldata', 'LE_info', 'obsSST_alldata', 'obsSST_info', 'obsSSH_alldata', 'obsSSH_info', ...
    'obsphotoC_TOT_zint_100m_alldata', 'obsphotoC_TOT_zint_100m_info', ...
    'obsHBLT_alldata', 'obsHBLT_info', 'obsU010_alldata', 'obsU010_info', 'obsV010_alldata', 'obsV010_info');


clim_SST=reshape(obsSST_data,[12, length(obsSST_data)/12]);
clim_mean_SST=mean(clim_SST,2);
clim_mean_SST_2d=repmat(clim_mean_SST, [1, length(obsSST_data)/12]);
yearly_SST=mean(clim_SST,1);
clim_mean_SST_rep=clim_mean_SST_2d(:);
mon_SST=obsSST_data-clim_mean_SST_rep;

clim_SSH=reshape(obsSSH_data(1:end-6),[12, length(obsSSH_data(1:end-6))/12]);
clim_mean_SSH=mean(clim_SSH,2);
clim_mean_SSH_2d=repmat(clim_mean_SSH, [1, length(obsSSH_data(1:end-6))/12]);
yearly_SSH=mean(clim_SSH,1);
clim_mean_SSH_rep=clim_mean_SSH_2d(:);
mon_SSH=obsSSH_data(1:end-6)-clim_mean_SSH_rep;

clim_photoC_TOT_zint_100m=reshape(obsphotoC_TOT_zint_100m_data,[12, length(obsphotoC_TOT_zint_100m_data)/12]);
clim_mean_photoC_TOT_zint_100m=mean(clim_photoC_TOT_zint_100m,2);
clim_mean_photoC_TOT_zint_100m_2d=repmat(clim_mean_photoC_TOT_zint_100m, [1, length(obsphotoC_TOT_zint_100m_data)/12]);
yearly_photoC_TOT_zint_100m=mean(clim_photoC_TOT_zint_100m,1);
clim_mean_photoC_TOT_zint_100m_rep=clim_mean_photoC_TOT_zint_100m_2d(:);
mon_photoC_TOT_zint_100m=obsphotoC_TOT_zint_100m_data-clim_mean_photoC_TOT_zint_100m_rep;

clim_HBLT=reshape(obsHBLT_data,[12, length(obsHBLT_data)/12]);
clim_mean_HBLT=mean(clim_HBLT,2);
clim_mean_HBLT_2d=repmat(clim_mean_HBLT, [1, length(obsHBLT_data)/12]);
yearly_HBLT=mean(clim_HBLT,1);
clim_mean_HBLT_rep=clim_mean_HBLT_2d(:);
mon_HBLT=obsHBLT_data-clim_mean_HBLT_rep;


plot(obsSST_alldata.year+obsSST_alldata.month/12, mon_SST/std(mon_SST), '-o', 'linewidth', 2);
hold on
plot(obsSSH_alldata.year(1:end-6)+obsSSH_alldata.month(1:end-6)/12, mon_SSH/std(mon_SSH), '-^', 'linewidth', 2);
plot(obsphotoC_TOT_zint_100m_alldata.year+obsphotoC_TOT_zint_100m_alldata.month/12, ...
    mon_photoC_TOT_zint_100m/std(mon_photoC_TOT_zint_100m), '-x', 'linewidth', 2);
plot(obsHBLT_alldata.year+obsHBLT_alldata.month/12, mon_HBLT/std(mon_HBLT), '-+', 'linewidth', 2);
hold off
legend({'SST', 'SSH', 'NPP', 'HBLT'});
grid minor;

corrcoef(mon_HBLT(end-323:end-12), mon_photoC_TOT_zint_100m(end-323:end-12))
corrcoef(mon_SSH(end-311:end), mon_photoC_TOT_zint_100m(end-323:end-12))
corrcoef(mon_SST(end-323:end-12), mon_photoC_TOT_zint_100m(end-323:end-12))