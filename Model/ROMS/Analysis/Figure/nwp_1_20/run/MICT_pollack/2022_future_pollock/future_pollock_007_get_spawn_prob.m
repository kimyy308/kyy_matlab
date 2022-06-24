% %  Updated 27-Apr-2021 by Yong-Yub Kim, 

close all; clear all;  clc;
warning off;

RCM_info.model_reana = 'nwp_1_10';

% RCM_info.name={ 'test2127'};
RCM_info.name={'test2117', 'test2118', 'test2119', 'test2120', 'test2121'};
% RCM_info.name={'test2127', 'test2128', 'test2129', 'test2130', 'test2131'};

% 
RCM_info.model = 'nwp_1_20';

RCM_info.dataroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, RCM_info.model, filesep];

RCM_info.saveroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, RCM_info.model, filesep];
RCM_info.phase = 'run';  % run or spinup
% RCM_info.region = {'EKB2'}; % NWP, AKP4, ES_KHOA, YS, ...
RCM_info.region = {'pollock_egg3'}; % NWP, AKP4, ES_KHOA, YS, ...

RCM_info.vars = {'SST'};

% D:\Data\Model\ROMS\nwp_1_20\test2127\pollock

tmp.pollock_depth_min=50;
tmp.pollock_depth_max=500;
tmp.pollock_lon_max=132;
tmp.pollock_temp_min=2;
tmp.pollock_temp_max=5;

% RCM_info.years = 1983:2021;  
% RCM_info.years = [2015:2050, 2081:2100];  
% RCM_info.years = [1983:1987];  
% RCM_info.years = [1988:1992];  
RCM_info.years = [1995:2014];  
% RCM_info.years = [1993:2021];  
% RCM_info.years = [2081:2100];  
% RCM_info.years = [2081:2100];  

% seasons_group={'February', 'January', 'JF-'};
seasons_group={'JF-'};

switch RCM_info.name{1}
    case 'test2127'
        tmp.ensname='ens2201';
    case 'test2117'
        tmp.ensname='ens2202';
end

% RCM_grid.dl = 1/20;
for seasons_groupi=1:length(seasons_group)
    for regionind2=1:length(RCM_info.region)
%% set path
        tmp.fs=filesep;  
        addpath(genpath('C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\MICT_pollack\2022_future_pollock\subroutine\'))

        tmp.dropboxpath = 'C:\Users\User\Dropbox';
        addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
        [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);
        addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Model' ...
            tmp.fs, 'ROMS', tmp.fs, 'Analysis', tmp.fs, 'Figure', tmp.fs, 'nwp_1_20', tmp.fs ...
            'run', tmp.fs, 'SSH', tmp.fs, '2phase_2nd', tmp.fs, 'subroutine']));

        RCM_info.season=seasons_group{seasons_groupi};
        [RCM_info.months,error_status] = Func_0019_get_month_from_season(RCM_info.season);

        tmp.testname=RCM_info.name{1};   % % need to change
        [tmp.abb, tmp.error_status] = Func_0018_abbname(tmp.testname);
        tmp.regionname=RCM_info.region{regionind2};
        [RCM_grid.refpolygon, RCM_grid.domain, tmp.error_status] = Func_0007_get_polygon_data_from_regionname(tmp.regionname);

%                 tmp.variable=RCM_info.vars{varind2};
% 
%                 dirs.figrawdir =strcat('D:\Research\Ph_D_course\2022_pollock_future\figure',filesep, 'all', filesep); % % where figure files will be saved            
%                 tmp.param_script =['C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\', RCM_info.model, '\run\fig_param\fig_param2_kyy_', tmp.regionname, '.m'];

%% set file info

%                 dirs.matdir = dirs.filedir;
%                 dirs.griddir = dirs.filedir; % % where grid data files are            

        for yearij = 1:length(RCM_info.years)
            tmp.tempyear=RCM_info.years(yearij);
            for monthij = 1:length(RCM_info.months)
                tic;
                tmp.tempmonth=RCM_info.months(monthij);
                dirs.filedir = ['D:', filesep, 'Data', filesep, 'Model', filesep, 'ROMS', ...
                    filesep, RCM_info.model, filesep, tmp.testname, filesep, 'pollock', filesep]; % % where data files are
                tmp.filename = [dirs.filedir, tmp.testname, '_', tmp.regionname, 'model_pollock_', ...
                         num2str(tmp.tempyear,'%04i'), '_', num2str(tmp.tempmonth, '%02i'), '.nc'];
                ncinfo(tmp.filename);
                if yearij==1 & monthij ==1
                    RCM_grid.lon_rho=ncread(tmp.filename, 'lon_rho');
                    tmp.lon_1d=RCM_grid.lon_rho(:,1);
                    tmp.lon_lim_ind=max(find(tmp.lon_1d<=132));

                    RCM_grid.h=ncread(tmp.filename, 'h');
                    RCM_grid.depthmask=NaN(size(RCM_grid.h));
                    RCM_grid.depthmask(RCM_grid.h>=50 & RCM_grid.h<=500)=1;
                    RCM_grid.depthmask(tmp.lon_lim_ind:end,:)=NaN;
%                             pcolor(RCM_grid.depthmask'); shading flat; colorbar;
                end
                tmp.tlen=length(ncread(tmp.filename, 'time'));
                tmp.temp_surf=NaN([size(RCM_grid.h), tmp.tlen, length(RCM_info.name)]);
                tmp.temp_surf(:,:,:,1)=ncread(tmp.filename, 'temp_surf');
%                         RCM_data.temp_surf(:,:,:,1)=RCM_data.temp_surf.*RCM_grid.depthmask;

                for testnameind2=2:length(RCM_info.name)
                    tmp.testname=RCM_info.name{testnameind2};   % % need to change
                    dirs.filedir = ['D:', filesep, 'Data', filesep, 'Model', filesep, 'ROMS', ...
                        filesep, RCM_info.model, filesep, tmp.testname, filesep, 'pollock', filesep]; % % where data files are
                    tmp.filename = [dirs.filedir, tmp.testname, '_', tmp.regionname, 'model_pollock_', ...
                         num2str(tmp.tempyear,'%04i'), '_', num2str(tmp.tempmonth, '%02i'), '.nc'];
                    tmp.temp_surf(:,:,:,testnameind2)=ncread(tmp.filename, 'temp_surf');
                end
                RCM_data.low_prob=NaN(size(RCM_grid.h,1),size(RCM_grid.h,2));
                RCM_data.high_prob=NaN(size(RCM_grid.h,1),size(RCM_grid.h,2));
                RCM_data.prob=NaN(size(RCM_grid.h,1),size(RCM_grid.h,2));
                for lonij=1:size(RCM_grid.h,1)
                    for latij=1:size(RCM_grid.h,2)
                        if RCM_grid.depthmask(lonij,latij)==1
                            for dayij=1:size(tmp.temp_surf,3)
                                tmp.sample=squeeze(tmp.temp_surf(lonij,latij,dayij,:));
                                if sum(isfinite(tmp.sample))~=0
                                    tmp.pd = fitdist(tmp.sample,'Normal');
                                    tmp.prob_range=cdf(tmp.pd,[tmp.pollock_temp_min,tmp.pollock_temp_max]).*100;
                                    RCM_data.low_prob(lonij,latij,dayij)=tmp.prob_range(1);
                                    RCM_data.high_prob(lonij,latij,dayij)=tmp.prob_range(2);
                                    RCM_data.prob(lonij,latij,dayij)=diff(tmp.prob_range);
                                end
                            end
                        end
                    end
                end
                toc;
                disp([num2str(tmp.tempyear), ', ', num2str(tmp.tempmonth)])
%                         pcolor(temp_surf(:,:,1)'); shading flat; colorbar;
                dirs.filedir = ['D:', filesep, 'Data', filesep, 'Model', filesep, 'ROMS', ...
                    filesep, RCM_info.model, filesep, tmp.ensname, filesep, 'pollock', filesep]; % % where data files are
                tmp.filename = [dirs.filedir, tmp.ensname, '_', tmp.regionname, 'model_pollock_', ...
                     num2str(tmp.tempyear,'%04i'), '_', num2str(tmp.tempmonth, '%02i'), '.nc'];
%                 ncinfo(tmp.filename)
                nccreate(tmp.filename, 'spawn_prob', ...
                    'Dimensions', {'lon', size(RCM_grid.h,1), 'lat', size(RCM_grid.h,2), 'time', tmp.tlen}, ...
                    'DeflateLevel', 1,'Shuffle', true);
                ncwrite(tmp.filename, 'spawn_prob', RCM_data.prob);
            end
        end

    end
end


% ncinfo(tmp.filename)



% % % % data=[2 4 6 7 10];
% % % % % % % pd = makedist('Normal','mu',mean(data),'sigma',std(data))
% % % % pd = fitdist(data','Normal')
% % % % y=cdf(pd,[2,5])
% % % % 
% % % % p=cdf(pd,data)
% % % % y=normpdf(data,mean(data),std(data))
% % % % plot(data,y)