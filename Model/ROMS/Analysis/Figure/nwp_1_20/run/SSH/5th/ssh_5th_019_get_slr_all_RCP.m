close all; clear all;  clc;
warning off;
% all_region2 ={'NWP','NWP2'}
% all_region2 ={'ES', 'YS', 'SS', 'NWP2'}
% all_testname2 = {'ens02','test53','test55','test56'};
% all_testname2 = {'test11', 'test12'};
% all_testname2 = {'test11', 'test12'};
% RCM_testnames = {'test61', 'test62', 'test63', 'test64'};
% RCM_testnames = {'test57', 'test58', 'test59', 'test60', 'ens'};
% RCM_testnames = {'test65', 'test66', 'test67', 'test68'};

GCM_testnames = {'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'NorESM1-M', 'MPI-ESM-LR', 'ens'};
% all_region2 ={'NWP','ES', 'SS', 'YS', 'AKP'}
% all_region2 ={'NWP','AKP2', 'YS', 'NES', 'SES'}

all_region2 ={'AKP4'};
% all_region2 ={'NWP'};

scennames={'rcp26', 'rcp45', 'rcp85'};
% % scenname='rcp26';
% scenname='rcp45';
% % scenname='rcp85';

close all;
% clearvars '*' -except regionind2 testnameind2 all_region2 all_testname2 testname testname2 sodatestname
% % % 
system_name=computer;

dropboxpath='C:\Users\User\Dropbox';
addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path

regionind2 = 1;

% for snu_desktopd
% testname=all_testname2{testnameind2}    % % need to change
inputyear = [2006:2030]; % % put year which you want to plot [year year ...]
inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]

varname ='zeta';
variable='SSH';
run('nwp_polygon_point.m');
regionname=all_region2{regionind2};



for folding=1:1
    switch(regionname)
        case('NWP') %% North western Pacific
            lonlat = [115, 164, 15, 52];  %% whole data area
            refpolygon(1,1)=lonlat(1);
            refpolygon(2,1)=lonlat(2);
            refpolygon(1,2)=lonlat(3);
            refpolygon(2,2)=lonlat(4);
        case('NWP2') %% North western Pacific
            lonlat = [115, 145, 25, 52];  %% whole data area
            refpolygon(1,1)=lonlat(1);
            refpolygon(2,1)=lonlat(2);
            refpolygon(1,2)=lonlat(3);
            refpolygon(2,2)=lonlat(4);
        case('ES') %% East Sea
            refpolygon=espolygon;
        case('NES') %% Northern East Sea
            refpolygon=nespolygon;
        case('SES') %% Southern East Sea
            refpolygon=sespolygon;
        case('SS') %% South Sea
            refpolygon=sspolygon;
        case('YS') %% Yellow Sea
            refpolygon=yspolygon;
        case('ECS') %% East China Sea
            refpolygon=ecspolygon;
        case('AKP') %% Around Korea Peninsula
            refpolygon=akppolygon;
        case('AKP2') %% Around Korea Peninsula
            refpolygon=akp2polygon;
        case('AKP3') %% Around Korea Peninsula
            refpolygon=akp3polygon;
        case('AKP4') %% Around Korea Peninsula
            refpolygon=akp4polygon;
        case('CA') %% Around Korea Peninsula
            refpolygon=capolygon;
        case('EKB') %% Around Korea Peninsula
            refpolygon=akp2polygon;
        otherwise
            ('?')
    end
    lonlat(1)=min(refpolygon(:,1));
    lonlat(2)=max(refpolygon(:,1));
    lonlat(3)=min(refpolygon(:,2));
    lonlat(4)=max(refpolygon(:,2));
end

for testind=1:length(GCM_testnames)
    for scenind=1:length(scennames)
    
    scenname = scennames{scenind};
    if (strcmp(scenname,'rcp26')==1)
        drivename='D';
    elseif (strcmp(scenname,'rcp45')==1)
        drivename='D';
    elseif (strcmp(scenname,'rcp85')==1)
        drivename='D';
    end
        
    GCM_testname=GCM_testnames{testind};
    if strcmp(GCM_testname, 'ens')==1
        switch scenname
            case 'rcp26'
                GCM_testname='ens09';
            case 'rcp45'
                GCM_testname='ens08';
            case 'rcp85'
                GCM_testname='ens10';
        end
        RCM_testname =GCM_testname;
    else
        if strcmp(GCM_testname, 'IPSL-CM5A-LR')==1 && strcmp(scenname, 'rcp45')==1
            RCM_testname='test57';
        elseif strcmp(GCM_testname, 'IPSL-CM5A-MR')==1 && strcmp(scenname, 'rcp45')==1
            RCM_testname='test58';
        elseif strcmp(GCM_testname, 'NorESM1-M')==1 && strcmp(scenname, 'rcp45')==1
            RCM_testname='test59';
        elseif strcmp(GCM_testname, 'MPI-ESM-LR')==1 && strcmp(scenname, 'rcp45')==1
            RCM_testname='test60';
        elseif strcmp(GCM_testname, 'IPSL-CM5A-LR')==1 && strcmp(scenname, 'rcp26')==1
            RCM_testname='test61';
        elseif strcmp(GCM_testname, 'IPSL-CM5A-MR')==1 && strcmp(scenname, 'rcp26')==1
            RCM_testname='test62';
        elseif strcmp(GCM_testname, 'NorESM1-M')==1 && strcmp(scenname, 'rcp26')==1
            RCM_testname='test63';
        elseif strcmp(GCM_testname, 'MPI-ESM-LR')==1 && strcmp(scenname, 'rcp26')==1
            RCM_testname='test64';
        elseif strcmp(GCM_testname, 'IPSL-CM5A-LR')==1 && strcmp(scenname, 'rcp85')==1
            RCM_testname='test65';
        elseif strcmp(GCM_testname, 'IPSL-CM5A-MR')==1 && strcmp(scenname, 'rcp85')==1
            RCM_testname='test66';
        elseif strcmp(GCM_testname, 'NorESM1-M')==1 && strcmp(scenname, 'rcp85')==1
            RCM_testname='test67';
        elseif strcmp(GCM_testname, 'MPI-ESM-LR')==1 && strcmp(scenname, 'rcp85')==1
            RCM_testname='test68';
        end      
    end

    filedir=strcat('D:\Data\Model\CMIP5\zos\', scenname, '\interp\', GCM_testname, '\');
    GCM_modelfilenames{testind} = strcat(filedir, GCM_testname,'_',regionname, '_ssh_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    GCM_cmemsfilenames{testind} = strcat(filedir, GCM_testname,'_',regionname, 'cmems_ssh_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    GCM_interpedfilenames{testind} = strcat(filedir, GCM_testname,'_',regionname, 'cmems_interped_ssh_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    GCM_movfilenames{testind} = strcat(filedir, GCM_testname,'_',regionname, 'moving_averaged_cmems_interped_ssh_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    
    sl=ncread(GCM_interpedfilenames{testind}, 'interped_sla_yearly_exp_fit');
    sl1=sl(:,:,1);
    sl2=sl(:,:,end);
    slr=sl2-sl1;
    meanslr=mean(slr(:),'omitnan');
    table_slr((scenind-1)*2+1, testind) = meanslr;
    
%     testname=RCM_testnames{testind};
    filedir=strcat(drivename, ':\Data\Model\ROMS\nwp_1_20\', RCM_testname, '\run\');
    RCM_modelfilenames{testind} = strcat(filedir, RCM_testname,'_',regionname, '_ssh_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    RCM_cmemsfilenames{testind} = strcat(filedir, RCM_testname,'_',regionname, 'cmems_ssh_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    RCM_interpedfilenames{testind} = strcat(filedir, RCM_testname,'_',regionname, 'cmems_interped_ssh_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    RCM_movfilenames{testind} = strcat(filedir, RCM_testname,'_',regionname, 'moving_averaged_cmems_interped_ssh_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    
    sl=ncread(RCM_interpedfilenames{testind}, 'interped_sla_yearly_exp_fit');
    sl1=sl(:,:,1);
    sl2=sl(:,:,end);
    slr=sl2-sl1;
    meanslr=mean(slr(:),'omitnan');
    table_slr((scenind-1)*2+2, testind) = meanslr;
    
    cmemsdir='Z:\내 드라이브\Data\Observation\CMEMS\';
    cmems_filename = strcat(filedir, RCM_testname,'_',regionname, 'cmems_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');

    end
end
