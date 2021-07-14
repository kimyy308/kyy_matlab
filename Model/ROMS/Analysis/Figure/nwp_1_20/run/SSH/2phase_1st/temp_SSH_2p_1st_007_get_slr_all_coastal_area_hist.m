close all; clear all;  clc;
warning off;

GCM_testnames= {'CNRM-ESM2-1', 'EC-Earth3-Veg', 'ACCESS-CM2', 'CNRM-CM6-1-HR', 'CMCC-ESM2'};


all_region2 ={'AKP4'};
% all_region2 ={'NWP'};
all_subregion2 ={'YS_KHOA', 'SS_KHOA', 'ES_KHOA'};

% scennames={'rcp26', 'rcp45', 'rcp85'};
scennames={'historical'};

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
tmp.fs=filesep;
if (strcmp(computer,'PCWIN64'))
    tmp.dropboxpath = 'C:\Users\User\Dropbox';
else
    tmp.dropboxpath = '/home/kimyy/Dropbox';
end
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));

[tmp.error_status, tmp.dropboxpath] = Func_0008_set_dropbox_path(computer);

regionind2 = 1;

% for snu_desktopd
% testname=all_testname2{testnameind2}    % % need to change
inputyear = [1993:2014]; % % put year which you want to plot [year year ...]
inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]

varname ='zeta';
variable='SSH';
run('nwp_polygon_point.m');
regionname=all_region2{regionind2};


for subregionind=1:length(all_subregion2)
    subregionname = all_subregion2{subregionind};
for testind=1:length(GCM_testnames)
    for scenind=1:length(scennames)
    
    scenname = scennames{scenind};
    if (strcmp(scenname,'rcp26')==1)
        drivename='E';
    elseif (strcmp(scenname,'historical')==1)
        drivename='D';
    elseif (strcmp(scenname,'rcp45')==1)
        drivename='E';
    elseif (strcmp(scenname,'rcp85')==1)
        drivename='E';
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
    [RCM_testname, error_status] = Func_0006_get_RCMname_from_GCM(GCM_testname, scenname)
    
    filedir=strcat(drivename, ':\Data\Model\ROMS\nwp_1_20\', RCM_testname, '\run\');
    RCM_modelfilenames{testind} = strcat(filedir, RCM_testname,'_',regionname, '_ssh_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    RCM_cmemsfilenames{testind} = strcat(filedir, RCM_testname,'_',regionname, 'cmems_ssh_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    RCM_interpedfilenames{testind} = strcat(filedir, RCM_testname,'_',regionname, 'cmems_interped_ssh_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    RCM_movfilenames{testind} = strcat(filedir, RCM_testname,'_',regionname, 'moving_averaged_cmems_interped_ssh_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    
%     ncinfo(RCM_interpedfilenames{testind})
    lon=ncread(RCM_interpedfilenames{testind}, 'lon_cmems');
    lat=ncread(RCM_interpedfilenames{testind}, 'lat_cmems');
    for folding=1:1
        switch(regionname)
            case('NWP') %% North western Pacific
                lonlat = [115, 164, 15, 52];  %% whole data area
                refpolygon(1,1)=lonlat(1);
                refpolygon(2,1)=lonlat(2);
                refpolygon(1,2)=lonlat(3);
                refpolygon(2,2)=lonlat(4);
            case('AKP4') %% Around Korea Peninsula
                refpolygon=akp4polygon;
            otherwise
                ('?')
        end
        lonlat(1)=min(refpolygon(:,1));
        lonlat(2)=max(refpolygon(:,1));
        lonlat(3)=min(refpolygon(:,2));
        lonlat(4)=max(refpolygon(:,2));

        switch(subregionname)
            case('ES_KHOA') %% East Sea
                subrefpolygon=es_khoapolygon;
            case('SS_KHOA') %% South Sea
                subrefpolygon=ss_khoapolygon;
            case('YS_KHOA') %% Yellow Sea
                subrefpolygon=ys_khoapolygon;
            otherwise
                ('?')
        end
        sublonlat(1)=min(subrefpolygon(:,1));
        sublonlat(2)=max(subrefpolygon(:,1));
        sublonlat(3)=min(subrefpolygon(:,2));
        sublonlat(4)=max(subrefpolygon(:,2));

        submask_model = double(inpolygon(lon,lat,subrefpolygon(:,1),subrefpolygon(:,2)));
        submask_model(submask_model==0)=NaN;
    end
    
    sl=ncread(RCM_interpedfilenames{testind}, 'interped_sla_yearly_poly1_fit');
    sl1=sl(:,:,1).* submask_model;
    sl2=sl(:,:,end).* submask_model;
    slr=sl2-sl1;
%     meanslr=mean(slr(:),'omitnan');
    [meanslr, error_status] = Func_0011_get_area_weighted_mean(slr, lon, lat)
    table_slr((scenind-1)*5+testind, subregionind) = meanslr;

%     pcolor(trend'); shading flat; colorbar
    cmemsdir='Z:\내 드라이브\Data\Observation\CMEMS\';
    cmems_filename = strcat(filedir, RCM_testname,'_',regionname, 'cmems_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');

    end
end
end
