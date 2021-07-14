close all; clear all;  clc;
warning off;
% all_region2 ={'NWP','NWP2'}
% all_region2 ={'ES', 'YS', 'SS', 'NWP2'}
% all_testname2 = {'ens02','test53','test55','test56'};
% all_testname2 = {'test11', 'test12'};
% all_testname2 = {'test11', 'test12'};
% RCM_testnames = {'test61', 'test62', 'test63', 'test64'};
RCM_testnames = {'test57', 'test58', 'test59', 'test60'};
% RCM_testnames = {'test65', 'test66', 'test67', 'test68'};
GCM_testnames = {'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'NorESM1-M', 'MPI-ESM-LR'};
% all_region2 ={'NWP','ES', 'SS', 'YS', 'AKP'}
% all_region2 ={'NWP','AKP2', 'YS', 'NES', 'SES'}

all_region2 ={'AKP4'};
% all_region2 ={'NWP'};


% scenname='rcp26';
scenname='rcp45';
% scenname='rcp85';

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

shadlev = [0 35];
rms_shadlev = [0 4];
%     trendlev = [-3 3];  %% trend lev
trendlev = [-10 10];  %% trend lev
abstrendlev =[2 7];
reltrendlev =[-5 5];
conlev  = 0:5:35;
meanplotlev =[-0.2 0.2];
meanplotlev2 =[-0.1 0.8] .*100.0;

trendplotlev = [0 7];
trenddifflev = [-10 10];
sshlev =[-0.3 0.3];
sshdifflev = [0 20];

% for snu_desktopd
% testname=all_testname2{testnameind2}    % % need to change
inputyear = [2006:2100]; % % put year which you want to plot [year year ...]
inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]

varname ='zeta'
variable='SSH'
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



%         load(['G:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
%             'ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);

% RCM_filedir = strcat('H:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
% GCM_filedir = strcat('E:\Data\Model\CMIP5\zos\', scenname, '\interp\', testname, '\'); % % where data files are

for testind=1:length(RCM_testnames)
    testname=RCM_testnames{testind};
    if (strcmp(scenname,'rcp26')==1)
        drivename='H';
    elseif (strcmp(scenname,'rcp45')==1)
        drivename='D';
    elseif (strcmp(scenname,'rcp85')==1)
        drivename='D';
    end
    filedir=strcat(drivename, ':\Data\Model\ROMS\nwp_1_20\', testname, '\run\');
    RCM_modelfilenames{testind} = strcat(filedir, testname,'_',regionname, '_ssh_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    RCM_cmemsfilenames{testind} = strcat(filedir, testname,'_',regionname, 'cmems_ssh_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    RCM_interpedfilenames{testind} = strcat(filedir, testname,'_',regionname, 'cmems_interped_ssh_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    RCM_movfilenames{testind} = strcat(filedir, testname,'_',regionname, 'moving_averaged_cmems_interped_ssh_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    
    testname=GCM_testnames{testind};
    filedir=strcat('D:\Data\Model\CMIP5\zos\', scenname, '\interp\', testname, '\');
    GCM_modelfilenames{testind} = strcat(filedir, testname,'_',regionname, '_ssh_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    GCM_cmemsfilenames{testind} = strcat(filedir, testname,'_',regionname, 'cmems_ssh_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    GCM_interpedfilenames{testind} = strcat(filedir, testname,'_',regionname, 'cmems_interped_ssh_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    GCM_movfilenames{testind} = strcat(filedir, testname,'_',regionname, 'moving_averaged_cmems_interped_ssh_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    
    cmemsdir='Z:\내 드라이브\Data\Observation\CMEMS\';
    cmems_filename = strcat(filedir, testname,'_',regionname, 'cmems_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');

end

% cmems_sla = ncread(cmems_filename, 'cmems_sla');

nyears=[1,2,5];
for testind=1:length(RCM_testnames)
    RCM_interped_sla(testind,:,:,:) = ncread(RCM_interpedfilenames{testind}, 'interped_sla');
    xlen=size(RCM_interped_sla,2); ylen = size(RCM_interped_sla,3); tlen = size(RCM_interped_sla,4);

    RCM_interped_sla_2d(testind,:,:)=reshape(RCM_interped_sla(testind,:,:,:), [xlen*ylen, tlen]);
    RCM_interped_sla_3d(testind,:,:)=reshape(RCM_interped_sla(testind,:,:),[xlen*ylen*12, tlen/12]);
    RCM_interped_sla_yearly_mean(testind,:)=mean(RCM_interped_sla_3d(testind,:,:),2,'omitnan') .* 100; % (m -> cm)

    RCM_interped_sla_mean(testind,:) = mean(RCM_interped_sla_2d(testind,:,:),2,'omitnan');

    RCM_interped_sla_4d(testind,:,:,:)=reshape(RCM_interped_sla(testind,:,:,:), [xlen*ylen, 12, tlen/12]);
    RCM_interped_sla_seasonal_mean=mean(mean(RCM_interped_sla_4d,2,'omitnan'),4,'omitnan');
    
    for t=1:tlen
    if mod(t,12)==0
        tt=12;
    else
        tt=mod(t,12);
    end
    RCM_interped_sla_seasonal_filtered(testind,t)=RCM_interped_sla_mean(testind,t)-RCM_interped_sla_seasonal_mean(testind,tt);
    end
    
    GCM_interped_sla(testind,:,:,:) = ncread(GCM_interpedfilenames{testind}, 'interped_sla');
    xlen=size(GCM_interped_sla,2); ylen = size(GCM_interped_sla,3); tlen = size(GCM_interped_sla,4);

    GCM_interped_sla_2d(testind,:,:)=reshape(GCM_interped_sla(testind,:,:,:), [xlen*ylen, tlen]);
    GCM_interped_sla_3d(testind,:,:)=reshape(GCM_interped_sla(testind,:,:),[xlen*ylen*12, tlen/12]);
    GCM_interped_sla_yearly_mean(testind,:)=mean(GCM_interped_sla_3d(testind,:,:),2,'omitnan') .* 100; % (m -> cm);
    
    GCM_interped_sla_mean(testind,:) = mean(GCM_interped_sla_2d(testind,:,:),2,'omitnan');

    GCM_interped_sla_4d(testind,:,:,:)=reshape(GCM_interped_sla(testind,:,:,:), [xlen*ylen, 12, tlen/12]);
    GCM_interped_sla_seasonal_mean=mean(mean(GCM_interped_sla_4d,2,'omitnan'),4,'omitnan');
    
    for t=1:tlen
    if mod(t,12)==0
        tt=12;
    else
        tt=mod(t,12);
    end
    GCM_interped_sla_seasonal_filtered(testind,t)=GCM_interped_sla_mean(testind,t)-GCM_interped_sla_seasonal_mean(testind,tt);
    end

    for nyeari=1:length(nyears)
        nyear=nyears(nyeari);
        nc_varname=['interped_', num2str(nyear), 'y_movmean'];
        eval([nc_varname, 's(testind,:,:,:)=squeeze(ncread(RCM_movfilenames{testind},', '''', nc_varname, '''', ',[1,1,1], [inf inf inf]));']);
    end
end


% % % % 
% % % % interped_2y_lowpass_sla = ncread(filename, 'interped_sla_2y_lowpass');
% % % % interped_5y_lowpass_sla = ncread(filename, 'interped_sla_5y_lowpass');
% % % % for t=1:tlen
% % % %     interped_sla(:,:,t)=interped_sla(:,:,t).*cmems_mask;
% % % %     interped_2y_lowpass_sla(:,:,t)=interped_2y_lowpass_sla(:,:,t).*cmems_mask;
% % % %     interped_5y_lowpass_sla(:,:,t)=interped_5y_lowpass_sla(:,:,t).*cmems_mask;
% % % % end
% % % % interped_sla_2d=reshape(interped_sla,[xlen*ylen, tlen]);
% % % % interped_sla_mean=mean(interped_sla_2d,1,'omitnan');
% % % % interped_2y_lowpass_sla_2d=reshape(interped_2y_lowpass_sla,[xlen*ylen, tlen]);
% % % % interped_2y_lowpass_sla_mean=mean(interped_2y_lowpass_sla_2d,1,'omitnan');
% % % % interped_5y_lowpass_sla_2d=reshape(interped_5y_lowpass_sla,[xlen*ylen, tlen]);
% % % % interped_5y_lowpass_sla_mean=mean(interped_5y_lowpass_sla_2d,1,'omitnan');
% % % % interped_sla_3d=reshape(interped_sla,[xlen*ylen*12, tlen/12]);
% % % % interped_sla_yearly_mean=mean(interped_sla_3d,1,'omitnan');


% interped_sla_4d=reshape(interped_sla,[xlen*ylen, 12, tlen/12]);
% interped_sla_seasonal_mean=mean(mean(interped_sla_4d,1,'omitnan'),3,'omitnan');
% for t=1:tlen
%     if mod(t,12)==0
%         tt=12;
%     else
%         tt=mod(t,12);
%     end
%     interped_sla_seasonal_filtered(t)=interped_sla_mean(t)-interped_sla_seasonal_mean(tt);
% end


% modelfilename = strcat(filedir, testname,'_',regionname, '_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
%         cmemsfilename = strcat(filedir, testname,'_',regionname, 'cmems_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
%         interpedfilename = strcat(filedir, testname,'_',regionname, 'cmems_interped_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
%         corrfilename = strcat(filedir, testname,'_',regionname, 'cmems_interped_ssh_corr_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
%         lpffilename = strcat(filedir, testname,'_',regionname, 'low_pass_filtered_cmems_interped_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
%         lpf_corrfilename = strcat(filedir, testname,'_',regionname, 'low_pass_filtered_cmems_interped_ssh_corr_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
%         movfilename = strcat(filedir, testname,'_',regionname, 'moving_averaged_cmems_interped_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
%         mov_corrfilename = strcat(filedir, testname,'_',regionname, 'moving_averaged_cmems_interped_ssh_corr_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
%         detrendfilename = strcat(filedir, testname,'_',regionname, 'detrended_cmems_interped_ssh_corr_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
% 
% filename = ['E:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
%             '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];
% filename2 = ['E:\Data\Model\ROMS\nwp_1_10\',testname2,'\run\',testname2,'_',regionname, ...
%             '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];
% 
%         
%         
% sodafilename =  ['E:\Data\Reanalysis\SODA\', sodatestname, '\',sodatestname,'_',regionname, ...
%             '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];
        
valnum=0;
run('C:\Users\User\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
wrmap = bwrmap(51:100,:);

if (strcmp(system_name,'PCWIN64'))
    % % for windows
    figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\5th_year\figure\all\',regionname,'\'); % % where figure files will be saved
%             param_script ='C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_EKB_RMS.m'
    param_script =['C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m']
%     filedir = strcat('E:\Data\Model\ROMS\nwp_1_10\', testname, '\run\'); % % where data files are
elseif (strcmp(system_name,'GLNXA64'))
end
%         std(mean(mean(comb_data_filtered,'omitnan'),'omitnan'))

run(param_script);

figdir=[figrawdir,'SSH\'];
if (exist(strcat(figdir) , 'dir') ~= 7)
    mkdir(strcat(figdir));
end 
outfile = strcat(figdir,regionname);


% cmems_trend=ncread(cmems_filename, 'cmems_trend');
% cmems_mask=ones(size(cmems_trend));
% cmems_mask(isnan(cmems_trend))=NaN;

% lon_rho=ncread(filename,'lon_rho');
% lat_rho=ncread(filename,'lat_rho');
% % % interped_sla = ncread(filename, 'interped_sla');
% % % interped_2y_lowpass_sla = ncread(filename, 'interped_sla_2y_lowpass');
% % % interped_5y_lowpass_sla = ncread(filename, 'interped_sla_5y_lowpass');
% % % xlen=size(interped_sla,1); ylen = size(interped_sla,2); tlen = size(interped_sla,3);
% % % for t=1:tlen
% % %     interped_sla(:,:,t)=interped_sla(:,:,t).*cmems_mask;
% % %     interped_2y_lowpass_sla(:,:,t)=interped_2y_lowpass_sla(:,:,t).*cmems_mask;
% % %     interped_5y_lowpass_sla(:,:,t)=interped_5y_lowpass_sla(:,:,t).*cmems_mask;
% % % end
% % % interped_sla_2d=reshape(interped_sla,[xlen*ylen, tlen]);
% % % interped_sla_mean=mean(interped_sla_2d,1,'omitnan');
% % % interped_2y_lowpass_sla_2d=reshape(interped_2y_lowpass_sla,[xlen*ylen, tlen]);
% % % interped_2y_lowpass_sla_mean=mean(interped_2y_lowpass_sla_2d,1,'omitnan');
% % % interped_5y_lowpass_sla_2d=reshape(interped_5y_lowpass_sla,[xlen*ylen, tlen]);
% % % interped_5y_lowpass_sla_mean=mean(interped_5y_lowpass_sla_2d,1,'omitnan');
% % % interped_sla_3d=reshape(interped_sla,[xlen*ylen*12, tlen/12]);
% % % interped_sla_yearly_mean=mean(interped_sla_3d,1,'omitnan');
% % % interped_sla_4d=reshape(interped_sla,[xlen*ylen, 12, tlen/12]);
% % % interped_sla_seasonal_mean=mean(mean(interped_sla_4d,1,'omitnan'),3,'omitnan');
% % % for t=1:tlen
% % %     if mod(t,12)==0
% % %         tt=12;
% % %     else
% % %         tt=mod(t,12);
% % %     end
% % %     interped_sla_seasonal_filtered(t)=interped_sla_mean(t)-interped_sla_seasonal_mean(tt);
% % % end
% % % 
% % % interped2_sla = ncread(filename2, 'corrected_interped_sla');
% % % interped2_2y_lowpass_sla = ncread(filename2, 'corrected_interped_sla_2y_lowpass');
% % % interped2_5y_lowpass_sla = ncread(filename2, 'corrected_interped_sla_5y_lowpass');
% % % xlen=size(interped2_sla,1); ylen = size(interped2_sla,2); tlen = size(interped2_sla,3);
% % % for t=1:tlen
% % %     interped2_sla(:,:,t)=interped2_sla(:,:,t).*cmems_mask;
% % %     interped2_2y_lowpass_sla(:,:,t)=interped2_2y_lowpass_sla(:,:,t).*cmems_mask;
% % %     interped2_5y_lowpass_sla(:,:,t)=interped2_5y_lowpass_sla(:,:,t).*cmems_mask;
% % % end
% % % interped2_sla_2d=reshape(interped2_sla,[xlen*ylen, tlen]);
% % % interped2_sla_mean=mean(interped2_sla_2d,1,'omitnan');
% % % interped2_2y_lowpass_sla_2d=reshape(interped2_2y_lowpass_sla,[xlen*ylen, tlen]);
% % % interped2_2y_lowpass_sla_mean=mean(interped2_2y_lowpass_sla_2d,1,'omitnan');
% % % interped2_5y_lowpass_sla_2d=reshape(interped2_5y_lowpass_sla,[xlen*ylen, tlen]);
% % % interped2_5y_lowpass_sla_mean=mean(interped2_5y_lowpass_sla_2d,1,'omitnan');
% % % interped2_sla_3d=reshape(interped2_sla,[xlen*ylen*12, tlen/12]);
% % % interped2_sla_yearly_mean=mean(interped2_sla_3d,1,'omitnan');
% % % interped2_sla_4d=reshape(interped2_sla,[xlen*ylen, 12, tlen/12]);
% % % interped2_sla_seasonal_mean=mean(mean(interped2_sla_4d,1,'omitnan'),3,'omitnan');
% % % for t=1:tlen
% % %     if mod(t,12)==0
% % %         tt=12;
% % %     else
% % %         tt=mod(t,12);
% % %     end
% % %     interped2_sla_seasonal_filtered(t)=interped2_sla_mean(t)-interped2_sla_seasonal_mean(tt);
% % % end
% % % 
% % % 
% % % 
% % % soda_sla = ncread(sodafilename, 'corrected_interped_sla');
% % % soda_2y_lowpass_sla = ncread(sodafilename, 'corrected_interped_sla_2y_lowpass');
% % % soda_5y_lowpass_sla = ncread(sodafilename, 'corrected_interped_sla_5y_lowpass');
% % % for t=1:tlen
% % %     soda_sla(:,:,t)=soda_sla(:,:,t).*cmems_mask;
% % %     soda_2y_lowpass_sla(:,:,t)=soda_2y_lowpass_sla(:,:,t).*cmems_mask;
% % %     soda_5y_lowpass_sla(:,:,t)=soda_5y_lowpass_sla(:,:,t).*cmems_mask;
% % % end
% % % soda_sla_2d=reshape(soda_sla,[xlen*ylen, tlen]);
% % % soda_sla_mean=mean(soda_sla_2d,1,'omitnan');
% % % soda_2y_lowpass_sla_2d=reshape(soda_2y_lowpass_sla,[xlen*ylen, tlen]);
% % % soda_2y_lowpass_sla_mean=mean(soda_2y_lowpass_sla_2d,1,'omitnan');
% % % soda_5y_lowpass_sla_2d=reshape(soda_5y_lowpass_sla,[xlen*ylen, tlen]);
% % % soda_5y_lowpass_sla_mean=mean(soda_5y_lowpass_sla_2d,1,'omitnan');
% % % soda_sla_3d=reshape(soda_sla,[xlen*ylen*12, tlen/12]);
% % % soda_sla_yearly_mean=mean(soda_sla_3d,1,'omitnan');
% % % soda_sla_4d=reshape(soda_sla,[xlen*ylen, 12, tlen/12]);
% % % soda_sla_seasonal_mean=mean(mean(soda_sla_4d,1,'omitnan'),3,'omitnan');
% % % for t=1:tlen
% % %     if mod(t,12)==0
% % %         tt=12;
% % %     else
% % %         tt=mod(t,12);
% % %     end
% % %     soda_sla_seasonal_filtered(t)=soda_sla_mean(t)-soda_sla_seasonal_mean(tt);
% % % end


% cmems_sla = ncread(RCM_cmemsfilenames{1}, 'cmems_sla');



% % for t=1:tlen
% %     cmems_sla(:,:,t)=cmems_sla(:,:,t).*cmems_mask;
% % 
% % end
% % cmems_sla_2d=reshape(cmems_sla,[xlen*ylen, tlen]);
% % cmems_sla_mean=mean(cmems_sla_2d,1,'omitnan');
% % 
% % cmems_sla_3d=reshape(cmems_sla,[xlen*ylen*12, tlen/12]);
% % cmems_sla_yearly_mean=mean(cmems_sla_3d,1,'omitnan');
% % cmems_sla_4d=reshape(cmems_sla,[xlen*ylen, 12, tlen/12]);
% % cmems_sla_seasonal_mean=mean(mean(cmems_sla_4d,1,'omitnan'),3,'omitnan');
% % for t=1:tlen
% %     if mod(t,12)==0
% %         tt=12;
% %     else
% %         tt=mod(t,12);
% %     end
% %     cmems_sla_seasonal_filtered(t)=cmems_sla_mean(t)-cmems_sla_seasonal_mean(tt);
% % end

% trend_filtered = ncread(filename,'trend_filtered');
% mean_trend_filtered = ncread(filename,'mean_trend_filtered');
% cmems_trend_filtered = ncread(filename,'cmems_trend_filtered');
% interped_trend_filtered = ncread(filename,'interped_trend_filtered_trend');



% start-------------------- make timedata for time series
        for i =1:length(inputyear) 
            tempyear=inputyear(i);
            for month=1:12
                xData((12*(i-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-01',]);
            end
            xData2(i) = datenum([num2str(tempyear),'-',num2str(6,'%02i'),'-30',]);
        end
% end-------------------- make timedata for time series  



% start-------------------- msl time series (yearly mean) (GCM)
for folding=1:1
        jpgname=strcat(outfile, '_', scenname, '_',regionname, '_GCM_msla_yearly_mean_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
                            
            hold on
            cmems_1y_mean=0;                        
            for testind=1:length(RCM_testnames)
                RCM_interped_1y_mean(testind)=mean(RCM_interped_sla_yearly_mean(testind,1:5));
                GCM_interped_1y_mean(testind)=mean(GCM_interped_sla_yearly_mean(testind,1:5));
                RCM_interped_sla_yearly_mean(testind,:)=RCM_interped_sla_yearly_mean(testind,:)+(cmems_1y_mean-RCM_interped_1y_mean(testind));
                GCM_interped_sla_yearly_mean(testind,:)=GCM_interped_sla_yearly_mean(testind,:)+(cmems_1y_mean-GCM_interped_1y_mean(testind));
%                 mslplot_all{testind}=plot(xData2,RCM_interped_sla_yearly_mean(testind,:),'r')
                mslplot2_all{testind}=plot(xData2,GCM_interped_sla_yearly_mean(testind,:),'b')
            end
%             mslplot=plot(xData2,cmems_sla_yearly_mean,'k')
            
            set(mslplot2_all{1},'Marker','*');
            set(mslplot2_all{2},'Marker','^');
            set(mslplot2_all{3},'Marker','o');
            set(mslplot2_all{4},'Marker','+');
            
            xlabel('Year')
            ylabel('Mean SSH (cm)')
            title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') '])
            datetick('x','yyyy','keepticks')
            axis tight;
            ylim(meanplotlev2)
%             set(mslplot_all{5},'LineWidth',2);
%             set(mslplot2_all{5},'LineWidth',2);
%             set(mslplot,'LineWidth',2);
            set(gca,'FontSize',15);
            grid on
            lgd=legend('GCM-IPSL-LR','GCM-IPSL-MR','GCM-Nor','GCM-MPI','SAT');
            set(lgd,'FontSize',15);
            set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
            set(lgd,'Orientation','horizontal');
            for tind=1:size(GCM_interped_sla_yearly_mean,2)
                model_std(tind)=std(GCM_interped_sla_yearly_mean(:,tind));
            end
            meanstd=mean(model_std);
            txt1=text(xData2(5), min(meanplotlev2)+diff(meanplotlev2)/16.0 ,['Mean std = ', num2str(round(meanstd,2))], 'FontSize', 20); 
%             constant_cor=corrcoef(GCM_interped_sla_yearly_mean(length(RCM_testnames),:),cmems_sla_yearly_mean);
%             txt1=text(xData2(10), min(meanplotlev2)+diff(meanplotlev2)/16.0 ,['R(GCM) = ', num2str(round(constant_cor(1,2),2))], 'FontSize', 20); 

            set(gcf,'PaperPosition', [0 0 36 12]) 
            hold off
            saveas(gcf,jpgname,'jpg');
            grid off
%         end
        close all;
end
% end-------------------- msl time series (yearly mean)(GCM)

% start-------------------- msl time series (yearly mean) (RCM)
for folding=1:1
        jpgname=strcat(outfile, '_', scenname, '_',regionname, '_RCM_msla_yearly_mean_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
%         if (exist(jpgname , 'file') ~= 2)
%             for varind=1:length(inputyear)*12
%                 msl_filt(varind)=squeeze(mean(mean(ncread(filename,'ssh_filtered',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
%             end
%     %         msl=squeeze(mean(mean(comb_data_filtered,1,'omitnan'),2,'omitnan'));
%             msl_filt=msl_filt-mean(msl_filt);    
%             interped_1y_mean=mean(interped_sla_seasonal_filtered(1:12*5));
%             interped2_1y_mean=mean(interped2_sla_seasonal_filtered(1:12*5));
%             soda_1y_mean=mean(soda_sla_seasonal_filtered(1:12*5));
%             cmems_1y_mean=mean(cmems_sla_seasonal_filtered(1:12*5));
                            
            hold on
            cmems_1y_mean=0;                        
            for testind=1:length(RCM_testnames)
                RCM_interped_1y_mean(testind)=mean(RCM_interped_sla_yearly_mean(testind,1:5));
                GCM_interped_1y_mean(testind)=mean(GCM_interped_sla_yearly_mean(testind,1:5));
                RCM_interped_sla_yearly_mean(testind,:)=RCM_interped_sla_yearly_mean(testind,:)+(cmems_1y_mean-RCM_interped_1y_mean(testind));
                GCM_interped_sla_yearly_mean(testind,:)=GCM_interped_sla_yearly_mean(testind,:)+(cmems_1y_mean-GCM_interped_1y_mean(testind));
                mslplot_all{testind}=plot(xData2,RCM_interped_sla_yearly_mean(testind,:),'r')
%                 mslplot2_all{testind}=plot(xData2,GCM_interped_sla_yearly_mean(testind,:),'b')
            end
%             mslplot=plot(xData2,cmems_sla_yearly_mean,'k')
            
            set(mslplot_all{1},'Marker','*');
            set(mslplot_all{2},'Marker','^');
            set(mslplot_all{3},'Marker','o');
            set(mslplot_all{4},'Marker','+');
            
            xlabel('Year')
            ylabel('Mean SSH (cm)')
            title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') '])
            datetick('x','yyyy','keepticks')
            axis tight;
            ylim(meanplotlev2)
%             set(mslplot_all{5},'LineWidth',2);
%             set(mslplot2_all{5},'LineWidth',2);
%             set(mslplot,'LineWidth',2);
            set(gca,'FontSize',15);
            grid on
            lgd=legend('RCM-IPSL-LR','RCM-IPSL-MR','RCM-Nor','RCM-MPI','SAT');
            set(lgd,'FontSize',15);
            set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
            set(lgd,'Orientation','horizontal');
            for tind=1:size(RCM_interped_sla_yearly_mean,2)
                model_std(tind)=std(RCM_interped_sla_yearly_mean(:,tind));
            end
            meanstd=mean(model_std);
            txt1=text(xData2(5), min(meanplotlev2)+diff(meanplotlev2)/16.0 ,['Mean std = ', num2str(round(meanstd,2))], 'FontSize', 20); 
%             constant_cor=corrcoef(GCM_interped_sla_yearly_mean(length(RCM_testnames),:),cmems_sla_yearly_mean);
%             txt1=text(xData2(10), min(meanplotlev2)+diff(meanplotlev2)/16.0 ,['R(GCM) = ', num2str(round(constant_cor(1,2),2))], 'FontSize', 20); 

            set(gcf,'PaperPosition', [0 0 36 12]) 
            hold off
            saveas(gcf,jpgname,'jpg');
            grid off
%         end
        close all;
end
% end-------------------- msl time series (yearly mean) (RCM)


