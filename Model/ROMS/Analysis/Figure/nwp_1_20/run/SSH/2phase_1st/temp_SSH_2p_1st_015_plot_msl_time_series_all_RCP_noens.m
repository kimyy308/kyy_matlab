close all; clear all;  clc;
warning off;

RCM_testnames = {'test2102', 'test2103', 'test2104', 'test2105', 'test2106'};
GCM_testnames = {'CNRM-ESM2-1', 'EC-Earth3-Veg', 'ACCESS-CM2', 'CNRM-CM6-1-HR', 'CMCC-ESM2'};
% all_region2 ={'NWP','ES', 'SS', 'YS', 'AKP'}
% all_region2 ={'NWP','AKP2', 'YS', 'NES', 'SES'}

all_region2 ={'AKP4'};
% all_region2 ={'NWP'};

scenname='historical';
% scenname='rcp26';
% scenname='rcp45';
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
meanplotlev2 =[-0.05 0.15] .*100.0;

trendplotlev = [0 7];
trenddifflev = [-10 10];
sshlev =[-0.3 0.3];
sshdifflev = [0 20];

% for snu_desktopd
% testname=all_testname2{testnameind2}    % % need to change
inputyear = [1993:2014]; % % put year which you want to plot [year year ...]
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



%         load(['G:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
%             'ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);

% RCM_filedir = strcat('H:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
% GCM_filedir = strcat('E:\Data\Model\CMIP5\zos\', scenname, '\interp\', testname, '\'); % % where data files are

for testind=1:length(RCM_testnames)
    testname=RCM_testnames{testind};
    if (strcmp(scenname,'rcp26')==1)
        drivename='D';
    elseif (strcmp(scenname,'rcp45')==1)
        drivename='D';
    elseif (strcmp(scenname,'rcp85')==1)
        drivename='D';
    elseif (strcmp(scenname,'historical')==1)
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
    filedir=strcat('D:\Data\Model\CMIP6\zos\', scenname, '\interp\', testname, '\');
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

lon_cmems= ncread(RCM_cmemsfilenames{1}, 'lon_cmems');
lat_cmems= ncread(RCM_cmemsfilenames{1}, 'lat_cmems');
nyears=[1,2,5];
for testind=1:length(RCM_testnames)

    RCM_interped_sla(testind,:,:,:) = ncread(RCM_interpedfilenames{testind}, 'interped_sla');
    xlen=size(RCM_interped_sla,2); ylen = size(RCM_interped_sla,3); tlen = size(RCM_interped_sla,4);

    RCM_interped_sla_2d(testind,:,:)=reshape(RCM_interped_sla(testind,:,:,:), [xlen*ylen, tlen]);
    RCM_interped_sla_3d(testind,:,:,:,:)=reshape(RCM_interped_sla(testind,:,:),[xlen, ylen, 12, tlen/12]);
    RCM_interped_sla_yearly(testind,:,:,:)=mean(RCM_interped_sla_3d(testind,:,:,:,:), 4);
    for ti=1:size(RCM_interped_sla_yearly,4)
        [RCM_interped_sla_yearly_mean(testind,ti), error_status] = ...
            Func_0011_get_area_weighted_mean(squeeze(RCM_interped_sla_yearly(testind,:,:,ti)), lon_cmems, lat_cmems); %(m -> cm)
        RCM_interped_sla_yearly_mean(testind,ti)= RCM_interped_sla_yearly_mean(testind,ti).*100;
    end
%     RCM_interped_sla_yearly_mean(testind,:)=mean(RCM_interped_sla_3d(testind,:,:),2,'omitnan') .* 100; % (m -> cm)

    RCM_interped_sla_mean(testind,:) = mean(RCM_interped_sla_2d(testind,:,:),2,'omitnan');

    RCM_interped_sla_4d(testind,:,:,:)=reshape(RCM_interped_sla(testind,:,:,:), [xlen*ylen, 12, tlen/12]);
    RCM_interped_sla_seasonal_mean=mean(mean(RCM_interped_sla_4d,2,'omitnan'),4,'omitnan');

    CMEMS_interped_sla(:,:,:) = ncread(RCM_cmemsfilenames{testind}, 'cmems_sla');
    CMEMS_interped_sla_2d(:,:)=reshape(CMEMS_interped_sla(:,:,:), [xlen*ylen, tlen]);
    CMEMS_interped_sla_3d(:,:,:,:)=reshape(CMEMS_interped_sla(:,:),[xlen ylen 12, tlen/12]);
    CMEMS_interped_sla_yearly(:,:,:)=mean(CMEMS_interped_sla_3d,3);
    for ti=1:size(RCM_interped_sla_yearly,4)
        [CMEMS_interped_sla_yearly_mean(ti), error_status] = ...
            Func_0011_get_area_weighted_mean(squeeze(CMEMS_interped_sla_yearly(:,:,ti)), lon_cmems, lat_cmems); %(m -> cm)
        CMEMS_interped_sla_yearly_mean(ti)= CMEMS_interped_sla_yearly_mean(ti).*100;
    end
%     CMEMS_interped_sla_yearly_mean(:)=mean(CMEMS_interped_sla_3d(:,:),1,'omitnan') .* 100; % (m -> cm)

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
%     GCM_interped_sla_3d(testind,:,:)=reshape(GCM_interped_sla(testind,:,:),[xlen*ylen*12, tlen/12]);
    GCM_interped_sla_3d(testind,:,:,:,:)=reshape(GCM_interped_sla(testind,:,:),[xlen, ylen, 12, tlen/12]);
    GCM_interped_sla_yearly(testind,:,:,:)=mean(GCM_interped_sla_3d(testind,:,:,:,:), 4);
    for ti=1:size(RCM_interped_sla_yearly,4)
        [GCM_interped_sla_yearly_mean(testind,ti), error_status] = ...
            Func_0011_get_area_weighted_mean(squeeze(GCM_interped_sla_yearly(testind,:,:,ti)), lon_cmems, lat_cmems); %(m -> cm)
        GCM_interped_sla_yearly_mean(testind,ti)= GCM_interped_sla_yearly_mean(testind,ti).*100;
    end
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

% % %     for nyeari=1:length(nyears)
% % %         nyear=nyears(nyeari);
% % %         nc_varname=['interped_', num2str(nyear), 'y_movmean'];
% % %         eval([nc_varname, 's(testind,:,:,:)=squeeze(ncread(RCM_movfilenames{testind},', '''', nc_varname, '''', ',[1,1,1], [inf inf inf]));']);
% % %     end
end
        
valnum=0;
run('C:\Users\User\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
wrmap = bwrmap(51:100,:);

if (strcmp(system_name,'PCWIN64'))
    % % for windows
    figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\6th_year\figure\all\',regionname,'\'); % % where figure files will be saved
%             param_script ='C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_EKB_RMS.m'
    param_script =['C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m'];
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

% start-------------------- make timedata for time series
        for i =1:length(inputyear) 
            tempyear=inputyear(i);
            for month=1:12
                xData((12*(i-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-01',]);
            end
            xData2(i) = datenum([num2str(tempyear),'-',num2str(6,'%02i'),'-30',]);
        end
% end-------------------- make timedata for time series  



correction_period=1; % 1 or 1:5
% start-------------------- msl time series (yearly mean) (GCM)
for folding=1:1
        jpgname=strcat(outfile, '_', scenname, '_',regionname, '_GCM_msla_yearly_mean_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
                            
            hold on
            cmems_1y_mean=0;       
            CMEMS_interped_1y_mean=mean(CMEMS_interped_sla_yearly_mean(correction_period));
            CMEMS_interped_sla_yearly_mean(:)=CMEMS_interped_sla_yearly_mean(:)-CMEMS_interped_1y_mean;
            for testind=1:length(RCM_testnames)
                RCM_interped_1y_mean(testind)=mean(RCM_interped_sla_yearly_mean(testind,correction_period));
                GCM_interped_1y_mean(testind)=mean(GCM_interped_sla_yearly_mean(testind,correction_period));
                RCM_interped_sla_yearly_mean(testind,:)=RCM_interped_sla_yearly_mean(testind,:)+(cmems_1y_mean-RCM_interped_1y_mean(testind));
                GCM_interped_sla_yearly_mean(testind,:)=GCM_interped_sla_yearly_mean(testind,:)+(cmems_1y_mean-GCM_interped_1y_mean(testind));
%                 mslplot_all{testind}=plot(xData2,RCM_interped_sla_yearly_mean(testind,:),'r')
                mslplot2_all{testind}=plot(xData2,GCM_interped_sla_yearly_mean(testind,:),'b')
                [coef_GCM(testind), pval_GCM(testind)] = corr(GCM_interped_sla_yearly_mean(testind,:)', CMEMS_interped_sla_yearly_mean(:));
                [coef_kendall_GCM(testind), pval_kendall_GCM(testind)] = corr(GCM_interped_sla_yearly_mean(testind,:)', CMEMS_interped_sla_yearly_mean(:), 'type', 'kendall');
            end
            
            mslplot2_all{testind+1}=plot(xData2,CMEMS_interped_sla_yearly_mean(:),'k')
            
%             mslplot=plot(xData2,cmems_sla_yearly_mean,'k')
            
            set(mslplot2_all{1},'Marker','*');
            set(mslplot2_all{2},'Marker','^');
            set(mslplot2_all{3},'Marker','o');
            set(mslplot2_all{4},'Marker','+');
            set(mslplot2_all{5},'Marker','square');
            set(mslplot2_all{6},'Marker','pentagram');
            set(mslplot2_all{6},'LineWidth',2);

            xlabel('Year')
            ylabel('Mean SSH (cm)')
            title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') '])
            datetick('x','yyyy','keepticks')
            axis tight;
            ylim(meanplotlev2)
%             set(mslplot_all{5},'LineWidth',2);
%             set(mslplot2_all{5},'LineWidth',2);
%             set(mslplot,'LineWidth',2);
            set(gca,'FontSize',20);
            grid on
            lgd=legend({'GCM-CNRM'; 'GCM-EC-Veg'; 'GCM-ACC'; 'GCM-CNRM-HR'; 'GCM-CMCC'; 'CMEMS'}, 'NumColumns',2);
%             lgd=legend({'GCM-CNRM'; 'GCM-EC-Veg'; 'GCM-ACC'; 'GCM-CNRM-HR', 'GCM-CMCC'});
            set(lgd,'FontSize',20);
            set(lgd,'Position',[0.03, 0.84, 0.775, 0.03]);
            set(lgd,'Orientation','horizontal');
            for tind=1:size(GCM_interped_sla_yearly_mean,2)
                model_std(tind)=std(GCM_interped_sla_yearly_mean(:,tind));
            end
            meanstd=mean(model_std);
%             txt1=text(xData2(5), min(meanplotlev2)+diff(meanplotlev2)/16.0 ,['Mean std = ', num2str(round(meanstd,2))], 'FontSize', 20); 

            %             constant_cor=corrcoef(GCM_interped_sla_yearly_mean(length(RCM_testnames),:),cmems_sla_yearly_mean);
%             txt1=text(xData2(10), min(meanplotlev2)+diff(meanplotlev2)/16.0 ,['R(GCM) = ', num2str(round(constant_cor(1,2),2))], 'FontSize', 20); 

            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver 22 17]) 
            hold off
            saveas(gcf,jpgname,'jpg'); RemoveWhiteSpace([], 'file', jpgname);

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
            CMEMS_interped_1y_mean=mean(CMEMS_interped_sla_yearly_mean(correction_period));
            CMEMS_interped_sla_yearly_mean(:)=CMEMS_interped_sla_yearly_mean(:)-CMEMS_interped_1y_mean;
            for testind=1:length(RCM_testnames)
                RCM_interped_1y_mean(testind)=mean(RCM_interped_sla_yearly_mean(testind,correction_period));
                GCM_interped_1y_mean(testind)=mean(GCM_interped_sla_yearly_mean(testind,correction_period));
                RCM_interped_sla_yearly_mean(testind,:)=RCM_interped_sla_yearly_mean(testind,:)+(cmems_1y_mean-RCM_interped_1y_mean(testind));
                GCM_interped_sla_yearly_mean(testind,:)=GCM_interped_sla_yearly_mean(testind,:)+(cmems_1y_mean-GCM_interped_1y_mean(testind));
                mslplot_all{testind}=plot(xData2,RCM_interped_sla_yearly_mean(testind,:),'r')
%                 mslplot2_all{testind}=plot(xData2,GCM_interped_sla_yearly_mean(testind,:),'b')
%                 [coef(testind), pval(testind)] = corr(GCM_interped_sla_yearly_mean(testind,:)', RCM_interped_sla_yearly_mean(testind,:)');
                [coef_RCM(testind), pval_RCM(testind)] = corr(RCM_interped_sla_yearly_mean(testind,:)', CMEMS_interped_sla_yearly_mean(:));
                [coef_kendall_RCM(testind), pval_kendall_RCM(testind)] = corr(RCM_interped_sla_yearly_mean(testind,:)', CMEMS_interped_sla_yearly_mean(:), 'type', 'kendall');
            end
            
            mslplot_all{testind+1}=plot(xData2,CMEMS_interped_sla_yearly_mean(:),'k')
            
            set(mslplot_all{1},'Marker','*');
            set(mslplot_all{2},'Marker','^');
            set(mslplot_all{3},'Marker','o');
            set(mslplot_all{4},'Marker','+');
            set(mslplot_all{5},'Marker','square');
            set(mslplot_all{6},'Marker','pentagram');
            set(mslplot_all{6},'LineWidth',2);
            
            xlabel('Year')
            ylabel('Mean SSH (cm)')
            title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') '])
            datetick('x','yyyy','keepticks')
            axis tight;
            ylim(meanplotlev2)
%             set(mslplot_all{5},'LineWidth',2);
%             set(mslplot2_all{5},'LineWidth',2);
%             set(mslplot,'LineWidth',2);
            set(gca,'FontSize',20);
            grid on
            lgd=legend({'RCM-CNRM'; 'RCM-EC-Veg'; 'RCM-ACC'; 'RCM-CNRM-HR'; 'RCM-CMCC'; 'CMEMS'}, 'NumColumns',2);
            set(lgd,'FontSize',20);
            set(lgd,'Position',[0.03, 0.84, 0.775, 0.03]);
            set(lgd,'Orientation','horizontal');
            for tind=1:size(RCM_interped_sla_yearly_mean,2)
                model_std(tind)=std(RCM_interped_sla_yearly_mean(:,tind));
            end
            meanstd=mean(model_std);
%             txt1=text(xData2(5), min(meanplotlev2)+diff(meanplotlev2)/16.0 ,['Mean std = ', num2str(round(meanstd,2))], 'FontSize', 20); 

%             constant_cor=corrcoef(GCM_interped_sla_yearly_mean(length(RCM_testnames),:),cmems_sla_yearly_mean);
%             txt1=text(xData2(10), min(meanplotlev2)+diff(meanplotlev2)/16.0 ,['R(GCM) = ', num2str(round(constant_cor(1,2),2))], 'FontSize', 20); 

            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver 22 17]) 
            hold off
            saveas(gcf,jpgname,'jpg'); RemoveWhiteSpace([], 'file', jpgname);

            grid off
%         end
        close all;
end
% end-------------------- msl time series (yearly mean) (RCM)


