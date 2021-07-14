close all; clear all;  clc;
warning off;
% all_region2 ={'NWP','NWP2'}
% all_region2 ={'ES', 'YS', 'SS', 'NWP2'}
% all_testname2 = {'ens02','test53','test55','test56'};
% all_testname2 = {'test11', 'test12'};
% all_testname2 = {'test11', 'test12'};
% RCM_testnames = {'test61', 'test62', 'test63', 'test64'};
RCM_testnames = {'test65', 'test66', 'test67', 'test68'};
GCM_testnames = {'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'NorESM1-M', 'MPI-ESM-LR'};
% all_region2 ={'NWP','ES', 'SS', 'YS', 'AKP'}
% all_region2 ={'NWP','AKP2', 'YS', 'NES', 'SES'}

all_region2 ={'AKP4'};
% all_region2 ={'NWP'};



for regionind2=1:length(all_region2)
% scenname='rcp26';
scenname='rcp85';

close all;
clearvars '*' -except regionind2 testnameind2 all_region2 scenname ...
    all_testname2 testname testname2 sodatestname RCM_testnames GCM_testnames
% % % 
system_name=computer;

dropboxpath='C:\Users\User\Dropbox';
addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path


shadlev = [0 35];
rms_shadlev = [0 4];
%     trendlev = [-3 3];  %% trend lev
trendlev = [-10 10];  %% trend lev
abstrendlev =[2 7];
reltrendlev =[-5 5];
conlev  = 0:5:35;
meanplotlev =[-0.2 0.2];
meanplotlev2 =[19 23];

trendplotlev = [0 7];
trenddifflev = [-10 10];
sstlev =[-0.3 0.3];
sstdifflev = [0 20];

% for snu_desktopd
% testname=all_testname2{testnameind2}    % % need to change
inputyear = [2006:2060]; % % put year which you want to plot [year year ...]
inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]

% varname ='zeta';
variable='SST';
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

% % %         flag configuration
for folding=1:1
    fig_flags{1,1}='GCM time series';
    fig_flags{2,1}='RCM time series';


    for flagi=1:10
        fig_flags{flagi,2}=0;
    end
    fig_flags{1,2}=1;
    fig_flags{2,2}=0;
end

%         load(['G:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
%             'sst_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);

% RCM_filedir = strcat('H:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
% GCM_filedir = strcat('E:\Data\Model\CMIP5\zos\', scenname, '\interp\', testname, '\'); % % where data files are

for testind=1:length(RCM_testnames)
    testname=RCM_testnames{testind};
    if (strcmp(scenname,'rcp26')==1)
        drivename='H';
    elseif (strcmp(scenname,'rcp85')==1)
        drivename='G';
    end
    filedir=strcat(drivename, ':\Data\Model\ROMS\nwp_1_20\', testname, '\run\');
    RCM_modelfilenames{testind} = strcat(filedir, testname,'_',regionname, '_sst_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    RCM_cmemsfilenames{testind} = strcat(filedir, testname,'_',regionname, 'cmems_sst_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    RCM_interpedfilenames{testind} = strcat(filedir, testname,'_',regionname, 'cmems_interped_sst_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
%     RCM_movfilenames{testind} = strcat(filedir, testname,'_',regionname, 'moving_averaged_cmems_interped_sst_trend_', ...
%         num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    
    testname=GCM_testnames{testind};
%     filedir=strcat('E:\Data\Model\CMIP5\zos\', scenname, '\interp\', testname, '\');
    filedir=strcat('E:\Data\Model\CMIP5\surface\',scenname,'_');

    GCM_modelfilenames{testind} = strcat(filedir, testname,'_',regionname, '_sst_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    GCM_cmemsfilenames{testind} = strcat(filedir, testname,'_',regionname, 'cmems_sst_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    GCM_interpedfilenames{testind} = strcat(filedir, testname,'_',regionname, 'cmems_interped_sst_trend_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
%     GCM_movfilenames{testind} = strcat(filedir, testname,'_',regionname, 'moving_averaged_cmems_interped_sst_trend_', ...
%         num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
    
    cmemsdir='E:\Data\Observation\CMEMS\';
    cmems_filename = strcat(filedir, testname,'_',regionname, 'cmems_sst_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
end

% msst1=ncread('E:\Data\Model\CMIP5\surface\rcp26_IPSL-CM5A-LR_NWP_CMIP5_sst_analysis_2006_2060.nc','raw_sst');
% msst2=ncread('E:\Data\Model\CMIP5\surface\rcp85_IPSL-CM5A-LR_NWP_CMIP5_sst_analysis_2006_2060.nc','raw_sst');
% plot(squeeze(mean(mean(msst1,1,'omitnan'),2,'omitnan')))
% hold on
% plot(squeeze(mean(mean(msst2,1,'omitnan'),2,'omitnan')))
% hold off

% cmems_sst = ncread(cmems_filename, 'cmems_sst');

nyears=[1,2,5];
for testind=1:length(RCM_testnames)
    RCM_interped_sst(testind,:,:,:) = ncread(RCM_interpedfilenames{testind}, 'interped_sst');
    xlen=size(RCM_interped_sst,2); ylen = size(RCM_interped_sst,3); tlen = size(RCM_interped_sst,4);

    RCM_interped_sst_2d(testind,:,:)=reshape(RCM_interped_sst(testind,:,:,:), [xlen*ylen, tlen]);
    RCM_interped_sst_3d(testind,:,:)=reshape(RCM_interped_sst(testind,:,:),[xlen*ylen*12, tlen/12]);
    RCM_interped_sst_yearly_mean(testind,:)=mean(RCM_interped_sst_3d(testind,:,:),2,'omitnan');

    RCM_interped_sst_mean(testind,:) = mean(RCM_interped_sst_2d(testind,:,:),2,'omitnan');

    RCM_interped_sst_4d(testind,:,:,:)=reshape(RCM_interped_sst(testind,:,:,:), [xlen*ylen, 12, tlen/12]);
    RCM_interped_sst_seasonal_mean=mean(mean(RCM_interped_sst_4d,2,'omitnan'),4,'omitnan');
    
    for t=1:tlen
    if mod(t,12)==0
        tt=12;
    else
        tt=mod(t,12);
    end
    RCM_interped_sst_seasonal_filtered(testind,t)=RCM_interped_sst_mean(testind,t)-RCM_interped_sst_seasonal_mean(testind,tt);
    end
    
    GCM_interped_sst(testind,:,:,:) = ncread(GCM_interpedfilenames{testind}, 'interped_sst') -273.15;
    xlen=size(GCM_interped_sst,2); ylen = size(GCM_interped_sst,3); tlen = size(GCM_interped_sst,4);

    GCM_interped_sst_2d(testind,:,:)=reshape(GCM_interped_sst(testind,:,:,:), [xlen*ylen, tlen]);
    GCM_interped_sst_3d(testind,:,:)=reshape(GCM_interped_sst(testind,:,:),[xlen*ylen*12, tlen/12]);
    GCM_interped_sst_yearly_mean(testind,:)=mean(GCM_interped_sst_3d(testind,:,:),2,'omitnan'); % (m -> cm);
    
    GCM_interped_sst_mean(testind,:) = mean(GCM_interped_sst_2d(testind,:,:),2,'omitnan');

    GCM_interped_sst_4d(testind,:,:,:)=reshape(GCM_interped_sst(testind,:,:,:), [xlen*ylen, 12, tlen/12]);
    GCM_interped_sst_seasonal_mean=mean(mean(GCM_interped_sst_4d,2,'omitnan'),4,'omitnan');
    
    for t=1:tlen
    if mod(t,12)==0
        tt=12;
    else
        tt=mod(t,12);
    end
    GCM_interped_sst_seasonal_filtered(testind,t)=GCM_interped_sst_mean(testind,t)-GCM_interped_sst_seasonal_mean(testind,tt);
    end

%     for nyeari=1:length(nyears)
%         nyear=nyears(nyeari);
%         nc_varname=['interped_', num2str(nyear), 'y_movmean'];
%         eval([nc_varname, 's(testind,:,:,:)=squeeze(ncread(RCM_movfilenames{testind},', '''', nc_varname, '''', ',[1,1,1], [inf inf inf]));']);
%     end
end

        
valnum=0;
run('C:\Users\User\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
wrmap = bwrmap(51:100,:);

if (strcmp(system_name,'PCWIN64'))
    % % for windows
    figrawdir =strcat('F:\OneDrive - 서울대학교\MEPL\project\SSH\5th_year\figure\all\',regionname,'\'); % % where figure files will be saved
%             param_script ='C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_EKB_RMS.m'
    param_script =['C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m']
%     filedir = strcat('E:\Data\Model\ROMS\nwp_1_10\', testname, '\run\'); % % where data files are
elseif (strcmp(system_name,'GLNXA64'))
end
%         std(mean(mean(comb_data_filtered,'omitnan'),'omitnan'))

run(param_script);

figdir=[figrawdir,'sst\'];
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



% start-------------------- msl time series (yearly mean) (GCM)
for folding=1:1
        jpgname=strcat(outfile, '_', scenname, '_',regionname, '_GCM_msst_yearly_mean_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
%         if (exist(jpgname , 'file') ~= 2)
%             for varind=1:length(inputyear)*12
%                 msl_filt(varind)=squeeze(mean(mean(ncread(filename,'sst_filtered',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
%             end
%     %         msl=squeeze(mean(mean(comb_data_filtered,1,'omitnan'),2,'omitnan'));
%             msl_filt=msl_filt-mean(msl_filt);    
%             interped_1y_mean=mean(interped_sst_seasonal_filtered(1:12*5));
%             interped2_1y_mean=mean(interped2_sst_seasonal_filtered(1:12*5));
%             soda_1y_mean=mean(soda_sst_seasonal_filtered(1:12*5));
%             cmems_1y_mean=mean(cmems_sst_seasonal_filtered(1:12*5));
                            
            hold on
            cmems_1y_mean=0;                        
            for testind=1:length(RCM_testnames)
%                 RCM_interped_1y_mean(testind)=mean(RCM_interped_sst_yearly_mean(testind,1:5));
%                 GCM_interped_1y_mean(testind)=mean(GCM_interped_sst_yearly_mean(testind,1:5));
%                 RCM_interped_sst_yearly_mean(testind,:)=RCM_interped_sst_yearly_mean(testind,:)+(cmems_1y_mean-RCM_interped_1y_mean(testind));
%                 GCM_interped_sst_yearly_mean(testind,:)=GCM_interped_sst_yearly_mean(testind,:)+(cmems_1y_mean-GCM_interped_1y_mean(testind));
%                 mslplot_all{testind}=plot(xData2,RCM_interped_sst_yearly_mean(testind,:),'r')
                mslplot2_all{testind}=plot(xData2,GCM_interped_sst_yearly_mean(testind,:),'b')
            end
%             mslplot=plot(xData2,cmems_sst_yearly_mean,'k')
            
            set(mslplot2_all{1},'Marker','*');
            set(mslplot2_all{2},'Marker','^');
            set(mslplot2_all{3},'Marker','o');
            set(mslplot2_all{4},'Marker','+');
            
            xlabel('Year')
            ylabel('Mean SST (^oC)')
            title([regionname, ', Mean sst(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') '])
            datetick('x','yyyy','keepticks')
            axis tight;
            ylim(colorbar_lev_meanplot)
%             set(mslplot_all{5},'LineWidth',2);
%             set(mslplot2_all{5},'LineWidth',2);
%             set(mslplot,'LineWidth',2);
            set(gca,'FontSize',15);
            grid on
            lgd=legend('GCM-IPSL-LR','GCM-IPSL-MR','GCM-Nor','GCM-MPI');
            set(lgd,'FontSize',15);
            set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
            set(lgd,'Orientation','horizontal');
            for tind=1:size(GCM_interped_sst_yearly_mean,2)
                model_std(tind)=std(GCM_interped_sst_yearly_mean(:,tind));
            end
            meanstd=mean(model_std);
            txt1=text(xData2(5), min(colorbar_lev_meanplot)+diff(colorbar_lev_meanplot)/16.0 ,['Mean std = ', num2str(round(meanstd,2))], 'FontSize', 20); 
%             constant_cor=corrcoef(GCM_interped_sst_yearly_mean(length(RCM_testnames),:),cmems_sst_yearly_mean);
%             txt1=text(xData2(10), min(meanplotlev2)+diff(meanplotlev2)/16.0 ,['R(GCM) = ', num2str(round(constant_cor(1,2),2))], 'FontSize', 20); 

            set(gcf,'PaperPosition', [0 0 36 12]) 
            hold off
            saveas(gcf,jpgname,'jpg');
            grid off
%         end
        close all;
end
% end-------------------- msl time series (yearly mean)

% start-------------------- msl time series (yearly mean) (RCM)
for folding=1:1
        jpgname=strcat(outfile, '_', scenname, '_',regionname, '_RCM_msst_yearly_mean_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
%         if (exist(jpgname , 'file') ~= 2)
% % % % %             for varind=1:length(inputyear)*12
% % % % %                 msl_filt(varind)=squeeze(mean(mean(ncread(filename,'sst_filtered',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
% % % % %             end
% % % % %     %         msl=squeeze(mean(mean(comb_data_filtered,1,'omitnan'),2,'omitnan'));
% % % % %             msl_filt=msl_filt-mean(msl_filt);    
% % % % %             interped_1y_mean=mean(interped_sst_seasonal_filtered(1:12*5));
% % % % %             interped2_1y_mean=mean(interped2_sst_seasonal_filtered(1:12*5));
% % % % %             soda_1y_mean=mean(soda_sst_seasonal_filtered(1:12*5));
% % % % %             cmems_1y_mean=mean(cmems_sst_seasonal_filtered(1:12*5));
                            
            hold on
            cmems_1y_mean=0;                        
            for testind=1:length(RCM_testnames)
%                 RCM_interped_1y_mean(testind)=mean(RCM_interped_sst_yearly_mean(testind,1:5));
%                 GCM_interped_1y_mean(testind)=mean(GCM_interped_sst_yearly_mean(testind,1:5));
                mslplot_all{testind}=plot(xData2,RCM_interped_sst_yearly_mean(testind,:),'r')
            end
%             mslplot=plot(xData2,cmems_sst_yearly_mean,'k')
            
            set(mslplot_all{1},'Marker','*');
            set(mslplot_all{2},'Marker','^');
            set(mslplot_all{3},'Marker','o');
            set(mslplot_all{4},'Marker','+');
            
            xlabel('Year')
            ylabel('Mean SST (^oC)')
            title([regionname, ', Mean sst(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') '])
            datetick('x','yyyy','keepticks')
            axis tight;
            ylim(colorbar_lev_meanplot)

%             set(mslplot_all{5},'LineWidth',2);
%             set(mslplot2_all{5},'LineWidth',2);
%             set(mslplot,'LineWidth',2);
            set(gca,'FontSize',15);
            grid on
            lgd=legend('RCM-IPSL-LR','RCM-IPSL-MR','RCM-Nor','RCM-MPI');
            set(lgd,'FontSize',15);
            set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
            set(lgd,'Orientation','horizontal');
            for tind=1:size(RCM_interped_sst_yearly_mean,2)
                model_std(tind)=std(RCM_interped_sst_yearly_mean(:,tind));
            end
            meanstd=mean(model_std);
            txt1=text(xData2(5), min(colorbar_lev_meanplot)+diff(colorbar_lev_meanplot)/16.0 ,['Mean std = ', num2str(round(meanstd,2))], 'FontSize', 20); 
%             constant_cor=corrcoef(GCM_interped_sst_yearly_mean(length(RCM_testnames),:),cmems_sst_yearly_mean);
%             txt1=text(xData2(10), min(meanplotlev2)+diff(meanplotlev2)/16.0 ,['R(GCM) = ', num2str(round(constant_cor(1,2),2))], 'FontSize', 20); 

            set(gcf,'PaperPosition', [0 0 36 12]) 
            hold off
            saveas(gcf,jpgname,'jpg');
            grid off
%         end
        close all;
end
% end-------------------- msl time series (yearly mean)

end
