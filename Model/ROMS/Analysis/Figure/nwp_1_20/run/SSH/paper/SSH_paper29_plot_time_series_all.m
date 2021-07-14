close all; clear all;  clc;
warning off;
% all_region2 ={'NWP','NWP2'}
% all_region2 ={'ES', 'YS', 'SS', 'NWP2'}
% all_testname2 = {'ens02','test53','test55','test56'};
% all_testname2 = {'test11', 'test12'};
% all_testname2 = {'test11', 'test12'};
testname = 'test11';
testname2 = 'test12';
sodatestname = 'SODA_3_4_2';
% all_region2 ={'NWP','ES', 'SS', 'YS', 'AKP'}
% all_region2 ={'NWP','AKP2', 'YS', 'NES', 'SES'}

all_region2 ={'AKP4'};
% all_region2 ={'NWP'};


close all;
clearvars '*' -except regionind2 testnameind2 all_region2 all_testname2 testname testname2 sodatestname
% % % 
% % % Read Model SST
% % % interp
% % % get RMS
% % % get BIAS
system_name=computer;
if (strcmp(system_name,'PCWIN64'))
    % % for windows
    dropboxpath='C:\Users\KYY\Dropbox';
    addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
    addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
    addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
    addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path
elseif (strcmp(system_name,'GLNXA64'))
    dropboxpath='/home/kimyy/Dropbox';
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
    addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
end

regionind2 = 1;

shadlev = [0 35];
rms_shadlev = [0 4];
%     trendlev = [-3 3];  %% trend lev
trendlev = [-10 10];  %% trend lev
abstrendlev =[2 7];
reltrendlev =[-5 5];
conlev  = 0:5:35;
meanplotlev =[-0.2 0.2];
trendplotlev = [0 7];
trenddifflev = [-10 10];
sshlev =[-0.3 0.3];
sshdifflev = [0 20];

% for snu_desktopd
% testname=all_testname2{testnameind2}    % % need to change
inputyear = [1994:2017]; % % put year which you want to plot [year year ...]
inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]

varname ='zeta'
variable='SSH'
run('nwp_polygon_point.m');
regionname=all_region2{regionind2};
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

% % % for EKB
% regionname='EKB';
% lonlat = [127, 129.5, 38, 40.5];

%         load(['G:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
%             'ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
filename = ['E:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
            '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];
filename2 = ['E:\Data\Model\ROMS\nwp_1_10\',testname2,'\run\',testname2,'_',regionname, ...
            '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];

        
        
sodafilename =  ['E:\Data\Reanalysis\SODA\', sodatestname, '\',sodatestname,'_',regionname, ...
            '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];


        
valnum=0;
run('C:\Users\kyy\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
wrmap = bwrmap(51:100,:);


if (strcmp(system_name,'PCWIN64'))
    % % for windows
    figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\4th_year\figure\all\',testname,'\',regionname,'\'); % % where figure files will be saved
%             param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_EKB_RMS.m'
    param_script =['C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_', regionname, '.m']
    filedir = strcat('E:\Data\Model\ROMS\nwp_1_10\', testname, '\run\'); % % where data files are
elseif (strcmp(system_name,'GLNXA64'))
end
%         std(mean(mean(comb_data_filtered,'omitnan'),'omitnan'))

run(param_script);

figdir=[figrawdir,'Trend\SSH\'];
if (exist(strcat(figdir) , 'dir') ~= 7)
    mkdir(strcat(figdir));
end 
outfile = strcat(figdir,regionname);

cmems_trend=ncread(filename, 'cmems_trend');
cmems_mask=ones(size(cmems_trend));
cmems_mask(isnan(cmems_trend))=NaN;

lon_rho=ncread(filename,'lon_rho');
lat_rho=ncread(filename,'lat_rho');
interped_ssh = ncread(filename, 'interped_ssh');
interped_2y_lowpass_ssh = ncread(filename, 'interped_2y_lowpass');
xlen=size(interped_ssh,1); ylen = size(interped_ssh,2); tlen = size(interped_ssh,3);
for t=1:tlen
    interped_ssh(:,:,t)=interped_ssh(:,:,t).*cmems_mask;
    interped_2y_lowpass_ssh(:,:,t)=interped_2y_lowpass_ssh(:,:,t).*cmems_mask;
end
interped_ssh_2d=reshape(interped_ssh,[xlen*ylen, tlen]);
interped_ssh_mean=mean(interped_ssh_2d,1,'omitnan');
interped_2y_lowpass_ssh_2d=reshape(interped_2y_lowpass_ssh,[xlen*ylen, tlen]);
interped_2y_lowpass_ssh_mean=mean(interped_2y_lowpass_ssh_2d,1,'omitnan');
interped_ssh_3d=reshape(interped_ssh,[xlen*ylen*12, tlen/12]);
interped_ssh_yearly_mean=mean(interped_ssh_3d,1,'omitnan');
interped_ssh_4d=reshape(interped_ssh,[xlen*ylen, 12, tlen/12]);
interped_ssh_seasonal_mean=mean(mean(interped_ssh_4d,1,'omitnan'),3,'omitnan');
for t=1:tlen
    if mod(t,12)==0
        tt=12;
    else
        tt=mod(t,12);
    end
    interped_ssh_seasonal_filtered(t)=interped_ssh_mean(t)-interped_ssh_seasonal_mean(tt);
end

interped2_ssh = ncread(filename2, 'interped_ssh');
interped2_2y_lowpass_ssh = ncread(filename2, 'interped_2y_lowpass');
xlen=size(interped2_ssh,1); ylen = size(interped2_ssh,2); tlen = size(interped2_ssh,3);
for t=1:tlen
    interped2_ssh(:,:,t)=interped2_ssh(:,:,t).*cmems_mask;
    interped2_2y_lowpass_ssh(:,:,t)=interped2_2y_lowpass_ssh(:,:,t).*cmems_mask;
end
interped2_ssh_2d=reshape(interped2_ssh,[xlen*ylen, tlen]);
interped2_ssh_mean=mean(interped2_ssh_2d,1,'omitnan');
interped2_2y_lowpass_ssh_2d=reshape(interped2_2y_lowpass_ssh,[xlen*ylen, tlen]);
interped2_2y_lowpass_ssh_mean=mean(interped2_2y_lowpass_ssh_2d,1,'omitnan');
interped2_ssh_3d=reshape(interped2_ssh,[xlen*ylen*12, tlen/12]);
interped2_ssh_yearly_mean=mean(interped2_ssh_3d,1,'omitnan');
interped2_ssh_4d=reshape(interped2_ssh,[xlen*ylen, 12, tlen/12]);
interped2_ssh_seasonal_mean=mean(mean(interped2_ssh_4d,1,'omitnan'),3,'omitnan');
for t=1:tlen
    if mod(t,12)==0
        tt=12;
    else
        tt=mod(t,12);
    end
    interped2_ssh_seasonal_filtered(t)=interped2_ssh_mean(t)-interped2_ssh_seasonal_mean(tt);
end



soda_ssh = ncread(sodafilename, 'interped_ssh');
soda_2y_lowpass_ssh = ncread(sodafilename, 'interped_2y_lowpass');
for t=1:tlen
    soda_ssh(:,:,t)=soda_ssh(:,:,t).*cmems_mask;
    soda_2y_lowpass_ssh(:,:,t)=soda_2y_lowpass_ssh(:,:,t).*cmems_mask;
end
soda_ssh_2d=reshape(soda_ssh,[xlen*ylen, tlen]);
soda_ssh_mean=mean(soda_ssh_2d,1,'omitnan');
soda_2y_lowpass_ssh_2d=reshape(soda_2y_lowpass_ssh,[xlen*ylen, tlen]);
soda_2y_lowpass_ssh_mean=mean(soda_2y_lowpass_ssh_2d,1,'omitnan');
soda_ssh_3d=reshape(soda_ssh,[xlen*ylen*12, tlen/12]);
soda_ssh_yearly_mean=mean(soda_ssh_3d,1,'omitnan');
soda_ssh_4d=reshape(soda_ssh,[xlen*ylen, 12, tlen/12]);
soda_ssh_seasonal_mean=mean(mean(soda_ssh_4d,1,'omitnan'),3,'omitnan');
for t=1:tlen
    if mod(t,12)==0
        tt=12;
    else
        tt=mod(t,12);
    end
    soda_ssh_seasonal_filtered(t)=soda_ssh_mean(t)-soda_ssh_seasonal_mean(tt);
end


cmems_sla = ncread(filename, 'cmems_sla');
cmems_2y_lowpass_sla = ncread(filename, 'cmems_2y_lowpass');
for t=1:tlen
    cmems_sla(:,:,t)=cmems_sla(:,:,t).*cmems_mask;
    cmems_2y_lowpass_sla(:,:,t)=cmems_2y_lowpass_sla(:,:,t).*cmems_mask;
end
cmems_sla_2d=reshape(cmems_sla,[xlen*ylen, tlen]);
cmems_sla_mean=mean(cmems_sla_2d,1,'omitnan');
cmems_2y_lowpass_sla_2d=reshape(cmems_2y_lowpass_sla,[xlen*ylen, tlen]);
cmems_2y_lowpass_sla_mean=mean(cmems_2y_lowpass_sla_2d,1,'omitnan');
cmems_sla_3d=reshape(cmems_sla,[xlen*ylen*12, tlen/12]);
cmems_sla_yearly_mean=mean(cmems_sla_3d,1,'omitnan');
cmems_sla_4d=reshape(cmems_sla,[xlen*ylen, 12, tlen/12]);
cmems_sla_seasonal_mean=mean(mean(cmems_sla_4d,1,'omitnan'),3,'omitnan');
for t=1:tlen
    if mod(t,12)==0
        tt=12;
    else
        tt=mod(t,12);
    end
    cmems_sla_seasonal_filtered(t)=cmems_sla_mean(t)-cmems_sla_seasonal_mean(tt);
end

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

% start-------------------- msl time series (seasonal filtered)
for folding=1:1
        jpgname=strcat(outfile, '_', testname, '_',regionname, '_msl_all_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
%         if (exist(jpgname , 'file') ~= 2)
%             for varind=1:length(inputyear)*12
%                 msl_filt(varind)=squeeze(mean(mean(ncread(filename,'ssh_filtered',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
%             end
%     %         msl=squeeze(mean(mean(comb_data_filtered,1,'omitnan'),2,'omitnan'));
%             msl_filt=msl_filt-mean(msl_filt);    
            interped_1y_mean=mean(interped_ssh_seasonal_filtered(1:12*5));
            interped2_1y_mean=mean(interped2_ssh_seasonal_filtered(1:12*5));
            soda_1y_mean=mean(soda_ssh_seasonal_filtered(1:12*5));
            cmems_1y_mean=mean(cmems_sla_seasonal_filtered(1:12*5));
            interped_ssh_seasonal_filtered=interped_ssh_seasonal_filtered+(cmems_1y_mean-interped_1y_mean);
            interped2_ssh_seasonal_filtered=interped2_ssh_seasonal_filtered+(cmems_1y_mean-interped2_1y_mean);
            soda_ssh_seasonal_filtered=soda_ssh_seasonal_filtered+(cmems_1y_mean-soda_1y_mean);
            
%             p=polyfit(xData,msl_filt,1);
%             msl2=xData*p(1)+p(2);
            mslplot=plot(xData,cmems_sla_seasonal_filtered,'k')
            hold on
            msl2plot=plot(xData,interped_ssh_seasonal_filtered,'r')
            msl3plot=plot(xData,soda_ssh_seasonal_filtered,'b')
            msl4plot=plot(xData,interped2_ssh_seasonal_filtered,'g')
            

            xlabel('year')
            ylabel('Mean SSH (m)')
            title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') '])
            datetick('x','yyyy','keepticks')
            axis tight;
            ylim(meanplotlev)
            set(mslplot,'LineWidth',2);
            set(msl2plot,'LineWidth',2);
            set(msl3plot,'LineWidth',2);
            set(msl4plot,'LineWidth',2);
            set(gca,'FontSize',15);
            grid on
            lgd=legend('CMEMS','IRCM','SODA','RCM');
            set(lgd,'FontSize',10);
            % set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
            set(gcf,'PaperPosition', [0 0 36 12]) 
            set(lgd,'Orientation','horizontal');
            hold off
            saveas(gcf,jpgname,'jpg');
            grid off
%         end
        close all;
end
% end-------------------- msl time series (seasonal filtered)

% start-------------------- msl time series
for folding=1:1

        jpgname=strcat(outfile, '_', testname, '_',regionname, '_msl_nonfilt_all_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
%         if (exist(jpgname , 'file') ~= 2)
            for varind=1:length(inputyear)*12
                msl(varind)=squeeze(mean(mean(ncread(filename,'raw_ssh',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
            end
            interped_1y_mean=mean(interped_ssh_mean(1:12*5));
            interped2_1y_mean=mean(interped2_ssh_mean(1:12*5));

            soda_1y_mean=mean(soda_ssh_mean(1:12*5));
            cmems_1y_mean=mean(cmems_sla_mean(1:12*5));
            interped_ssh_mean=interped_ssh_mean+(cmems_1y_mean-interped_1y_mean);
            interped2_ssh_mean=interped2_ssh_mean+(cmems_1y_mean-interped_1y_mean);
            soda_ssh_mean=soda_ssh_mean+(cmems_1y_mean-soda_1y_mean);

            
%             msl=msl-mean(msl);     
%             
%             p=polyfit(xData,msl,1);
%             msl2=xData*p(1)+p(2);
%             mslplot=plot(xData,msl,'k')
%             hold on
%             mslplot2=plot(xData,msl2,'Color','r')

            mslplot=plot(xData,cmems_sla_mean,'k')
            hold on
            msl2plot=plot(xData,interped_ssh_mean,'r')
            msl3plot=plot(xData,soda_ssh_mean,'b')
            msl4plot=plot(xData,interped2_ssh_mean,'g')

            xlabel('year')
            ylabel('Mean SSH (m)')
            mean_trend=ncread(filename,'mean_trend');
            title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') '])
            datetick('x','yyyy','keepticks')
            axis tight;
            ylim(meanplotlev)
            set(mslplot,'LineWidth',2);
            set(msl2plot,'LineWidth',2);
            set(msl3plot,'LineWidth',2);
            set(msl4plot,'LineWidth',2);            
            set(gca,'FontSize',15);
            grid on
            lgd=legend('CMEMS','IRCM','SODA','RCM');
            set(lgd,'FontSize',10);
            set(gcf,'PaperPosition', [0 0 36 12]) 
            set(lgd,'Orientation','horizontal');
            hold off
            saveas(gcf,jpgname,'jpg');
            grid off
            close all;
%         end
end
% end-------------------- msl time series

% start-------------------- yearly msl time series
for folding=1:1
        jpgname=strcat(outfile, '_', testname, '_',regionname, '_msl_yearly_all_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
%         if (exist(jpgname , 'file') ~= 2)
            interped_1y_yearly_mean=mean(interped_ssh_yearly_mean(1:5));
            interped2_1y_yearly_mean=mean(interped2_ssh_yearly_mean(1:5));            
            soda_1y_yearly_mean=mean(soda_ssh_yearly_mean(1:5));
            cmems_1y_yearly_mean=mean(cmems_sla_yearly_mean(1:5));
            interped_ssh_yearly_mean=interped_ssh_yearly_mean+(cmems_1y_yearly_mean-interped_1y_yearly_mean);
            interped2_ssh_yearly_mean=interped2_ssh_yearly_mean+(cmems_1y_yearly_mean-interped2_1y_yearly_mean);
            soda_ssh_yearly_mean=soda_ssh_yearly_mean+(cmems_1y_yearly_mean-soda_1y_yearly_mean);

            
%             msl=msl-mean(msl);     
%             
%             p=polyfit(xData,msl,1);
%             msl2=xData*p(1)+p(2);
%             mslplot=plot(xData,msl,'k')
%             hold on
%             mslplot2=plot(xData,msl2,'Color','r')

            mslplot=plot(xData2,cmems_sla_yearly_mean,'k')
            hold on
            msl2plot=plot(xData2,interped_ssh_yearly_mean,'r')
            msl3plot=plot(xData2,soda_ssh_yearly_mean,'b')
            msl4plot=plot(xData2,interped2_ssh_yearly_mean,'g')
            
            xlabel('year')
            ylabel('Mean SSH (m)')
            mean_trend=ncread(filename,'mean_trend');
            title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') '])
            datetick('x','yyyy','keepticks')
            axis tight;
            ylim(meanplotlev)
            set(mslplot,'LineWidth',2);
            set(msl2plot,'LineWidth',2);
            set(msl3plot,'LineWidth',2);
            set(msl4plot,'LineWidth',2);
            set(gca,'FontSize',15);
            grid on
            lgd=legend('CMEMS','IRCM','SODA','RCM');
            set(lgd,'FontSize',10);
            set(gcf,'PaperPosition', [0 0 36 12]) 
            set(lgd,'Orientation','horizontal');
            hold off
            saveas(gcf,jpgname,'jpg');
            grid off
            close all;
%         end
end
% end-------------------- yearly msl time series

% start-------------------- 2y lowpass msl time series
for folding=1:1
        jpgname=strcat(outfile, '_', testname, '_',regionname, '_msl_2y_lowpass_all_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
%         if (exist(jpgname , 'file') ~= 2)
            interped_1y_mean=mean(interped_2y_lowpass_ssh_mean(1:12*5),'omitnan');
            interped2_1y_mean=mean(interped2_2y_lowpass_ssh_mean(1:12*5),'omitnan');
            soda_1y_mean=mean(soda_2y_lowpass_ssh_mean(1:12*5),'omitnan');
            cmems_1y_mean=mean(cmems_2y_lowpass_sla_mean(1:12*5),'omitnan');
            
            interped_2y_lowpass_ssh_mean=interped_2y_lowpass_ssh_mean+(cmems_1y_mean-interped_1y_mean);
            interped2_2y_lowpass_ssh_mean=interped2_2y_lowpass_ssh_mean+(cmems_1y_mean-interped2_1y_mean);
            soda_2y_lowpass_ssh_mean=soda_2y_lowpass_ssh_mean+(cmems_1y_mean-soda_1y_mean);

            mslplot=plot(xData,cmems_2y_lowpass_sla_mean,'k')
            hold on
            msl2plot=plot(xData,interped_2y_lowpass_ssh_mean,'r')
            msl3plot=plot(xData,soda_2y_lowpass_ssh_mean,'b')
            msl4plot=plot(xData,interped2_2y_lowpass_ssh_mean,'g')
            
            xlabel('year')
            ylabel('Mean SSH (m)')
            mean_trend=ncread(filename,'mean_trend');
            title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') '])
            datetick('x','yyyy','keepticks')
%             axis tight;
            ylim(meanplotlev)
            xlim([min(xData) max(xData)])
            set(mslplot,'LineWidth',2);
            set(msl2plot,'LineWidth',2);
            set(msl3plot,'LineWidth',2);
            set(msl4plot,'LineWidth',2);
            set(gca,'FontSize',15);
            grid on
            lgd=legend('CMEMS','IRCM','SODA','RCM');
            set(lgd,'FontSize',10);
            set(gcf,'PaperPosition', [0 0 36 12]) 
            set(lgd,'Orientation','horizontal');
            hold off
            saveas(gcf,jpgname,'jpg');
            grid off
            close all;
%         end
end
% end-------------------- 2y lowpass msl time series



% start-------------------- seasonal msl
for folding=1:1

        climdir = [figdir,'\CLIM\'];
        if (exist(strcat(climdir) , 'dir') ~= 7)
            mkdir(strcat(climdir));
        end 
        climoutfile = strcat(climdir,regionname);
        jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_seasonal_msl_all_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),  '.jpg'); %% ~_year_month.jpg
        
%         clim_ssh=ncread(filename,'clim_ssh');
%         for ijij=1:12
%             temp_clim_ssh=clim_ssh(:,:,ijij);
%             mean_clim_ssh(ijij)=mean(temp_clim_ssh(:), 'omitnan');
%         end
%         mean_clim_ssh=mean_clim_ssh-mean(mean_clim_ssh);
        
%         if (exist(jpgname , 'file') ~= 2)
            comb_interped_data=ncread(filename,'interped_ssh');
            lon_cmems=ncread(filename, 'lon_cmems');
            lat_cmems=ncread(filename, 'lat_cmems');
            len_lon=length(lon_cmems);
            len_lat=length(lat_cmems);
            
            interped_1y_seasonal_mean=mean(interped_ssh_seasonal_mean(1));
            interped2_1y_seasonal_mean=mean(interped2_ssh_seasonal_mean(1));            
            soda_1y_seasonal_mean=mean(soda_ssh_seasonal_mean(1));
            cmems_1y_seasonal_mean=mean(cmems_sla_seasonal_mean(1));
            interped_ssh_seasonal_mean=interped_ssh_seasonal_mean+(cmems_1y_seasonal_mean-interped_1y_seasonal_mean);
            interped2_ssh_seasonal_mean=interped2_ssh_seasonal_mean+(cmems_1y_seasonal_mean-interped2_1y_seasonal_mean);            
            soda_ssh_seasonal_mean=soda_ssh_seasonal_mean+(cmems_1y_seasonal_mean-soda_1y_seasonal_mean);

            mslplot=plot(1:12,cmems_sla_seasonal_mean,'k')
            hold on
            msl2plot=plot(1:12,interped_ssh_seasonal_mean,'r')
            msl3plot=plot(1:12,soda_ssh_seasonal_mean,'b')
            msl4plot=plot(1:12,interped2_ssh_seasonal_mean,'g')
            
            xlabel('month')
            ylabel('Mean SSH (m)')
            title([regionname, ', seasonal msl(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
            axis tight;
            ylim(meanplotlev)
            set(mslplot,'LineWidth',2);
            set(msl2plot,'LineWidth',2);
            set(msl3plot,'LineWidth',2);
            set(msl4plot,'LineWidth',2);
            set(gca,'FontSize',20);
            grid on
            lgd=legend('CMEMS','IRCM','SODA','RCM');
            set(lgd,'FontSize',10);
            set(gcf,'PaperPosition', [0 0 36 12]) 
            set(lgd,'Orientation','horizontal');
            hold off
            saveas(gcf,jpgname,'jpg');
            grid off
            close all;
%         end
end
% end-------------------- seasonal msl