close all; clear all;  clc;
warning off;
% all_region2 ={'NWP','NWP2'}
% all_region2 ={'ES', 'YS', 'SS', 'NWP2'}
% all_testname2 = {'ens02','test53','test55','test56'};
all_testname2 = {'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'NorESM1-M', 'MPI-ESM-LR'};
all_testname3 = {'test57', 'test58', 'test59', 'test60'};

% all_region2 ={'NWP','ES', 'SS', 'YS', 'AKP'}
all_region2 ={'NWP', 'AKP2'}

% all_region2 ={'NWP'};

% all_region2 ={'AKP2'}
for regionind2=1:length(all_region2)
    close all;
    clearvars '*' -except regionind2 testnameind2 all_region2 all_testname2 all_testname3
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

    shadlev = [0 35];
    rms_shadlev = [0 4];
%     trendlev = [-3 3];  %% trend lev
    trendlev = [-10 10];  %% trend lev
    abstrendlev =[4 7];
    reltrendlev =[-5 5];
    conlev  = 0:5:35;
    meanplotlev =[-0.3 0.3];
    trendplotlev = [4 6.5];
    sshlev =[-0.7 1.3];
    sshdifflev = [40 70];

    % for snu_desktopd
%     testname=all_testname2{testnameind2}    % % need to change
    inputyear = [2006:2100]; % % put year which you want to plot [year year ...]
    inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]

    varname ='zeta'
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

%         load(['G:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',testname,'_',regionname, ...
%             'ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
    
    allname = ['G:\Data\Model\ROMS\nwp_1_20\','rcp_45_all','_',regionname, ...
                'ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat'];
%     if (exist(allname , 'file') ~= 2)
        for testind2=1:length(all_testname2)
        testname=all_testname2{testind2};
            load(['G:\Data\Model\CMIP5\zos\rcp45\Omon\',testname,'\',testname,'_',regionname, ...
                    'ssh_trend_rcp45_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
            trend_all{testind2}=trend;
            trend_filtered_all{testind2}=double(trend_filtered);
            trend_clim_all{testind2}=double(trend_clim);
            mean_clim_trend_all{testind2}=squeeze(mean(mean(trend_clim,1,'omitnan'),2,'omitnan'));
            lon_all{testind2}=lon;
            lat_all{testind2}=lat;
            msl_all{testind2}=squeeze(mean(mean(comb_data,1,'omitnan'),2,'omitnan'));
            msl_filt_all{testind2}=squeeze(mean(mean(comb_data_filtered,1,'omitnan'),2,'omitnan'));
        end
        for testind3=1:length(all_testname3)
            testname=all_testname3{testind3};
            filename = ['G:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',testname,'_',regionname, ...
                        '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];
            trend = ncread(filename,'trend');
            trend_filtered = double(ncread(filename,'trend_filtered'));
            lon_rho = ncread(filename, 'lon_rho');
            lat_rho = ncread(filename, 'lat_rho');
            trend_clim = ncread(filename,'clim_ssh_trend')*1000.0;
            trend_all{testind3+length(all_testname2)}=trend;
            trend_filtered_all{testind3+length(all_testname2)}=double(trend_filtered);
            trend_clim_all{testind3+length(all_testname2)}=trend_clim;
            mean_clim_trend_all{testind3+length(all_testname2)}=squeeze(mean(mean(trend_clim,1,'omitnan'),2,'omitnan'));
            lon_all{testind3+length(all_testname2)}=lon_rho;
            lat_all{testind3+length(all_testname2)}=lat_rho;
            msl_all{testind3+length(all_testname2)}=squeeze(mean(mean(comb_data,1,'omitnan'),2,'omitnan'));
            msl_filt_all{testind3+length(all_testname2)}=squeeze(mean(mean(comb_data_filtered,1,'omitnan'),2,'omitnan'));
        end
        save(allname)
%     else
%         load(allname)
%     end

    valnum=0;
    wrmap = bwrmap(51:100,:);
% %     valid cell number
%      for vi=1:size(comb_spatial_meanressh,1)
%          for vj=1:size(comb_spatial_meanressh,2)
%              if (isnan(comb_spatial_meanressh(vi,vj,1))~=1)
%                 valnum=valnum+1;
%              end
%          end
%      end

%     plot(squeeze(mean(mean(comb_data_filtered,1,'omitnan'),2,'omitnan')))
%     isize = size(comb_data_filtered,1)
%     jsize = size(comb_data_filtered,2)
%     lsize = size(comb_data_filtered,3)
%     comb_yearly_data_filtered=reshape(comb_data_filtered,[isize, jsize, 12, lsize/12]);
%     mean_yearly_data_filtered=squeeze(mean(mean(mean(comb_yearly_data_filtered,1,'omitnan'),2,'omitnan'),3,'omitnan'));
%     trendtime=14:29;
%     p=polyfit(trendtime,mean_yearly_data_filtered(14:29)',1);
%     yearly_interped_trend=p(1);
%     yearly_interped_trend = yearly_interped_trend * 1000.0; %% m/y -> mm/y


%         if (strcmp(system_name,'PCWIN64'))
        % % for windows
        figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\4th_year\figure\nwp_1_20\all\',regionname,'\'); % % where figure files will be saved
%             param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m'
        param_script =['C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_CMIP5_', regionname, '.m']
        filedir = strcat('G:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
%         elseif (strcmp(system_name,'GLNXA64'))
%         end
%         std(mean(mean(comb_data_filtered,'omitnan'),'omitnan'))

    var='SSH';
    run(param_script);

    figdir=[figrawdir,'Trend\SSH\'];
    if (exist(strcat(figdir) , 'dir') ~= 7)
        mkdir(strcat(figdir));
    end 
    outfile = strcat(figdir,regionname);

    lon_rho=ncread(filename,'lon_rho');
    lat_rho=ncread(filename,'lat_rho');
    trend_filtered = ncread(filename,'trend_filtered');
    mean_trend_filtered = ncread(filename,'mean_trend_filtered');


% start-------------------- make timedata for time series
    for i =1:length(inputyear) 
        tempyear=inputyear(i);
        for month=1:12
            xData((12*(i-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-01',]);
        end
    end
% end-------------------- make timedata for time series  

% % % % start-------------------- msl time series (seasonal filtered)
% % % 
% % %     figdir=[figrawdir,'Trend\SSH\'];
% % %     if (exist(strcat(figdir) , 'dir') ~= 7)
% % %         mkdir(strcat(figdir));
% % %     end 
% % %     outfile = strcat(figdir,regionname);
% % % 
% % %     jpgname=strcat(outfile, '_', testname, '_',regionname, '_msl_all_', ...
% % %         num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
% % %     if (exist(jpgname , 'file') ~= 2)
% % %         for mslind=1:length(msl_filt_all)
% % %             msl_filt_all{mslind}=msl_filt_all{mslind}-mean(msl_filt_all{mslind});
% % %             msl_filt_all_p{mslind}=polyfit(xData,msl_filt_all{mslind},1);
% % %             msl_filt_all2{mslind}=xData*msl_filt_all_p{mslind}(1)+msl_filt_all_p{mslind}(2);
% % %         end
% % %         --------------
% % %         mslplot=plot(xData,msl_filt,'k')
% % %         hold on
% % %         mslplot2=plot(xData,msl2,'Color','r')
% % %         xlabel('year')
% % %         ylabel('Mean SSH (m)')
% % %         title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(mean_trend_filtered,2)), ' mm/y'])
% % %         datetick('x','yyyy','keepticks')
% % %         axis tight;
% % %         ylim(meanplotlev)
% % %         set(mslplot,'LineWidth',2);
% % %         set(gca,'FontSize',20);
% % %         grid on
% % %         hold off
% % %         saveas(gcf,jpgname,'jpg');
% % %         grid off
% % %     end
% % %     close all;
% % % % end-------------------- msl time series (seasonal filtered)

% % % % start-------------------- msl time series
% % %     jpgname=strcat(outfile, '_', testname, '_',regionname, '_msl_nonfilt_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
% % %     if (exist(jpgname , 'file') ~= 2)
% % %         for varind=1:length(inputyear)*12
% % %             msl(varind)=squeeze(mean(mean(ncread(filename,'raw_ssh',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
% % %         end
% % %         msl=msl-mean(msl);    
% % %         p=polyfit(xData,msl,1);
% % %         msl2=xData*p(1)+p(2);
% % %         mslplot=plot(xData,msl,'k')
% % %         hold on
% % %         mslplot2=plot(xData,msl2,'Color','r')
% % %         xlabel('year')
% % %         ylabel('Mean SSH (m)')
% % %         mean_trend=ncread(filename,'mean_trend');
% % %         title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(mean_trend,2)), ' mm/y'])
% % %         datetick('x','yyyy','keepticks')
% % %         axis tight;
% % %         ylim(meanplotlev)
% % %         set(mslplot,'LineWidth',2);
% % %         set(gca,'FontSize',20);
% % %         grid on
% % %         hold off
% % %         saveas(gcf,jpgname,'jpg');
% % %         grid off
% % %         close all;
% % %     end
% % % % end-------------------- msl time series

% % % % start-------------------- climatological msl time series
% % %     climdir = [figdir,'\CLIM\'];
% % %     if (exist(strcat(climdir) , 'dir') ~= 7)
% % %         mkdir(strcat(climdir));
% % %     end 
% % %     climoutfile = strcat(climdir,regionname);
% % %     for monthij=1:12
% % %         jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_clim_msl_', ...
% % %             num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',  num2str(monthij, '%02i'), '.jpg'); %% ~_year_month.jpg
% % %         if (exist(jpgname , 'file') ~= 2)
% % %             if (exist('msl')==0)
% % %                 for varind=1:length(inputyear)*12
% % %                     msl(varind)=squeeze(mean(mean(ncread(filename,'raw_ssh',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
% % %                 end
% % %                 climmsl=reshape(msl,[12,length(inputyear)]);
% % %             else
% % %                 climmsl=reshape(msl,[12,length(inputyear)]);
% % %             end
% % %             tempmsl=squeeze(climmsl(monthij,:));
% % %             tempmsl=tempmsl-mean(tempmsl);
% % %             p=polyfit(inputyear,tempmsl,1);
% % %             msl2=inputyear*p(1)+p(2);
% % %             mslplot=plot(inputyear,tempmsl,'k')
% % %             hold on
% % %             mslplot2=plot(inputyear,msl2,'Color','r')
% % %             xlabel('year')
% % %             ylabel('Mean SSH (m)')
% % %             title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), ...
% % %                 ', ', calendarname{monthij}(1:3), '), ',num2str(round(mean_clim_trend(monthij),2)), ' mm/y'])
% % % %                 datetick('x','yyyy','keepticks')
% % %             axis tight;
% % %             ylim(meanplotlev)
% % %             set(mslplot,'LineWidth',2);
% % %             set(gca,'FontSize',20);
% % %             grid on
% % %             hold off
% % %             saveas(gcf,jpgname,'jpg');
% % %             grid off
% % %             close all;
% % %         end
% % %     end
% % % % end-------------------- climatological msl time series

% start-------------------- climatological msl trend (all)
    climdir = [figdir,'\CLIM\'];
    if (exist(strcat(climdir) , 'dir') ~= 7)
        mkdir(strcat(climdir));
    end 
    climoutfile = strcat(climdir,regionname);
    jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_clim_msl_trend_all_', ...
        num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',  num2str(monthij, '%02i'), '.jpg'); %% ~_year_month.jpg
%     if (exist(jpgname , 'file') ~= 2)
        for climind=1:length(mean_clim_trend_all)
            mean_clim_trend_all2(climind,:)=mean_clim_trend_all{climind};
        end
        mslplot1=plot(1:12,mean(mean_clim_trend_all2(1:4,:),1),'k');
        hold on
        mslplot2=plot(1:12,mean(mean_clim_trend_all2(5:8,:),1),'r');
        xlabel('month')
        ylabel('trend (mm/yr)')
%         title([regionname, ', climatological trend(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
        axis tight;
        ylim(trendplotlev)
        set(mslplot1,'LineWidth',2);
        set(mslplot2,'LineWidth',2);
        set(gca,'FontSize',20);
        grid on
        hold off
        lgd=legend('Glo.Ens', 'Reg.Ens');
        set(lgd,'FontSize',20);
        % set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
        set(lgd,'Orientation','horizontal');
        saveas(gcf,jpgname,'jpg');
        grid off
        close all;
%     end
% end-------------------- climatological msl trend (all)

% start-------------------- climatological msl trend (each model)
    for cmipind=1:length(all_testname2)
        testname=all_testname2{cmipind};
        climdir = [figdir,'\CLIM\'];
        if (exist(strcat(climdir) , 'dir') ~= 7)
            mkdir(strcat(climdir));
        end 
        climoutfile = strcat(climdir,regionname);
        jpgname=strcat(climoutfile, '_', testname, '_',regionname, '_clim_msl_trend_IPSL_LR_', ...
            num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '_',  num2str(monthij, '%02i'), '.jpg'); %% ~_year_month.jpg
%         if (exist(jpgname , 'file') ~= 2)
            for climind=1:length(mean_clim_trend_all)
                mean_clim_trend_all2(climind,:)=mean_clim_trend_all{climind};
            end
            mslplot1=plot(1:12,mean(mean_clim_trend_all2(cmipind,:),1),'k');
            hold on
            mslplot2=plot(1:12,mean(mean_clim_trend_all2(cmipind+4,:),1),'r');
            xlabel('month')
            ylabel('trend (mm/yr)')
    %         title([regionname, ', climatological trend(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
            axis tight;
            ylim(trendplotlev)
            set(mslplot1,'LineWidth',2);
            set(mslplot2,'LineWidth',2);
            set(gca,'FontSize',20);
            grid on
            hold off
            lgd=legend('Global', 'Regional');
            set(lgd,'FontSize',20);
            % set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
            set(lgd,'Orientation','horizontal', 'Location', 'south');
            saveas(gcf,jpgname,'jpg');
            grid off
            close all;
%         end
    end
% end-------------------- climatological msl trend (each model)

% if (strcmp(regionname,'AKP2'))
    for modelind=1:4
        mask_akp2{modelind} = double(inpolygon(lon_all{modelind},lat_all{modelind},akp2polygon(:,1),akp2polygon(:,2)));

        mask_ss{modelind} = double(inpolygon(lon_all{modelind},lat_all{modelind},akp2polygon(:,1),akp2polygon(:,2)));

        mask_es{modelind} = double(inpolygon(lon_all{modelind},lat_all{modelind},espolygon(:,1),espolygon(:,2)));
        mask_ys{modelind} = double(inpolygon(lon_all{modelind},lat_all{modelind},yspolygon(:,1),yspolygon(:,2)));
        
        %         mask_ss{modelind} = double(inpolygon(lon_all{modelind},lat_all{modelind},sspolygon(:,1),sspolygon(:,2)));
       
        mask_akp2{4+modelind} = double(inpolygon(lon_rho,lat_rho,akp2polygon(:,1),akp2polygon(:,2)));
        mask_ss{4+modelind} = double(inpolygon(lon_rho,lat_rho,akp2polygon(:,1),akp2polygon(:,2)));
        mask_es{4+modelind} = double(inpolygon(lon_rho,lat_rho,espolygon(:,1),espolygon(:,2)));
        mask_ys{4+modelind} = double(inpolygon(lon_rho,lat_rho,yspolygon(:,1),yspolygon(:,2)));
        
        es_trend_filtered_all{modelind}=trend_filtered_all{modelind}.*mask_es{modelind};
        es_trend_filtered_all{modelind+4}=double(trend_filtered_all{modelind+4}.*mask_es{modelind+4});
        ys_trend_filtered_all{modelind}=trend_filtered_all{modelind}.*mask_ys{modelind};
        ys_trend_filtered_all{modelind+4}=double(trend_filtered_all{modelind+4}.*mask_ys{modelind+4});
%         trend_filtered_all{modelind}=trend_filtered_all{modelind}.*mask_akp2{modelind};
%         trend_filtered_all{modelind+4}=trend_filtered_all{modelind+4}.*mask_akp2{modelind+4};

%         mask_ss{modelind}(isfinite(es_trend_filtered_all{modelind}))=0;
%         mask_ss{modelind}(isfinite(ys_trend_filtered_all{modelind}))=0;
        mask_ss{modelind}=mask_ss{modelind}-mask_es{modelind}-mask_ys{modelind};
        mask_ss{modelind}(mask_ss{modelind}<0.01)=0;
        mask_ss{modelind+4}=mask_ss{modelind+4}-mask_es{modelind+4}-mask_ys{modelind+4};
        mask_ss{modelind+4}(mask_ss{modelind+4}<0.01)=0;
        
        ss_trend_filtered_all{modelind}=trend_filtered_all{modelind}.*mask_ss{modelind};
        ss_trend_filtered_all{modelind+4}=double(trend_filtered_all{modelind+4}.*mask_ss{modelind+4});
    end
    for modelind=1:8
%         es_trend_filtered_all{modelind}(es_trend_filtered_all{modelind}==0)=NaN;
%         ys_trend_filtered_all{modelind}(ys_trend_filtered_all{modelind}==0)=NaN;
%         ss_trend_filtered_all{modelind}(ss_trend_filtered_all{modelind}==0)=NaN;
        es_trend_filtered_all{modelind}(mask_es{modelind}==0)=NaN;
        ys_trend_filtered_all{modelind}(mask_ys{modelind}==0)=NaN;
        ss_trend_filtered_all{modelind}(mask_ss{modelind}==0)=NaN;
%         trend_filtered_all{modelind}(mask_akp2{modelind}==0)=NaN;
        es_valid_cell{modelind}=sum(sum(isfinite(es_trend_filtered_all{modelind})));
        ys_valid_cell{modelind}=sum(sum(isfinite(ys_trend_filtered_all{modelind})));
        ss_valid_cell{modelind}=sum(sum(isfinite(ss_trend_filtered_all{modelind})));
        akp2_valid_cell{modelind}=sum(sum(isfinite(trend_filtered_all{modelind})));
%         es_mean_trend_filtered(modelind)=mean(mean(es_trend_filtered_all{modelind},'omitnan','double'),'omitnan','double')
        es_mean_trend_filtered(modelind)=mean(mean(es_trend_filtered_all{modelind},'omitnan'),'omitnan')
        ys_mean_trend_filtered(modelind)=mean(mean(ys_trend_filtered_all{modelind},'omitnan'),'omitnan')
        ss_mean_trend_filtered(modelind)=mean(mean(ss_trend_filtered_all{modelind},'omitnan'),'omitnan')
        es_mean_trend_filtered2(modelind)=mean(mean(double(es_trend_filtered_all{modelind}(~isnan(es_trend_filtered_all{modelind})))))
        ys_mean_trend_filtered2(modelind)=mean(mean(double(ys_trend_filtered_all{modelind}(~isnan(ys_trend_filtered_all{modelind})))))
        ss_mean_trend_filtered2(modelind)=mean(mean(double(ss_trend_filtered_all{modelind}(~isnan(ss_trend_filtered_all{modelind})))))
        es_mean_trend_filtered3(modelind)=sum(sum(es_trend_filtered_all{modelind},'omitnan'),'omitnan') / es_valid_cell{modelind};
        ys_mean_trend_filtered3(modelind)=sum(sum(ys_trend_filtered_all{modelind},'omitnan'),'omitnan') / ys_valid_cell{modelind};
        ss_mean_trend_filtered3(modelind)=sum(sum(ss_trend_filtered_all{modelind},'omitnan'),'omitnan') / ss_valid_cell{modelind};
        
        mean_trend_filtered0(modelind)=nanmean(nanmean(double(trend_filtered_all{modelind})))
        mean_trend_filtered(modelind)=mean(mean(double(trend_filtered_all{modelind}),'omitnan'),'omitnan')
        mean_trend_filtered2(modelind)=mean(mean(double(trend_filtered_all{modelind}(~isnan(trend_filtered_all{modelind})))))
        mean_trend_filtered3(modelind)=sum(sum(trend_filtered_all{modelind},'omitnan'),'omitnan') / akp2_valid_cell{modelind};
        mean_trend_filtered4(modelind)= (es_mean_trend_filtered2(modelind).*es_valid_cell{modelind} + ys_mean_trend_filtered2(modelind).*ys_valid_cell{modelind} + ss_mean_trend_filtered2(modelind).*ss_valid_cell{modelind}) ...
            / double(es_valid_cell{modelind} + ys_valid_cell{modelind} + ss_valid_cell{modelind})
        
        rise_mean(modelind)=mean_trend_filtered2(modelind)*94.0/10.0
    end
    
    pcolor(es_trend_filtered_all{5}'+ss_trend_filtered_all{5}'+ys_trend_filtered_all{5}')
    pcolor(es_trend_filtered_all{5}'+ys_trend_filtered_all{5}')
    pcolor(ss_trend_filtered_all{5}')
  
    mean(mean_trend_filtered3(5))
    mean(es_mean_trend_filtered3(5))
    mean(ys_mean_trend_filtered3(5))
    mean(ss_mean_trend_filtered3(5))

    
    
    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
    m_pcolor(lon_all{5}',lat_all{5}',es_trend_filtered_all{5}')
    shading flat
    hold on
    m_pcolor(lon_all{5}',lat_all{5}',ss_trend_filtered_all{5}')
    m_pcolor(lon_all{5}',lat_all{5}',ys_trend_filtered_all{5}')
    shading flat
    hold off
    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
      
      abc = trend_filtered_all{5}+trend_filtered_all{6}+trend_filtered_all{7}+trend_filtered_all{8};
      abc=abc/4.0;
      mean(mean(abc,'omitnan'),'omitnan')

      abc = es_trend_filtered_all{5}+es_trend_filtered_all{6}+es_trend_filtered_all{7}+es_trend_filtered_all{8};
      abc=abc/4.0;
      mean(mean(abc,'omitnan'),'omitnan')
      
      abc = ss_trend_filtered_all{5}+ss_trend_filtered_all{6}+ss_trend_filtered_all{7}+ss_trend_filtered_all{8};
      abc=abc/4.0;
      mean(mean(abc,'omitnan'),'omitnan')
      
      abc = ys_trend_filtered_all{5}+ys_trend_filtered_all{6}+ys_trend_filtered_all{7}+ys_trend_filtered_all{8};
      abc=abc/4.0;
      mean(mean(abc,'omitnan'),'omitnan')
      
      iiii=5
      abc=mean(mean(trend_filtered_all{iiii},'omitnan'),'omitnan')
      abcd=mean(mean(es_trend_filtered_all{iiii},'omitnan'),'omitnan')
      abcde=mean(mean(ys_trend_filtered_all{iiii},'omitnan'),'omitnan')
      abcdef=mean(mean(ss_trend_filtered_all{iiii},'omitnan'),'omitnan')
      finabc=(abcd*es_valid_cell{iiii}+abcde*ys_valid_cell{iiii}+abcdef*ss_valid_cell{iiii})/(es_valid_cell{iiii}+ys_valid_cell{iiii}+ss_valid_cell{iiii})
      
      
      sum(sum(~isnan(es_trend_filtered_all{5}))) + sum(sum(~isnan(ys_trend_filtered_all{5}))) + sum(sum(~isnan(ss_trend_filtered_all{5})))
      sum(sum(~isnan(trend_filtered_all{5})))
      
      
%       f_trend_filtered_all{iiii}=zeros(size(trend_filtered_all{iiii}));
%       f_es_trend_filtered_all{iiii}(~isnan(es_trend_filtered_all{iiii}))=es_trend_filtered_all{iiii}(~isnan(es_trend_filtered_all{iiii}));
%       f_ys_trend_filtered_all{iiii}(~isnan(ys_trend_filtered_all{iiii}))=ys_trend_filtered_all{iiii}(~isnan(ys_trend_filtered_all{iiii}));
%       f_ss_trend_filtered_all{iiii}(~isnan(ss_trend_filtered_all{iiii}))=ss_trend_filtered_all{iiii}(~isnan(ss_trend_filtered_all{iiii}));
%       f_trend_filtered_all{iiii}=f_es_trend_filtered_all{iiii} + f_ys_trend_filtered_all{iiii} + f_ss_trend_filtered_all{iiii};
%       sum(sum(f_trend_filtered_all{iiii}))
%       sum(sum(trend_filtered_all{iiii},'omitnan'),'omitnan')
      
      mean(es_trend_filtered_all{iiii}(~isnan(es_trend_filtered_all{iiii})))
      mean(ys_trend_filtered_all{iiii}(~isnan(ys_trend_filtered_all{iiii})))
      mean(ss_trend_filtered_all{iiii}(~isnan(ss_trend_filtered_all{iiii})))
      mean(mean(es_trend_filtered_all{iiii},'omitnan'),'omitnan')
% end
end
