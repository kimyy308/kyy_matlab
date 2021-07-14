close all; clear all;  clc;
warning off;
% all_region2 ={'NWP','NWP2'}
% all_region2 ={'ES', 'YS', 'SS', 'NWP2'}
% all_testname2 = {'ens02','test53','test55','test56'};
% all_testname2 = {'test11', 'test12'};
% all_testname2 = {'test53', 'test54', 'test55', 'test56', 'ens03'};
all_testname2 = {'test53', 'test54', 'test55', 'test56'};
% all_testname2 = {'test56'};

% all_region2 ={'NWP','ES', 'SS', 'YS', 'AKP'}
% all_region2 ={'NWP','AKP4', 'YS', 'YSECS', 'ECS2', 'NES', 'SES', 'ES'}

% all_region2 ={'NWP', 'AKP4'};
all_region2 ={'AKP4'};

% all_region2 ={'TEST'}
for testnameind2=1:length(all_testname2)
    for regionind2=1:length(all_region2)
        close all;
        clearvars '*' -except regionind2 testnameind2 all_region2 all_testname2 all_ysecs_corr all_es_corr ...
            all_ysecs_u_corr all_es_u_corr all_ysecs_v_corr all_es_v_corr
        % % % 
        % % % Read Model SST
        % % % interp
        % % % get RMS
        % % % get BIAS
        system_name=computer;
        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            dropboxpath='C:\Users\user\Dropbox';
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
        abstrendlev =[2 7];
        reltrendlev =[-5 5];
        conlev  = 0:5:35;
        meanplotlev =[-0.15 0.15];
        trendplotlev = [0 7];
        trenddifflev = [-10 10];
        sshlev =[-0.3 0.3];
        sshdifflev = [0 20];

        % for snu_desktopd
        testname=all_testname2{testnameind2}    % % need to change
        inputyear = [1993:2005]; % % put year which you want to plot [year year ...]
%         inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
        inputmonth = [1:12]; % % put month which you want to plot [month month ...]
        gcmtestname = Func_0004_get_GCMname_from_RCM(testname);
        scenname= Func_0003_RCM_CMIP5_scenname(testname);
        varname ='zeta';
        variable='SSH';
        run('nwp_polygon_point.m');
        regionname=all_region2{regionind2};
        
% % %         switch region
        [refpolygon, lonlat] = Func_0007_get_polygon_data_from_regionname(regionname);
        
% % %         flag configuration
        fig_flags{43,1}='pattern_correlation';

        for flagi=1:100
            fig_flags{flagi,2}=0;
        end
        
        fig_flags{43,2}=2;
        fig_flags{71,2}=2;
        fig_flags{72,2}=2;
        
%         load(['G:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
%             'ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
        
%         filename = ['E:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
%                     '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];

        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\5th_year\figure\nwp_1_20\',testname,'\',regionname,'\'); % % where figure files will be saved
%             param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_EKB_RMS.m'
            param_script =['C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname, '.m']
            filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
            cmip5dir = strcat('D:\Data\Model\CMIP5\'); % % where data files are
            gcmfiledir = strcat('D:\Data\Model\CMIP5\zos\', scenname, '\interp\', gcmtestname, '\'); % % where data files are
        elseif (strcmp(system_name,'GLNXA64'))
        end
%         std(mean(mean(comb_data_filtered,'omitnan'),'omitnan'))

        modelfilename = strcat(filedir, testname,'_',regionname, '_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
        cmemsfilename = strcat(filedir, testname,'_',regionname, 'cmems_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
        interpedfilename = strcat(filedir, testname,'_',regionname, 'cmems_interped_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
        corrfilename = strcat(filedir, testname,'_',regionname, 'cmems_interped_ssh_corr_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
        lpffilename = strcat(filedir, testname,'_',regionname, 'low_pass_filtered_cmems_interped_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
        lpf_corrfilename = strcat(filedir, testname,'_',regionname, 'low_pass_filtered_cmems_interped_ssh_corr_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
        movfilename = strcat(filedir, testname,'_',regionname, 'moving_averaged_cmems_interped_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
        mov_corrfilename = strcat(filedir, testname,'_',regionname, 'moving_averaged_cmems_interped_ssh_corr_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
        detrendfilename = strcat(filedir, testname,'_',regionname, 'detrended_cmems_interped_ssh_corr_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
        cmemsvelfilename = strcat(filedir, testname,'_',regionname, 'cmems_vel_rms_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
        interpedvelfilename = strcat(filedir, testname,'_',regionname, 'cmems_interped_vel_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
        avhrrfilename = strcat(filedir, testname, '_',regionname,'_rms_clim_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
        gcminterpedfilename = strcat(gcmfiledir, gcmtestname,'_',regionname, 'cmems_interped_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
        gcmcmemsvelfilename=strcat(cmip5dir,'\surface\', gcmtestname,'_',regionname, 'cmems_vel_rms_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
        
        valnum=0;
        run('C:\Users\user\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
        wrmap = bwrmap(51:100,:);
        
        

        run(param_script);
        
        figdir=[figrawdir,'Trend\SSH\'];
        if (exist(strcat(figdir) , 'dir') ~= 7)
            mkdir(strcat(figdir));
        end 
        
        outfile = strcat(figdir,regionname);
        trendoutfile = strcat(figdir,regionname);
        
        cmems_trend=ncread(cmemsfilename, 'cmems_trend');
        cmems_trend_filtered=ncread(cmemsfilename, 'cmems_trend_filtered');

        cmems_mask=ones(size(cmems_trend));
        cmems_mask(isnan(cmems_trend))=NaN;
        
        lon_rho=ncread(modelfilename,'lon_rho');
        lat_rho=ncread(modelfilename,'lat_rho');
        trend_filtered = ncread(modelfilename,'trend_filtered');
%         mean_trend_filtered = ncread(filename,'mean_trend_filtered');
        cmems_trend_filtered = ncread(cmemsfilename,'cmems_trend_filtered');
        interped_trend_filtered = ncread(interpedfilename,'interped_trend_filtered');

% start-------------------- make timedata(xData) for time series
        for i =1:length(inputyear) 
            tempyear=inputyear(i);
            for month=1:12
                xData((12*(i-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-01',]);
            end
        end
        

% start-------------------- model pattern correlation    
        fig_flag=fig_flags{43,2};
        while (fig_flag)
            cmems_adt=ncread(cmemsfilename,'cmems_adt');
            interped_ssh=ncread(interpedfilename,'interped_ssh');
            cmems_sla=ncread(cmemsfilename,'cmems_sla');
            interped_sla=ncread(interpedfilename,'interped_sla');
            gcm_interped_ssh = ncread(gcminterpedfilename, 'interped_ssh');
            m_gcm_interped_ssh = mean(gcm_interped_ssh,3,'omitnan');
            cut_lon_rho = ncread(cmemsfilename, 'lon_cmems');
            cut_lat_rho = ncread(cmemsfilename, 'lat_cmems');

            
%             month_season = [6,7,8];
%             for meani = 1:length(month_season)
%                 tempind{meani}= (inputyear-1982)*12 + month_season(meani);
%             end
%             ind_season = sort([tempind{:}]);
%             season_interped_sst = comb_interped_sst(:,:,ind_season);
%             comb_avhrr_sst = ncread(avhrrfilename, 'avhrr_sst');
%             season_avhrr_sst = comb_avhrr_sst(:,:,ind_season);
%             mean_avhrr_sst = mean(season_avhrr_sst, 3);
            
            
            gcm_mask = m_gcm_interped_ssh;
            gcm_mask(isfinite(gcm_mask))=1;
            m_cmems_adt=mean(cmems_adt,3,'omitnan');
%             m_interped_ssh=mean(interped_ssh,3,'omitnan').*gcm_mask-0.25;
            m_interped_ssh=mean(interped_ssh,3,'omitnan');
            m_cmems_sla=mean(cmems_sla,3,'omitnan');
            m_interped_sla=mean(interped_sla,3,'omitnan');
            m_interped_trend_filtered=mean(interped_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered))));
            m_cmems_trend_filtered=mean(cmems_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered))));
            
            ysecsmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, ysecs_khoapolygon);
            ysmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, yspolygon);
            esmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, es_khoapolygon);
            
            ysecs_m_interped_ssh = m_interped_ssh.*ysecsmask;
            ysecs_m_cmems_adt = m_cmems_adt.*ysecsmask;
            ysecs_m_gcm_interped_ssh = m_gcm_interped_ssh.*ysecsmask;

            es_m_interped_ssh = m_interped_ssh.*esmask;
            es_m_gcm_interped_ssh = m_gcm_interped_ssh.*esmask;
            es_m_cmems_adt = m_cmems_adt.*esmask;

% % % % % %             check 2d area for correlation
% % % % % %             m_interped_ssh_mask=m_interped_ssh;
% % % % % %             m_interped_ssh_mask(isfinite(m_interped_ssh_mask))=1;
% % % % % %             m_cmems_adt_mask=m_cmems_adt;
% % % % % %             m_cmems_adt_mask(isfinite(m_cmems_adt_mask))=1;
% % % % % %             pcolor((m_interped_ssh'-mean(m_interped_ssh', 'omitnan')).*m_interped_ssh_mask'.*m_cmems_adt_mask'); colorbar; caxis([-0.2 0.2]);
% % % % % %             
% % % % % %             m_gcm_interped_ssh_mask=m_gcm_interped_ssh;
% % % % % %             m_gcm_interped_ssh_mask(isfinite(m_gcm_interped_ssh_mask))=1;
% % % % % %             m_gcm_cmems_adt_mask=m_cmems_adt;
% % % % % %             m_gcm_cmems_adt_mask(isfinite(m_gcm_cmems_adt_mask))=1;
% % % % % %             pcolor((m_gcm_interped_ssh'-mean(m_gcm_interped_ssh', 'omitnan')).*m_gcm_interped_ssh_mask'.*m_gcm_cmems_adt_mask'); colorbar; caxis([-0.2 0.2]);
% % % % % %             
% % % % % %             ysecs_m_interped_ssh_mask=ysecs_m_interped_ssh;
% % % % % %             ysecs_m_interped_ssh_mask(isfinite(ysecs_m_interped_ssh_mask))=1;
% % % % % %             ysecs_m_cmems_adt_mask=ysecs_m_cmems_adt;
% % % % % %             ysecs_m_cmems_adt_mask(isfinite(ysecs_m_cmems_adt_mask))=1;
% % % % % %             pcolor((ysecs_m_interped_ssh'-mean(ysecs_m_interped_ssh', 'omitnan')).*ysecs_m_interped_ssh_mask'.*ysecs_m_cmems_adt_mask'); colorbar; caxis([-0.2 0.2]);
% % % % % %             
% % % % % %             ysecs_m_gcm_interped_ssh_mask=ysecs_m_gcm_interped_ssh;
% % % % % %             ysecs_m_gcm_interped_ssh_mask(isfinite(ysecs_m_gcm_interped_ssh_mask))=1;
% % % % % %             ysecs_m_gcm_cmems_adt_mask=ysecs_m_cmems_adt;
% % % % % %             ysecs_m_gcm_cmems_adt_mask(isfinite(ysecs_m_gcm_cmems_adt_mask))=1;
% % % % % %             pcolor((ysecs_m_gcm_interped_ssh'-mean(ysecs_m_gcm_interped_ssh', 'omitnan')).*ysecs_m_gcm_interped_ssh_mask'.*ysecs_m_gcm_cmems_adt_mask'); colorbar; caxis([-0.2 0.2]);
% % % % % %             
% % % % % %             pcolor((ysecs_m_cmems_adt'-mean(ysecs_m_cmems_adt', 'omitnan')).*ysecs_m_gcm_interped_ssh_mask'.*ysecs_m_gcm_cmems_adt_mask'); colorbar; caxis([-0.2 0.2]);
% % % % % % 
% % % % % %             es_m_interped_ssh_mask=es_m_interped_ssh;
% % % % % %             es_m_interped_ssh_mask(isfinite(es_m_interped_ssh_mask))=1;
% % % % % %             es_m_cmems_adt_mask=es_m_cmems_adt;
% % % % % %             es_m_cmems_adt_mask(isfinite(es_m_cmems_adt_mask))=1;
% % % % % %             pcolor((es_m_interped_ssh'-mean(es_m_interped_ssh', 'omitnan')).*es_m_interped_ssh_mask'.*es_m_cmems_adt_mask'); colorbar; caxis([-0.2 0.2]);
% % % % % %             
% % % % % %             es_m_gcm_interped_ssh_mask=es_m_gcm_interped_ssh;
% % % % % %             es_m_gcm_interped_ssh_mask(isfinite(es_m_gcm_interped_ssh_mask))=1;
% % % % % %             es_m_gcm_cmems_adt_mask=es_m_cmems_adt;
% % % % % %             es_m_gcm_cmems_adt_mask(isfinite(es_m_gcm_cmems_adt_mask))=1;
% % % % % %             pcolor((es_m_gcm_interped_ssh'-mean(es_m_gcm_interped_ssh', 'omitnan')).*es_m_gcm_interped_ssh_mask'.*es_m_gcm_cmems_adt_mask'); colorbar; caxis([-0.2 0.2]);
% % % % % %             
% % % % % % 
% % % % % %             es_m_gcm_cmems_adt_mask=es_m_cmems_adt;
% % % % % %             es_m_gcm_cmems_adt_mask(isfinite(es_m_gcm_cmems_adt_mask))=1;
% % % % % %             pcolor((es_m_cmems_adt'-mean(es_m_cmems_adt', 'omitnan')).*es_m_gcm_interped_ssh_mask'.*es_m_gcm_cmems_adt_mask'); colorbar; caxis([-0.2 0.2]);


            
            
%              pcolor(m_interped_ssh'-0.25); shading flat; colorbar; caxis([-0.2 1.0]); colormap jet
             aaa= movmean(m_interped_ssh(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_adt))),1);
             bbb= movmean(m_cmems_adt(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_adt))),1);
             corrcoef(aaa,bbb);
             ccc= movmean(m_gcm_interped_ssh(logical(~isnan(m_gcm_interped_ssh).*~isnan(m_cmems_adt))),1);
             ddd= movmean(m_cmems_adt(logical(~isnan(m_gcm_interped_ssh).*~isnan(m_cmems_adt))),1);
             corrcoef(ccc,ddd);
             
             eee= ysecs_m_interped_ssh(logical(~isnan(ysecs_m_interped_ssh).*~isnan(ysecs_m_cmems_adt)));
             fff= ysecs_m_cmems_adt(logical(~isnan(ysecs_m_interped_ssh).*~isnan(ysecs_m_cmems_adt)));
             ysecscorr=corrcoef(eee,fff)
             
             eeeg= ysecs_m_interped_ssh(logical(~isnan(ysecs_m_interped_ssh).*~isnan(ysecs_m_cmems_adt).*~isnan(ysecs_m_gcm_interped_ssh)));
             fffg= ysecs_m_cmems_adt(logical(~isnan(ysecs_m_interped_ssh).*~isnan(ysecs_m_cmems_adt).*~isnan(ysecs_m_gcm_interped_ssh)));
             rcmg_ysecscorr=corrcoef(eeeg,fffg)
             
             
             gg= es_m_interped_ssh(logical(~isnan(es_m_interped_ssh).*~isnan(es_m_cmems_adt)));
             hh= es_m_cmems_adt(logical(~isnan(es_m_interped_ssh).*~isnan(es_m_cmems_adt)));
             escorr=corrcoef(gg,hh)
%              ysecsescorr=corrcoef([eee; gg], [fff; hh])
                
             ggg= es_m_interped_ssh(logical(~isnan(es_m_interped_ssh).*~isnan(es_m_cmems_adt).*~isnan(es_m_gcm_interped_ssh)));
             hhg= es_m_cmems_adt(logical(~isnan(es_m_interped_ssh).*~isnan(es_m_cmems_adt).*~isnan(es_m_gcm_interped_ssh)));
             rcmg_escorr=corrcoef(ggg,hhg)
             
             ii= ysecs_m_gcm_interped_ssh(logical(~isnan(ysecs_m_gcm_interped_ssh).*~isnan(ysecs_m_cmems_adt)));
             jj= ysecs_m_cmems_adt(logical(~isnan(ysecs_m_gcm_interped_ssh).*~isnan(ysecs_m_cmems_adt)));
             gcmysecscorr=corrcoef(ii,jj)
             
             kk= es_m_gcm_interped_ssh(logical(~isnan(es_m_gcm_interped_ssh).*~isnan(es_m_cmems_adt)));
             ll= es_m_cmems_adt(logical(~isnan(es_m_gcm_interped_ssh).*~isnan(es_m_cmems_adt)));
             gcmescorr=corrcoef(kk,ll)
%              gcmysecsescorr=corrcoef([ii; kk], [jj; ll])
            
             all_ysecs_corr(testnameind2,1)=gcmysecscorr(1,2);
             all_ysecs_corr(testnameind2,2)=rcmg_ysecscorr(1,2);
             all_es_corr(testnameind2,1)=gcmescorr(1,2);
             all_es_corr(testnameind2,2)=rcmg_escorr(1,2);
             
%             plot(aaa-mean(aaa)); hold on; plot(bbb-mean(bbb)); plot(ccc-mean(ccc)-0.5); plot(ddd-mean(ddd)-0.5); legend('rcm', 'sat', 'gcm', 'sat-gcm'); hold off
%             plot(eee-mean(eee)); hold on; plot(fff-mean(fff)); plot(ii-mean(ii)-0.5); plot(jj-mean(jj)-0.5); legend('ysecs-rcm', 'ysecs-sat', 'ysecs-gcm', 'ysecs-sat-gcm'); hold off
%             plot(gg-mean(gg)); hold on; plot(hh-mean(hh)); plot(kk-mean(kk)-0.5); plot(ll-mean(ll)-0.5); legend('es-rcm', 'es-sat', 'es-gcm', 'es-sat-gcm'); hold off
            

%             corr2(m_interped_ssh(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_adt))), m_cmems_adt(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_adt))))
            corr2(m_interped_ssh(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_adt))), m_cmems_adt(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_adt))))
            corr2(interped_sla(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_sla))), m_cmems_sla(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_sla))))
%             corr2(interped_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered))) - m_interped_trend_filtered, ...
%             cmems_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered)))- m_cmems_trend_filtered)

        % end-------------------- model pattern correlation


%         % start-------------------- model pattern correlation (yearly)
%             pngname=strcat(outfile, '_', testname,'_',regionname, '_', nc_varname,'_','yearly_pattern_corr_',num2str(min(inputyear),'%04i'), ...
%                                 '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
%             if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
%                 for yearij=1:length(inputyear)
%                     tempyear=inputyear(yearij);
% 
%                     yearly_cmems_adt=ncread(filename,'cmems_adt',[1 1 (yearij-1)*12+1], [inf inf 12]);
%                     yearly_interped_ssh=ncread(filename,'interped_ssh',[1 1 (yearij-1)*12+1], [inf inf 12]);
%                     yearly_cmems_sla=ncread(filename,'cmems_sla',[1 1 (yearij-1)*12+1], [inf inf 12]);
%                     yearly_interped_sla=ncread(filename,'interped_sla',[1 1 (yearij-1)*12+1], [inf inf 12]);
% 
%                     m_cmems_adt=mean(yearly_cmems_adt,3,'omitnan');
%                     m_interped_ssh=mean(yearly_interped_ssh,3,'omitnan');
%                     m_cmems_sla=mean(yearly_cmems_sla,3,'omitnan');
%                     m_interped_sla=mean(yearly_interped_sla,3,'omitnan');
%         %             m_interped_trend_filtered=mean(interped_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered))));
%         %             m_cmems_trend_filtered=mean(cmems_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered))));
% 
%                     ssh_corr(yearij) = corr2(m_interped_ssh(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_adt))), m_cmems_adt(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_adt))));
%                     sla_corr(yearij) = corr2(interped_sla(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_sla))), m_cmems_sla(logical(~isnan(m_interped_ssh).*~isnan(m_cmems_sla))));
%         %             corr2(interped_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered))) - m_interped_trend_filtered, ...
%         %                 cmems_trend_filtered(logical(~isnan(interped_trend_filtered).*~isnan(cmems_trend_filtered)))- m_cmems_trend_filtered)
%                 end
%                 plot(ssh_corr);
%                 hold on
%                 plot(sla_corr);
%                 hold off
%                 close all
%             end
            fig_flag=0;
        end

% start-------------------- u_rms
        fig_flag=fig_flags{71,2};
        while (fig_flag)
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
             pngname=strcat(outfile, '_', testname,'_',regionname, '_cmems_', 'u_','rms_',num2str(min(inputyear),'%04i'), ...
            '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                run(param_script);
                if (exist('lon_min' , 'var') ~= 1)
                    lon_cmems=ncread(interpedfilename, 'lon_cmems');
                    lat_cmems=ncread(interpedfilename, 'lat_cmems');
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                    cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end
                rms_u_rho=ncread(cmemsvelfilename,'rms_u_rho');
                comb_interped_u_rho = ncread(cmemsvelfilename, 'interped_u_rho');
                mean_interped_u_rho = mean(comb_interped_u_rho, 3);
                comb_cmems_u_rho = ncread(cmemsvelfilename, 'cmems_u_rho');
                mean_cmems_u_rho = mean(comb_cmems_u_rho, 3);
                
                
                gcm_rms_u_rho = ncread(gcmcmemsvelfilename,'rms_u_rho');
                gcm_rms_mask = gcm_rms_u_rho;
                gcm_rms_mask(isfinite(gcm_rms_u_rho))=1;
                
                u_rho_corr = corrcoef(mean_interped_u_rho(logical(~isnan(mean_interped_u_rho).*~isnan(mean_cmems_u_rho))), mean_cmems_u_rho(logical(~isnan(mean_interped_u_rho).*~isnan(mean_cmems_u_rho))))
                u_rho_gcmcorr = corrcoef(mean_interped_u_rho(logical(~isnan(mean_interped_u_rho).*~isnan(mean_cmems_u_rho).*~isnan(gcm_rms_mask))), mean_cmems_u_rho(logical(~isnan(mean_interped_u_rho).*~isnan(mean_cmems_u_rho).*~isnan(gcm_rms_mask))))
                
                
                m_interped_u_rho = mean_interped_u_rho;
                m_cmems_u_rho=mean_cmems_u_rho;
                comb_gcm_interped_u_rho = ncread(gcmcmemsvelfilename, 'interped_u_rho');
                m_gcm_interped_u_rho = mean(comb_gcm_interped_u_rho,3);
                
                ysecsmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, ysecs_khoapolygon);
                ysmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, yspolygon);
                esmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, es_khoapolygon);

                ysecs_m_interped_u_rho = m_interped_u_rho.*ysecsmask;
                ysecs_m_cmems_u_rho = m_cmems_u_rho.*ysecsmask;
                ysecs_m_gcm_interped_u_rho = m_gcm_interped_u_rho.*ysecsmask;

                es_m_interped_u_rho = m_interped_u_rho.*esmask;
                es_m_gcm_interped_u_rho = m_gcm_interped_u_rho.*esmask;
                es_m_cmems_u_rho = m_cmems_u_rho.*esmask;
                
                eeeg= ysecs_m_interped_u_rho(logical(~isnan(ysecs_m_interped_u_rho).*~isnan(ysecs_m_cmems_u_rho).*~isnan(ysecs_m_gcm_interped_u_rho)));
                fffg= ysecs_m_cmems_u_rho(logical(~isnan(ysecs_m_interped_u_rho).*~isnan(ysecs_m_cmems_u_rho).*~isnan(ysecs_m_gcm_interped_u_rho)));
                rcmg_ysecscorr=corrcoef(eeeg,fffg)

                ggg= es_m_interped_u_rho(logical(~isnan(es_m_interped_u_rho).*~isnan(es_m_cmems_u_rho).*~isnan(es_m_gcm_interped_u_rho)));
                hhg= es_m_cmems_u_rho(logical(~isnan(es_m_interped_u_rho).*~isnan(es_m_cmems_u_rho).*~isnan(es_m_gcm_interped_u_rho)));
                rcmg_escorr=corrcoef(ggg,hhg)

                ii= ysecs_m_gcm_interped_u_rho(logical(~isnan(ysecs_m_gcm_interped_u_rho).*~isnan(ysecs_m_cmems_u_rho)));
                jj= ysecs_m_cmems_u_rho(logical(~isnan(ysecs_m_gcm_interped_u_rho).*~isnan(ysecs_m_cmems_u_rho)));
                gcmysecscorr=corrcoef(ii,jj)

                kk= es_m_gcm_interped_u_rho(logical(~isnan(es_m_gcm_interped_u_rho).*~isnan(es_m_cmems_u_rho)));
                ll= es_m_cmems_u_rho(logical(~isnan(es_m_gcm_interped_u_rho).*~isnan(es_m_cmems_u_rho)));
                gcmescorr=corrcoef(kk,ll)
                
                all_ysecs_u_corr(testnameind2,1)=gcmysecscorr(1,2);
                all_ysecs_u_corr(testnameind2,2)=rcmg_ysecscorr(1,2);
                all_es_u_corr(testnameind2,1)=gcmescorr(1,2);
                all_es_u_corr(testnameind2,2)=rcmg_escorr(1,2);
 
                clear lon_rho rms_u_rho
            end

            fig_flag=0;
        end

% start-------------------- v_rms
        fig_flag=fig_flags{72,2};
        while (fig_flag)
            figdir2=[figrawdir,'CLIM\'];
            if (exist(strcat(figdir2) , 'dir') ~= 7)
                mkdir(strcat(figdir2));
            end 
            outfile = strcat(figdir2,regionname);
             pngname=strcat(outfile, '_', testname,'_',regionname, '_cmems_', 'v_','rms_',num2str(min(inputyear),'%04i'), ...
            '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
            if (exist(pngname , 'file') ~= 2 || fig_flag==2)      
                run(param_script);
                if (exist('lon_min' , 'var') ~= 1)
                    lon_cmems=ncread(interpedfilename, 'lon_cmems');
                    lat_cmems=ncread(interpedfilename, 'lat_cmems');
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_cmems, lat_cmems);
                    cut_lon_rho = lon_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_cmems(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end
                rms_v_rho=ncread(cmemsvelfilename,'rms_v_rho');
                comb_interped_v_rho = ncread(cmemsvelfilename, 'interped_v_rho');
                mean_interped_v_rho = mean(comb_interped_v_rho, 3);
                comb_cmems_v_rho = ncread(cmemsvelfilename, 'cmems_v_rho');
                mean_cmems_v_rho = mean(comb_cmems_v_rho, 3);
                
                gcm_rms_v_rho = ncread(gcmcmemsvelfilename,'rms_v_rho');
                gcm_rms_mask = gcm_rms_v_rho;
                gcm_rms_mask(isfinite(gcm_rms_v_rho))=1;
                
                v_rho_corr = corrcoef(mean_interped_v_rho(logical(~isnan(mean_interped_v_rho).*~isnan(mean_cmems_v_rho))), mean_cmems_v_rho(logical(~isnan(mean_interped_v_rho).*~isnan(mean_cmems_v_rho))));
                v_rho_gcmcorr = corrcoef(mean_interped_v_rho(logical(~isnan(mean_interped_v_rho).*~isnan(mean_cmems_v_rho).*~isnan(gcm_rms_mask))), mean_cmems_v_rho(logical(~isnan(mean_interped_v_rho).*~isnan(mean_cmems_v_rho).*~isnan(gcm_rms_mask))))

                
                m_interped_v_rho = mean_interped_v_rho;
                m_cmems_v_rho=mean_cmems_v_rho;
                comb_gcm_interped_v_rho = ncread(gcmcmemsvelfilename, 'interped_v_rho');
                m_gcm_interped_v_rho = mean(comb_gcm_interped_v_rho,3);
                
                ysecsmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, ysecs_khoapolygon);
                ysmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, yspolygon);
                esmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, es_khoapolygon);

                ysecs_m_interped_v_rho = m_interped_v_rho.*ysecsmask;
                ysecs_m_cmems_v_rho = m_cmems_v_rho.*ysecsmask;
                ysecs_m_gcm_interped_v_rho = m_gcm_interped_v_rho.*ysecsmask;

                es_m_interped_v_rho = m_interped_v_rho.*esmask;
                es_m_gcm_interped_v_rho = m_gcm_interped_v_rho.*esmask;
                es_m_cmems_v_rho = m_cmems_v_rho.*esmask;
                
                eeeg= ysecs_m_interped_v_rho(logical(~isnan(ysecs_m_interped_v_rho).*~isnan(ysecs_m_cmems_v_rho).*~isnan(ysecs_m_gcm_interped_v_rho)));
                fffg= ysecs_m_cmems_v_rho(logical(~isnan(ysecs_m_interped_v_rho).*~isnan(ysecs_m_cmems_v_rho).*~isnan(ysecs_m_gcm_interped_v_rho)));
                rcmg_ysecscorr=corrcoef(eeeg,fffg)

                ggg= es_m_interped_v_rho(logical(~isnan(es_m_interped_v_rho).*~isnan(es_m_cmems_v_rho).*~isnan(es_m_gcm_interped_v_rho)));
                hhg= es_m_cmems_v_rho(logical(~isnan(es_m_interped_v_rho).*~isnan(es_m_cmems_v_rho).*~isnan(es_m_gcm_interped_v_rho)));
                rcmg_escorr=corrcoef(ggg,hhg)

                ii= ysecs_m_gcm_interped_v_rho(logical(~isnan(ysecs_m_gcm_interped_v_rho).*~isnan(ysecs_m_cmems_v_rho)));
                jj= ysecs_m_cmems_v_rho(logical(~isnan(ysecs_m_gcm_interped_v_rho).*~isnan(ysecs_m_cmems_v_rho)));
                gcmysecscorr=corrcoef(ii,jj)

                kk= es_m_gcm_interped_v_rho(logical(~isnan(es_m_gcm_interped_v_rho).*~isnan(es_m_cmems_v_rho)));
                ll= es_m_cmems_v_rho(logical(~isnan(es_m_gcm_interped_v_rho).*~isnan(es_m_cmems_v_rho)));
                gcmescorr=corrcoef(kk,ll)
                
                all_ysecs_v_corr(testnameind2,1)=gcmysecscorr(1,2);
                all_ysecs_v_corr(testnameind2,2)=rcmg_ysecscorr(1,2);
                all_es_v_corr(testnameind2,1)=gcmescorr(1,2);
                all_es_v_corr(testnameind2,2)=rcmg_escorr(1,2);
                
% 
%                 for loni=1:size(model_sla,1)
%                     for lati=1:size(model_sla,2)
%                         model_sla_std(loni,lati)=std(model_sla(loni,lati,:));
%                     end
%                 end
%                 cmems_land=ones(size(model_sla));
%                 cmems_land(isnan(model_sla))=1;
%                 cmems_land(isfinite(model_sla))=NaN;
%                 save([filedir,regionname, '_', testname, '_cmems_std_land','.mat'], 'cmems_land');
                

%                 mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
%                 mask_model(mask_model==0)=NaN;
%                 rms_u_rho=rms_u_rho .* mask_model;

%                 yearstr_min=num2str(inputyear(1));
%                 yearstr_max=num2str(inputyear(end));
%                 save([filedir,regionname, '_cmems_linear_', 'det_std','_',yearstr_min, '_', yearstr_max, '.mat'], 'cut_lon_rho', 'cut_lat_rho', 'model_sla_std');
                

                close all;
                clear lon_rho rms_v_rho
            end

            fig_flag=0;
        end
        
    end
end