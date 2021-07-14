close all; clear all;  clc;
warning off;

GCM_testnames  = {'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'NorESM1-M', 'MPI-ESM-LR'};
RCM_testnames  = {'test57', 'test58', 'test59', 'test60'};
scenname='rcp45';

all_region2 ={'AKP4'};

varname='zeta';
variable = 'SSH';
meanplotlev2 =[-0.1 0.8] .*100.0;


% scenname='rcp45';
for regionind2=1:length(all_region2)
        close all;

        inputyear = [2006:2100]; % % put year which you want to plot [year year ...]
        yearstr_min=num2str(inputyear(1));
        yearstr_max=num2str(inputyear(end));
        inputmonth = [1:12]; % % put month which you want to plot [month month ...]
        system_name=computer;
        for folding=1:1
            if (strcmp(system_name,'PCWIN64'))
                % % for windows
                dropboxpath='C:\Users\User\Dropbox';
                addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
                addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
                addpath(genpath([dropboxpath '\source\matlab\Common\cptcmap']));
                addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
                addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
                addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Roms_tools\Preprocessing_tools']));
                addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path
            elseif (strcmp(system_name,'GLNXA64'))
                dropboxpath='/home/kimyy/Dropbox';
                addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
                addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
                addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
                addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
            end
            
            run('nwp_polygon_point.m');
            regionname=all_region2{regionind2};
            switch(regionname)
            case('NWP')
                refpolygon=nwppolygon;
            case('AKP4') %% Around Korea Peninsula
                refpolygon=akp4polygon;
            end
            lonlat(1)=min(refpolygon(:,1));
            lonlat(2)=max(refpolygon(:,1));
            lonlat(3)=min(refpolygon(:,2));
            lonlat(4)=max(refpolygon(:,2));
        end
       
        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\5th_year\figure\paper2\fin_fig\'); % % where figure files will be saved
        elseif (strcmp(system_name,'GLNXA64'))
        end
        param_script ='C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m';
        run(param_script);

        correction_right_fig=[-0.1600,0,0,0]; % right, up, width, height
        correction_upper_fig=[0,0,0,0];
        correction_large_fig=[0,0,0,0.020];
        hold on;
        tifname=strcat(figrawdir, 'fig08','.tif'); %% ~_year_month.jpg

%         f1=figure(1);
        byrmap = customcolormap_preset('red-yellow-blue');
        yrmap = byrmap(129:256,:);
        

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
            RCM_interpedfilenames{testind} = strcat(filedir, testname,'_',regionname, 'cmems_interped_ssh_trend_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');

            testname=GCM_testnames{testind};
            filedir=strcat('D:\Data\Model\CMIP5\zos\', scenname, '\interp\', testname, '\');
            GCM_modelfilenames{testind} = strcat(filedir, testname,'_',regionname, '_ssh_trend_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
            GCM_interpedfilenames{testind} = strcat(filedir, testname,'_',regionname, 'cmems_interped_ssh_trend_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');

            cmemsdir='Z:\내 드라이브\Data\Observation\CMEMS\';
            cmems_filename = strcat(filedir, testname,'_',regionname, 'cmems_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');

        end
        
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
        end
        
% start-------------------- make timedata for time series
        for i =1:length(inputyear) 
            tempyear=inputyear(i);
            for month=1:12
                xData((12*(i-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-01',]);
            end
            xData2(i) = datenum([num2str(tempyear),'-',num2str(6,'%02i'),'-30',]);
        end
% end-------------------- make timedata for time series  
        
% start-------------------- GCM rcp 4.5 time series
        testnameind=1;
        sb{testnameind}=subplot(4,4,[1 2 5 6]);  % Auto-fitted to the figure.
        pos_sb{testnameind}=get(sb{testnameind}, 'pos'); % Get the position.delete(sb1); % Delete the subplot axes
        delete(sb{testnameind}); % Delete the subplot axes
        ax{testnameind}=axes;
        set(ax{testnameind},'pos',pos_sb{testnameind});
        
         hold on
        cmems_1y_mean=0;                        
        for testind=1:length(RCM_testnames)
            GCM_interped_1y_mean(testind)=mean(GCM_interped_sla_yearly_mean(testind,1:5));
            GCM_interped_sla_yearly_mean(testind,:)=GCM_interped_sla_yearly_mean(testind,:)+(cmems_1y_mean-GCM_interped_1y_mean(testind));
            mslplot_all{testind}=plot(xData2,GCM_interped_sla_yearly_mean(testind,:),'b', 'parent', ax{testnameind}, 'LineWidth', 2);
        end

        set(mslplot_all{1},'Marker','*');
        set(mslplot_all{2},'Marker','^');
        set(mslplot_all{3},'Marker','o');
        set(mslplot_all{4},'Marker','+');

        xlabel('Year')
        ylabel('Mean SSH (cm)')
        datetick('x','yyyy','keepticks')
        axis tight;
        ylim(meanplotlev2)
        set(ax{testnameind},'FontSize',m_grid_fontsize);
        grid on
        lgd{testnameind}=legend(ax{testnameind}, 'GCM-IPSL-LR','GCM-IPSL-MR','GCM-Nor','GCM-MPI', 'Location','northwest');
        set(lgd{testnameind},'FontSize',m_grid_fontsize);
%         set(lgd{testnameind},'Position',[0.07 0.83, 0.3, 0.03]);
%         set(lgd,'Orientation','horizontal');  
        
        text(xData2(3), 50, '(A) GCM, RCP 4.5', 'FontSize', m_grid_fontsize)
        
        
% start-------------------- RCM rcp4.5 time series
        
        testnameind=2;
        sb{testnameind}=subplot(4,4,[9 10 13 14]);  % Auto-fitted to the figure.
        pos_sb{testnameind}=get(sb{testnameind}, 'pos'); % Get the position.delete(sb1); % Delete the subplot axes
        delete(sb{testnameind}); % Delete the subplot axes
        ax{testnameind}=axes;
        set(ax{testnameind},'pos',pos_sb{testnameind});
        
         hold on
        cmems_1y_mean=0;                        
        for testind=1:length(RCM_testnames)
            RCM_interped_1y_mean(testind)=mean(RCM_interped_sla_yearly_mean(testind,1:5));
            RCM_interped_sla_yearly_mean(testind,:)=RCM_interped_sla_yearly_mean(testind,:)+(cmems_1y_mean-RCM_interped_1y_mean(testind));
            mslplot_all{testind}=plot(xData2,RCM_interped_sla_yearly_mean(testind,:),'r', 'parent', ax{testnameind}, 'LineWidth', 2);
        end

        set(mslplot_all{1},'Marker','*');
        set(mslplot_all{2},'Marker','^');
        set(mslplot_all{3},'Marker','o');
        set(mslplot_all{4},'Marker','+');

        xlabel('Year')
        ylabel('Mean SSH (cm)')
        datetick('x','yyyy','keepticks')
        axis tight;
        ylim(meanplotlev2)
        set(ax{testnameind},'FontSize',m_grid_fontsize);
        grid on
        lgd{testnameind}=legend(ax{testnameind}, 'RCM-IPSL-LR','RCM-IPSL-MR','RCM-Nor','RCM-MPI', 'Location','northwest');
        set(lgd{testnameind},'FontSize',m_grid_fontsize);
%         set(lgd{testnameind},'Position',[0.07 0.83, 0.3, 0.03]);
%         set(lgd,'Orientation','horizontal');  
        text(xData2(3), 50, '(C) RCM, RCP 4.5', 'FontSize', m_grid_fontsize)
      
                
        
        
        scenname ='rcp85';
        RCM_testnames  = {'test65', 'test66', 'test67', 'test68'};

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
            RCM_interpedfilenames{testind} = strcat(filedir, testname,'_',regionname, 'cmems_interped_ssh_trend_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');

            testname=GCM_testnames{testind};
            filedir=strcat('D:\Data\Model\CMIP5\zos\', scenname, '\interp\', testname, '\');
            GCM_modelfilenames{testind} = strcat(filedir, testname,'_',regionname, '_ssh_trend_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
            GCM_interpedfilenames{testind} = strcat(filedir, testname,'_',regionname, 'cmems_interped_ssh_trend_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');

            cmemsdir='Z:\내 드라이브\Data\Observation\CMEMS\';
            cmems_filename = strcat(filedir, testname,'_',regionname, 'cmems_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');

            end

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
        end
        
        
        
        % start-------------------- GCM rcp 8.5 time series
        testnameind=3;
        sb{testnameind}=subplot(4,4,[3 4 7 8]);  % Auto-fitted to the figure.
        pos_sb{testnameind}=get(sb{testnameind}, 'pos'); % Get the position.delete(sb1); % Delete the subplot axes
        delete(sb{testnameind}); % Delete the subplot axes
        ax{testnameind}=axes;
        set(ax{testnameind},'pos',pos_sb{testnameind});
        
         hold on
        cmems_1y_mean=0;                        
        for testind=1:length(RCM_testnames)
            GCM_interped_1y_mean(testind)=mean(GCM_interped_sla_yearly_mean(testind,1:5));
            GCM_interped_sla_yearly_mean(testind,:)=GCM_interped_sla_yearly_mean(testind,:)+(cmems_1y_mean-GCM_interped_1y_mean(testind));
            mslplot_all{testind}=plot(xData2,GCM_interped_sla_yearly_mean(testind,:),'b', 'parent', ax{testnameind}, 'LineWidth', 2);
        end

        set(mslplot_all{1},'Marker','*');
        set(mslplot_all{2},'Marker','^');
        set(mslplot_all{3},'Marker','o');
        set(mslplot_all{4},'Marker','+');

        xlabel('Year')
        ylabel('Mean SSH (cm)')
        datetick('x','yyyy','keepticks')
        axis tight;
        ylim(meanplotlev2)
        set(ax{testnameind},'FontSize',m_grid_fontsize);
        set(ax{testnameind},'YAxisLocation','right');
        grid on
        lgd{testnameind}=legend(ax{testnameind}, 'GCM-IPSL-LR','GCM-IPSL-MR','GCM-Nor','GCM-MPI', 'Location','northwest');
        set(lgd{testnameind},'FontSize',m_grid_fontsize);
%         set(lgd{testnameind},'Position',[0.07 0.83, 0.3, 0.03]);
%         set(lgd,'Orientation','horizontal');  
        
        text(xData2(3), 50, '(B) GCM, RCP 8.5', 'FontSize', m_grid_fontsize)
  
        
             % start-------------------- RCM rcp 8.5 time series
        testnameind=4;
        sb{testnameind}=subplot(4,4,[11 12 15 16]);  % Auto-fitted to the figure.
        pos_sb{testnameind}=get(sb{testnameind}, 'pos'); % Get the position.delete(sb1); % Delete the subplot axes
        delete(sb{testnameind}); % Delete the subplot axes
        ax{testnameind}=axes;
        set(ax{testnameind},'pos',pos_sb{testnameind});
        
         hold on
        cmems_1y_mean=0;                        
        for testind=1:length(RCM_testnames)
            RCM_interped_1y_mean(testind)=mean(RCM_interped_sla_yearly_mean(testind,1:5));
            RCM_interped_sla_yearly_mean(testind,:)=RCM_interped_sla_yearly_mean(testind,:)+(cmems_1y_mean-RCM_interped_1y_mean(testind));
            mslplot_all{testind}=plot(xData2,RCM_interped_sla_yearly_mean(testind,:),'r', 'parent', ax{testnameind}, 'LineWidth', 2);
        end

        set(mslplot_all{1},'Marker','*');
        set(mslplot_all{2},'Marker','^');
        set(mslplot_all{3},'Marker','o');
        set(mslplot_all{4},'Marker','+');

        xlabel('Year')
        ylabel('Mean SSH (cm)')
        datetick('x','yyyy','keepticks')
        axis tight;
        ylim(meanplotlev2)
        set(ax{testnameind},'FontSize',m_grid_fontsize);
        set(ax{testnameind},'YAxisLocation','right');
        grid on
        lgd{testnameind}=legend(ax{testnameind}, 'RCM-IPSL-LR','RCM-IPSL-MR','RCM-Nor','RCM-MPI', 'Location','northwest');
        set(lgd{testnameind},'FontSize',m_grid_fontsize);
%         set(lgd{testnameind},'Position',[0.07 0.83, 0.3, 0.03]);
%         set(lgd,'Orientation','horizontal');  
        
        text(xData2(3), 50, '(D) RCM, RCP 8.5', 'FontSize', m_grid_fontsize)
        
        
        
        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width*2.2 paper_position_height*2.4]) 

%         saveas(gcf,tifname,'tif');
        print('-dtiff','-r500',tifname)

        hold off;
        close all;
        
end