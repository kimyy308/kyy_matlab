close all; clear all;  clc;
warning off;

% all_region2 ={'NWP','ES', 'SS', 'YS', 'AKP'}
% all_region2 ={'AKP2'};
all_region2 ={'NES', 'SES', 'YS'}
for regionind2=1:length(all_region2)
        close all;
        clearvars '*' -except regionind2 all_region2
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
        conlev  = 0:5:35;
        meanplotlev =[-0.12 0.12];

        % for snu_desktopd
%         testname=all_testname2{testnameind2}    % % need to change
        inputyear = [1993:2017]; % % put year which you want to plot [year year ...]
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
    
        
%         filename = ['E:\Data\Observation\CMEMS\',regionname, ...
%             'cmems_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];
%         cmems_msl=ncread(filename,'cmems_ssh_filtered');
            load(['E:\Data\Observation\CMEMS\',regionname, ...
        'cmems_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);

        testname='test11';
        filename = ['E:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
            '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];
        msldata{1}=ncread(filename,'ssh_filtered');
        tr_model(1)=ncread(filename,'mean_trend_filtered');
        
        
        valnum=0;

        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\4th_year\figure\nwp_1_10\','all','\',regionname,'\'); % % where figure files will be saved
            param_script =['C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_', regionname, '.m']
%             filedir = strcat('H:\Data\Model\ROMS\nwp_1_10\', testname, '\run\'); % % where data files are
        elseif (strcmp(system_name,'GLNXA64'))
        end
%         std(mean(mean(comb_interped_data_filtered,'omitnan'),'omitnan'))
%         std(mean(mean(comb_cmems_data_filtered,'omitnan'),'omitnan'))


        figdir=[figrawdir,'Trend\'];
        if (exist(strcat(figdir) , 'dir') ~= 7)
            mkdir(strcat(figdir));
        end 
        outfile = strcat(figdir,regionname);
    

        for i =1:length(inputyear) 
            tempyear=inputyear(i);
            for month=1:12
                xData((12*(i-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-01',]);
            end
        end
        
        for i=1:length(msldata)
            msl{i}=  squeeze(mean(mean(msldata{i},1,'omitnan'),2,'omitnan'));
            msl{i}=msl{i}-mean(msl{i});
            ppp{i}=polyfit(xData,msl{i}',1);
            msl2{i}=xData*ppp{i}(1)+ppp{i}(2);
        end
        msl{1}=msl{1}-msl2{1}';
        
        cmems_msl=squeeze(mean(mean(comb_cmems_data_filtered,1,'omitnan'),2,'omitnan'));
        cmems_msl=cmems_msl-mean(cmems_msl);        
        cmems_p=polyfit(xData,cmems_msl',1);
        cmems_msl2=xData*cmems_p(1)+cmems_p(2);
        
        cmems_msl=cmems_msl-cmems_msl2';
%         msl=squeeze(mean(mean(comb_interped_data_filtered,1,'omitnan'),2,'omitnan'));
%         msl=msl-mean(msl);    
%         p=polyfit(xData,msl',1);
%         msl2=xData*p(1)+p(2);
        
        hold on
%         mslplot{5}=plot(xData,msl{5},'k');
%         cmems_mslplot=plot(xData,cmems_msl,'r')
%         lgd=legend('Model Ensemble','Reconstructed SSH');
        for i=1:length(msldata)
            mslplot{i}=plot(xData,msl{i},'k');
%             mslplot2{i}=plot(xData,msl2{i},'Color','r');
%             set(mslplot{i}, 'Color', [0.8 0.8 0.8])
            set(mslplot{i}, 'Color', 'k')

%             set(mslplot2{i}, 'Color', 'k')
            set(get(get(mslplot{i},'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
%             set(get(get(mslplot2{i},'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
%         set(get(get(mslplot{5},'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
        cmems_mslplot=plot(xData,cmems_msl,'r')
%         cmems_mslplot2=plot(xData,cmems_msl2,'Color','r')
        set(get(get(cmems_mslplot,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
%         set(get(get(cmems_mslplot2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

        
        jpgname=strcat(outfile, '_', testname, '_',regionname, '_detrended_msl_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
        xlabel('year')
        ylabel('Mean SSH (m)')
        title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') '])
% %         num2str(round(mean_trend_filtered,2)), ' mm/y'
        datetick('x','yyyy','keepticks')
        axis tight;
        ylim(meanplotlev)
%         set(mslplot{5}, 'Color', 'k')
%         set(mslplot2{5}, 'Color', 'k')
%         set(cmems_mslplot2, 'Color', 'r')
%         set(mslplot{5},'LineWidth',2);
        set(cmems_mslplot,'LineWidth',2);
%         tr_ens=round(ppp{5}(1)*365*1000,2);
        tr_cmems=round(cmems_p(1)*365*1000,2);
        lgd=legend(['SODA+CMEMS bndy Model : ', num2str(round(tr_model(1),2)),'mm/year'], ['CMEMS SSH : ', num2str(tr_cmems),'mm/year']);
%         lgd=legend(mslplot{5},'Model Ensemble','Location', 'northwest');
%         lgd2=legappend(cmems_mslplot,'Reconstructed SSH','Location','northeast');
        set(gca,'FontSize',20);
        grid on
        hold off
        constant_cor=corrcoef(msl{1},cmems_msl);
        txt1=text(xData(5), min(meanplotlev)+diff(meanplotlev)/16.0 ,['R = ', num2str(round(constant_cor(1,2),2))], 'FontSize', 20); 
        set(gcf,'PaperPosition', [0 0 40 18]) 
        saveas(gcf,jpgname,'jpg');
        grid off

%         cmems_msl=squeeze(mean(mean(comb_cmems_data_filtered,1,'omitnan'),2,'omitnan'));
%         cmems_msl=cmems_msl-mean(cmems_msl);        
%         p=polyfit(xData,cmems_msl',1);
%         cmems_msl2=xData*p(1)+p(2);
%         cmems_mslplot=plot(xData,cmems_msl,'k')
%         hold on
%         cmems_mslplot2=plot(xData,cmems_msl2,'Color','r')
%         jpgname=strcat(outfile, '_', testname, '_',regionname,'_cmems_msl_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
%         xlabel('year')
%         ylabel('Mean SSH (m)')
%         title([regionname ', Mean cmems SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(mean_cmems_trend_filtered,2)), ' mm/y'])
%         datetick('x','yymmm','keepticks')
%         axis tight;
%         ylim(meanplotlev)
%         set(cmems_mslplot,'LineWidth',2);
%         set(gca,'FontSize',15);
%         grid on
%         hold off
%         saveas(gcf,jpgname,'jpg');
%         grid off
%         close all;



%         for i =1:length(inputyear) 
%             tempyear=inputyear(i);
%     %         for month=1:12
%                 xData2(i) = datenum([num2str(6,'%02i'),'-30-',num2str(tempyear)]);
%     %         end
%         end
%         %     plot(squeeze(mean(mean(comb_interped_data_filtered,1,'omitnan'),2,'omitnan')))
%         isize = size(comb_interped_data_filtered,1)
%         jsize = size(comb_interped_data_filtered,2)
%         lsize = size(comb_interped_data_filtered,3)
%         comb_yearly_interped_data_filtered=reshape(comb_interped_data_filtered,[isize, jsize, 12, lsize/12]);
%         mean_yearly_interped_data_filtered=squeeze(mean(mean(mean(comb_yearly_interped_data_filtered,1,'omitnan'),2,'omitnan'),3,'omitnan'));
%         trendtime=1:length(xData2);
%         p=polyfit(trendtime,mean_yearly_interped_data_filtered(1:length(xData2))',1);
%         yearly_interped_trend=p(1);
%         yearly_interped_trend = yearly_interped_trend * 1000.0; %% m/y -> mm/y
% 
%         yearly_msl=squeeze(mean(mean(mean(comb_yearly_interped_data_filtered,1,'omitnan'),2,'omitnan'),3,'omitnan'));
%         yearly_msl=yearly_msl-mean(yearly_msl);    
%         p=polyfit(xData2,yearly_msl',1);
%         yearly_msl2=xData2*p(1)+p(2);
%         yearly_mslplot=plot(xData2,yearly_msl,'k','Marker','o','MarkerFaceColor','k','MarkerSize',4)
%         hold on
%         yearly_mslplot2=plot(xData2,yearly_msl2,'Color','r')
%         jpgname=strcat(outfile, '_', testname, '_',regionname, '_yearly_msl_',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
%         xlabel('Year')
%         ylabel('Mean SSH (m)')
%         title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(yearly_interped_trend,2)), ' mm/y'])
%         ylim(meanplotlev)
%         datetick('x','yymmm','keepticks')
%         axis tight;
%         ylim(meanplotlev)
%         set(yearly_mslplot,'LineWidth',2);
%         set(gca,'FontSize',15);
%         grid on
%         hold off
%         saveas(gcf,jpgname,'jpg');
%         grid off
%         close all;
% 
% 
%         yearly_corrected_msl=squeeze(mean(mean(mean(comb_yearly_interped_data_filtered,1,'omitnan'),2,'omitnan'),3,'omitnan'));
%         for i=1:length(xData2)
%             yearly_corrected_msl(i)=yearly_corrected_msl(i)+(i-1)*0.0017;
%         end
%         yearly_corrected_msl=yearly_corrected_msl-mean(yearly_corrected_msl);    
%         p=polyfit(xData2,yearly_corrected_msl',1);
%         yearly_corrected_msl2=xData2*p(1)+p(2);
%         yearly_corrected_mslplot=plot(xData2,yearly_corrected_msl,'k','Marker','o','MarkerFaceColor','k','MarkerSize',4)
%         hold on
%         yearly_corrected_mslplot2=plot(xData2,yearly_corrected_msl2,'Color','r')
%         jpgname=strcat(outfile, '_', testname, '_',regionname, '_yearly_corrected_msl_',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
%         xlabel('Year')
%         ylabel('Mean SSH (m)')
%         title([regionname, ',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), ...
%             '), ',num2str(round(yearly_interped_trend+1.7,2)),' (',num2str(round(yearly_interped_trend,2)),'+',num2str(round(1.7,2)),')', ' mm/y'])
%         ylim(meanplotlev)
%         datetick('x','yymmm','keepticks')
%         axis tight;
%         ylim(meanplotlev)
%         set(yearly_corrected_mslplot,'LineWidth',2);
%         set(gca,'FontSize',15);
%         grid on
%         hold off
%         saveas(gcf,jpgname,'jpg');
%         grid off
%         close all;



%         comb_yearly_cmems_data_filtered=reshape(comb_cmems_data_filtered,[isize, jsize, 12, lsize/12]);
%         mean_yearly_cmems_data_filtered=squeeze(mean(mean(mean(comb_yearly_cmems_data_filtered,1,'omitnan'),2,'omitnan'),3,'omitnan'));
%         trendtime=1:length(xData2);
%         p=polyfit(trendtime,mean_yearly_cmems_data_filtered(1:length(xData2))',1);
%         yearly_cmems_trend=p(1);
%         yearly_cmems_trend = yearly_cmems_trend * 1000.0; %% m/y -> mm/y
% 
%         yearly_cmems_msl=squeeze(mean(mean(mean(comb_yearly_cmems_data_filtered,1,'omitnan'),2,'omitnan'),3,'omitnan'));
%         yearly_cmems_msl=yearly_cmems_msl-mean(yearly_cmems_msl);        
%         p=polyfit(xData2,yearly_cmems_msl',1);
%         yearly_cmems_msl2=xData2*p(1)+p(2);
%         yearly_cmems_mslplot=plot(xData2,yearly_cmems_msl,'k','Marker','o','MarkerFaceColor','k','MarkerSize',4)
%         hold on
%         yearly_cmems_mslplot2=plot(xData2,yearly_cmems_msl2,'Color','r')
%         jpgname=strcat(outfile, '_', testname, '_',regionname,'_yearly_cmems_msl_',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
%         xlabel('Year')
%         ylabel('Mean SSH (m)')
%         title([regionname ', Mean cmems SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(yearly_cmems_trend,2)), ' mm/y'])
%         ylim(meanplotlev)
%         datetick('x','yymmm','keepticks')
%         axis tight;
%         ylim(meanplotlev)
%         set(yearly_cmems_mslplot,'LineWidth',2);
%         set(gca,'FontSize',15);
%         grid on
%         hold off
%         saveas(gcf,jpgname,'jpg');
%         grid off
        close all;
    end
