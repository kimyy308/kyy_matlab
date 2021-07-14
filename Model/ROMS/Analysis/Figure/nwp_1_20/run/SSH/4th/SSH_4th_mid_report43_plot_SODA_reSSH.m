close all; clear all;  clc;
warning off;
% all_region2 ={'NWP','NWP2'}
% all_region2 ={'ES', 'YS', 'SS', 'NWP2'}
% all_testname2 = {'ens02','test53','test55','test56'};
% all_testname2 = {'test54'};

% all_region2 ={'NWP','ES', 'SS', 'YS', 'AKP'}
all_region2 ={'AKP2'};

% all_region2 ={'NWP'}
for regionind2=1:length(all_region2)

% for testnameind2=1:length(all_testname2)
%     
% end
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
        meanplotlev =[-0.3 0.3];

        % for snu_desktopd
%         testname=all_testname2{testnameind2}    % % need to change
        inputyear = [1980:2008]; % % put year which you want to plot [year year ...]
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
    
        testname='test42';
        load(['H:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',testname,'_',regionname, ...
            'ssh_recon_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
        msldata{1}=comb_interped_data_filtered;
        
        testname='test49';
        load(['H:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',testname,'_',regionname, ...
            'ssh_recon_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
        msldata{2}=comb_interped_data_filtered;

        valnum=0;

        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\4th_year\figure\nwp_1_20\','all','\',regionname,'\'); % % where figure files will be saved
            param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m'
%             filedir = strcat('H:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
            recondir='E:\Data\Observation\OISST\monthly\';
        elseif (strcmp(system_name,'GLNXA64'))
        end
        std(mean(mean(comb_interped_data_filtered,'omitnan'),'omitnan'))
        std(mean(mean(comb_recon_data_filtered,'omitnan'),'omitnan'))


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
        
        recon_msl=squeeze(mean(mean(comb_recon_data_filtered,1,'omitnan'),2,'omitnan'));
        recon_msl=recon_msl-mean(recon_msl);        
        recon_p=polyfit(xData,recon_msl',1);
        recon_msl2=xData*recon_p(1)+recon_p(2);
        
%         msl=squeeze(mean(mean(comb_interped_data_filtered,1,'omitnan'),2,'omitnan'));
%         msl=msl-mean(msl);    
%         p=polyfit(xData,msl',1);
%         msl2=xData*p(1)+p(2);
        
        hold on
%         mslplot{5}=plot(xData,msl{5},'k');
%         recon_mslplot=plot(xData,recon_msl,'r')
%         lgd=legend('Model Ensemble','Reconstructed SSH');
        for i=1:length(msldata)
            mslplot{i}=plot(xData,msl{i},'k');
            mslplot2{i}=plot(xData,msl2{i},'Color','r');
            set(mslplot{i}, 'Color', 'k')
            set(mslplot2{i}, 'Color', 'k')
            set(get(get(mslplot{i},'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
            set(get(get(mslplot2{i},'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
%         set(get(get(mslplot{5},'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
        recon_mslplot=plot(xData,recon_msl,'r')
        recon_mslplot2=plot(xData,recon_msl2,'Color','r')
        set(get(get(recon_mslplot,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
        set(get(get(recon_mslplot2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

        
        jpgname=strcat(outfile, '_', 'all', '_',regionname, '_msl_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
        xlabel('year')
        ylabel('Mean SSH (m)')
        title([regionname, ', Mean SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') '])
% %         num2str(round(mean_trend_filtered,2)), ' mm/y'
        datetick('x','yyyy','keepticks')
        axis tight;
        ylim(meanplotlev)
        set(mslplot{5}, 'Color', 'k')
        set(mslplot2{5}, 'Color', 'k')
        set(recon_mslplot2, 'Color', 'r')
        set(mslplot{5},'LineWidth',2);
        set(recon_mslplot,'LineWidth',2);
%         tr_ens=round(ppp{5}(1)*365*1000,2);
        tr_recon=round(recon_p(1)*365*1000,2);
        lgd=legend(['Model : ', num2str(tr_ens),'mm/year'], ['Reconstructed SSH : ', num2str(tr_recon),'mm/year']);
%         lgd=legend(mslplot{5},'Model Ensemble','Location', 'northwest');
%         lgd2=legappend(recon_mslplot,'Reconstructed SSH','Location','northeast');
        set(gca,'FontSize',20);
        grid on
        hold off
        constant_cor=corrcoef(msl{5},recon_msl);
        txt1=text(xData(5), min(meanplotlev)+diff(meanplotlev)/16.0 ,['R = ', num2str(round(constant_cor(1,2),2))], 'FontSize', 20); 
        set(gcf,'PaperPosition', [0 0 40 18]) 
        saveas(gcf,jpgname,'jpg');
        grid off

%         recon_msl=squeeze(mean(mean(comb_recon_data_filtered,1,'omitnan'),2,'omitnan'));
%         recon_msl=recon_msl-mean(recon_msl);        
%         p=polyfit(xData,recon_msl',1);
%         recon_msl2=xData*p(1)+p(2);
%         recon_mslplot=plot(xData,recon_msl,'k')
%         hold on
%         recon_mslplot2=plot(xData,recon_msl2,'Color','r')
%         jpgname=strcat(outfile, '_', testname, '_',regionname,'_recon_msl_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
%         xlabel('year')
%         ylabel('Mean SSH (m)')
%         title([regionname ', Mean recon SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(mean_recon_trend_filtered,2)), ' mm/y'])
%         datetick('x','yymmm','keepticks')
%         axis tight;
%         ylim(meanplotlev)
%         set(recon_mslplot,'LineWidth',2);
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



%         comb_yearly_recon_data_filtered=reshape(comb_recon_data_filtered,[isize, jsize, 12, lsize/12]);
%         mean_yearly_recon_data_filtered=squeeze(mean(mean(mean(comb_yearly_recon_data_filtered,1,'omitnan'),2,'omitnan'),3,'omitnan'));
%         trendtime=1:length(xData2);
%         p=polyfit(trendtime,mean_yearly_recon_data_filtered(1:length(xData2))',1);
%         yearly_recon_trend=p(1);
%         yearly_recon_trend = yearly_recon_trend * 1000.0; %% m/y -> mm/y
% 
%         yearly_recon_msl=squeeze(mean(mean(mean(comb_yearly_recon_data_filtered,1,'omitnan'),2,'omitnan'),3,'omitnan'));
%         yearly_recon_msl=yearly_recon_msl-mean(yearly_recon_msl);        
%         p=polyfit(xData2,yearly_recon_msl',1);
%         yearly_recon_msl2=xData2*p(1)+p(2);
%         yearly_recon_mslplot=plot(xData2,yearly_recon_msl,'k','Marker','o','MarkerFaceColor','k','MarkerSize',4)
%         hold on
%         yearly_recon_mslplot2=plot(xData2,yearly_recon_msl2,'Color','r')
%         jpgname=strcat(outfile, '_', testname, '_',regionname,'_yearly_recon_msl_',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
%         xlabel('Year')
%         ylabel('Mean SSH (m)')
%         title([regionname ', Mean recon SSH(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(yearly_recon_trend,2)), ' mm/y'])
%         ylim(meanplotlev)
%         datetick('x','yymmm','keepticks')
%         axis tight;
%         ylim(meanplotlev)
%         set(yearly_recon_mslplot,'LineWidth',2);
%         set(gca,'FontSize',15);
%         grid on
%         hold off
%         saveas(gcf,jpgname,'jpg');
%         grid off
        close all;
    end
