close all; clear all;  clc;
% % horizontal SST trend plot (connected with SSH_mid_report1)
warning off;
% all_region2 ={'AKP','NWP','ES', 'SS', 'YS'}
% all_region2 ={'AKP'};

% all_testname2 = {'test54'};

all_region2 ={'AKP2'}
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
    elseif (strcmp(system_name,'GLNXA64'))
        dropboxpath='/home/kimyy/Dropbox';
        addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
        addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
        addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
        addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
    end

    shadlev = [0 35];
    rms_shadlev = [0 4];
    trendlev = [-0.1 0.1];
    conlev  = 0:5:35;
    meanplotlev =[-5 35];


    % for snu_desktop
%     testname=all_testname2{testnameind2}
    inputyear = [1982:2005]; % % put year which you want to plot [year year ...]
    inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]

        varname ='temp'
        run('nwp_polygon_point.m');
        regionname=all_region2{regionind2};
        switch(regionname)
        case('NWP') %% North western Pacific
            lonlat = [115, 164, 15, 52];  %% whole data area
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
        
        testname='test53';
        load(['G:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',regionname,'sst_rms_and_bias_',num2str(inputyear(1),'%04i'),'_',num2str(inputyear(end)),'.mat']);
        msstdata{1}=comb_interped_data;
        
        testname='test54';
        load(['G:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',regionname,'sst_rms_and_bias_',num2str(inputyear(1),'%04i'),'_',num2str(inputyear(end)),'.mat']);
        msstdata{2}=comb_interped_data;
        
        testname='test55';
        load(['G:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',regionname,'sst_rms_and_bias_',num2str(inputyear(1),'%04i'),'_',num2str(inputyear(end)),'.mat']);
        msstdata{3}=comb_interped_data;
        
        testname='test56';
        load(['G:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',regionname,'sst_rms_and_bias_',num2str(inputyear(1),'%04i'),'_',num2str(inputyear(end)),'.mat']);
        msstdata{4}=comb_interped_data;
        
        testname='ens03';
        load(['G:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',regionname,'sst_rms_and_bias_',num2str(inputyear(1),'%04i'),'_',num2str(inputyear(end)),'.mat']);
        msstdata{5}=comb_interped_data;
        
            valnum=0;
    % %     valid cell number
%          for vi=1:size(comb_spatial_meanavhrr,1)
%              for vj=1:size(comb_spatial_meanavhrr,2)
%                  if (isnan(comb_spatial_meanavhrr(vi,vj,1))~=1)
%                     valnum=valnum+1;
%                  end
%              end
%          end

        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\4th_year\figure\nwp_1_20\','all','\',regionname,'\'); % % where figure files will be saved
            param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m'
            filedir = strcat('G:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
            avhrrdir='E:\Data\Observation\OISST\monthly\';
        elseif (strcmp(system_name,'GLNXA64'))
        end


        trendlev = [-0.1 0.1];

        figdir=[figrawdir,'Trend\'];
        if (exist(strcat(figdir) , 'dir') ~= 7)
            mkdir(strcat(figdir));
        end 
        outfile = strcat(figdir,regionname);

        for i =1:length(inputyear) 
            tempyear=inputyear(i);
            for month=1:12
                xData((12*(i-1))+month) = datenum([num2str(month,'%02i'),'-01-',num2str(tempyear)]);
            end
        end
        
        for i=1:length(msstdata)
            msst{i}=  squeeze(mean(mean(msstdata{i},1,'omitnan'),2,'omitnan'));
%             msst{i}=msst{i}-mean(msst{i});
            ppp{i}=polyfit(xData,msst{i}',1);
            msst2{i}=xData*ppp{i}(1)+ppp{i}(2);
        end
        
        avhrr_msst=squeeze(mean(mean(comb_avhrr_data,1,'omitnan'),2,'omitnan'));
%         avhrr_msst=avhrr_msst-mean(avhrr_msst);        
        avhrr_p=polyfit(xData,avhrr_msst',1);
        avhrr_msst2=xData*avhrr_p(1)+avhrr_p(2);
        
        hold on
        
        for i=1:length(msstdata)
            msstplot{i}=plot(xData,msst{i},'k');
            msstplot2{i}=plot(xData,msst2{i},'Color','r');
            set(msstplot{i}, 'Color', [0.8 0.8 0.8])
            set(msstplot2{i}, 'Color', [0.8 0.8 0.8])
            set(get(get(msstplot{i},'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            set(get(get(msstplot2{i},'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
        set(get(get(msstplot{5},'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
        avhrr_msstplot=plot(xData,avhrr_msst,'r')
        avhrr_msstplot2=plot(xData,avhrr_msst2,'Color','r')
        set(get(get(avhrr_msstplot,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
        set(get(get(avhrr_msstplot2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        
        


%         msst=squeeze(mean(mean(comb_interped_data,1,'omitnan'),2,'omitnan'));
        % msst=msst-mean(msst);    
%         p=polyfit(xData,msst',1);
%         msst2=xData*p(1)+p(2);
%         msstplot=plot(xData,msst,'k')
%         hold on
%         msstplot2=plot(xData,msst2,'Color','r')
        jpgname=strcat(outfile, '_', 'all','_',regionname, '_msst_seasonal', '.jpg'); %% ~_year_month.jpg
        xlabel('month')
        ylabel('Mean SSTA (^oC)')
        title([regionname, ', Mean SSTA(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), '])
%         num2str(round(mean_trend_filtered,3)), ' ^oC/y'
%         set(msstplot,'LineWidth',2);
        set(msstplot{5}, 'Color', 'k')
        set(msstplot2{5}, 'Color', 'k')
        set(avhrr_msstplot2, 'Color', 'r')
        set(msstplot{5},'LineWidth',2);
        set(avhrr_msstplot,'LineWidth',2);
        tr_ens=round(ppp{5}(1)*365,3);
        tr_avhrr=round(avhrr_p(1)*365,3);
        lgd=legend(['Model Ensemble : ', num2str(tr_ens),'^oC/year'], ['AVHRR : ', num2str(tr_avhrr),'^oC/year']);

        datetick('x','yyyy','keepticks')     
        set(gca,'FontSize',20);
        axis tight;
        ylim(meanplotlev);
        grid on
        hold off
        constant_cor=corrcoef(msst{5},avhrr_msst);
        txt1=text(xData(5), min(meanplotlev)+diff(meanplotlev)/16.0 ,['R = ', num2str(round(constant_cor(1,2),2))], 'FontSize', 20); 
        set(gcf,'PaperPosition', [0 0 40 18]) 
        saveas(gcf,jpgname,'jpg');
        grid off
        close all;

%         avhrr_msst=squeeze(mean(mean(comb_avhrr_data,1,'omitnan'),2,'omitnan'));
%         % avhrr_msst=avhrr_msst-mean(avhrr_msst);        
%         p=polyfit(xData,avhrr_msst',1);
%         avhrr_msst2=xData*p(1)+p(2);
%         avhrr_msstplot=plot(xData,avhrr_msst,'k')
%         hold on
%         avhrr_msstplot2=plot(xData,avhrr_msst2,'Color','r')
%         jpgname=strcat(outfile, '_', testname,'_',regionname, '_avhrr_msst', '.jpg'); %% ~_year_month.jpg
%         xlabel('month')
%         ylabel('Mean SSTA (^oC)')
%         title([regionname,', Mean AVHRR SSTA(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(mean_avhrr_trend_filtered,3)), '^oC/y'])
%         set(avhrr_msstplot,'LineWidth',2);
%         datetick('x','yyyy','keepticks')
%         set(gca,'FontSize',15);
%         axis tight
%         ylim(meanplotlev);
%         grid on
%         hold off
%         saveas(gcf,jpgname,'jpg');
%         grid off
        close all;
    end


% SSH_4th_mid_report6