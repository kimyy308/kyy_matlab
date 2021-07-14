close all; clear all;  clc;
warning off;

all_var ={'zeta'};
for varind=1:length(all_var)
    close all;
    clearvars '*' -except varind all_var regionind2 all_region2
    % % % 
    testname='test49';
    inputyear = [1993:2008]; % % put year which you want to plot [year year ...]
    inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
    varname =all_var{varind};  %% reference variable -> temperature
    Vtransform=2; Vstretching=4; theta_s=10; theta_b=1; N=40; hc=250;
    var_unit='m';
    
    system_name=computer;
    if (strcmp(system_name,'PCWIN64'))
        % % for windows
        dropboxpath='C:\Users\KYY\Dropbox';
        addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
        addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
        addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
        addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
        addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path
    end

% % %     filename=['E:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',num2str(inputyear(1),'%04i'),'\test42_monthly_',num2str(inputyear(1),'%04i'),'_01.nc'];
% % %     model_depth = ncread(filename,'h',[1 1], [1 920]);
    
% % %     % zeta = ncread(filename,'zeta',[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
% % %     model_depth_z=zlevs(Vtransform,Vstretching,model_depth,0,theta_s,theta_b,hc,N,'w');
% % %     model_depth_z2=squeeze(diff(model_depth_z,1));
% % %     model_depth_z3=repmat(model_depth_z2',1,1,12*length(inputyear));
    
% % % %     load ('E:\Data\Model\ROMS\nwp_1_20\test42\run\NWP_u_all_time_1980_2008.mat')
% % % %     discharge_u=squeeze(mean(mean(comb_west_data,1,'omitnan'),2,'omitnan'));
% % % %     comb_west_data_u=comb_west_data;
% % %     
% % % %     load (['E:\Data\Model\ROMS\nwp_1_20\',testname,'\run\NWP_v_all_time_1980_2008.mat'])
% % % %     % discharge_v=squeeze(mean(mean(comb_west_data,1,'omitnan'),2,'omitnan'));
% % % %     comb_west_data_v=comb_west_data;
% % % 
% % %     load (['E:\Data\Model\ROMS\nwp_1_20\',testname,'\run\NWP_',varname,'_all_time_', ...
% % %         num2str(inputyear(1),'%04i'),'_',num2str(inputyear(end),'%04i'),'.mat'])
% % % 
% % %     if (strcmp(system_name,'PCWIN64'))
% % %         % % for windows
% % %         figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\3rd_year\figure\nwp_1_20\',testname,'\NIFS\'); % % where figure files will be saved
% % %         param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m'
% % %         filedir = strcat('E:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
% % %         run(param_script);
% % %     elseif (strcmp(system_name,'GLNXA64'))
% % %     end
% % % 
% % %     figdir=[figrawdir,'Trend\'];
% % % 
% % %     if (exist(strcat(figdir) , 'dir') ~= 7)
% % %         mkdir(strcat(figdir));
% % %     end 
% % %     outfile = strcat(figdir);
% % % 
% % %     % plot(squeeze(mean(mean(comb_west_data,1,'omitnan'),2,'omitnan')))
% % % 
% % %     hor_dist=m_lldist([115 115], [15 15.0483])*1000;
% % % 
% % %     % discharge=squeeze(mean(mean(comb_west_data.*comb_west_data_v.*model_depth_z3.*hor_dist,1,'omitnan'),2,'omitnan'));
% % % 
% % %     % discharge=squeeze(sum(sum(comb_west_data.*comb_west_data_u.*model_depth_z3.*hor_dist,1,'omitnan'),2,'omitnan'));
% % %     discharge=squeeze(sum(sum(comb_west_data.*model_depth_z3.*hor_dist,1,'omitnan'),2,'omitnan'));
% % % 
% % %     % % ...   ./squeeze(sum(sum(model_depth_z3*hor_dist,1),2));
% % % 
% % % 
% % % 
% % % 
% % %     reshap_discharge=reshape(discharge,[12,length(inputyear)]);
% % % 
% % %     seasonal_mean=mean(reshap_discharge,2);
% % % 
% % %     for i=1:length(inputyear) 
% % %         reshap_discharge2(:,i)=reshap_discharge(:,i)-seasonal_mean;
% % %     end
% % %     discharge2=reshape(reshap_discharge2,[1 12*length(inputyear)])';
% % % 
% % %     for i =1:length(inputyear) 
% % %         tempyear=inputyear(i);
% % %         for month=1:12
% % %             xData((12*(i-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-15',]);
% % %         end
% % %     end
% % % 
% % %     % plot(xData,discharge_v.*discharge)
% % %     % plot(xData,discharge_v)
% % % 
% % % 
% % % 
% % %     p=polyfit(xData,(discharge)',1);
% % %     line_trend2=xData*p(1)+p(2);
% % %     msaltplot=plot(xData,(discharge)','color','k');
% % %     hold on
% % %     msaltplot2=plot(xData,line_trend2,'Color','r');
% % % 
% % %     hold off
% % %     jpgname=strcat(outfile,testname, '_salt_west_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
% % %     xlabel('year')
% % %     ylabel('sum of salt at the boundary (meter^2)')
% % %     title(['western boundary salinity, ',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i')])
% % %     datetick('x','yymmm','keepticks')
% % %     axis tight;
% % %     % ylim(meanplotlev)
% % %     set(msaltplot,'LineWidth',2);
% % %     set(msaltplot2,'LineWidth',2);
% % %     set(gca,'FontSize',20);
% % %     % set(gca,'YTick',min(meanplotlev):1:max(meanplotlev))
% % %     % txt1=text(xData(3), min(meanplotlev)+1 ,['R = ', num2str(round(constant_cor(1,2),2)), ', '], 'FontSize', m_quiver_ref_text_fontsize); 
% % %     lgd=legend(['salt transport']);
% % % 
% % %     set(gcf, 'PaperUnits', 'points');
% % %             set(gcf, 'PaperSize', [1000, 400]);
% % %             set(gcf,'PaperPosition', [0 0 1000 400])
% % % 
% % %     saveas(gcf,jpgname,'jpg');                              

    load (['E:\Data\Model\ROMS\nwp_1_20\input\',testname,'\NWP_',varname,'_all_time_', ...
        num2str(inputyear(1),'%04i'),'_',num2str(inputyear(end),'%04i'),'.mat'])
    
    for i =1:length(inputyear) 
        tempyear=inputyear(i);
        for month=1:12
            xData((12*(i-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-15',]);
        end
    end
    
    mdata_west=mean(comb_west_data_bndy,1);
    p=polyfit(xData,(mdata_west),1);
    line_trend2=xData*p(1)+p(2);
    mdataplot=plot(xData,(mdata_west),'color','k');
    hold on
    mdataplot2=plot(xData,line_trend2,'Color','r');

    hold off
    jpgname=strcat(outfile,testname, '_',varname,'_west_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
    xlabel('year')
    ylabel(['mean ', varname,'(',var_unit,')'])
    title(['western boundary ', varname, ', ', ...
        num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), ...
        ',',num2str(round(p(1)*1000*365,2)), 'mm/yr'])
    datetick('x','yyyy','keepticks')
    axis tight;
    % ylim(meanplotlev)
    set(mdataplot,'LineWidth',2);
%     set(msaltplot2,'LineWidth',2);
    set(gca,'FontSize',20);
    % set(gca,'YTick',min(meanplotlev):1:max(meanplotlev))
    % txt1=text(xData(3), min(meanplotlev)+1 ,['R = ', num2str(round(constant_cor(1,2),2)), ', '], 'FontSize', m_quiver_ref_text_fontsize); 
    lgd=legend(varname);

    set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [1000, 400]);
            set(gcf,'PaperPosition', [0 0 1000 400])

    saveas(gcf,jpgname,'jpg');        
    
    
    
    mdata_east=mean(comb_east_data_bndy,1);
    p=polyfit(xData,(mdata_east),1);
    line_trend2=xData*p(1)+p(2);
    mdataplot=plot(xData,(mdata_east),'color','k');
    hold on
    mdataplot2=plot(xData,line_trend2,'Color','r');

    hold off
    jpgname=strcat(outfile,testname, '_',varname,'_east_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
    xlabel('year')
    ylabel(['mean ', varname,'(',var_unit,')'])
    title(['eastern boundary ', varname, ', ', ...
        num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), ...
        ',',num2str(round(p(1)*1000*365,2)), 'mm/yr'])    
    datetick('x','yyyy','keepticks')
    axis tight;
    % ylim(meanplotlev)
    set(mdataplot,'LineWidth',2);
%     set(msaltplot2,'LineWidth',2);
    set(gca,'FontSize',20);
    % set(gca,'YTick',min(meanplotlev):1:max(meanplotlev))
    % txt1=text(xData(3), min(meanplotlev)+1 ,['R = ', num2str(round(constant_cor(1,2),2)), ', '], 'FontSize', m_quiver_ref_text_fontsize); 
    lgd=legend(varname);

    set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [1000, 400]);
            set(gcf,'PaperPosition', [0 0 1000 400])

    saveas(gcf,jpgname,'jpg');    
    
    
    mdata_south=mean(comb_south_data_bndy,1);
    p=polyfit(xData,(mdata_south),1);
    line_trend2=xData*p(1)+p(2);
    mdataplot=plot(xData,(mdata_south),'color','k');
    hold on
    mdataplot2=plot(xData,line_trend2,'Color','r');

    hold off
    jpgname=strcat(outfile,testname, '_',varname,'_south_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
    xlabel('year')
    ylabel(['mean ', varname,'(',var_unit,')'])
    title(['southern boundary ', varname, ', ', ...
        num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), ...
        ',',num2str(round(p(1)*1000*365,2)), 'mm/yr'])    
    datetick('x','yyyy','keepticks')
    axis tight;
    % ylim(meanplotlev)
    set(mdataplot,'LineWidth',2);
%     set(msaltplot2,'LineWidth',2);
    set(gca,'FontSize',20);
    % set(gca,'YTick',min(meanplotlev):1:max(meanplotlev))
    % txt1=text(xData(3), min(meanplotlev)+1 ,['R = ', num2str(round(constant_cor(1,2),2)), ', '], 'FontSize', m_quiver_ref_text_fontsize); 
    lgd=legend(varname);

    set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [1000, 400]);
            set(gcf,'PaperPosition', [0 0 1000 400])

    saveas(gcf,jpgname,'jpg');     
    
    
    
    mdata_north=mean(comb_north_data_bndy,1);
    p=polyfit(xData,(mdata_north),1);
    line_trend2=xData*p(1)+p(2);
    mdataplot=plot(xData,(mdata_north),'color','k');
    hold on
    mdataplot2=plot(xData,line_trend2,'Color','r');

    hold off
    jpgname=strcat(outfile,testname, '_',varname,'_north_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
    xlabel('year')
    ylabel(['mean ', varname,'(',var_unit,')'])
    title(['northern boundary ', varname, ', ', ...
        num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'), ...
        ',',num2str(round(p(1)*1000*365,2)), 'mm/yr'])    
    datetick('x','yyyy','keepticks')
    axis tight;
    % ylim(meanplotlev)
    set(mdataplot,'LineWidth',2);
%     set(msaltplot2,'LineWidth',2);
    set(gca,'FontSize',20);
    % set(gca,'YTick',min(meanplotlev):1:max(meanplotlev))
    % txt1=text(xData(3), min(meanplotlev)+1 ,['R = ', num2str(round(constant_cor(1,2),2)), ', '], 'FontSize', m_quiver_ref_text_fontsize); 
    lgd=legend(varname);

    set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [1000, 400]);
            set(gcf,'PaperPosition', [0 0 1000 400])

    saveas(gcf,jpgname,'jpg');    
end