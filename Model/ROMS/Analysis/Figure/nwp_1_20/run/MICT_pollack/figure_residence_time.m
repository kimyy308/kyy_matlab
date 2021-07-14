close all; clear all; clc;

system_name=computer;
if (strcmp(system_name,'PCWIN64'))
    % % for windows
    dropboxpath='C:\Users\KYY\Dropbox'; %% SNU_desktop, kyy_laptop
    addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
    addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
    addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run']));

elseif (strcmp(system_name,'GLNXA64'))
    dropboxpath='/home/kimyy/Dropbox'; %% DAMO
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
    addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_20/run']));
end

warning off;
% % set calendar name
calendarname=cell(1,12); calendarname{1} = 'January'; calendarname{2} = 'February'; calendarname{3} = 'March'; calendarname{4} = 'April'; calendarname{5} = 'May'; calendarname{6} = 'June';
calendarname{7} = 'July'; calendarname{8} = 'August'; calendarname{9} = 'September'; calendarname{10} = 'October'; calendarname{11} = 'November'; calendarname{12} = 'December';

% workdir='/data1/kimyy/Model/LTRANS/LTRANSv2b/LTRANSv2b/output/seo_ens_mean/';

totmonth=1:3;
xlimday=90;
inputyear=1982:199;
totmonth=1:3; 
tind=1;
for yearij = 1:length(inputyear)
    tempyear = inputyear(yearij);
    ftime(tind) = datenum(tempyear,totmonth(1),15);
%         ftime(tind) = datenum(tempyear,totmonth(1),15) - datenum(1900,12,31);
    tind=tind+1;
end
scalename={'e_folding', 'half', '2_3', 'tenth'};
titlename={'37%(e-folding)','50%(half)','67%(2/3)','10%(1/10)'};
testname='nwp_1_10_EnOI';
 switch(testname)
        case('seo_nwp') 
            dirname='avg_ens_10km_mean';
        case('nwp_1_10_EnOI') 
            dirname='nwp_1_10_EnOI';
        otherwise
            ('?')
 end

figdir=['D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\3rd_year\figure\',dirname,'\LTRANS\residence_time\'];
if (exist(strcat(figdir) , 'dir') ~= 7)
    mkdir(strcat(figdir));
end 

for scaleind=1:length(scalename)
    filename=[testname,'_comb_res_time_',scalename{scaleind},'.mat']
    load(filename);
    for nmonth=1:length(totmonth)
        month=totmonth(nmonth);
        figname=[figdir, scalename{scaleind},'_', num2str(month,'%04i'),'.jpg'];
        resplot=plot(ftime, squeeze(comb_res_time(:,month)) ,'color','k','Marker','o','MarkerFaceColor','k','MarkerSize',4)
        xlabel('year')
        ylabel('Residence time (day)')
        title(['(Residence time','-',titlename{scaleind},'-',calendarname{month},')'])
        set(gca,'YTick',(0:10:100));
        set(gca,'XTick',(ftime(1):366.0:ftime(end)));
        datetick('x','yy','keepticks')
        axis tight;
        ylim([0 100])
        set(resplot,'LineWidth',2);
        grid on
        set(gcf,'PaperPosition', [0 0 36 12]) 
        set(gca,'FontSize',13);
        saveas(gcf,figname,'jpg');
        grid off
    end
        figname=[figdir, scalename{scaleind},'_mean','.jpg'];
        resplot=plot(ftime, squeeze(mean(comb_res_time(:,:),2)) ,'color','k','Marker','o','MarkerFaceColor','k','MarkerSize',4)
        xlabel('year')
        ylabel('Residence time (day)')
        title(['(Residence time','-',titlename{scaleind},'-','Mean(Jan-Mar)',')'])
        set(gca,'YTick',(0:10:100));
        set(gca,'XTick',(ftime(1):366.0:ftime(end)));
        datetick('x','yy','keepticks')
        axis tight;
        ylim([0 100])
        set(resplot,'LineWidth',2);
        grid on
        set(gcf,'PaperPosition', [0 0 36 12]) 
        set(gca,'FontSize',13);
        saveas(gcf,figname,'jpg');
        grid off
end


