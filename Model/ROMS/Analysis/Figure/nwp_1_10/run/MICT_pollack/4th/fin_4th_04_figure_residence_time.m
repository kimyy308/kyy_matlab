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

testname='nwp_1_10_EnOI';
% testname='avg_ens_10km_mean';

if (strcmp(testname,'avg_ens_10km_mean')==1)
    workdir=['E:\Data\Model\LTRANS\output\test06_DA_3rd\'];
else
    workdir=['/data1/kimyy/Model/LTRANS/LTRANSv2b/LTRANSv2b/output/',testname,'/'];
end

totmonth=1:3;
xlimday=90;
inputyear=1983:2019;
totmonth=1:3;
tind=1;
for yearij = 1:length(inputyear)
    tempyear = inputyear(yearij);
    ftime(tind) = datenum(tempyear,totmonth(1),15);
%         ftime(tind) = datenum(tempyear,totmonth(1),15) - datenum(1900,12,31);
    tind=tind+1;
end

% scalename={'e_folding', 'half', '2_3', 'tenth'};
titlename={'37%(e-folding)','50%(half)','67%(2/3)','10%(1/10)'};

scalename={'e_folding'};
% titlename={'37%(e-folding)'};

juv_year = [1975, 1976, 1977, 1978, 1979, ...
1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, ...
1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997];

juv_catch = [55293, 82512, 104212, 92974, 68120, 68236, ...
    115493, 99190, 56237, 66737, 38029, 32466, 13550, 2890, 8427, ...
16518, 10111, 5038, 7563, 3141, 2261, 3823, 910];
diff_juv_catch=diff(juv_catch);

for scaleind=1:length(scalename)
%     if (strcmp(testname,'avg_ens_10km_mean')==1)
        filename=[testname,'_comb_res_time_',scalename{scaleind},'_',num2str(min(inputyear)),'_',num2str(max(inputyear)),'.mat']
%     else
%         filename=[testname,'_comb_res_time_',scalename{scaleind},'.mat']
%     end
    load(filename);
    for nmonth=1:length(totmonth)
        month=totmonth(nmonth);
        figname=['D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\4th_year\Figure\',testname,'\LTRANS\residence_time\', ...
            scalename{scaleind},'_',num2str(inputyear(1),'%04i'),'_',num2str(inputyear(end),'%04i'),'_', num2str(month,'%04i'),'.jpg'];
        resplot=plot(ftime, squeeze(comb_res_time(inputyear(1)-1982:inputyear(end)-1982,month)) ,'color','k','Marker','o','MarkerFaceColor','k','MarkerSize',4)
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
        set(gcf,'PaperPosition', [0 0 1.2*length(inputyear) 12]) 
        set(gca,'FontSize',13);
        saveas(gcf,figname,'jpg');
        grid off
    end
        figname=['D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\4th_year\Figure\',testname,'\LTRANS\residence_time\', ...
            scalename{scaleind},'_',num2str(inputyear(1),'%04i'),'_',num2str(inputyear(end),'%04i'),'_mean_ylim_30','.jpg'];
        
%         yyaxis left
        resplot=plot(ftime, squeeze(mean(comb_res_time(inputyear(1)-1982:inputyear(end)-1982,:),2)) ,'color','k','Marker','o','MarkerFaceColor','k','MarkerSize',4)
        xlabel('year')
        title(['(Residence time','-',titlename{scaleind},'-','Mean(Jan-Mar)',')'])
%         set(gca,'YTick',(0:10:100));
        set(gca,'YTick',(0:10:90));
        set(gca,'XTick',(ftime(1):366.0:ftime(end)));
        datetick('x','yy','keepticks')
        axis tight;
        ylim([0 90])
        ylabel('Residence time (day)')
        set(resplot,'LineWidth',2);

% %         yyaxis right  
% %         comb_res_catch(1:length(inputyear))=NaN;
% %         for yearij = 1:length(inputyear)
% %             tempyear = inputyear(yearij);
% %             catch_ind= find(juv_year==tempyear);
% %             if (catch_ind>0)
% %                 comb_res_catch(yearij)= diff_juv_catch(catch_ind-1);
% %             end
% %         end
% % %         if (inputyear(end)<=1997)
% % %             juv_catch(find(juv_year==inputyear(1))-1:find(juv_year==inputyear(end))-1);
% % %         else
% % %             juv_catch(find(juv_year==inputyear(1))-1:end);
% % %         end
% %         resplot2=plot(ftime, comb_res_catch ,'color','r','Marker','o','MarkerFaceColor','r','MarkerSize',4);
% %         ylim([-50000 150000])
% %         ylabel('Juvenile catch change,compared to last year (t)')
% %         set(resplot2,'LineWidth',1);
% %         hold on   
% %         resplot3=plot(ftime,zeros(length(inputyear),1),'r','linewidth',1,'MarkerFaceColor','r','linestyle','--');
% %         hold off
        

        grid on
%         set(gcf,'PaperPosition', [0 0 2.5*length(inputyear) 12]) 
        set(gca,'FontSize',13);
        saveas(gcf,figname,'jpg');
        grid off
end

% mean(mean(comb_res_time(16:20,:)))
