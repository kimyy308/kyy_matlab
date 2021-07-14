% % Updated 04-Jun-2018 by Yong-Yub Kim (get observed monthly transport)


clc;close all;clear all;
warning off;

linux=0; windows=1;

if (windows ==1)
    % % for windows
    dropboxpath='C:\Users\KYY\Dropbox'; %% SNU_desktop
    addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
    addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
    addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run']));
elseif (linux==1)
    % % for linux
    dropboxpath='/home01/kimyy/Dropbox'; %% ROMS server
%     dropboxpath='/home/kimyy/Dropbox'; %% DAMO server
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
    addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_20/run']));
end

testname = 'test42';
refyear = 1980;
year = [1980:2005];
month = [1:12]; % % put month which you want to plot [month month ...]
% inputdir = strcat(['E:\Data\Model\ROMS\nwp_1_20\test37\input\']);
% outputdir = strcat(['E:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',num2str(year,'%04i'),'\']); % % where data files are

for i=1:length(year)
    tempyear=num2str(year(i),'%04i');
    filename = ['E:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',tempyear,'\nwp_1_20_monthly_',testname,'_',tempyear,'.txt'];
    startRow = 2;
    formatSpec = '%8f%9f%9f%9f%9f%9f%s%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    dataArray{7} = strtrim(dataArray{7});
    fclose(fileID);
    t_korea = dataArray{:, 1};
    t_tsugaru = dataArray{:, 2};
    t_soyat = dataArray{:, 3};
    t_aiwankur = dataArray{:, 4};
    t_o_intruy = dataArray{:, 5};
    t_ellowsea = dataArray{:, 6};
%     runyear=num2str(year(i)-refyear+1,'%04i');
    for month=1:12
% %     get model transport and set time
        xData((12*(i-1))+month) = datenum([num2str(month,'%02i'),'-01-',tempyear]);
        korea((12*(i-1))+month) = t_korea(month);
        tsugaru((12*(i-1))+month) = t_tsugaru(month);
        soya((12*(i-1))+month) = t_soyat(month);
        taiwan((12*(i-1))+month) = t_aiwankur(month);
        onshore((12*(i-1))+month) = t_o_intruy(month);
        yellow((12*(i-1))+month) = t_ellowsea(month);
    end
end
ES_diff=korea-tsugaru-soya;
ECS_diff=korea-taiwan+onshore-yellow;
% % get observation data (from .xls)

[~, ~, raw] = xlsread('E:\Data\Observation\Transport_Korea\obs_tsushima_197001_201209_ORIGIN.xls','TWC','A2:F517');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % 숫자형이 아닌 셀 찾기
raw(R) = {NaN}; % 숫자형이 아닌 셀 바꾸기

obs_ks = reshape([raw{:}],size(raw));

clearvars raw R;

year_obs_30y=obs_ks(find(obs_ks(:,1)<2011 & obs_ks(:,1)>1979),1);
month_obs_30y=obs_ks(find(obs_ks(:,1)<2011 & obs_ks(:,1)>1979),2);
tr_obs_30y=obs_ks(find(obs_ks(:,1)<2011 & obs_ks(:,1)>1979),5);

korea_obs=tr_obs_30y(find(year_obs_30y >= year(1) & year_obs_30y <= year(length(year))))';    %% get observation transport at the korea strait

korea_obs_run=korea_obs;
% for i=1:length(year)
%     korea_obs_run((12*(i-1))+1:(12*(i-1))+12)=tr_obs_30y(find(year_obs_30y == refyear))';    %% get observation transport at the korea strait
% end




% % 1970, Jan ~ 2012, sep  korea strait transport observation climatology data
clim =[1.853 2.057 2.383 2.627 2.726 2.693 2.821 3.039 3.002 3.047 2.761 2.257];
% % Fukudome et al, 2010 climatology mean feb 1997 to feb 2007 
fuku_clim =[2.01 2.28 2.52 2.60 2.66 2.70 2.68 3.05 2.93 3.10 2.84 2.36];
clim_std = [0.274404 0.274396 0.224372 0.218135 0.195140 0.211314 0.219307 0.288883 0.262485 0.274527 0.297884 0.31302];

clearvars filename startRow formatSpec fileID dataArray ans;




% % % East Sea transport (run)

es4plot=plot(xData, korea_obs_run,'b');
hold on
es1plot=plot(xData, korea,'k');
es2plot=plot(xData, tsugaru,'r');
es3plot=plot(xData, soya,'y');
es5plot=plot(xData, ES_diff,'g');
datetick('x',12,'keeplimits')
% datetick('x',12,'keeplimits')

% diff(1:12*year)=korea(1:12*year)-soyat(1:12*year)-tsugaru(1:12*year);

% es4plot=plot(diff(1:12*year),'y');

% es1plot=plot(korea(1:72),'k');
% hold on
% es2plot=plot(tsugaru(1:72),'r');
% es3plot=plot(soyat(1:72),'b');
% diff(1:72)=korea(1:72)-soyat(1:72)-tsugaru(1:72);
% es4plot=plot(diff(1:72),'y');
set(gca,'YLim',[0 4.5]);
% set(gca,'XLim',[1 12*8]);
% set(gca,'XLim',[1 72]);
set(es1plot,'LineWidth',2);
set(es2plot,'LineWidth',2);
set(es3plot,'LineWidth',2);
set(es4plot,'LineWidth',2);
title(['Model monthly mean EJS transports(',testname,')'],'fontsize',17);
xlabel('Time (month)','color','k','FontSize',17,'fontweight','bold');
ylabel('Volume Transport (sv)','color','k','FontSize',17,'fontweight','bold');
lgd=legend('Korea(Obs)','Korea(Model)','Tsugaru', 'Soya');
set(lgd,'FontSize',10);
set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
set(lgd,'Orientation','horizontal');
set(gca,'FontSize',13);
grid on;

for i=1:length(year)
    figdir =strcat('E:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',tempyear,'\figures\transport\'); % % where figure files will be saved
    if (exist(strcat(figdir) , 'dir') ~= 7)
        mkdir(strcat(figdir));
    end 
    saveas(gcf,[figdir, 'ES_transp_',num2str(year(1)),'_',num2str(year(length(year))),'.png'],'png');
end
hold off;









% % % East China Sea transport (run)

ecs4plot=plot(xData, korea_obs_run,'b');
hold on
ecs1plot=plot(xData, korea,'k');
ecs2plot=plot(xData, taiwan,'g');
ecs3plot=plot(xData, -onshore,'m');
ecs5plot=plot(xData, ECS_diff,'g');
datetick('x',12,'keeplimits')
% datetick('x',12,'keeplimits')

% diff(1:12*year)=korea(1:12*year)-soyat(1:12*year)-tsugaru(1:12*year);

% es4plot=plot(diff(1:12*year),'y');

% es1plot=plot(korea(1:72),'k');
% hold on
% es2plot=plot(tsugaru(1:72),'r');
% es3plot=plot(soyat(1:72),'b');
% diff(1:72)=korea(1:72)-soyat(1:72)-tsugaru(1:72);
% es4plot=plot(diff(1:72),'y');
set(gca,'YLim',[-1 4.5]);
% set(gca,'XLim',[1 12*8]);
% set(gca,'XLim',[1 72]);
set(ecs1plot,'LineWidth',2);
set(ecs2plot,'LineWidth',2);
set(ecs3plot,'LineWidth',2);
set(ecs4plot,'LineWidth',2);
title(['Model monthly mean ECS transports(',testname,')'],'fontsize',17);
xlabel('Time (month)','color','k','FontSize',17,'fontweight','bold');
ylabel('Volume Transport (sv)','color','k','FontSize',17,'fontweight','bold');
lgd=legend('Korea(Obs)','Korea(Model)','Taiwan', 'onshore(out(-))');
set(lgd,'FontSize',10);
set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
set(lgd,'Orientation','horizontal');
set(gca,'FontSize',13);
grid on;

for i=1:length(year)
    figdir =strcat('E:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',tempyear,'\figures\transport\'); % % where figure files will be saved
    if (exist(strcat(figdir) , 'dir') ~= 7)
        mkdir(strcat(figdir));
    end 
    saveas(gcf,[figdir, 'ECS_transp_',num2str(year(1)),'_',num2str(year(length(year))),'.png'],'png');
end
hold off;















% 
% % % % east sea transport difference
% % es4plot=plot(diff(1:60),'-o','MarkerSize',5,'Color','k');
% es4plot=plot(xData, diff(1:12*year),'Color','k');
% datetick('x','yyyy')
% hold on
% set(gca,'YLim',[-0.1 0.1]);
% % set(gca,'XLim',[1 12*year]);
% set(es4plot,'LineWidth',2);
% title('transport difference','fontsize',17);
% xlabel('Time (month)','color','k','FontSize',17,'fontweight','bold');
% set(gca,'FontSize',13);
% grid on;
% saveas(gcf,[figdir, 'trans_model_diff_',num2str(year),'.png'],'png');
% hold off;
% 
% % % % % east sea transport difference (5y)
% % % es4plot=plot(diff(1:60),'-o','MarkerSize',5,'Color','k');
% % es4plot=plot(diff(49:60),'Color','k');
% % hold on
% % set(gca,'YLim',[-0.1 0.1]);
% % set(gca,'XLim',[1 12]);
% % set(es4plot,'LineWidth',2);
% % title('transport difference (5y)','fontsize',17);
% % xlabel('Time (month)','color','k','FontSize',17,'fontweight','bold');
% % set(gca,'FontSize',15);
% % grid on;
% % saveas(gcf,[figdir, 'trans_model_diff_5y.png'],'png');
% % hold off;
% 
% 
% % % % % % east sea transport and volume difference (all)
% % % % es4plot=plot(diff(1:60),'-o','MarkerSize',5,'Color','k');
% % % for i=1:9
% % %     for j=1:12
% % %         day(12*i+j)=day(j);
% % %     end
% % % end
% % % for i=1:60
% % %     total_diff(i)=sum(diff(1:i) * 86400 .* day(1:i))
% % % end
% % % 
% % % vol_east=ncread('D:\MEPL\project\SSH\1st_year\data\test37\vol_diff_test37.nc','VOL_EAST');
% % % es4plot=plot(diff(1:60) * 86400 .* day(1:60),'Color','k');
% % % hold on
% % % es4plot2=plot(vol_east(1:60)/1000000.0,'Color','r');
% % % % set(gca,'YLim',[-0.1 0.1]);
% % % % set(gca,'XLim',[1 12]);
% % % % es4plot3=plot(total_diff(1:60),'Color','b');
% % % 
% % % set(es4plot,'LineWidth',2);
% % % set(es4plot2,'LineWidth',2);
% % % % set(es4plot3,'LineWidth',2);
% % % title('transport difference ()','fontsize',17);
% % % xlabel('Time (month)','color','k','FontSize',17,'fontweight','bold');
% % % ylabel('Volume/10^6 (m^3))','color','k','FontSize',17,'fontweight','bold');
% % % set(gca,'FontSize',15);
% % % grid on;
% % % lgd=legend('trans-diff','Volume')
% % % set(lgd,'FontSize',15);
% % % set(lgd,'Position',[0.13 0.87, 0.775, 0.03]);
% % % set(lgd,'Orientation','horizontal');
% % % set(gca,'FontSize',15);
% % % saveas(gcf,'D:\MEPL\project\SSH\1st_year\figure\test37\trans_vol_model_diff.png','png');
% % % hold off;
% % % 
% % % % es4plot3=plot(vol_east(1:120)/1000000.0,'Color','r');
% % % 
% % % % % % east sea vol(from transport) and volume difference (all)
% % % % es4plot=plot(diff(1:60),'-o','MarkerSize',5,'Color','k');
% % % day = [31 29 31 30 31 30 31 31 30 31 30 31];
% % % % day = [30 31 30 31 30 31 30 31 30 31 30 31];
% % % for i=1:9
% % %     for j=1:12
% % %         day(12*i+j)=day(j);
% % %     end
% % % end
% % % for i=1:60
% % %     total_diff(i)=sum(diff(1:i) * 86400 .* day(1:i))
% % % end
% % % 
% % % vol_east=ncread('D:\MEPL\project\SSH\1st_year\data\test37\vol_diff_test37.nc','VOL_EAST');
% % % % es4plot=plot(diff(1:60) * 86400 .* day(1:60),'Color','k');
% % % 
% % % es4plot2=plot(vol_east(1:60)/1000000.0,'Color','r');
% % % hold on
% % % % set(gca,'YLim',[-0.1 0.1]);
% % % % set(gca,'XLim',[1 12]);
% % % es4plot3=plot(total_diff(1:60),'Color','b');
% % % 
% % % % set(es4plot,'LineWidth',2);
% % % set(es4plot2,'LineWidth',2);
% % % set(es4plot3,'LineWidth',2);
% % % title('transport difference ()','fontsize',17);
% % % xlabel('Time (month)','color','k','FontSize',17,'fontweight','bold');
% % % ylabel('Volume/10^6 (m^3))','color','k','FontSize',17,'fontweight','bold');
% % % set(gca,'FontSize',15);
% % % grid on;
% % % lgd=legend('Volume','vol from trans')
% % % 
% % % set(lgd,'FontSize',15);
% % % set(lgd,'Position',[0.155 0.75 0.266 0.135 ]);
% % % set(gca,'FontSize',15);
% % % saveas(gcf,'D:\MEPL\project\SSH\1st_year\figure\test37\trans_vol_model_vol.png','png');
% % % hold off;
% % % 
% % % % es4plot3=plot(vol_east(1:120)/1000000.0,'Color','r');
% % % 
% % % 
% % % 
% % % % % % east sea vol(from transport) and volume difference (5y)
% % % % es4plot=plot(diff(1:60),'-o','MarkerSize',5,'Color','k');
% % % % day = [31 29 31 30 31 30 31 31 30 31 30 31];
% % % for i=1:9
% % %     for j=1:12
% % %         day(12*i+j)=day(j);
% % %     end
% % % end
% % % for i=1:60
% % %     total_diff(i)=sum(diff(1:i) * 86400 .* day(1:i))
% % % end
% % % 
% % % vol_east=ncread('D:\MEPL\project\SSH\1st_year\data\test37\vol_diff_test37.nc','VOL_EAST');
% % % % es4plot=plot(diff(1:60) * 86400 .* day(1:60),'Color','k');
% % % 
% % % es4plot2=plot(vol_east(49:60)/1000000.0,'Color','r');
% % % hold on
% % % % set(gca,'YLim',[-0.1 0.1]);
% % % set(gca,'XLim',[1 12]);
% % % es4plot3=plot(total_diff(49:60),'Color','b');
% % % 
% % % % set(es4plot,'LineWidth',2);
% % % set(es4plot2,'LineWidth',2);
% % % set(es4plot3,'LineWidth',2);
% % % title('transport difference (5y)','fontsize',17);
% % % xlabel('Time (month)','color','k','FontSize',17,'fontweight','bold');
% % % ylabel('Volume/10^6 (m^3))','color','k','FontSize',17,'fontweight','bold');
% % % set(gca,'FontSize',15);
% % % grid on;
% % % lgd=legend('Volume','vol from trans')
% % % set(lgd,'FontSize',15);
% % % set(lgd,'Position',[0.155 0.75 0.266 0.135 ]);
% % % set(gca,'FontSize',15);
% % % saveas(gcf,'D:\MEPL\project\SSH\1st_year\figure\test37\trans_vol_model_vol_5y.png','png');
% % % hold off;
% % % 
% % % % es4plot3=plot(vol_east(1:120)/1000000.0,'Color','r');
% 
% 
% % % % taiwan strait, onshore transport
% es1plot=plot(xData, aiwankur(1:12*year),'k');
% hold on
% es2plot=plot(xData, (-o_intruy(1:12*year)),'r');
% datetick('x','yyyy')
% % es3plot=plot(soyat(1:60),'b');
% set(gca,'YLim',[-4 4]);
% % set(gca,'XLim',[1 60]);
% set(es1plot,'LineWidth',2);
% set(es2plot,'LineWidth',2);
% % set(es3plot,'LineWidth',2);
% title('Monthly mean transports of the model results','fontsize',17);
% xlabel('Time (month)','color','k','FontSize',17,'fontweight','bold');
% ylabel('Volume Transport (sv)','color','k','FontSize',17,'fontweight','bold');
% lgd=legend('Taiwan','Onshore');
% set(lgd,'FontSize',15);
% set(lgd,'Position',[0.13 0.87, 0.775, 0.03]);
% set(lgd,'Orientation','horizontal');
% set(gca,'FontSize',13);
% grid on;
% saveas(gcf,[figdir, 'trans_taiwan_model.png'],'png');
% hold off;
% 
% 
% 
% % % % % east sea transport (cshwa)
% % load('D:\MEPL\project\SSH\1st_year\data\cshwa_transport\transp_08_8&3_dig.mat')
% % tp(1:12,:)=transp;
% % load('D:\MEPL\project\SSH\1st_year\data\cshwa_transport\transp_08_8&3_dig_4yr_1.mat')
% % tp(13:24,:)=transp;
% % es1plot=plot(tp(1:24,1),'k');
% % hold on
% % es2plot=plot(tp(1:24,2),'r');
% % es3plot=plot(tp(1:24,3),'b');
% % diffcshwa(1:24)=tp(1:24,1)-tp(1:24,3)-tp(1:24,2);
% % es4plot=plot(diff(1:60),'y');
% % set(gca,'YLim',[-0.5 4]);
% % set(gca,'XLim',[1 24]);
% % set(es1plot,'LineWidth',2);
% % set(es2plot,'LineWidth',2);
% % set(es3plot,'LineWidth',2);
% % title('Monthly mean transports of the model results','fontsize',17);
% % xlabel('Time (month)','color','k','FontSize',17,'fontweight','bold');
% % ylabel('Volume Transport (sv)','color','k','FontSize',17,'fontweight','bold');
% % lgd=legend('Korea','Tsugaru', 'Soya');
% % set(lgd,'FontSize',15);
% % set(lgd,'Position',[0.13 0.87, 0.775, 0.03]);
% % set(lgd,'Orientation','horizontal');
% % set(gca,'FontSize',15);
% % grid on;
% % saveas(gcf,[figdir, 'trans_model_cshwa.png'],'png');
% % hold off;
