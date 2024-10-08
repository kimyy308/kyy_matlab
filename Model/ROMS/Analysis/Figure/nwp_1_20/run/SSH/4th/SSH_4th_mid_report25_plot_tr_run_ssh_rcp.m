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

testname = 'test57';
refyear = 2006;
year = [2006:2090];
month = [1:12]; % % put month which you want to plot [month month ...]
% inputdir = strcat(['G:\Data\Model\ROMS\nwp_1_20\test37\input\']);
% outputdir = strcat(['G:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',num2str(year,'%04i'),'\']); % % where data files are

for i=1:length(year)
    tempyear=num2str(year(i),'%04i');
    filename = ['G:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',tempyear,'\nwp_1_20_monthly_',testname,'_',tempyear,'.txt'];
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
ECS_diff=korea-taiwan+onshore;
AKP_diff=-(soya+tsugaru-taiwan+onshore);
AKP_diff_yr=squeeze(mean(reshape(AKP_diff(1:348),[12,29]),1));
% % get observation data (from .xls)

% taiwan_res=reshape(taiwan,[29, 12]);
% onshore_res=reshape(onshore,[29, 12]);

% [~, ~, raw] = xlsread('G:\Data\Observation\Transport_Korea\obs_tsushima_197001_201209_ORIGIN.xls','TWC','A2:H553');
% raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
% 
% R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % 숫자형이 아닌 셀 찾기
% raw(R) = {NaN}; % 숫자형이 아닌 셀 바꾸기
% 
% obs_ks = reshape([raw{:}],size(raw));
% 
% clearvars raw R;

% year_obs_30y=obs_ks(find(obs_ks(:,1)<2010 & obs_ks(:,1)>1979),1);
% month_obs_30y=obs_ks(find(obs_ks(:,1)<2010 & obs_ks(:,1)>1979),2);
% tr_obs_30y=obs_ks(find(obs_ks(:,1)<2010 & obs_ks(:,1)>1979),5);

% korea_obs=tr_obs_30y(find(year_obs_30y >= year(1) & year_obs_30y <= year(length(year))))';    %% get observation transport at the korea strait

% korea_obs_run=korea_obs;
% for i=1:length(year)
%     korea_obs_run((12*(i-1))+1:(12*(i-1))+12)=tr_obs_30y(find(year_obs_30y == refyear))';    %% get observation transport at the korea strait
% end


% % Takikawa et al, 2005 sea level difference Jan 1965 to Dec 2001
% % Takikawa et al, 2005 ADCP mean feb 1997 to Aug 2002 
% % Fukudome et al, 2010 ADCP Mar 1997 to feb 2007 

% SLD(206:360)=NaN;
% ADCP=obs_ks(find(obs_ks(:,1)<=year(end) & obs_ks(:,1)>=year(1)),5);
% SLD=obs_ks(find(obs_ks(:,1)<=year(end) & obs_ks(:,1)>=year(1)),6);
% HYCOM=obs_ks(find(obs_ks(:,1)<=year(end) & obs_ks(:,1)>=year(1)),7);
% Fuk=obs_ks(find(obs_ks(:,1)<=year(end) & obs_ks(:,1)>=year(1)),8);
% Fuk(1:205)=NaN; 
% Fuk(327:360)=NaN;




clearvars filename startRow formatSpec fileID dataArray ans;




% % % East Sea transport (run)

% es4plot=plot(xData, korea_obs_run,'b');
hold on
% es1plot=plot(xData(1:36), korea(1:36),'k');
% es2plot=plot(xData(1:36), tsugaru(1:36),'r');
% es3plot=plot(xData(1:36), soya(1:36),'b');
% es5plot=plot(xData(1:36), ES_diff(1:36),'g');
es1plot=plot(xData, korea,'k');
es2plot=plot(xData, tsugaru,'r');
es3plot=plot(xData, soya,'b');
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
% set(es4plot,'LineWidth',2);
title(['Model monthly mean EJS transports(',testname,')'],'fontsize',17);
xlabel('Time (month)','color','k','FontSize',17,'fontweight','bold');
ylabel('Volume Transport (sv)','color','k','FontSize',17,'fontweight','bold');
lgd=legend('Korea(Model)','Tsugaru', 'Soya');
set(lgd,'FontSize',10);
% set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
set(gcf,'PaperPosition', [0 0 36 12]) 
set(lgd,'Orientation','horizontal');
set(gca,'FontSize',13);
grid on;
grid minor;

for i=1:length(year)
    figdir =strcat('G:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',tempyear,'\figures\transport\'); % % where figure files will be saved
    if (exist(strcat(figdir) , 'dir') ~= 7)
        mkdir(strcat(figdir));
    end 
    saveas(gcf,[figdir, 'ES_transp_',num2str(year(1)),'_',num2str(year(3)),'.png'],'png');
end
hold off;
close all;


% % % % Korea Strait transport (run)
% % korea2=korea(1:327);
% % SLD2=SLD(1:327);
% % Fuk2=Fuk(1:327);
% % hycomtrans4=hycomtrans3(1:327);
% % xData2=xData(1:327);
% 
% % es1plot=plot(xData, korea, 'color', 'b');
% % hold on
% % es2plot=plot(xData, SLD,'r');
% % es3plot=plot(xData, Fuk,'color',[0.91 0.41 0.17]);
% % es4plot=plot(xData, hycomtrans3,'g');
% % es1plot=plot(xData2, korea2, 'color', 'b');
% es1plot=plot(xData, korea, 'color', [0/255, 114/255, 189/255]);
% 
% hold on
% % es2plot=plot(xData2, SLD2,'r');
% es2plot=plot(xData, SLD,'color',[217/255, 83/255, 25/255]);
% 
% % es3plot=plot(xData2, Fuk2,'color',[0.91 0.41 0.17]);
% es3plot=plot(xData, Fuk,'color',[237/255 117/255 32/255]);
% 
% es4plot=plot(xData, HYCOM,'color',[0.8 0.8 0.8]);
% 
% % datetick('x',12,'keeplimits')
% datetick('x','yy','keeplimits')
% 
% axis tight;
% set(gca,'YLim',[0 5]);
% set(es1plot,'LineWidth',2);
% set(es2plot,'LineWidth',1.5);
% set(es3plot,'LineWidth',1.5);
% set(es4plot,'LineWidth',1.5);
% title(['Model monthly mean EJS transports(',testname,')'],'fontsize',17);
% xlabel('Time (year)','color','k','FontSize',17,'fontweight','bold');
% ylabel('Volume Transport (sv)','color','k','FontSize',17,'fontweight','bold');
% lgd=legend('Korea(Model)','Sea lev dif', 'ADCP(97.2-07.2)', 'hycom');
% set(lgd,'FontSize',13);
% % set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
% %         set(gcf,'PaperPosition', [0 0 36 18])   %% normal
%                 set(gcf,'PaperPosition', [0 0 54 18]) 
% 
% set(lgd,'Orientation','horizontal');
% set(gca,'FontSize',15);
% grid on;
% grid minor;
% 
% % for i=1:length(year)
%     figdir =strcat('G:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',tempyear,'\figures\transport\'); % % where figure files will be saved
%     if (exist(strcat(figdir) , 'dir') ~= 7)
%         mkdir(strcat(figdir));
%     end 
%     saveas(gcf,[figdir, 'KS_transp_',num2str(year(1)),'_',num2str(year(length(year))),'.png'],'png');
% % end
% hold off;




% % % East China Sea transport (run)

% ecs4plot=plot(xData(1:36), korea_obs_run(1:36),'b');
hold on
% ecs1plot=plot(xData(1:36), korea(1:36),'k');
% ecs2plot=plot(xData(1:36), taiwan(1:36),'g');
% ecs3plot=plot(xData(1:36), -onshore(1:36),'m');
% ecs5plot=plot(xData(1:36), ECS_diff(1:36),'r');
ecs1plot=plot(xData, korea,'k');
ecs2plot=plot(xData, taiwan,'g');
ecs3plot=plot(xData, -onshore,'m');
ecs5plot=plot(xData, ECS_diff,'r');
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
% set(ecs4plot,'LineWidth',2);
title(['Model monthly mean ECS transports(',testname,')'],'fontsize',17);
xlabel('Time (month)','color','k','FontSize',17,'fontweight','bold');
ylabel('Volume Transport (sv)','color','k','FontSize',17,'fontweight','bold');
lgd=legend('Korea(Model)','Taiwan', 'onshore(in)', 'diff');
set(lgd,'FontSize',10);
set(gcf,'PaperPosition', [0 0 36 12]) 
% set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
set(lgd,'Orientation','horizontal');
set(gca,'FontSize',13);
grid on;
grid minor;
for i=1:length(year)
    figdir =strcat('G:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',tempyear,'\figures\transport\'); % % where figure files will be saved
    if (exist(strcat(figdir) , 'dir') ~= 7)
        mkdir(strcat(figdir));
    end 
    saveas(gcf,[figdir, 'ECS_transp_',num2str(year(1)),'_',num2str(year(3)),'.png'],'png');
end
hold off;
close all;


% % % Around Korea Peninsula transport (run)
% akp4plot=plot(xData(1:36), tsugaru(1:36),'b');
% hold on
% akp1plot=plot(xData(1:36), soya(1:36),'k');
% akp2plot=plot(xData(1:36), taiwan(1:36),'g');
% akp3plot=plot(xData(1:36), -onshore(1:36),'m');
% akp5plot=plot(xData(1:36), AKP_diff(1:36),'r');
akp4plot=plot(xData, tsugaru,'b');
hold on
akp1plot=plot(xData, soya,'k');
akp2plot=plot(xData, taiwan,'g');
akp3plot=plot(xData, -onshore,'m');
akp5plot=plot(xData, AKP_diff,'r');
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
set(akp1plot,'LineWidth',2);
set(akp2plot,'LineWidth',2);
set(akp3plot,'LineWidth',2);
set(akp4plot,'LineWidth',2);
title(['Model monthly mean AKP transports(',testname,')'],'fontsize',17);
xlabel('Time (month)','color','k','FontSize',17,'fontweight','bold');
ylabel('Volume Transport (sv)','color','k','FontSize',17,'fontweight','bold');
lgd=legend('tsugaru', 'soya', 'taiwan', 'onshore(in)', 'diff');
set(lgd,'FontSize',10);
set(gcf,'PaperPosition', [0 0 36 12]) 
% set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
set(lgd,'Orientation','horizontal');
set(gca,'FontSize',13);
grid on;
grid minor;
for i=1:length(year)
    figdir =strcat('G:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',tempyear,'\figures\transport\'); % % where figure files will be saved
    if (exist(strcat(figdir) , 'dir') ~= 7)
        mkdir(strcat(figdir));
    end 
    saveas(gcf,[figdir, 'AKP_transp_',num2str(year(1)),'_',num2str(year(3)),'.png'],'png');
end
hold off;
close all;