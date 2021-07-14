clear all; clc; close all;

filename = 'D:\MEPL\project\SSH\2nd_year\data\test37\transport\roms_tr_test37_2018_01_22.txt';
startRow = 2;

formatSpec = '%8f%9f%9f%9f%9f%9f%s%[^\n\r]';

fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

dataArray{7} = strtrim(dataArray{7});

fclose(fileID);

korea = dataArray{:, 1};
tsugaru = dataArray{:, 2};
soyat = dataArray{:, 3};
aiwankur = dataArray{:, 4};
o_intruy = dataArray{:, 5};
ellowsea = dataArray{:, 6};
% VarName7 = dataArray{:, 7};

clearvars filename startRow formatSpec fileID dataArray ans;

% % 
% %  get observation data (Jan, 2001 - Sep, 2012) 
% % 

filename = 'D:\MEPL\project\SSH\2nd_year\data\test37\transport\obs_tsushima_original\2001jan_2012sep.txt';
endRow = 141;
formatSpec = '%5s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, endRow, 'Delimiter', '', 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));
rawData = dataArray{1};
for row=1:size(rawData, 1);
    % 숫자형이 아닌 접두사 및 접미사를 검색하고 제거하는 정규 표현식을 만듭니다.
    regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
    try
        result = regexp(rawData{row}, regexstr, 'names');
        numbers = result.numbers;
        
        % 천 단위가 아닌 위치에서 쉼표를 검색했습니다.
        invalidThousandsSeparator = false;
        if any(numbers==',');
            thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
            if isempty(regexp(numbers, thousandsRegExp, 'once'));
                numbers = NaN;
                invalidThousandsSeparator = true;
            end
        end
        % 숫자형 텍스트를 숫자로 변환합니다.
        if ~invalidThousandsSeparator;
            numbers = textscan(strrep(numbers, ',', ''), '%f');
            numericData(row, 1) = numbers{1};
            raw{row, 1} = numbers{1};
        end
    catch me
    end
end
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % 숫자형이 아닌 셀 찾기
raw(R) = {NaN}; % 숫자형이 아닌 셀 바꾸기
Obs2001_2012 = cell2mat(raw(:, 1));
clearvars filename endRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me R;

figdir = 'D:\MEPL\project\SSH\2nd_year\figure\test37\Transport\';



% filename = 'D:\MEPL\project\SSH\1st_year\data\test37\roms_tr_test37_2017_12_12_daily.txt';
% startRow = 2;
% 
% formatSpec = '%8f%9f%9f%9f%9f%9f%s%[^\n\r]';
% 
% fileID = fopen(filename,'r');
% 
% dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
% 
% dataArray{7} = strtrim(dataArray{7});
% 
% fclose(fileID);
% 
% daily_korea = dataArray{:, 1};
% daily_tsugaru = dataArray{:, 2};
% daily_soyat = dataArray{:, 3};
% daily_aiwankur = dataArray{:, 4};
% daily_o_intruy = dataArray{:, 5};
% daily_ellowsea = dataArray{:, 6};
% VarName7 = dataArray{:, 7};
% 
% clearvars filename startRow formatSpec fileID dataArray ans;
% 
% day = [31 29 31 30 31 30 31 31 30 31 30 31];
% for i= 1: 12
%     if (i==1)
%         m_daily_korea(i)=mean(daily_korea(1:day(1)));
%     else
%         m_daily_korea(i)=mean(daily_korea(sum(day(1:i-1))+1:sum(day(1:i))))
%     end
% end
% % 1:31
% % 31+1:32+29-1
% % 31+29+1:31+29+31-1


refyear=2001;
year=10;


plot(korea)
hold on
plot(tsugaru)
plot(soyat)
% plot(aiwankur)
% plot(o_intruy)
% plot(ellowsea)
hold off

% % 1970, Jan ~ 2012, sep  korea strait transport observation climatology data
clim =[1.853 2.057 2.383 2.627 2.726 2.693 2.821 3.039 3.002 3.047 2.761 2.257];
% % Fukudome et al, 2010 climatology mean feb 1997 to feb 2007 
fuku_clim =[2.01 2.28 2.52 2.60 2.66 2.70 2.68 3.05 2.93 3.10 2.84 2.36];
std = [0.274404 0.274396 0.224372 0.218135 0.195140 0.211314 0.219307 0.288883 0.262485 0.274527 0.297884 0.31302];
n = [43 43 43 42 43 43 43 43 43 42 42 42];
err= 1.96* std./sqrt(n);  %%(95% confidence interval)
krplot=plot(mean(reshape(korea',[year,12])),'k');
set(krplot,'LineWidth',4);
set(gca,'YLim',[0 4]);
set(gca,'XLim',[1 12]);
hold on
% obplot=plot(obs);
% obplot=plot(clim,'r');

obplot=errorbar(1:12,clim,std,'r');
set(obplot,'LineWidth',4);

% fukuplot=plot(fuku_clim,'b');
% set(fukuplot,'LineWidth',4);

% diff=abs(korea(49:60)-clim');
% diffplot=plot(diff,'y')


% set(diffplot,'LineWidth',4);
hold off
title('Monthly mean transports through the Korea Strait','fontsize',13);
xlabel('Time (month)','color','k','FontSize',17,'fontweight','bold');
ylabel('Volume Transport (sv)','color','k','FontSize',17,'fontweight','bold');
% lgd=legend('Model','Climatology','Difference')
lgd=legend('Model(5th Year)','Observation')
% lgd=legend('Model','Fukudome et al(2010)')
% lgd=legend('Model','Climatology','Fukudome et al(2010)')
set(lgd,'FontSize',15);
set(lgd,'Position',[0.555 0.25 0.266 0.135 ]);
set(gca,'FontSize',15);
saveas(gcf,[figdir, 'trans_clim.png'],'png');


% % % 2011, oct ~ 2012, sep  korea strait transport observation data
% obs=[2.004 2.264 2.569 2.431 2.528 2.699 2.858 3.178 2.866 3.436 3.286 2.257 ];

% % 2001, jan ~ dec  korea strait transport observation data
obs=[2.059 2.362 2.567 2.733 2.632 2.634 2.858 3.031 2.941 3.002 2.730 2.047];

% krplot=plot(korea(49:60),'k'); %%5y
% krplot=plot(korea(61:72),'k'); %%6y
krplot=plot(korea(1:12),'k'); %%7y
set(krplot,'LineWidth',4);
set(gca,'YLim',[0 4]);
set(gca,'XLim',[1 12]);
hold on
obplot=plot(obs,'b');
% obplot=plot(clim);
% diff=abs(korea(49:60)-obs');
% diffplot=plot(diff,'y')
set(obplot,'LineWidth',4);
% set(diffplot,'LineWidth',4);
hold off
title('Monthly mean transports through the Korea Strait (7y)','fontsize',17);
xlabel('Time (month)','color','k','FontSize',17,'fontweight','bold');
ylabel('Volume Transport (sv)','color','k','FontSize',17,'fontweight','bold');
% lgd=legend('Model','Observation','Difference')
lgd=legend('Model','Observation')
set(lgd,'FontSize',15);
set(gca,'FontSize',15);
saveas(gcf,[figdir, 'trans_obs_7y.png'],'png');


% % % east sea transport



startDate = datenum(['01-01-',num2str(refyear)]);
endDate = datenum(['12-31-',num2str(refyear+year-1)]);
xData = linspace(startDate,endDate,year*12);


es4plot=plot(xData, Obs2001_2012(1:12*year),'y');
hold on
es1plot=plot(xData, korea(1:12*year),'k');
es2plot=plot(xData, tsugaru(1:12*year),'r');
es3plot=plot(xData, soyat(1:12*year),'b');
datetick('x','yyyy')

diff(1:12*year)=korea(1:12*year)-soyat(1:12*year)-tsugaru(1:12*year);

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
set(es1plot,'LineWidth',4);
set(es2plot,'LineWidth',4);
set(es3plot,'LineWidth',4);
set(es4plot,'LineWidth',4);
title('Monthly mean transports of the model results','fontsize',17);
xlabel('Time (month)','color','k','FontSize',17,'fontweight','bold');
ylabel('Volume Transport (sv)','color','k','FontSize',17,'fontweight','bold');
lgd=legend('Korea(Obs)','Korea(Model)','Tsugaru', 'Soya');
set(lgd,'FontSize',10);
set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
set(lgd,'Orientation','horizontal');
set(gca,'FontSize',13);
grid on;
saveas(gcf,[figdir, 'trans_model_',num2str(year),'y.png'],'png');
hold off;





% % % east sea transport difference
% es4plot=plot(diff(1:60),'-o','MarkerSize',5,'Color','k');
es4plot=plot(xData, diff(1:12*year),'Color','k');
datetick('x','yyyy')
hold on
set(gca,'YLim',[-0.1 0.1]);
% set(gca,'XLim',[1 12*year]);
set(es4plot,'LineWidth',4);
title('transport difference','fontsize',17);
xlabel('Time (month)','color','k','FontSize',17,'fontweight','bold');
set(gca,'FontSize',13);
grid on;
saveas(gcf,[figdir, 'trans_model_diff_',num2str(year),'.png'],'png');
hold off;

% % % % east sea transport difference (5y)
% % es4plot=plot(diff(1:60),'-o','MarkerSize',5,'Color','k');
% es4plot=plot(diff(49:60),'Color','k');
% hold on
% set(gca,'YLim',[-0.1 0.1]);
% set(gca,'XLim',[1 12]);
% set(es4plot,'LineWidth',4);
% title('transport difference (5y)','fontsize',17);
% xlabel('Time (month)','color','k','FontSize',17,'fontweight','bold');
% set(gca,'FontSize',15);
% grid on;
% saveas(gcf,[figdir, 'trans_model_diff_5y.png'],'png');
% hold off;


% % % % % east sea transport and volume difference (all)
% % % es4plot=plot(diff(1:60),'-o','MarkerSize',5,'Color','k');
% % for i=1:9
% %     for j=1:12
% %         day(12*i+j)=day(j);
% %     end
% % end
% % for i=1:60
% %     total_diff(i)=sum(diff(1:i) * 86400 .* day(1:i))
% % end
% % 
% % vol_east=ncread('D:\MEPL\project\SSH\1st_year\data\test37\vol_diff_test37.nc','VOL_EAST');
% % es4plot=plot(diff(1:60) * 86400 .* day(1:60),'Color','k');
% % hold on
% % es4plot2=plot(vol_east(1:60)/1000000.0,'Color','r');
% % % set(gca,'YLim',[-0.1 0.1]);
% % % set(gca,'XLim',[1 12]);
% % % es4plot3=plot(total_diff(1:60),'Color','b');
% % 
% % set(es4plot,'LineWidth',4);
% % set(es4plot2,'LineWidth',4);
% % % set(es4plot3,'LineWidth',4);
% % title('transport difference ()','fontsize',17);
% % xlabel('Time (month)','color','k','FontSize',17,'fontweight','bold');
% % ylabel('Volume/10^6 (m^3))','color','k','FontSize',17,'fontweight','bold');
% % set(gca,'FontSize',15);
% % grid on;
% % lgd=legend('trans-diff','Volume')
% % set(lgd,'FontSize',15);
% % set(lgd,'Position',[0.13 0.87, 0.775, 0.03]);
% % set(lgd,'Orientation','horizontal');
% % set(gca,'FontSize',15);
% % saveas(gcf,'D:\MEPL\project\SSH\1st_year\figure\test37\trans_vol_model_diff.png','png');
% % hold off;
% % 
% % % es4plot3=plot(vol_east(1:120)/1000000.0,'Color','r');
% % 
% % % % % east sea vol(from transport) and volume difference (all)
% % % es4plot=plot(diff(1:60),'-o','MarkerSize',5,'Color','k');
% % day = [31 29 31 30 31 30 31 31 30 31 30 31];
% % % day = [30 31 30 31 30 31 30 31 30 31 30 31];
% % for i=1:9
% %     for j=1:12
% %         day(12*i+j)=day(j);
% %     end
% % end
% % for i=1:60
% %     total_diff(i)=sum(diff(1:i) * 86400 .* day(1:i))
% % end
% % 
% % vol_east=ncread('D:\MEPL\project\SSH\1st_year\data\test37\vol_diff_test37.nc','VOL_EAST');
% % % es4plot=plot(diff(1:60) * 86400 .* day(1:60),'Color','k');
% % 
% % es4plot2=plot(vol_east(1:60)/1000000.0,'Color','r');
% % hold on
% % % set(gca,'YLim',[-0.1 0.1]);
% % % set(gca,'XLim',[1 12]);
% % es4plot3=plot(total_diff(1:60),'Color','b');
% % 
% % % set(es4plot,'LineWidth',4);
% % set(es4plot2,'LineWidth',4);
% % set(es4plot3,'LineWidth',4);
% % title('transport difference ()','fontsize',17);
% % xlabel('Time (month)','color','k','FontSize',17,'fontweight','bold');
% % ylabel('Volume/10^6 (m^3))','color','k','FontSize',17,'fontweight','bold');
% % set(gca,'FontSize',15);
% % grid on;
% % lgd=legend('Volume','vol from trans')
% % 
% % set(lgd,'FontSize',15);
% % set(lgd,'Position',[0.155 0.75 0.266 0.135 ]);
% % set(gca,'FontSize',15);
% % saveas(gcf,'D:\MEPL\project\SSH\1st_year\figure\test37\trans_vol_model_vol.png','png');
% % hold off;
% % 
% % % es4plot3=plot(vol_east(1:120)/1000000.0,'Color','r');
% % 
% % 
% % 
% % % % % east sea vol(from transport) and volume difference (5y)
% % % es4plot=plot(diff(1:60),'-o','MarkerSize',5,'Color','k');
% % % day = [31 29 31 30 31 30 31 31 30 31 30 31];
% % for i=1:9
% %     for j=1:12
% %         day(12*i+j)=day(j);
% %     end
% % end
% % for i=1:60
% %     total_diff(i)=sum(diff(1:i) * 86400 .* day(1:i))
% % end
% % 
% % vol_east=ncread('D:\MEPL\project\SSH\1st_year\data\test37\vol_diff_test37.nc','VOL_EAST');
% % % es4plot=plot(diff(1:60) * 86400 .* day(1:60),'Color','k');
% % 
% % es4plot2=plot(vol_east(49:60)/1000000.0,'Color','r');
% % hold on
% % % set(gca,'YLim',[-0.1 0.1]);
% % set(gca,'XLim',[1 12]);
% % es4plot3=plot(total_diff(49:60),'Color','b');
% % 
% % % set(es4plot,'LineWidth',4);
% % set(es4plot2,'LineWidth',4);
% % set(es4plot3,'LineWidth',4);
% % title('transport difference (5y)','fontsize',17);
% % xlabel('Time (month)','color','k','FontSize',17,'fontweight','bold');
% % ylabel('Volume/10^6 (m^3))','color','k','FontSize',17,'fontweight','bold');
% % set(gca,'FontSize',15);
% % grid on;
% % lgd=legend('Volume','vol from trans')
% % set(lgd,'FontSize',15);
% % set(lgd,'Position',[0.155 0.75 0.266 0.135 ]);
% % set(gca,'FontSize',15);
% % saveas(gcf,'D:\MEPL\project\SSH\1st_year\figure\test37\trans_vol_model_vol_5y.png','png');
% % hold off;
% % 
% % % es4plot3=plot(vol_east(1:120)/1000000.0,'Color','r');


% % % taiwan strait, onshore transport
es1plot=plot(xData, aiwankur(1:12*year),'k');
hold on
es2plot=plot(xData, (-o_intruy(1:12*year)),'r');
datetick('x','yyyy')
% es3plot=plot(soyat(1:60),'b');
set(gca,'YLim',[-4 4]);
% set(gca,'XLim',[1 60]);
set(es1plot,'LineWidth',4);
set(es2plot,'LineWidth',4);
% set(es3plot,'LineWidth',4);
title('Monthly mean transports of the model results','fontsize',17);
xlabel('Time (month)','color','k','FontSize',17,'fontweight','bold');
ylabel('Volume Transport (sv)','color','k','FontSize',17,'fontweight','bold');
lgd=legend('Taiwan','Onshore');
set(lgd,'FontSize',15);
set(lgd,'Position',[0.13 0.87, 0.775, 0.03]);
set(lgd,'Orientation','horizontal');
set(gca,'FontSize',13);
grid on;
saveas(gcf,[figdir, 'trans_taiwan_model.png'],'png');
hold off;



% % % % east sea transport (cshwa)
% load('D:\MEPL\project\SSH\1st_year\data\cshwa_transport\transp_08_8&3_dig.mat')
% tp(1:12,:)=transp;
% load('D:\MEPL\project\SSH\1st_year\data\cshwa_transport\transp_08_8&3_dig_4yr_1.mat')
% tp(13:24,:)=transp;
% es1plot=plot(tp(1:24,1),'k');
% hold on
% es2plot=plot(tp(1:24,2),'r');
% es3plot=plot(tp(1:24,3),'b');
% diffcshwa(1:24)=tp(1:24,1)-tp(1:24,3)-tp(1:24,2);
% es4plot=plot(diff(1:60),'y');
% set(gca,'YLim',[-0.5 4]);
% set(gca,'XLim',[1 24]);
% set(es1plot,'LineWidth',4);
% set(es2plot,'LineWidth',4);
% set(es3plot,'LineWidth',4);
% title('Monthly mean transports of the model results','fontsize',17);
% xlabel('Time (month)','color','k','FontSize',17,'fontweight','bold');
% ylabel('Volume Transport (sv)','color','k','FontSize',17,'fontweight','bold');
% lgd=legend('Korea','Tsugaru', 'Soya');
% set(lgd,'FontSize',15);
% set(lgd,'Position',[0.13 0.87, 0.775, 0.03]);
% set(lgd,'Orientation','horizontal');
% set(gca,'FontSize',15);
% grid on;
% saveas(gcf,[figdir, 'trans_model_cshwa.png'],'png');
% hold off;
