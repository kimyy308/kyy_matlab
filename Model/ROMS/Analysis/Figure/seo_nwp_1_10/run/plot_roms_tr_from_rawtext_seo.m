clear all; clc; close all;

filename = 'E:\Data\Reanalysis\nwp_1_10_seo\transport\ghseo_rawfile\transport_enkf_8006_monthly.dat';
startRow = 2;

formatSpec = '%16f%16f%16f%16f%f%[^\n\r]';

fileID = fopen(filename,'r');

textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'ReturnOnError', false);

fclose(fileID);

% transportenkf8006monthly = [dataArray{1:end-1}];

korea = [dataArray{1}];
tsugaru = [dataArray{2}];
soya = [dataArray{3}];
diff=korea-tsugaru-soya;


clearvars filename startRow formatSpec fileID dataArray ans;




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



% % % east sea transport

figdir = 'D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\3rd year\figure\avg_ens_10km_mean\transport\';

nyear=length(korea)/12;


refyear=1980
startDate = datenum(['01-01-',num2str(refyear)]);
endDate = datenum(['12-31-',num2str(refyear+nyear-1)]);
xData = linspace(startDate,endDate,nyear*12);


es1plot=plot(xData, korea(1:12*nyear),'k');
hold on
es2plot=plot(xData, tsugaru(1:12*nyear),'r');
es3plot=plot(xData, soya(1:12*nyear),'b'); 
%     es4plot=plot(xData, tr_obs_30y(1:12*nyear),'y');
es5plot=plot(xData, diff(1:12*nyear),'g');
datetick('x','yy')

set(es1plot,'LineWidth',2);
set(es2plot,'LineWidth',2);
set(es3plot,'LineWidth',2);
%     set(es4plot,'LineWidth',2);
set(es5plot,'LineWidth',2);


% diff(1:12*year)=korea(1:12*year)-soya(1:12*year)-tsugaru(1:12*year);

% es4plot=plot(diff(1:12*year),'y');

% es1plot=plot(korea(1:72),'k');
% hold on
% es2plot=plot(tsugaru(1:72),'r');
% es3plot=plot(soyat(1:72),'b');
% diff(1:72)=korea(1:72)-soyat(1:72)-tsugaru(1:72);
% es4plot=plot(diff(1:72),'y');
set(gca,'YLim',[-4.5 4.5]);
% set(gca,'XLim',[1 12*8]);
% set(gca,'XLim',[1 72]);

title(['Monthly mean transports of the reanalysis results-textfile '],'fontsize',17);
xlabel('Time (month)','color','k','FontSize',17,'fontweight','bold');
ylabel('Volume Transport (sv)','color','k','FontSize',17,'fontweight','bold');
lgd=legend('Korea(Model)','Tsugaru', 'Soya','ks-ts-ss');
set(lgd,'FontSize',10);
set(lgd,'Position',[0.13 0.18, 0.775, 0.03]);
set(lgd,'Orientation','horizontal');
set(gca,'FontSize',13);
grid on;
hold off;
saveas(gcf,[figdir, 'trans_seo_',num2str(nyear),'y_','text','.png'],'png');










% % % east sea transport

figdir = 'D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\3rd year\figure\avg_ens_10km_mean\transport\';

refyear=1985
nyear=4;
startDate = datenum(['01-01-',num2str(refyear)]);
endDate = datenum(['12-31-',num2str(refyear+nyear-1)]);
xData = linspace(startDate,endDate,nyear*12);

es1plot=plot(xData, korea(49:48+12*nyear),'k');
hold on
es2plot=plot(xData, tsugaru(49:48+12*nyear),'r');
es3plot=plot(xData, soya(49:48+12*nyear),'b'); 
%     es4plot=plot(xData, tr_obs_30y(1:12*nyear),'y');
es5plot=plot(xData, diff(49:48+12*nyear),'g');
datetick('x','yy')

set(es1plot,'LineWidth',2);
set(es2plot,'LineWidth',2);
set(es3plot,'LineWidth',2);
%     set(es4plot,'LineWidth',2);
set(es5plot,'LineWidth',2);


% diff(1:12*year)=korea(1:12*year)-soya(1:12*year)-tsugaru(1:12*year);

% es4plot=plot(diff(1:12*year),'y');

% es1plot=plot(korea(1:72),'k');
% hold on
% es2plot=plot(tsugaru(1:72),'r');
% es3plot=plot(soyat(1:72),'b');
% diff(1:72)=korea(1:72)-soyat(1:72)-tsugaru(1:72);
% es4plot=plot(diff(1:72),'y');
set(gca,'YLim',[-4.5 4.5]);
% set(gca,'XLim',[1 12*8]);
% set(gca,'XLim',[1 72]);

title(['Monthly mean transports of the reanalysis results-textfile '],'fontsize',17);
xlabel('Time (month)','color','k','FontSize',17,'fontweight','bold');
ylabel('Volume Transport (sv)','color','k','FontSize',17,'fontweight','bold');
lgd=legend('Korea(Model)','Tsugaru', 'Soya','ks-ts-ss');
set(lgd,'FontSize',10);
set(lgd,'Position',[0.13 0.18, 0.775, 0.03]);
set(lgd,'Orientation','horizontal');
set(gca,'FontSize',13);
grid on;
hold off;
saveas(gcf,[figdir, 'trans_seo_85_88','y_','text','.png'],'png');




for i=1:4
    korea_8588(i)=mean(korea(48+12*(i-1):50+12*(i-1)));
    tsugaru_8588(i)=mean(tsugaru(48+12*(i-1):50+12*(i-1)));
    soya_8588(i)=mean(soya(48+12*(i-1):50+12*(i-1)));
    diff_8588(i)=mean(diff(48+12*(i-1):50+12*(i-1)));
end

startDate = datenum(['01-01-',num2str(refyear)]);
endDate = datenum(['12-31-',num2str(refyear+nyear-1)]);
xData = linspace(startDate,endDate,nyear);

es1plot=plot(xData, korea_8588,'k');
hold on
es2plot=plot(xData, tsugaru_8588,'r');
es3plot=plot(xData, soya_8588,'b'); 
%     es4plot=plot(xData, tr_obs_30y(1:12*nyear),'y');
es5plot=plot(xData, diff_8588,'g');
datetick('x','yy')


yeardiff=mean(korea_8588(3:4))-mean(korea_8588(1:2));


set(es1plot,'LineWidth',2);
set(es2plot,'LineWidth',2);
set(es3plot,'LineWidth',2);
%     set(es4plot,'LineWidth',2);
set(es5plot,'LineWidth',2);


% diff(1:12*year)=korea(1:12*year)-soya(1:12*year)-tsugaru(1:12*year);

% es4plot=plot(diff(1:12*year),'y');

% es1plot=plot(korea(1:72),'k');
% hold on
% es2plot=plot(tsugaru(1:72),'r');
% es3plot=plot(soyat(1:72),'b');
% diff(1:72)=korea(1:72)-soyat(1:72)-tsugaru(1:72);
% es4plot=plot(diff(1:72),'y');
set(gca,'YLim',[0 4.0]);
% set(gca,'XLim',[1 12*8]);
% set(gca,'XLim',[1 72]);

title(['Monthly mean transports of the reanalysis results-textfile '],'fontsize',17);
xlabel('Time (month)','color','k','FontSize',17,'fontweight','bold');
ylabel('Volume Transport (sv)','color','k','FontSize',17,'fontweight','bold');
lgd=legend('Korea(Model)','Tsugaru', 'Soya','ks-ts-ss');
set(lgd,'FontSize',10);
set(lgd,'Position',[0.13 0.78, 0.775, 0.03]);
set(lgd,'Orientation','horizontal');
set(gca,'FontSize',13);
grid on;
hold off;
saveas(gcf,[figdir, 'trans_seo_85_88_winter','y_','text','.png'],'png');