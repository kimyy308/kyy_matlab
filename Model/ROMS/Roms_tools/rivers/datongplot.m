close all; clear all;  clc;
warning off;

close all;
clearvars '*' -except regionind2 all_region2
% % % 

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


%% �ؽ�Ʈ ���Ͽ��� �����͸� �����ɴϴ�.
% ���� �ؽ�Ʈ ���Ͽ��� �����͸� �������� ���� ��ũ��Ʈ:
%
%    C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Roms_tools\rivers\datong_1961_2013.dat
%
% ������ �ٸ� �����ͳ� �ؽ�Ʈ ���Ϸ� �ڵ带 Ȯ���Ϸ��� ��ũ��Ʈ ��� �Լ��� �����Ͻʽÿ�.

% MATLAB���� ���� ��¥�� �ڵ� ������: 2018/11/19 10:12:43

%% ������ �ʱ�ȭ�մϴ�.
filename = 'C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Roms_tools\rivers\datong_1961_2013.dat';
startRow = 230;
endRow = 577;

%% �� �ؽ�Ʈ ������ ����:
%   ��1: double (%f)
%	��2: double (%f)
% �ڼ��� ������ ���� �������� TEXTSCAN�� �����Ͻʽÿ�.
formatSpec = '%7f%11f%[^\n\r]';

%% �ؽ�Ʈ ������ ���ϴ�.
fileID = fopen(filename,'r');

%% ���Ŀ� ���� ������ ���� �н��ϴ�.
% �� ȣ���� �� �ڵ带 �����ϴ� �� ���Ǵ� ������ ����ü�� ������� �մϴ�. �ٸ� ���Ͽ� ���� ������ �߻��ϴ� ��� �������� ������
% �ڵ带 �ٽ� �����Ͻʽÿ�.
dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% �ؽ�Ʈ ������ �ݽ��ϴ�.
fclose(fileID);

%% ������ �� ���� �����Ϳ� ���� ���� ó�� ���Դϴ�.
% �������� �������� ������ �� ���� �����Ϳ� ��Ģ�� ������� �ʾ����Ƿ� ���� ó�� �ڵ尡 ���Ե��� �ʾҽ��ϴ�. ������ �� ����
% �����Ϳ� ����� �ڵ带 �����Ϸ��� ���Ͽ��� ������ �� ���� ���� �����ϰ� ��ũ��Ʈ�� �ٽ� �����Ͻʽÿ�.

%% ��� ���� �����
datong1 = table(dataArray{1:end-1}, 'VariableNames', {'YYMM','CRD'});

%% �ӽ� ���� �����
clearvars filename startRow endRow formatSpec fileID dataArray ans;

testname='test42'
if (strcmp(system_name,'PCWIN64'))
    % % for windows
    figrawdir =strcat('D:\OneDrive - ������б�\MEPL\project\SSH\3rd_year\figure\nwp_1_20\',testname,'\NIFS\'); % % where figure files will be saved
    param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m'
    filedir = strcat('E:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
    NIFSdir='E:\Data\Observation\OISST\monthly\';
    run(param_script);
elseif (strcmp(system_name,'GLNXA64'))
end

figdir=[figrawdir,'Trend\'];

if (exist(strcat(figdir) , 'dir') ~= 7)
    mkdir(strcat(figdir));
end 
outfile = strcat(figdir);



inputyear=1980:2008;
discharge=table2array(datong1(:,2));

reshap_discharge=reshape(discharge,[12,29]);

seasonal_mean=mean(reshap_discharge,2);

for i=1:length(inputyear) 
    reshap_discharge2(:,i)=reshap_discharge(:,i)-seasonal_mean;
end
discharge2=reshape(reshap_discharge2,[1 348])';

for i =1:length(inputyear) 
    tempyear=inputyear(i);
    for month=1:12
        xData((12*(i-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-15',]);
    end
end

plot(xData,discharge2)

p=polyfit(xData,discharge2',1);
line_trend2=xData*p(1)+p(2);
msaltplot=plot(xData,discharge2','color','k');
hold on
msaltplot2=plot(xData,line_trend2,'Color','r')

hold off
jpgname=strcat(outfile,testname, '_yangtze_river_discharge_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
xlabel('year')
ylabel('river discharge (meter^3 / sec)')
title(['Yangtze river discharge, ',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i')])
datetick('x','yymmm','keepticks')
axis tight;
% ylim(meanplotlev)
set(msaltplot,'LineWidth',2);
set(msaltplot2,'LineWidth',2);
set(gca,'FontSize',20);
% set(gca,'YTick',min(meanplotlev):1:max(meanplotlev))
% txt1=text(xData(3), min(meanplotlev)+1 ,['R = ', num2str(round(constant_cor(1,2),2)), ', '], 'FontSize', m_quiver_ref_text_fontsize); 
lgd=legend(['Yangtze river discharge']);

set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [1000, 400]);
        set(gcf,'PaperPosition', [0 0 1000 400])

saveas(gcf,jpgname,'jpg');
