close all; clear all;  clc;
warning off;
% all_region2 ={'NWP','NWP2'}
% all_region2 ={'ES', 'YS', 'SS', 'NWP2'}
all_testname2 = {'IPSL-CM5A-LR','IPSL-CM5A-MR','NorESM1-M','MPI-ESM-LR'};
all_testname3 = {'test53', 'test54', 'test55', 'test56'};


close all;
clearvars '*' -except regionind2 all_region2 all_testname2 all_testname3
% % % 
% % % Read Model SST
% % % interp
% % % get RMS
% % % get BIAS
system_name=computer;
if (strcmp(system_name,'PCWIN64'))
    % % for windows
    dropboxpath='C:\Users\USER\Dropbox';
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
meanplotlev2 =[0 4.5];
% for snu_desktopd
%         testname=all_testname2{testnameind2}    % % need to change
inputyear = [1997:2005]; % % put year which you want to plot [year year ...]
inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
regionname = 'KS';
for testind=1:length(all_testname2)
% %             get CMIP5 transport            
%     filename = ['D:\Data\Model\CMIP5\transport\', all_testname2{testind}, '_korea_tr', num2str(min(inputyear),'%04i'), '_', num2str(max(inputyear),'%04i'), '.txt'];
    filename = ['D:\Data\Model\CMIP5\transport\', all_testname2{testind}, '_korea_tr', num2str(1976,'%04i'), '_', num2str(2005,'%04i'), '.txt'];
    delimiter = ' ';
    formatSpec = '%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string',  'ReturnOnError', false);
    fclose(fileID);
    trdata3{testind} = dataArray{1:end-1};
    clearvars filename delimiter formatSpec fileID dataArray ans;
% %             get RCM transport
    for i=1:length(inputyear)
        tempyear=num2str(inputyear(i),'%04i');
        filename2 = ['D:\Data\Model\ROMS\nwp_1_20\',all_testname3{testind},'\run\transport','\nwp_1_20_monthly_',all_testname3{testind},'_',tempyear,'.txt'];
        startRow = 2;
        formatSpec = '%8f%9f%9f%9f%9f%9f%s%[^\n\r]';
        fileID = fopen(filename2,'r');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
        dataArray{7} = strtrim(dataArray{7});
        fclose(fileID);
        t_korea = dataArray{:, 1};
        for month=1:12
            korea((12*(i-1))+month) = t_korea(month);
        end
    end
    trdata2{testind}=korea';
end
trdata_temp3=cell2mat(trdata3);
trdata_temp(1:(max(inputyear)-min(inputyear))*12+12,:)=trdata_temp3((min(inputyear)-1976)*12+1 : (max(inputyear)-1976)*12+12,:);
for testind=1:length(all_testname2)
    trdata{testind}=trdata_temp(:,testind);
end
trdata{testind+1}=mean(trdata_temp,2);  % ens
trdata2_temp=cell2mat(trdata2);
trdata2{testind+1}=mean(trdata2_temp,2);     % ens   

% % get observation transport
[~, ~, raw] = xlsread('Z:\내 드라이브\Data\Observation\Transport_Korea\obs_tsushima_197001_201209_ORIGIN_new.xls','TWC','A2:H553');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % 숫자형이 아닌 셀 찾기
raw(R) = {NaN}; % 숫자형이 아닌 셀 바꾸기

obs_ks = reshape([raw{:}],size(raw));

clearvars raw R;

ADCP=obs_ks(find(obs_ks(:,1)<=inputyear(end) & obs_ks(:,1)>=inputyear(1)),5);
SLD=obs_ks(find(obs_ks(:,1)<=inputyear(end) & obs_ks(:,1)>=inputyear(1)),6);
HYCOM=obs_ks(find(obs_ks(:,1)<=inputyear(end) & obs_ks(:,1)>=inputyear(1)),7);
Fuk=obs_ks(find(obs_ks(:,1)<=inputyear(end) & obs_ks(:,1)>=inputyear(1)),8);
SLD(isfinite(Fuk))=NaN;

corrcoef(trdata{5}(isfinite(ADCP)), ADCP(isfinite(ADCP)))
corrcoef(trdata2{5}(isfinite(ADCP)), ADCP(isfinite(ADCP)))

if (strcmp(system_name,'PCWIN64'))
    % % for windows
    figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\5th_year\figure\nwp_1_20\','all','\'); % % where figure files will be saved
    param_script ='C:\Users\USER\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m'
elseif (strcmp(system_name,'GLNXA64'))
end



figdir=[figrawdir,'Transport\'];
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
% % % % % % 
% % % % % % hold on
% % % % % % 
% % % % % % 
% % % % % % for i=1:length(trdata)
% % % % % %     trplot{i}=plot(xData,trdata{i},'b');
% % % % % %     if i<length(trdata)
% % % % % %         trplot{i}.Color(4) = .1;
% % % % % %         set(get(get(trplot{i},'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% % % % % %     else
% % % % % %     end
% % % % % %     set(trplot{i},'LineWidth',2);
% % % % % % 
% % % % % % end
% % % % % % 
% % % % % % for i=1:length(trdata2)
% % % % % %     trplot2{i}=plot(xData,trdata2{i},'r');
% % % % % %     if i<length(trdata2)
% % % % % %         trplot2{i}.Color(4) = .1;
% % % % % %         set(get(get(trplot2{i},'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% % % % % %     else
% % % % % %     end
% % % % % %     set(trplot2{i},'LineWidth',2);
% % % % % % end
% % % % % % 
% % % % % % trplot3{1}=plot(xData,SLD,'k');
% % % % % % trplot3{2}=plot(xData,Fuk,'k');
% % % % % % set(get(get(trplot3{2},'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% % % % % % set(trplot3{1},'LineWidth',2);
% % % % % % set(trplot3{2},'LineWidth',2);
% % % % % % 
% % % % % % jpgname=strcat(outfile, '_', 'all', '_tr_hist_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
% % % % % % xlabel('Year')
% % % % % % ylabel('Transport (Sv)')
% % % % % % title([regionname, ', Transport(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') '])
% % % % % % datetick('x','yyyy','keepticks')
% % % % % % axis tight;
% % % % % % lgd=legend(['GCM Model Ensemble'], ['RCM Model Ensemble'], ['Observation']);
% % % % % % set(gca,'FontSize',20);
% % % % % % grid on
% % % % % % hold off
% % % % % % %         constant_cor=corrcoef(tr{5},recon_tr);
% % % % % % %         txt1=text(xData(5), min(meanplotlev)+diff(meanplotlev)/16.0 ,['R = ', num2str(round(constant_cor(1,2),2))], 'FontSize', 20); 
% % % % % % set(gcf,'PaperPosition', [0 0 40 18]) 
% % % % % % saveas(gcf,jpgname,'jpg');
% % % % % % grid off
% % % % % % 
% % % % % % 
% % % % % % close all;


% cmap = get(groot,'defaultaxescolororder');
% cmap_b = rgb2hsv(cmap);
% cmap_b(:,2) = cmap_b(:,2)*.3;
% cmap_b(:,3) = cmap_b(:,3)*.3+.7;
% cmap_b = hsv2rgb(cmap_b);


% start-------------------- tr time series(GCM)
for folding=1:1
      jpgname=strcat(outfile, '_', 'GCM', '_tr_hist_', ...
          num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
            hold on
            for testind=1:length(all_testname2)+1
                mslplot2_all{testind}=plot(xData(isfinite(Fuk)),trdata{testind}(isfinite(Fuk)),'b')
                constant_cor(:,:,testind)=corrcoef(trdata{testind}(isfinite(Fuk)),Fuk(isfinite(Fuk)));
            end
            mslplot=plot(xData(isfinite(Fuk)),Fuk(isfinite(Fuk)),'k')
            
            set(mslplot2_all{1},'Marker','*');
            set(mslplot2_all{2},'Marker','^');
            set(mslplot2_all{3},'Marker','o');
            set(mslplot2_all{4},'Marker','+');
%             set(mslplot2_all{5},'Marker','+');
            
            xlabel('Year')
            ylabel('Transport (Sv)')
            title([regionname, ', Transport(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') '])
%             datetick('x','yyyymm','keepticks')
            datetick('x','yyyy')

            axis tight;
            ylim(meanplotlev2)
%             set(mslplot_all{5},'LineWidth',2);
            set(mslplot2_all{5},'LineWidth',2);
            set(mslplot,'LineWidth',2);
            set(gca,'FontSize',15);
            grid on
            lgd=legend('GCM-IPSL-LR','GCM-IPSL-MR','GCM-Nor','GCM-MPI','Ens','ADCP');
            set(lgd,'FontSize',15);
            set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
            set(lgd,'Orientation','horizontal');
%             constant_cor=corrcoef(RCM_interped_sla_yearly_mean(length(RCM_testnames),:),cmems_sla_yearly_mean);
%             txt1=text(xData2(5), min(meanplotlev2)+diff(meanplotlev2)/16.0 ,['R(RCM) = ', num2str(round(constant_cor(1,2),2))], 'FontSize', 20); 
%             constant_cor=corrcoef(GCM_interped_sla_yearly_mean(length(RCM_testnames),:),cmems_sla_yearly_mean);
%             txt1=text(xData2(10), min(meanplotlev2)+diff(meanplotlev2)/16.0 ,['R(GCM) = ', num2str(round(constant_cor(1,2),2))], 'FontSize', 20); 

%             for tind=1:size(GCM_interped_sla_yearly_mean,2)
%                 model_std(tind)=std(GCM_interped_sla_yearly_mean(:,tind));
%             end
%             meanstd=mean(model_std);
            txt1=text(xData(5), min(meanplotlev2)+diff(meanplotlev2)/16.0 , ...
                ['Corr = ', num2str(round(constant_cor(1,2,1),2)), ', '  ...
                    num2str(round(constant_cor(1,2,2),2)), ', '  ...
                    num2str(round(constant_cor(1,2,3),2)), ', '  ...
                    num2str(round(constant_cor(1,2,4),2)), ', '  ...
                    num2str(round(constant_cor(1,2,5),2)),   ...
                    ], 'FontSize', 20); 

            set(gcf,'PaperPosition', [0 0 36 12]) 
            hold off
            saveas(gcf,jpgname,'jpg');
            grid off
%         end
        close all;
end
% end-------------------- tr time series (GCM)


% start-------------------- tr time series(RCM)
for folding=1:1
      jpgname=strcat(outfile, '_', 'RCM', '_tr_hist_', ...
          num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
            hold on
            for testind=1:length(all_testname2)+1
                mslplot2_all{testind}=plot(xData(isfinite(Fuk)),trdata2{testind}(isfinite(Fuk)),'r')
                constant_cor(:,:,testind)=corrcoef(trdata2{testind}(isfinite(Fuk)),Fuk(isfinite(Fuk)));
            end
            mslplot=plot(xData(isfinite(Fuk)),Fuk(isfinite(Fuk)),'k')
            
            set(mslplot2_all{1},'Marker','*');
            set(mslplot2_all{2},'Marker','^');
            set(mslplot2_all{3},'Marker','o');
            set(mslplot2_all{4},'Marker','+');
%             set(mslplot2_all{5},'Marker','+');
            
            xlabel('Year')
            ylabel('Transport (Sv)')
            title([regionname, ', Transport(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') '])
%             datetick('x','yyyymm','keepticks')
            datetick('x','yyyy')

            axis tight;
            ylim(meanplotlev2)
%             set(mslplot_all{5},'LineWidth',2);
            set(mslplot2_all{5},'LineWidth',2);
            set(mslplot,'LineWidth',2);
            set(gca,'FontSize',15);
            grid on
            lgd=legend('GCM-IPSL-LR','GCM-IPSL-MR','GCM-Nor','GCM-MPI','Ens','ADCP');
            set(lgd,'FontSize',15);
            set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
            set(lgd,'Orientation','horizontal');
%             constant_cor=corrcoef(RCM_interped_sla_yearly_mean(length(RCM_testnames),:),cmems_sla_yearly_mean);
%             txt1=text(xData2(5), min(meanplotlev2)+diff(meanplotlev2)/16.0 ,['R(RCM) = ', num2str(round(constant_cor(1,2),2))], 'FontSize', 20); 
%             constant_cor=corrcoef(GCM_interped_sla_yearly_mean(length(RCM_testnames),:),cmems_sla_yearly_mean);
%             txt1=text(xData2(10), min(meanplotlev2)+diff(meanplotlev2)/16.0 ,['R(GCM) = ', num2str(round(constant_cor(1,2),2))], 'FontSize', 20); 

%             for tind=1:size(GCM_interped_sla_yearly_mean,2)
%                 model_std(tind)=std(GCM_interped_sla_yearly_mean(:,tind));
%             end
%             meanstd=mean(model_std);
            txt1=text(xData(5), min(meanplotlev2)+diff(meanplotlev2)/16.0 , ...
                ['Corr = ', num2str(round(constant_cor(1,2,1),2)), ', '  ...
                    num2str(round(constant_cor(1,2,2),2)), ', '  ...
                    num2str(round(constant_cor(1,2,3),2)), ', '  ...
                    num2str(round(constant_cor(1,2,4),2)), ', '  ...
                    num2str(round(constant_cor(1,2,5),2)),   ...
                    ], 'FontSize', 20); 

            set(gcf,'PaperPosition', [0 0 36 12]) 
            hold off
            saveas(gcf,jpgname,'jpg');
            grid off
%         end
        close all;
end
% end-------------------- tr time series (RCM)