close all; clear all;  clc;
warning off;
% all_region2 ={'NWP','NWP2'}
% all_region2 ={'ES', 'YS', 'SS', 'NWP2'}
all_testname2 = {'IPSL-CM5A-LR','IPSL-CM5A-MR','NorESM1-M','MPI-ESM-LR'};
all_testname3 = {'test53', 'test54', 'test55', 'test56'};
% all_testname3 = {'test65', 'test66', 'test67', 'test68'};

scenname = 'historical';
% scenname = 'rcp85';


close all;
clearvars '*' -except regionind2 all_region2 all_testname2 all_testname3 scenname
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
meanplotlev2 =[0 5];
% for snu_desktopd
%         testname=all_testname2{testnameind2}    % % need to change
inputyear = [1997:2005]; % % put year which you want to plot [year year ...]
inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
regionname = 'KS';
cmip5dir='D:\Data\Model\CMIP5\';


[~, ~, raw] = xlsread('Z:\내 드라이브\Data\Observation\Transport_Korea\obs_tsushima_197001_201209_ORIGIN_new.xls','TWC','A2:H553');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % 숫자형이 아닌 셀 찾기
raw(R) = {NaN}; % 숫자형이 아닌 셀 바꾸기

obs_ks = reshape([raw{:}],size(raw));

clearvars raw R;

% % Takikawa et al, 2005 sea level difference Jan 1965 to Dec 2001
% % Takikawa et al, 2005 ADCP mean feb 1997 to Aug 2002 
% % Fukudome et al, 2010 ADCP Mar 1997 to feb 2007 

ADCP=obs_ks(find(obs_ks(:,1)<=inputyear(end) & obs_ks(:,1)>=inputyear(1)),5);
SLD=obs_ks(find(obs_ks(:,1)<=inputyear(end) & obs_ks(:,1)>=inputyear(1)),6);
HYCOM=obs_ks(find(obs_ks(:,1)<=inputyear(end) & obs_ks(:,1)>=inputyear(1)),7);
Fuk=obs_ks(find(obs_ks(:,1)<=inputyear(end) & obs_ks(:,1)>=inputyear(1)),8);


for testind=1:length(all_testname2)
% %             get CMIP5 transport   
    all_tr_filename = ['D:\Data\Model\CMIP5\transport\', all_testname2{testind}, '_korea_tr1976_2005.txt'];
    gcm_all_tr=textread(all_tr_filename);
    for yearind=1:length(inputyear)
        tempyear=num2str(inputyear(yearind),'%04i');
        t_gcm_korea = gcm_all_tr((inputyear(yearind) - 1976)*12 +1 : (inputyear(yearind) - 1976)*12 +12,1)
        gcm_korea((yearind-1)*12+1:(yearind-1)*12+12)=t_gcm_korea;
        gcm_korea_yearly(yearind)=mean(t_gcm_korea);
    end
    trdata{testind}=gcm_korea';
    trdata3(testind,:)=gcm_korea';
    trdata_yearly{testind}=gcm_korea_yearly';
    trdata3_yearly(testind,:)=gcm_korea_yearly';
%     (Fuk - trdata{testind})    % Errors
%     (Fuk - trdata{testind}).^2   % Squared Error
%     mean((Fuk - trdata{testind}).^2, 'omitnan')   % Mean Squared Error
    gcm_rmse(testind) = sqrt(mean((Fuk - trdata{testind}).^2, 'omitnan'));  % Root Mean Squared Error
    tempcorr=corrcoef(Fuk(isfinite(Fuk)), trdata{testind}(isfinite(Fuk)));
    gcm_corr(testind) = tempcorr(1,2);
    
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
        korea_yearly(i) = mean(t_korea);
    end
    trdata2{testind}=korea';
    trdata2_yearly{testind}=korea_yearly';
    
    %     (Fuk - trdata2{testind})    % Errors
%     (Fuk - trdata2{testind}).^2   % Squared Error
%     mean((Fuk - trdata2{testind}).^2, 'omitnan')   % Mean Squared Error
    rcm_rmse(testind) = sqrt(mean((Fuk - trdata2{testind}).^2, 'omitnan'));  % Root Mean Squared Error
    tempcorr=corrcoef(Fuk(isfinite(Fuk)), trdata2{testind}(isfinite(Fuk)));
    rcm_corr(testind) = tempcorr(1,2);
end

trdata{testind+1}=mean(trdata3,1);  % ens
trdata_yearly{testind+1}=mean(trdata3_yearly,1);
trdata2_temp=cell2mat(trdata2);
trdata2_yearly_temp=cell2mat(trdata2_yearly);
trdata2{testind+1}=mean(trdata2_temp,2);     % ens   
trdata2_yearly{testind+1}=mean(trdata2_yearly_temp,2);     % ens   

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

for i =1:length(inputyear) 
    tempyear=inputyear(i);
    xData_yearly(i) = datenum([num2str(tempyear),'-',num2str(6,'%02i'),'-30']);
end


% start-------------------- tr time series(GCM)
for folding=1:1
      jpgname=strcat(outfile, '_', 'GCM', '_tr_', scenname, '_', ...
          num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
            hold on
            for testind=1:length(all_testname2)
                mslplot2_all{testind}=plot(xData,trdata{testind},'b');
            end
%             mslplot=plot(xData),Fuk(isfinite(Fuk)),'k')
            mslplot2_all{testind+1} = plot(xData,Fuk,'k');
            
            set(mslplot2_all{1},'Marker','*');
            set(mslplot2_all{2},'Marker','^');
            set(mslplot2_all{3},'Marker','o');
            set(mslplot2_all{4},'Marker','+');
%             set(mslplot2_all{5},'Marker','+');
            
            xlabel('Year')
            ylabel('Transport (Sv)')
%             title([regionname, ', Transport(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') '])
%             datetick('x','yyyymm','keepticks')
            datetick('x','yyyy')

            axis tight;
            ylim(meanplotlev2)
%             set(mslplot_all{5},'LineWidth',2);
            set(mslplot2_all{5},'LineWidth',2);
%             set(mslplot,'LineWidth',2);
            set(gca,'FontSize',25);
            grid on
            lgd=legend('GCM-IPSL-LR','GCM-IPSL-MR','GCM-Nor','GCM-MPI','Obs');
            set(lgd,'FontSize',15);
%             set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
            set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);

            set(lgd,'Orientation','horizontal');

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
      jpgname=strcat(outfile, '_', 'RCM', '_tr_', scenname, '_', ...
          num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
            hold on
            for testind=1:length(all_testname2)
                mslplot2_all{testind}=plot(xData,trdata2{testind},'r');
            end
            mslplot2_all{testind+1} = plot(xData,Fuk,'k');

            set(mslplot2_all{1},'Marker','*');
            set(mslplot2_all{2},'Marker','^');
            set(mslplot2_all{3},'Marker','o');
            set(mslplot2_all{4},'Marker','+');
%             set(mslplot2_all{5},'Marker','+');
            
            xlabel('Year')
            ylabel('Transport (Sv)')
%             title([regionname, ', Transport(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') '])
%             datetick('x','yyyymm','keepticks')
            datetick('x','yyyy')

            axis tight;
            ylim(meanplotlev2)
%             set(mslplot_all{5},'LineWidth',2);
            set(mslplot2_all{5},'LineWidth',2);
%             set(mslplot,'LineWidth',2);
            set(gca,'FontSize',25);
            grid on
            lgd=legend('RCM-IPSL-LR','RCM-IPSL-MR','RCM-Nor','RCM-MPI','Obs');
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


            set(gcf,'PaperPosition', [0 0 36 12])   
            hold off
            saveas(gcf,jpgname,'jpg');
            grid off
%         end
        close all;
end
% end-------------------- tr time series (RCM)

% 
% % start-------------------- yearly tr time series(GCM)
% for folding=1:1
%       jpgname=strcat(outfile, '_', 'GCM', '_yearly_tr_', scenname, '_', ...
%           num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
%             hold on
%             for testind=1:length(all_testname2)+1
%                 mslplot2_all{testind}=plot(xData_yearly,trdata_yearly{testind},'b');
%             end
% %             mslplot=plot(xData),Fuk(isfinite(Fuk)),'k')
%             
%             set(mslplot2_all{1},'Marker','*');
%             set(mslplot2_all{2},'Marker','^');
%             set(mslplot2_all{3},'Marker','o');
%             set(mslplot2_all{4},'Marker','+');
% %             set(mslplot2_all{5},'Marker','+');
%             
%             xlabel('Year')
%             ylabel('Transport (Sv)')
%             title([regionname, ', Transport(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') '])
% %             datetick('x','yyyymm','keepticks')
%             datetick('x','yyyy')
% 
%             axis tight;
%             ylim(meanplotlev2)
% %             set(mslplot_all{5},'LineWidth',2);
%             set(mslplot2_all{5},'LineWidth',2);
% %             set(mslplot,'LineWidth',2);
%             set(gca,'FontSize',25);
%             grid on
%             lgd=legend('GCM-IPSL-LR','GCM-IPSL-MR','GCM-Nor','GCM-MPI','Ens');
%             set(lgd,'FontSize',15);
%             set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
%             set(lgd,'Orientation','horizontal');
% 
%                         set(gcf,'PaperPosition', [0 0 26 20])   
%             hold off
%             saveas(gcf,jpgname,'jpg');
%             save(['D:\Data\Model\CMIP5\transport\all\', scenname, '_GCM_all.mat'], 'xData_yearly', 'trdata_yearly');
%             grid off
% %         end
%         close all;
% end
% % end-------------------- yearly tr time series (GCM)
% 
% % start-------------------- yearly tr time series(RCM)
% for folding=1:1
%       jpgname=strcat(outfile, '_', 'RCM', '_yearly_tr_', scenname, '_', ...
%           num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
%             hold on
%             for testind=1:length(all_testname2)+1
%                 mslplot2_all{testind}=plot(xData_yearly,trdata2_yearly{testind},'r')
%             end
%             
%             set(mslplot2_all{1},'Marker','*');
%             set(mslplot2_all{2},'Marker','^');
%             set(mslplot2_all{3},'Marker','o');
%             set(mslplot2_all{4},'Marker','+');
% %             set(mslplot2_all{5},'Marker','+');
%             
%             xlabel('Year')
%             ylabel('Transport (Sv)')
%             title([regionname, ', Transport(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') '])
% %             datetick('x','yyyymm','keepticks')
%             datetick('x','yyyy')
% 
%             axis tight;
%             ylim(meanplotlev2)
% %             set(mslplot_all{5},'LineWidth',2);
%             set(mslplot2_all{5},'LineWidth',2);
% %             set(mslplot,'LineWidth',2);
%             set(gca,'FontSize',25);
%             grid on
%             lgd=legend('GCM-IPSL-LR','GCM-IPSL-MR','GCM-Nor','GCM-MPI','Ens');
%             set(lgd,'FontSize',15);
%             set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
%             set(lgd,'Orientation','horizontal');
% %             constant_cor=corrcoef(RCM_interped_sla_yearly_mean(length(RCM_testnames),:),cmems_sla_yearly_mean);
% %             txt1=text(xData2(5), min(meanplotlev2)+diff(meanplotlev2)/16.0 ,['R(RCM) = ', num2str(round(constant_cor(1,2),2))], 'FontSize', 20); 
% %             constant_cor=corrcoef(GCM_interped_sla_yearly_mean(length(RCM_testnames),:),cmems_sla_yearly_mean);
% %             txt1=text(xData2(10), min(meanplotlev2)+diff(meanplotlev2)/16.0 ,['R(GCM) = ', num2str(round(constant_cor(1,2),2))], 'FontSize', 20); 
% 
% %             for tind=1:size(GCM_interped_sla_yearly_mean,2)
% %                 model_std(tind)=std(GCM_interped_sla_yearly_mean(:,tind));
% %             end
% %             meanstd=mean(model_std);
% 
% 
%                         set(gcf,'PaperPosition', [0 0 26 20])   
%             hold off
%             saveas(gcf,jpgname,'jpg');
%             save(['D:\Data\Model\ROMS\nwp_1_20\transport_all\', scenname, '_RCM_all.mat'], 'xData_yearly', 'trdata2_yearly');
%             grid off
% %         end
%         close all;
% end
% % end-------------------- yearly tr time series (RCM)
% 












