close all; clear all;  clc;
warning off;

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
meanplotlev2 =[2 3];
% for snu_desktopd
%         testname=all_testname2{testnameind2}    % % need to change
inputyear = [2006:2100]; % % put year which you want to plot [year year ...]
inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
regionname = 'KS';
cmip5dir='D:\Data\Model\CMIP5\';


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

% start-------------------- yearly tr time series(GCM)
for folding=1:1
      jpgname=strcat(outfile, '_', 'ALL_scen_ens', '_yearly_tr', '_', ...
          num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
            hold on
            load(['D:\Data\Model\CMIP5\transport\all\', 'rcp26', '_GCM_all.mat']);
            gcm_rcp26_trdata_yearly=trdata_yearly;
            load(['D:\Data\Model\CMIP5\transport\all\', 'rcp85', '_GCM_all.mat']);
            gcm_rcp85_trdata_yearly=trdata_yearly;
            load(['D:\Data\Model\ROMS\nwp_1_20\transport_all\', 'rcp26', '_RCM_all.mat']);
            rcm_rcp26_trdata_yearly=trdata2_yearly;
            load(['D:\Data\Model\ROMS\nwp_1_20\transport_all\', 'rcp85', '_RCM_all.mat']);
            rcm_rcp85_trdata_yearly=trdata2_yearly;

%             for testind=1:length(all_testname2)+1
%                 mslplot2_all{testind}=plot(xData_yearly,trdata_yearly{testind},'b');
%             end
%             mslplot=plot(xData),Fuk(isfinite(Fuk)),'k')
            mslplot2_all{1}=plot(xData_yearly,gcm_rcp26_trdata_yearly{5},'b--');
            mslplot2_all{2}=plot(xData_yearly,gcm_rcp85_trdata_yearly{5},'b');
            mslplot2_all{3}=plot(xData_yearly,rcm_rcp26_trdata_yearly{5},'r--');
            mslplot2_all{4}=plot(xData_yearly,rcm_rcp85_trdata_yearly{5},'r');
            
            
%             set(mslplot2_all{1},'Marker','*');
%             set(mslplot2_all{2},'Marker','^');
%             set(mslplot2_all{3},'Marker','o');
%             set(mslplot2_all{4},'Marker','+');

%             set(mslplot2_all{5},'Marker','+');
            
            xlabel('Year')
            ylabel('Transport (Sv)')
            title([regionname, ', Transport(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),') '])
%             datetick('x','yyyymm','keepticks')
            datetick('x','yyyy')

            axis tight;
            ylim(meanplotlev2)
            for i=1:4
                set(mslplot2_all{i},'LineWidth',2);
            end
%             set(mslplot_all{5},'LineWidth',2);
%             set(mslplot2_all{5},'LineWidth',2);
%             set(mslplot,'LineWidth',2);
            set(gca,'FontSize',25);
            grid on
            lgd=legend('GCM-RCP2.6','GCM-RCP8.5','RCM-RCP2.6','RCM-RCP8.5');
            set(lgd,'FontSize',15);
            set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
            set(lgd,'Orientation','horizontal');

                        set(gcf,'PaperPosition', [0 0 26 20])   
            hold off
            saveas(gcf,jpgname,'jpg');
            grid off
%         end
        close all;
end
% end-------------------- yearly tr time series (GCM)
