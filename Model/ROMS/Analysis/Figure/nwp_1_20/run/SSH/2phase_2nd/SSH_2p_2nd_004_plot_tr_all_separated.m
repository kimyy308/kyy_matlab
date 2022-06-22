close all; clear all;  clc;
warning off;

tmp.fs=filesep;
if (strcmp(computer,'PCWIN64'))
    tmp.dropboxpath = 'C:\Users\User\Dropbox';
else
    tmp.dropboxpath = '/home/kimyy/Dropbox';
end
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
[tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Model' ...
    tmp.fs, 'ROMS', tmp.fs, 'Analysis', tmp.fs, 'Figure', tmp.fs, 'nwp_1_20', tmp.fs ...
    'run', tmp.fs, 'SSH', tmp.fs, '2phase_2nd', tmp.fs, 'subroutine']));



% RCM_info.name = {'test2117', 'test2118', 'test2119', 'test2120', 'test2121'};
RCM_info.name = {'test2127', 'test2128', 'test2129', 'test2130', 'test2131'};

% RCM_info.name = {'CNRM-ESM2-1', 'EC-Earth3-Veg', 'ACCESS-CM2'};


% close all;
% clearvars '*' -except regionind2 all_region2 RCM_info.name RCM_info.name scenname
% % % 



% shadlev = [0 35];
% rms_shadlev = [0 4];
% trendlev = [-10 10];  %% trend lev
% conlev  = 0:5:35;
% meanplotlev =[-0.3 0.3];
figs.meanplotlev2 =[0 4];


% RCM_info.years = [1985:2014]; % % put year which you want to plot [year year ...]
RCM_info.years = [2015:2100]; % % put year which you want to plot [year year ...]

RCM_info.months = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
tmp.regionname = 'KS';
% cmip5dir='D:\Data\Model\CMIP5\';

for testind=1:length(RCM_info.name)
        tmp.testname=RCM_info.name{testind};
        [RCM_info.scenname, tmp.error_status] = Func_0013_RCM_CMIP6_scenname(tmp.testname);

% %             get RCM transport
    for i=1:length(RCM_info.years)
        tmp.tempyear=num2str(RCM_info.years(i),'%04i');
        tmp.filename_model = ['D:\Data\Model\ROMS\nwp_1_20\transport_barot\',tmp.testname,'\','nwp_1_20_monthly_', tmp.testname,'_',tmp.tempyear,'.txt'];
%         nwp_1_20_monthly_test2107_2015
        tmp.startRow = 2;
        tmp.formatSpec = '%8f%9f%9f%9f%9f%9f%s%[^\n\r]';
        tmp.fileID = fopen(tmp.filename_model,'r');
        tmp.dataArray = textscan(tmp.fileID, tmp.formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines' ,tmp.startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
        tmp.dataArray{7} = strtrim(tmp.dataArray{7});
        fclose(tmp.fileID);
        tmp.t_korea = tmp.dataArray{:, 1};
        tmp.t_tsugaru = tmp.dataArray{:, 2};
        tmp.t_soya = tmp.dataArray{:, 3};
        for month=1:12
            tmp.korea_all((12*(i-1))+month) = tmp.t_korea(month);
            tmp.tsugaru_all((12*(i-1))+month) = tmp.t_tsugaru(month);
            tmp.soya_all((12*(i-1))+month) = tmp.t_soya(month);
        end
        tmp.korea_yearly_all(i) = mean(tmp.t_korea);
        tmp.tsugaru_yearly_all(i) = mean(tmp.t_tsugaru);
        tmp.soya_yearly_all(i) = mean(tmp.t_soya);
    end
    RCM_data.tr{testind}=tmp.korea_all';
    RCM_data.tr_tsugaru{testind}=tmp.tsugaru_all';
    RCM_data.tr_soya{testind}=tmp.soya_all';
    RCM_data.tr_yearly{testind}=tmp.korea_yearly_all';
    RCM_data.tr_yearly_tsugaru{testind}=tmp.tsugaru_yearly_all';
    RCM_data.tr_yearly_soya{testind}=tmp.soya_yearly_all';
end
% GCM transport
% trdata{testind+1}=mean(trdata3,1);  % ens
% trdata_yearly{testind+1}=mean(trdata3_yearly,1);

RCM_data.tr_temp=cell2mat(RCM_data.tr);
RCM_data.tr_yearly_temp=cell2mat(RCM_data.tr_yearly);
RCM_data.tr{testind+1}=mean(RCM_data.tr_temp,2);     % ens   
RCM_data.tr_yearly{testind+1}=mean(RCM_data.tr_yearly_temp,2);     % ens   

% % get observation transport
[~, ~, raw] = xlsread('D:\Data\Observation\Transport_Korea\obs_tsushima_197001_201209_ORIGIN_new.xls','TWC','A2:H553');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % 숫자형이 아닌 셀 찾기
raw(R) = {NaN}; % 숫자형이 아닌 셀 바꾸기

OBS_data.obs_ks_all = reshape([raw{:}],size(raw));

clearvars raw R;

OBS_data.ADCP=OBS_data.obs_ks_all(find(OBS_data.obs_ks_all(:,1)<=RCM_info.years(end) & OBS_data.obs_ks_all(:,1)>=RCM_info.years(1)),5);
OBS_data.SLD=OBS_data.obs_ks_all(find(OBS_data.obs_ks_all(:,1)<=RCM_info.years(end) & OBS_data.obs_ks_all(:,1)>=RCM_info.years(1)),6);
OBS_data.HYCOM=OBS_data.obs_ks_all(find(OBS_data.obs_ks_all(:,1)<=RCM_info.years(end) & OBS_data.obs_ks_all(:,1)>=RCM_info.years(1)),7);
OBS_data.Fuk=OBS_data.obs_ks_all(find(OBS_data.obs_ks_all(:,1)<=RCM_info.years(end) & OBS_data.obs_ks_all(:,1)>=RCM_info.years(1)),8);
OBS_data.SLD(isfinite(OBS_data.Fuk))=NaN;

OBS_data.OBS_data.ADCP_res=reshape(OBS_data.ADCP,[12, length(OBS_data.ADCP)/12]);
OBS_data.OBS_data.ADCP_yearly=mean(OBS_data.OBS_data.ADCP_res,1);
% 
% corrcoef(trdata{5}(isfinite(OBS_data.ADCP)), OBS_data.ADCP(isfinite(OBS_data.ADCP)))
% corrcoef(RCM_data.tr{5}(isfinite(OBS_data.ADCP)), OBS_data.ADCP(isfinite(OBS_data.ADCP)))

if (strcmp(computer,'PCWIN64'))
    % % for windows
    dirs.figrawdir =strcat('D:\MEPL\project\SSH\7th_year(2022)\figure\nwp_1_20\','all','\'); % % where figure files will be saved
    param_script ='C:\Users\USER\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m';
elseif (strcmp(system_name,'GLNXA64'))
end

dirs.figdir=[dirs.figrawdir,'Transport\'];
if (exist(strcat(dirs.figdir) , 'dir') ~= 7)
    mkdir(strcat(dirs.figdir));
end 
tmp.outfile = strcat(dirs.figdir,tmp.regionname);

for i =1:length(RCM_info.years) 
    tmp.tempyear=RCM_info.years(i);
    for month=1:12
        xData((12*(i-1))+month) = datenum([num2str(tmp.tempyear),'-',num2str(month,'%02i'),'-01',]);
    end
end

for i =1:length(RCM_info.years) 
    tmp.tempyear=RCM_info.years(i);
    xData_yearly(i) = datenum([num2str(tmp.tempyear),'-',num2str(6,'%02i'),'-30']);
end

% start-------------------- tr time series(RCM)
      jpgname=strcat(tmp.outfile, '_', 'RCM', '_tr_', RCM_info.scenname, '_', ...
          num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), '.jpg'); %% ~_year_month.jpg
            hold on
            for testind=1:length(RCM_info.name)
                figs.mslplot2_all{testind}=plot(xData,RCM_data.tr{testind},'r');
            end
            figs.mslplot2_all{testind+1}=plot(xData,OBS_data.ADCP, 'k');            
            set(figs.mslplot2_all{1},'Marker','*');
            set(figs.mslplot2_all{2},'Marker','^');
            set(figs.mslplot2_all{3},'Marker','o');
            set(figs.mslplot2_all{4},'Marker','+');
            set(figs.mslplot2_all{5},'Marker','square');
            set(figs.mslplot2_all{6},'Marker','pentagram');
            set(figs.mslplot2_all{6},'LineWidth',2);
            
            xlabel('Year')
            ylabel('Transport (Sv)')
            title([tmp.regionname, ', Transport(',num2str(min(RCM_info.years),'%04i'),'-',num2str(max(RCM_info.years),'%04i'),') '])
%             datetick('x','yyyymm','keepticks')
            datetick('x','yyyy')

            axis tight;
            ylim(figs.meanplotlev2)
%             set(mslplot_all{5},'LineWidth',2);
%             set(figs.mslplot2_all{5},'LineWidth',2);
%             set(mslplot,'LineWidth',2);
            set(gca,'FontSize',25);
            grid on
            figs.lgd=legend({'RCM-CNRM'; 'RCM-EC-Veg'; 'RCM-ACC'; 'RCM-CNRM-HR'; 'RCM-CMCC'; 'OBS_data.ADCP'}, 'NumColumns',3);
            set(figs.lgd,'FontSize',15);
%             set(figs.lgd,'Position',[0.13 0.88, 0.775, 0.03]);  %top cneter
            set(figs.lgd,'Orientation','horizontal');
            set(figs.lgd,'location','south'); 

%             constant_cor=corrcoef(RCM_interped_sla_yearly_mean(length(RCM_testnames),:),cmems_sla_yearly_mean);
%             txt1=text(xData2(5), min(figs.meanplotlev2)+diff(figs.meanplotlev2)/16.0 ,['R(RCM) = ', num2str(round(constant_cor(1,2),2))], 'FontSize', 20); 
%             constant_cor=corrcoef(GCM_interped_sla_yearly_mean(length(RCM_testnames),:),cmems_sla_yearly_mean);
%             txt1=text(xData2(10), min(figs.meanplotlev2)+diff(figs.meanplotlev2)/16.0 ,['R(GCM) = ', num2str(round(constant_cor(1,2),2))], 'FontSize', 20); 

%             for tind=1:size(GCM_interped_sla_yearly_mean,2)
%                 model_std(tind)=std(GCM_interped_sla_yearly_mean(:,tind));
%             end
%             meanstd=mean(model_std);


                        set(gcf,'PaperPosition', [0 0 26 20])   
            hold off
            saveas(gcf,jpgname,'jpg'); RemoveWhiteSpace([], 'file', jpgname);
            grid off
%         end
        close all;

% end-------------------- tr time series (RCM)


% % start-------------------- yearly tr time series(GCM)
% for folding=1:1
%       jpgname=strcat(tmp.outfile, '_', 'GCM', '_yearly_tr_', scenname, '_', ...
%           num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), '.jpg'); %% ~_year_month.jpg
%             hold on
%             for testind=1:length(RCM_info.name)+1
%                 figs.mslplot2_all{testind}=plot(xData_yearly,trdata_yearly{testind},'b');
%             end
% %             mslplot=plot(xData),OBS_data.Fuk(isfinite(OBS_data.Fuk)),'k')
%             
%             set(figs.mslplot2_all{1},'Marker','*');
%             set(figs.mslplot2_all{2},'Marker','^');
%             set(figs.mslplot2_all{3},'Marker','o');
%             set(figs.mslplot2_all{4},'Marker','+');
% %             set(figs.mslplot2_all{5},'Marker','+');
%             
%             xlabel('Year')
%             ylabel('Transport (Sv)')
%             title([tmp.regionname, ', Transport(',num2str(min(RCM_info.years),'%04i'),'-',num2str(max(RCM_info.years),'%04i'),') '])
% %             datetick('x','yyyymm','keepticks')
%             datetick('x','yyyy')
% 
%             axis tight;
%             ylim(figs.meanplotlev2)
% %             set(mslplot_all{5},'LineWidth',2);
%             set(figs.mslplot2_all{5},'LineWidth',2);
% %             set(mslplot,'LineWidth',2);
%             set(gca,'FontSize',25);
%             grid on
%             figs.lgd=legend('GCM-IPSL-LR','GCM-IPSL-MR','GCM-Nor','GCM-MPI','Ens');
%             set(figs.lgd,'FontSize',15);
%             set(figs.lgd,'Position',[0.13 0.88, 0.775, 0.03]);
%             set(figs.lgd,'Orientation','horizontal');
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

% start-------------------- yearly tr time series(RCM)
      jpgname=strcat(tmp.outfile, '_', 'RCM', '_yearly_tr_', RCM_info.scenname, '_', ...
          num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), '.jpg'); %% ~_year_month.jpg
            hold on
            for testind=1:length(RCM_info.name)
                figs.mslplot2_all{testind}=plot(xData_yearly,RCM_data.tr_yearly{testind},'r')
            end
            figs.mslplot2_all{testind+1}=plot(xData_yearly,OBS_data.OBS_data.ADCP_yearly, 'k');

            set(figs.mslplot2_all{1},'Marker','*');
            set(figs.mslplot2_all{2},'Marker','^');
            set(figs.mslplot2_all{3},'Marker','o');
            set(figs.mslplot2_all{4},'Marker','+');
            set(figs.mslplot2_all{5},'Marker','square');
            set(figs.mslplot2_all{6},'Marker','pentagram');
            set(figs.mslplot2_all{6},'LineWidth',2);
            
            xlabel('Year')
            ylabel('Transport (Sv)')
            title([tmp.regionname, ', Transport(',num2str(min(RCM_info.years),'%04i'),'-',num2str(max(RCM_info.years),'%04i'),') '])
%             datetick('x','yyyymm','keepticks')
            datetick('x','yyyy')

            axis tight;
            ylim(figs.meanplotlev2)
%             set(mslplot_all{5},'LineWidth',2);
%             set(figs.mslplot2_all{5},'LineWidth',2);
%             set(mslplot,'LineWidth',2);
            set(gca,'FontSize',25);
            grid on
            figs.lgd=legend({'RCM-CNRM'; 'RCM-EC-Veg'; 'RCM-ACC'; 'RCM-CNRM-HR'; 'RCM-CMCC'; 'ADCP'}, 'NumColumns',3);
            set(figs.lgd,'FontSize',15);
%             set(figs.lgd,'Position',[0.13 0.88, 0.775, 0.03]);  top center

            set(figs.lgd,'Orientation','horizontal');
            set(figs.lgd,'location','south'); 

%             constant_cor=corrcoef(RCM_interped_sla_yearly_mean(length(RCM_testnames),:),cmems_sla_yearly_mean);
%             txt1=text(xData2(5), min(figs.meanplotlev2)+diff(figs.meanplotlev2)/16.0 ,['R(RCM) = ', num2str(round(constant_cor(1,2),2))], 'FontSize', 20); 
%             constant_cor=corrcoef(GCM_interped_sla_yearly_mean(length(RCM_testnames),:),cmems_sla_yearly_mean);
%             txt1=text(xData2(10), min(figs.meanplotlev2)+diff(figs.meanplotlev2)/16.0 ,['R(GCM) = ', num2str(round(constant_cor(1,2),2))], 'FontSize', 20); 

%             for tind=1:size(GCM_interped_sla_yearly_mean,2)
%                 model_std(tind)=std(GCM_interped_sla_yearly_mean(:,tind));
%             end
%             meanstd=mean(model_std);


            set(gcf,'PaperPosition', [0 0 26 20])   
            hold off
            saveas(gcf,jpgname,'jpg'); RemoveWhiteSpace([], 'file', jpgname);
            save(['D:\Data\Model\ROMS\nwp_1_20\transport_all\', RCM_info.scenname, '_RCM_all.mat'], 'xData', 'xData_yearly', 'RCM_info', 'RCM_data', 'OBS_data');
            grid off
%         end
        close all;

% end-------------------- yearly tr time series (RCM)



% start-------------------- yearly soya tr time series(RCM)
      jpgname=strcat(tmp.outfile, '_', 'RCM', '_yearly_soya_tr_', RCM_info.scenname, '_', ...
          num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), '.jpg'); %% ~_year_month.jpg
            hold on
            for testind=1:length(RCM_info.name)
                figs.mslplot2_all{testind}=plot(xData_yearly,RCM_data.tr_yearly_soya{testind},'r')
            end
%             figs.mslplot2_all{testind+1}=plot(xData_yearly,OBS_data.OBS_data.ADCP_yearly, 'k');

            set(figs.mslplot2_all{1},'Marker','*');
            set(figs.mslplot2_all{2},'Marker','^');
            set(figs.mslplot2_all{3},'Marker','o');
            set(figs.mslplot2_all{4},'Marker','+');
            set(figs.mslplot2_all{5},'Marker','square');
%             set(figs.mslplot2_all{6},'Marker','pentagram');
%             set(figs.mslplot2_all{6},'LineWidth',2);
            
            xlabel('Year')
            ylabel('Transport (Sv)')
            title([tmp.regionname, ', Transport(',num2str(min(RCM_info.years),'%04i'),'-',num2str(max(RCM_info.years),'%04i'),') '])
%             datetick('x','yyyymm','keepticks')
            datetick('x','yyyy')

            axis tight;
            ylim([0 1.5])
%             set(mslplot_all{5},'LineWidth',2);
%             set(figs.mslplot2_all{5},'LineWidth',2);
%             set(mslplot,'LineWidth',2);
            set(gca,'FontSize',25);
            grid on
            figs.lgd=legend({'RCM-CNE'; 'RCM-ECV'; 'RCM-ACC'; 'RCM-CNH'; 'RCM-CMC'}, 'NumColumns',3);
            set(figs.lgd,'FontSize',15);
%             set(figs.lgd,'Position',[0.13 0.88, 0.775, 0.03]);  top center

            set(figs.lgd,'Orientation','horizontal');
            set(figs.lgd,'location','north'); 

%             constant_cor=corrcoef(RCM_interped_sla_yearly_mean(length(RCM_testnames),:),cmems_sla_yearly_mean);
%             txt1=text(xData2(5), min(figs.meanplotlev2)+diff(figs.meanplotlev2)/16.0 ,['R(RCM) = ', num2str(round(constant_cor(1,2),2))], 'FontSize', 20); 
%             constant_cor=corrcoef(GCM_interped_sla_yearly_mean(length(RCM_testnames),:),cmems_sla_yearly_mean);
%             txt1=text(xData2(10), min(figs.meanplotlev2)+diff(figs.meanplotlev2)/16.0 ,['R(GCM) = ', num2str(round(constant_cor(1,2),2))], 'FontSize', 20); 

%             for tind=1:size(GCM_interped_sla_yearly_mean,2)
%                 model_std(tind)=std(GCM_interped_sla_yearly_mean(:,tind));
%             end
%             meanstd=mean(model_std);


            set(gcf,'PaperPosition', [0 0 26 20])   
            hold off
            saveas(gcf,jpgname,'jpg'); RemoveWhiteSpace([], 'file', jpgname);
            save(['D:\Data\Model\ROMS\nwp_1_20\transport_all\', RCM_info.scenname, '_RCM_all.mat'], 'xData', 'xData_yearly', 'RCM_info', 'RCM_data', 'OBS_data');
            grid off
%         end
        close all;

% end-------------------- yearly tr time series (RCM)