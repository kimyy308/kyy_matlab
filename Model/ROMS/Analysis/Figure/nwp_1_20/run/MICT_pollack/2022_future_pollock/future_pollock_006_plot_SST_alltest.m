% %  Updated 27-Apr-2021 by Yong-Yub Kim, 

close all; clear all;  clc;
warning off;

tmp.testname_reana='test06';
RCM_info.model_reana = 'nwp_1_10';

% RCM_info.name={ 'test2127'};
% RCM_info.name={'test2117', 'test2118', 'test2119', 'test2120', 'test2121'};
RCM_info.name={'test2127', 'test2128', 'test2129', 'test2130', 'test2131'};
% 
RCM_info.model = 'nwp_1_20';

RCM_info.dataroot_reana = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, RCM_info.model_reana, filesep, 'backup_surf', filesep];
RCM_info.dataroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, RCM_info.model, filesep, 'backup_surf', filesep];
RCM_info.saveroot_reana = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, RCM_info.model_reana, filesep, 'backup_surf', filesep];
RCM_info.saveroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, RCM_info.model, filesep, 'backup_surf', filesep];
RCM_info.phase = 'run';  % run or spinup
RCM_info.region = {'EKB2'}; % NWP, AKP4, ES_KHOA, YS, ...
% RCM_info.region = {'pollock_egg3'}; % NWP, AKP4, ES_KHOA, YS, ...

RCM_info.vars = {'SST'};

% RCM_info.years = 1983:2021;  
% RCM_info.years = [2015:2050, 2081:2100];  
% RCM_info.years = [1983:1987];  
% RCM_info.years = [1988:1992];  
RCM_info.years = [1995:2014];  
% RCM_info.years = [1993:2021];  
% RCM_info.years = [2081:2100];  
% RCM_info.years = [2015:2100];  

% seasons_group={'February', 'January', 'JF-'};
seasons_group={'JF-'};

% RCM_grid.dl = 1/20;
for seasons_groupi=1:length(seasons_group)
    for testnameind2=1:length(RCM_info.name)
        for regionind2=1:length(RCM_info.region)
            for varind2=1:length(RCM_info.vars)

                tmp.fs=filesep;  
                addpath(genpath('C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\MICT_pollack\2022_future_pollock\subroutine\'))

                tmp.dropboxpath = 'C:\Users\User\Dropbox';
                addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
                [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);
                addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Model' ...
                    tmp.fs, 'ROMS', tmp.fs, 'Analysis', tmp.fs, 'Figure', tmp.fs, 'nwp_1_20', tmp.fs ...
                    'run', tmp.fs, 'SSH', tmp.fs, '2phase_2nd', tmp.fs, 'subroutine']));

                RCM_info.season=seasons_group{seasons_groupi};
                [RCM_info.months,error_status] = Func_0019_get_month_from_season(RCM_info.season);

                tmp.testname=RCM_info.name{testnameind2};   % % need to change
                [tmp.abb, tmp.error_status] = Func_0018_abbname(tmp.testname);
                tmp.regionname=RCM_info.region{regionind2};
                [RCM_grid.refpolygon, RCM_grid.domain, tmp.error_status] = Func_0007_get_polygon_data_from_regionname(tmp.regionname);

    %             tmp.variable ='temp';
                tmp.variable=RCM_info.vars{varind2};

                dirs.figrawdir =strcat('D:\Research\Ph_D_course\2022_pollock_future\figure',filesep, 'all', filesep); % % where figure files will be saved            
                tmp.param_script =['C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\', RCM_info.model, '\run\fig_param\fig_param2_kyy_', tmp.regionname, '.m'];
                dirs.filedir = strcat('D:\Data\Model\ROMS\', RCM_info.model, '\backup_surf\', tmp.testname, '\run\'); % % where data files are          
                dirs.matdir = strcat('D:\Data\Model\ROMS\', RCM_info.model, '\', tmp.testname, '\run\mean\');
                dirs.matdir_reana = strcat('D:\Data\Model\ROMS\', RCM_info.model_reana, '\', tmp.testname_reana, '\run\mean\');
                dirs.griddir = strcat('D:\Data\Model\ROMS\', RCM_info.model, '\backup_surf\'); % % where grid data files are            

                tmp.matname = [dirs.matdir, tmp.testname, '_', tmp.regionname, '_', tmp.variable,...
                        '_mean_', num2str(min(RCM_info.years),'%04i'), '-', num2str(max(RCM_info.years),'%04i'),...
                    '_', RCM_info.season, '.mat'];
                load(tmp.matname, 'RCM_data');
                RCM_all_data.(tmp.testname)=RCM_data;
                
                if max(RCM_info.years)<=2021
                    tmp.matname_reana = [dirs.matdir_reana, tmp.testname_reana, '_', tmp.regionname, '_', tmp.variable,...
                            '_mean_', num2str(min(RCM_info.years),'%04i'), '-', num2str(max(RCM_info.years),'%04i'),...
                        '_', RCM_info.season, '.mat'];
                    load(tmp.matname_reana, 'RCM_data');
                    RCM_all_data.reana=RCM_data;
                end
                RCM_all_data.abb{testnameind2}=tmp.abb;
            end
        end
    end
end
if max(RCM_info.years)<=2021
    RCM_all_data.abb{end+1}='reanalysis';
end
% start-------------------- yearly tr time series(RCM)
dirs.figdir=[dirs.figrawdir,'ts', tmp.fs, tmp.regionname, tmp.fs, tmp.variable, tmp.fs, ...
        num2str(min(RCM_info.years)), '_', num2str(max(RCM_info.years)), tmp.fs];
    if (exist(strcat(dirs.figdir) , 'dir') ~= 7)
        mkdir(strcat(dirs.figdir));
    end 

%% set figure file name
    tmp.tifname=strcat(dirs.figdir, 'all', '_yearly_ts_', tmp.variable, '_', ...
        num2str(min(RCM_info.years),'%04i'), '_',num2str(max(RCM_info.years),'%04i'), ...
        '_', RCM_info.season,'.tif'); %% ~_year_month.jpg
    
            hold on
            RCM_all_data.ens.spa_yearly_mean=zeros(1, length(RCM_all_data.(RCM_info.name{1}).spa_yearly_mean));
            %% each test
            for testind=1:length(RCM_info.name)
                mslplot2_all{testind}=plot(RCM_info.years,RCM_all_data.(RCM_info.name{testind}).spa_yearly_mean,'r');
                RCM_all_data.ens.spa_yearly_mean=RCM_all_data.ens.spa_yearly_mean + RCM_all_data.(RCM_info.name{testind}).spa_yearly_mean / length(RCM_info.name);
            end
            %% reanalysis
            if max(RCM_info.years)<=2021
                mslplot2_all{testind+1}=plot(RCM_info.years,RCM_all_data.('reana').spa_yearly_mean, 'k');
                set(mslplot2_all{end},'LineWidth',2);
            end
            set(mslplot2_all{1},'Marker','*');
            set(mslplot2_all{2},'Marker','^');
            set(mslplot2_all{3},'Marker','o');
            set(mslplot2_all{4},'Marker','+');
            set(mslplot2_all{5},'Marker','square');
%             set(mslplot2_all{6},'Marker','pentagram');
            %% ens
            mslplot2_all{testind+2}=plot(RCM_info.years,RCM_all_data.ens.spa_yearly_mean,'r');
            set(mslplot2_all{end},'LineWidth',3);
            RCM_all_data.abb{end+1}='RCM-ENS';

            
            xlabel('Year')
            ylabel('Temperature (^oC)')
%             title([tmp.regionname, ', temp(',num2str(min(RCM_info.years),'%04i'),'-',num2str(max(RCM_info.years),'%04i'),') '])
%             datetick('x','yyyy')

            axis tight;
%             ylim([-2 13]) %EKB2, 95-14
%             ylim([2 12]) %pollock_egg3, 95-14
%             ylim([2 18]) %EKB2, 2081-2100
            ylim([-3 18]) %EKB2, 2015-2100

%             set(mslplot_all{5},'LineWidth',2);
%             set(mslplot2_all{5},'LineWidth',2);
%             set(mslplot,'LineWidth',2);
            set(gca,'FontSize',25);
            grid on
%             lgd=legend({'RCM-CNRM'; 'RCM-EC-Veg'; 'RCM-ACC'; 'RCM-CNRM-HR'; 'RCM-CMCC'; 'ADCP'}, 'NumColumns',3);
            lgd=legend(RCM_all_data.abb, 'NumColumns',3);

            set(lgd,'FontSize',15);
%             set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);  top center

            set(lgd,'Orientation','horizontal');
            set(lgd,'location','south'); 

%             constant_cor=corrcoef(RCM_interped_sla_yearly_mean(length(RCM_testnames),:),cmems_sla_yearly_mean);
%             txt1=text(xData2(5), min(meanplotlev2)+diff(meanplotlev2)/16.0 ,['R(RCM) = ', num2str(round(constant_cor(1,2),2))], 'FontSize', 20); 
%             constant_cor=corrcoef(GCM_interped_sla_yearly_mean(length(RCM_testnames),:),cmems_sla_yearly_mean);
%             txt1=text(xData2(10), min(meanplotlev2)+diff(meanplotlev2)/16.0 ,['R(GCM) = ', num2str(round(constant_cor(1,2),2))], 'FontSize', 20); 

%             for tind=1:size(GCM_interped_sla_yearly_mean,2)
%                 model_std(tind)=std(GCM_interped_sla_yearly_mean(:,tind));
%             end
%             meanstd=mean(model_std);
% 
%             set(gcf,'PaperPosition', [0 0 48 12])   

            

            set(gcf,'PaperPosition', [0 0 24 12])   
            hold off
            saveas(gcf,tmp.tifname,'jpg'); RemoveWhiteSpace([], 'file', tmp.tifname);
%             save(['D:\Data\Model\ROMS\nwp_1_20\transport_all\', scenname, '_RCM_all.mat'], 'xData_yearly', 'trdata2_yearly');
            grid off
%         end
        close all;
% end-------------------- yearly tr time series (RCM)