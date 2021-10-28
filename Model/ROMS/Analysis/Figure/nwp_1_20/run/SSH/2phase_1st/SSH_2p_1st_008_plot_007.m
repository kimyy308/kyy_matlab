close all; clear all;  clc;
% %  Updated 05-Jul-2021 by Yong-Yub Kim, structure

% % % configuration of RCM
% RCM_info_all.name={'test2102', 'test2103', 'test2104', 'test2105', 'test2106'};
RCM_info_all.name={'test57', 'test58', 'test59', 'test60', ...
    'test61', 'test62', 'test63', 'test64', ...
    'test65', 'test66', 'test67', 'test68'};
% RCM_info_all.name={'test2102', 'test2103', 'test2104'};
RCM_info_all.abbs = {'RCM-IPSL-L', 'RCM-IPSL-M', 'RCM-Nor', 'RCM-MPI', ...
    'RCM-IPSL-L', 'RCM-IPSL-M', 'RCM-Nor', 'RCM-MPI', ...
    'RCM-IPSL-L', 'RCM-IPSL-M', 'RCM-Nor', 'RCM-MPI'};
% RCM_info_all.name={  'test2107'};
% RCM_info_all.abbs = {  'RCM-CNRM'};

RCM_info_all.model = 'nwp_1_20';
RCM_info_all.dataroot = ['G:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, 'nwp_1_20', filesep];
RCM_info_all.saveroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, 'nwp_1_20', filesep];
RCM_info_all.phase = 'run';
% RCM_info_all.region = {'NWP', 'AKP4'};
RCM_info_all.region = {'YS_KHOA'};
% RCM_info_all.years = 1985:2014;
RCM_info_all.years = 2006:2100;
RCM_info_all.months = 1:12;
RCM_grid.dl = 1/20;

% % % configuration of flags (make file)
flags.make_std_mat = 1;
flags.make_std_nc = 1;


% %  working
for testnameind=1:length(RCM_info_all.name)
    for regionind=1:length(RCM_info_all.region)
        clearvars '*' -except regionind testnameind flags flg_lev param ...
            RCM_info_all RCM_grid GCM_info GCM_grid CMEMS_info CMEMS_grid ...
            RCM_mean_data RCM_mean_trend GCM_mean_data GCM_mean_trend RCM_info_all GCM_003_info ...
            RCM_info_all RCM_data_all
        
        tmp.variable = 'temp';
        tmp.variable_GCM = 'thetao';
        tmp.fs = filesep; % file separator win = '\', linux = '/'

        % %     set dropbox path
        tmp.dropboxpath = 'C:\Users\User\Dropbox';
        addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));

% %     set temporary variables (testname, regionname, filesep, ...)
        RCM_info_all.testname = RCM_info_all.name{testnameind};
        RCM_info_all.regionname = RCM_info_all.region{regionind};
        RCM_info_all.abb = RCM_info_all.abbs{testnameind};
        [RCM_info_all.scenario, tmp.error_status] =  Func_0003_RCM_CMIP5_scenname(RCM_info_all.testname);
        
%         GCM_info.testname = GCM_003_info.name{testnameind};
%         GCM_info.regionname = RCM_info_all.regionname;
%         GCM_info.abb = GCM_info.abbs{testnameind};
%         GCM_info.scenario = RCM_info_all.scenario;
                
        [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);
        addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Model' ...
            tmp.fs, 'ROMS', tmp.fs, 'Analysis', tmp.fs, 'Figure', tmp.fs, 'nwp_1_20', tmp.fs ...
            'run', tmp.fs, 'SSH', tmp.fs, '2phase_1st', tmp.fs, 'subroutine']));
        
        [RCM_grid.refpolygon, RCM_grid.domain, tmp.error_status] = Func_0007_get_polygon_data_from_regionname(RCM_info_all.regionname);

        tmp.param_script = [tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', ...
            tmp.fs, 'Model', tmp.fs, 'ROMS', tmp.fs, 'Analysis', tmp.fs, 'Figure', ...
            tmp.fs, 'nwp_1_20', tmp.fs, 'run', tmp.fs, 'fig_param', tmp.fs, ...
            'fig_param2_kyy_', RCM_info_all.regionname, '.m'];
        RCM_info_all.filedir = [RCM_info_all.dataroot, RCM_info_all.testname, tmp.fs, 'run', tmp.fs];
        RCM_info_all.savedir = [RCM_info_all.saveroot, RCM_info_all.testname, tmp.fs, 'run', tmp.fs];
%         GCM_info.filedir = [GCM_info.dataroot, tmp.variable_GCM, tmp.fs, GCM_info.scenario, tmp.fs, ...
%             'Omon', tmp.fs, GCM_info.testname, tmp.fs];
        
        RCM_info_all.matname=[RCM_info_all.savedir,RCM_info_all.testname,'_',RCM_info_all.regionname, '_RCM_YSBCW_', ...
        num2str(min(RCM_info_all.years),'%04i'),'_',num2str(max(RCM_info_all.years),'%04i'),'.mat'];

%         RCM_info_all.matname_mean_data_trend = [RCM_info_all.savedir,'ENS3','_',RCM_info_all.regionname, '_RCM_ssh_mean_data_trend_', ...
%             num2str(min(RCM_info_all.years),'%04i'),'_',num2str(max(RCM_info_all.years),'%04i'),'.mat'];
%         GCM_003_info.matname_mean_data_trend = [RCM_info_all.savedir,'ENS3','_',RCM_info_all.regionname, '_GCM_ssh_mean_data_trend_', ...
%             num2str(min(RCM_info_all.years),'%04i'),'_',num2str(max(RCM_info_all.years),'%04i'),'.mat'];
        
% % %         time set
        for folding=1:1
            tmp.tind=1;
            for yearij = 1:length(RCM_info_all.years)
                for month=1:length(RCM_info_all.months) 
                    tmp.year = RCM_info_all.years(yearij);
%                     tmp.month = RCM_info_all.months(monthij);
                    RCM_time.ftime(tmp.tind) = datenum(tmp.year,month,15) - datenum(1900,12,31);
                    tmp.tind=tmp.tind+1;
                end
            end
            for month=1:length(RCM_info_all.months)  
                tmp.year = RCM_info_all.years(yearij);
                RCM_time.climtime(month) = datenum(1950,month,15) - datenum(1900,12,31);
            end

            for i =1:length(RCM_info_all.years) 
                tmp.year = RCM_info_all.years(yearij);
                for month=1:length(RCM_info_all.months)  
                    RCM_time.xData((12*(i-1))+month) = datenum([num2str(tmp.year),'-',num2str(month,'%02i'),'-01',]); 
                end
            end

            RCM_time.trendtime=RCM_info_all.years(1):1/length(RCM_info_all.months) : RCM_info_all.years(end)+1-1/length(RCM_info_all.months) ;
            RCM_time.trendtime_yearly=RCM_info_all.years(1) : RCM_info_all.years(end);
        end     
        
        RCM_info_all.matname = [RCM_info_all.savedir,RCM_info_all.testname,'_',RCM_info_all.regionname, '_RCM_YSBCW_', ...
        num2str(min(RCM_info_all.years),'%04i'),'_',num2str(max(RCM_info_all.years),'%04i'),'.mat'];
%         GCM_info.matname = [RCM_info_all.savedir,RCM_info_all.testname,'_',RCM_info_all.regionname, '_GCM_ssh_', ...
%             num2str(min(RCM_info_all.years),'%04i'),'_',num2str(max(RCM_info_all.years),'%04i'),'.mat'];

        if (exist(RCM_info_all.matname , 'file') == 2)   
            load(RCM_info_all.matname); 
%             load(GCM_info.matname);
%             load(GCM_info.matname_interped);
        end
        RCM_data_all.(RCM_info_all.testname)=RCM_data;
    end
end

% for testnameind=9:length(RCM_info_all.name)
for testnameind=9:12
    RCM_info_all.testname = RCM_info_all.name{testnameind};
    plot(RCM_info_all.years, RCM_data_all.(RCM_info_all.testname).all_volume)
    hold on
end
hold off
lgd=legend({'RCM-IPSL-L', 'RCM-IPSL-M', 'RCM-Nor', 'RCM-MPI'});


for testnameind=1:length(RCM_info_all.name)
    RCM_info_all.testname = RCM_info_all.name{testnameind};
    plot(RCM_info_all.years, RCM_data_all.(RCM_info_all.testname).all_southern_limit)
    hold on
end
hold off

RCM_info_all.scenarios = {'rcp26', 'rcp45', 'rcp85'};
RCM_info_all.rcp26_names={'test61', 'test62', 'test63', 'test64'};
RCM_info_all.rcp45_names={'test57', 'test58', 'test59', 'test60'};
RCM_info_all.rcp85_names={'test65', 'test66', 'test67', 'test68'};

for scenind=1:length(RCM_info_all.scenarios)
    RCM_info_all.scenario = RCM_info_all.scenarios{scenind};
    for testnameind=1:length(RCM_info_all.([RCM_info_all.scenario,'_names']))
        RCM_info_all.testname = RCM_info_all.([RCM_info_all.scenario,'_names']){testnameind};
        if testnameind==1
            RCM_data_all.(['ens_', RCM_info_all.scenario]).all_volume = ...
                RCM_data_all.(RCM_info_all.testname).all_volume / length(RCM_info_all.([RCM_info_all.scenario,'_names']));
            RCM_data_all.(['ens_', RCM_info_all.scenario]).all_southern_limit = ...
                RCM_data_all.(RCM_info_all.testname).all_southern_limit / length(RCM_info_all.([RCM_info_all.scenario,'_names']));
        else
            RCM_data_all.(['ens_', RCM_info_all.scenario]).all_volume = RCM_data_all.(['ens_', RCM_info_all.scenario]).all_volume ...
                + RCM_data_all.(RCM_info_all.testname).all_volume / length(RCM_info_all.([RCM_info_all.scenario,'_names']));
            RCM_data_all.(['ens_', RCM_info_all.scenario]).all_southern_limit = RCM_data_all.(['ens_', RCM_info_all.scenario]).all_southern_limit ...
                + RCM_data_all.(RCM_info_all.testname).all_southern_limit / length(RCM_info_all.([RCM_info_all.scenario,'_names']));
        end
    end
end

val_transparent = 0.3;
RCM_info_all.cmap_rcp26 = [0,0,1]; % blue
RCM_info_all.cmap_rcp26_b = rgb2hsv(RCM_info_all.cmap_rcp26);
RCM_info_all.cmap_rcp26_b(:,2) =  val_transparent;
RCM_info_all.cmap_rcp26_b = hsv2rgb(RCM_info_all.cmap_rcp26_b);

RCM_info_all.cmap_rcp45 = [1, 165/255, 0];
RCM_info_all.cmap_rcp45_b = rgb2hsv(RCM_info_all.cmap_rcp45);
RCM_info_all.cmap_rcp45_b(:,2) =  val_transparent;
RCM_info_all.cmap_rcp45_b = hsv2rgb(RCM_info_all.cmap_rcp45_b);

RCM_info_all.cmap_rcp85 = [1, 0, 0];
RCM_info_all.cmap_rcp85_b = rgb2hsv(RCM_info_all.cmap_rcp85);
RCM_info_all.cmap_rcp85_b(:,2) =  val_transparent;
RCM_info_all.cmap_rcp85_b = hsv2rgb(RCM_info_all.cmap_rcp85_b);

for scenind=1:length(RCM_info_all.scenarios)
    RCM_info_all.scenario = RCM_info_all.scenarios{scenind};    
    plot(RCM_info_all.years, RCM_data_all.(['ens_', RCM_info_all.scenario]).all_volume, ...
        'Linewidth', 2, 'color', RCM_info_all.(['cmap_', RCM_info_all.scenario]))
    hold on
end
hold off
lgd=legend({'RCP 2.6', 'RCP 4.5', 'RCP 8.5'});
set(lgd,'FontSize',15);
set(gca, 'FontSize', 15);
ylabel('Volume (m^3)');
xlabel('Years');
xlim([2006 2100])
set(gcf,'PaperPosition', [0 0 22 17]) 
hold off
tifname= ['Z:\내 드라이브\MEPL\project\SSH\6th_year\figure\nwp_1_20\all\ts\CMIP5_RCM_YSBCW_ts_vol_', ...
    num2str(min(RCM_info_all.years),'%04i'),'_',num2str(max(RCM_info_all.years),'%04i'), '.tif'];
saveas(gcf,tifname,'tif'); RemoveWhiteSpace([], 'file', tifname);


for scenind=1:length(RCM_info_all.scenarios)
    RCM_info_all.scenario = RCM_info_all.scenarios{scenind};    
    plot(RCM_info_all.years, RCM_data_all.(['ens_', RCM_info_all.scenario]).all_southern_limit, ...
        'Linewidth', 2, 'color', RCM_info_all.(['cmap_', RCM_info_all.scenario]))
    hold on
end
hold off
lgd=legend({'RCP 2.6', 'RCP 4.5', 'RCP 8.5'});
set(lgd,'FontSize',15);
set(gca, 'FontSize', 15);
ylabel('Latitude (^oN)');
xlabel('Years');
xlim([2006 2100])
set(gcf,'PaperPosition', [0 0 22 17]) 
hold off
tifname= ['Z:\내 드라이브\MEPL\project\SSH\6th_year\figure\nwp_1_20\all\ts\CMIP5_RCM_YSBCW_ts_south_lim_', ...
    num2str(min(RCM_info_all.years),'%04i'),'_',num2str(max(RCM_info_all.years),'%04i'), '.tif'];
saveas(gcf,tifname,'tif'); RemoveWhiteSpace([], 'file', tifname);
