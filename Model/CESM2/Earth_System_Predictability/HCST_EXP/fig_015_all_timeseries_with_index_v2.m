% %  Created 12-Apr-2023 by Yong-Yub Kim
clc; clear all; close all;
warning off;

%% set path
[error_status, tmp.hostname] = system('hostname');
tmp.hostname=tmp.hostname(1:end-1);
switch tmp.hostname
    case 'Yong-Yubs-iMac-Pro.local'
        tmp.dropboxpath = '/Volumes/kyy_raid/kimyy/Dropbox';
        tmp.kimyypath = '/Volumes/kyy_raid/kimyy';        
    case {'da1', 'da2', 'da3', 'da4'}
        tmp.dropboxpath = '/mnt/lustre/proj/kimyy/Dropbox';
end

tmp.fs=filesep;
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);

%% model configuration

% dirs.hcstroot=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/HCST_EXP/archive/atm/', cfg.var];
% dirs.obsroot=['/Volumes/kyy_raid/kimyy/Observation/GPCP/monthly_reg_cam'];
% dirs.figroot=['/Volumes/kyy_raid/kimyy/Figure/CESM2/ESP/HCST_EXP/archive/atm/', cfg.var];


% cfg.len_t_y = length(cfg.iyears);
cfg_ts.vars={'SSH', 'photoC_TOT_zint_100m', 'NO3'};
cfg_ts.vars={'GPP', 'PRECT', 'PSL', 'TREFHT', 'TWS'};

% cfg_ts.vars={'SSH', 'photoC_TOT_zint_100m'};
% cfg_ts.vars={'NO3'};
% cfg_ts.vars={'SST'};
% cfg_ts.vars={'DIC_ALT_CO2', 'DIC', 'TEMP', 'SALT'};
% cfg_ts.vars={'SHF'};
% cfg_ts.vars={'SALT'};
% cfg_ts.vars={'HMXL', 'HBLT'};
% cfg_ts.vars={'DIC_ALT_CO2'};
% cfg_ts.vars={'TEMP'};
% cfg_ts.vars={'SALT'};
% cfg_ts.vars={'U', 'V'};
% cfg_ts.var='NO3';

cfg_ts.vars={'SSH'};
% cfg_ts.vars={'photoC_TOT_zint_100m', 'NO3'};
cfg_ts.vars={'SST'};
cfg_ts.vars={'photoC_TOT_zint_100m'};
cfg_ts.vars={'TREFHT'};


% grid.regions=[45 80 -30 -20];
% grid.regions=[135 200 10 20];
grid.regions=[240 290 -40 -20];
grid.regions=[210 240 30 60];
grid.regions=[160 160 20 20];

% sta_lonlat = {[160 160 20 20]};
sta_lonlat = {[300 300 25 25]};
sta_lonlat = {[60 60 -25 -25]};
sta_lonlat = {[10 10 -10 -10], [120 120 -60 -60], [160 160 -10 -10], ...
    [160 160 15 15], [160 160 20 20], [180 180 80 80], [200 200 -30 -30], ...
    [270 270 -30 -30], [300 300 40 40], [330 330 -30 -30], [340 340 25 25], ...
    [345 345 -10 -10], [60 60 -25 -25], [90 90 10 10], ...
    [135 135 15 15], [165 165 40 40], [180 180 20 20], [320 320 -25 -25], ...
    [340 340 -10 -10], [60 60 -30 -30], [330 330 -25 -25]};

% sta_lonlat = {[135 135 15 15], [165 165 40 40], [180 180 20 20], [320 320 -25 -25], ...
%     [340 340 -10 -10], [60 60 -30 -30]};

sta_lonlat = {[330 330 -25 -25]};
sta_lonlat = {[10 10 10 10], [120 120 13 13], [15 15 0 0], [260 260 30 30], [30 30 -26 -26], [322 322 -8 -8], ...
    [60 60 30 30], [80 80 25 25], [90 90 33 33]}; % for land variables

sta_lonlat = {[140 140 5 5], [180 180 5 5], [200 200 -15 -15], [340 340 60 60], [350 350 60 60], [280 280 25 25], ...
    [300 300 25 25], [350 350 -5 -5], [5 5 3 3], [50 50 0 0], [55 55 -20 -20], [90 90 0 0]};


sta_lonlat = {[140 155 25 35], [142 155 30 40], [150 150 20 50], [140 140 25 25], [145 145 25 25], ...
    [150 150 25 25], [155 155 25 25], [140 140 30 30], [150 150 30 30], [155 155 30 30], ...
    [140 140 35 35], [145 145 35 35], [150 150 35 35], [155 155 35 35]};
sta_lonlat ={[142 155 30 40], [140 155 25 35]}; 
% sta_lonlat ={[140 155 25 35]};
sta_lonlat ={[114 164 -10 19], [340 350 50 60], ...
    [264 280 20 30], [264 301 9 20], [160 180 -60 -33], ...
    [42 52 -26 -10], [1 9 0 6], [67 95 -13 23]}; %% Alexia
% sta_lonlat ={[114 164 -10 19]}; 
sta_lonlat ={[50 70 -35 -20], [160 190 35 50], [120 200 10 20], [150 210 -30 -10], [170 240 -60 -40], ...
    [240 300 -50 -20], [315 335 -30 -10], [300 360 -50 -40], [280 310 35 45], ...
    [280 360 25 45], [315 345 15 25], [70 140 -60 -50], [180 300 70 90], [30 180 70 90]}; 
sta_lonlat = {[0 25 35 45], [110 120 -40 -12], [122 133 23 30], [199 206 18 23], [268 272 -2 0]};
sta_lonlat = {[210 270 -20 10]};
% sta_lonlat = {[200 200 20 20], [210 270 -5 5], [230 250 30 40], [320 360 -20 0], [340 360 30 40], [90 120, -10 10], [90 120 -20 0], [90 120 -40 0]};
sta_lonlat = {[275 324 -18 17], [260 281 26 33], ...
    [197 229 58 66], [297 339 61 81], [25 55 50 57], ...
    [342 9 -6 12], [30 59 14 44], [70 89 6 32], ...
    [97 152 -9 18], [113 179 -46 -12]};

%% 145 30 should be rechecked

cfg.vlayer=1; % surface, vertical slice
% cfg.vlayer=24; % surface, vertical slice
% cfg.vlayer=32; % atm surface

% cfg.vlayer=1:10; % 10layer. don't put more than 15
cfg.vlayer_1st=min(cfg.vlayer);
cfg.vlayer_cnt=max(cfg.vlayer)-cfg.vlayer_1st+1;

vstr=['v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i')];

cfg.iyears=1960:2020;



% monthss={[12:14], [3:5], [6:8], [9:11]};
monthss={[1:12]};

for mmmi=1:length(monthss)
    months=monthss{mmmi};
    % months = 1:12;
    % months = 1:3;
    % months = 12:14;
    % months = 3:5;
    % months = 6:8;
    % months = 9:11;
    
    month1=months(months<=12);
    month2=months(months>12);
    
    if length(months)==12
        str_prepos='annual';
    else
        str_prepos='seasonal';
    end
    
    once=1;
    
    for vari=1:length(cfg_ts.vars)
        cfg_ts.var=cfg_ts.vars{vari};
        cfg.comp=Func_0025_CESM2_cmpname_var(cfg_ts.var);
        dirs.figroot=[tmp.kimyypath, '/Figure/CESM2/ESP/HCST_EXP/archive/', cfg.comp,'/', cfg_ts.var, filesep, vstr];
    for stai=1:length(sta_lonlat)
        grid.regions=sta_lonlat{stai};
    
        cfg_ts.obs=f_obs_name(cfg_ts.var);
        
    %     str_regions = ['_lon',num2str(grid.regions(1), '%02i'), '_', num2str(grid.regions(2), '%02i'), ...
    %         '__lat',num2str(grid.regions(3),'%02i'),'_',num2str(grid.regions(4),'%02i')];
        str_regions = ['_lon',num2str(grid.regions(1)), '_', num2str(grid.regions(2)), ...
            '__lat',num2str(grid.regions(3)),'_',num2str(grid.regions(4))];
        
        %% ts_data read
        load(['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/ts_data/', ...
            'ts_data_all_', cfg_ts.var,str_regions,'_', vstr, '_obs_', cfg_ts.obs, '.mat'])
        
        if strcmp(cfg_ts.var, 'SST')
            data_obs.SST(data_obs.SST<-900)=NaN;
            data_obs.SST_ym(data_obs.SST_ym<-900)=NaN;
            m_data_obs.SST=Func_0011_get_area_weighted_mean( data_obs.SST, grid.cut_tlong,grid.cut_tlat);
            m_data_obs.SST_ym=Func_0011_get_area_weighted_mean( data_obs.SST_ym, grid.cut_tlong,grid.cut_tlat);

%             m_data_assm.SST=m_data_assm.SST -273.15;
%             m_data_assm.SST_ym=m_data_assm.SST_ym -273.15;
%             for li=1:5
%                 str_li=['ly', num2str(li)];
%                 m_data_hcst.SST.(str_li)=m_data_hcst.SST.(str_li)-273.15;
%                 m_data_hcst.SST_ym.(str_li)=m_data_hcst.SST_ym.(str_li) -273.15;
%             end
%             m_data_lens2.SST=m_data_lens2.SST -273.15;
%             m_data_lens2.SST_ym=m_data_lens2.SST_ym -273.15;
        end
    
    
        len_m=size(m_data_assm.(cfg_ts.var),2);
        len_mem_oda=cfg_assm.len_mem;
        len_mem_lens2=cfg_lens2.len_mem;
        
        %% without obs
        if length(months)==12
            %% annual graph
            for ly=1:5
                fig_h = figure('name','ts','visible','off');
                fig_h.Position= [0,0,1000,500];
        
                tmp.lyear_str=num2str(ly);
                hold on
                for mi=1:50
                    fig_ts.plots=plot(cfg.iyears, m_data_lens2.([cfg_ts.var, '_ym'])(mi,:), 'g');
                    set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                end
                plot(cfg.iyears, mean(m_data_lens2.([cfg_ts.var, '_ym']),1), 'g', 'linewidth', 3);
                
                for mi=1:20
                    fig_ts.plots=plot(cfg.iyears+ly-1, m_data_hcst.([cfg_ts.var, '_ym']).(['ly',tmp.lyear_str])(mi,:), 'b');
                    set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                    
                    fig_ts.plots=plot(cfg.iyears, m_data_assm.([cfg_ts.var, '_ym'])(mi,:), 'r');
                    set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');        
                end
                plot(cfg.iyears+ly-1, mean(m_data_hcst.([cfg_ts.var, '_ym']).(['ly',tmp.lyear_str]),1), 'b', 'linewidth', 3);
                
                plot(cfg.iyears, mean(m_data_assm.([cfg_ts.var, '_ym']),1), 'r', 'linewidth', 3);
                
                legend ('LE', 'ODA', 'HIND', 'Location', 'Southoutside', 'Orientation', 'Horizontal')
                title([cfg_ts.var, ', ', num2str(sta_lonlat{stai}(1)),'E, ', ...
                    num2str(sta_lonlat{stai}(2)),'E, ', ...
                    num2str(sta_lonlat{stai}(3)), 'N, ', ...
                    num2str(sta_lonlat{stai}(4)), 'N'])
            
                grid minor
                xlim([1960 2025])
                set(gca, 'fontsize', 20)    
            
                dirs.figdir= [dirs.figroot, filesep, 'ens_all', filesep, cfg_ts.var, '_time_series_indv', filesep, 'l',tmp.lyear_str];
                if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
                cfg.figname=[dirs.figdir, filesep, 'ts_all_l',tmp.lyear_str, '_', num2str(sta_lonlat{stai}(1)), 'E_', ...
                    num2str(sta_lonlat{stai}(2)), 'E_', num2str(sta_lonlat{stai}(3)), 'N_',num2str(sta_lonlat{stai}(4)), 'N_', cfg_ts.var, '.tif'];
                print(fig_h, cfg.figname, '-dpng');
                RemoveWhiteSpace([], 'file', cfg.figname);
                close all;
            end
        
            fig_h = figure('name','ts','visible','off');
            fig_h.Position= [0,0,1000,500];
        
            tmp.lyear_str=num2str(ly);
            hold on
            cfg.iyears_4ym=movmean(cfg.iyears, 4, 'Endpoints', 'discard');
            for mi=1:50
                fig_ts.plots=plot(cfg.iyears_4ym, movmean(m_data_lens2.([cfg_ts.var, '_ym'])(mi,:),4,'Endpoints', 'discard'), 'g');
                set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
            plot(cfg.iyears_4ym, movmean(mean(m_data_lens2.([cfg_ts.var, '_ym']),1),4,'Endpoints', 'discard'), 'g', 'linewidth', 3);
            for mi=1:20
        %         fig_ts.plots=plot(cfg.iyears+ly-1, m_data_hcst.([cfg_ts.var, '_ym']).(['ly',tmp.lyear_str])(mi,:), 'b');
        %         set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                
                tmp.data(mi,:)= m_data_hcst.([cfg_ts.var, '_ym']).ly2(mi,:) + ...
                    m_data_hcst.([cfg_ts.var, '_ym']).ly3(mi,:) + ...
                    m_data_hcst.([cfg_ts.var, '_ym']).ly4(mi,:) + ...
                    m_data_hcst.([cfg_ts.var, '_ym']).ly5(mi,:);
                tmp.data(mi,:)=tmp.data(mi,:)/4;
        
                fig_ts.plots=plot(cfg.iyears+2.5, tmp.data(mi,:), 'b');
                set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                
                fig_ts.plots=plot(cfg.iyears_4ym, movmean(m_data_assm.([cfg_ts.var, '_ym'])(mi,:),4,'Endpoints', 'discard'), 'r');
                set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');        
            end
        %     plot(cfg.iyears+ly-1, mean(m_data_hcst.([cfg_ts.var, '_ym']).(['ly',tmp.lyear_str]),1), 'b', 'linewidth', 3);
            plot(cfg.iyears+2.5, mean(tmp.data,1), 'b', 'linewidth', 3);
            plot(cfg.iyears_4ym, movmean(mean(m_data_assm.([cfg_ts.var, '_ym']),1),4,'Endpoints', 'discard'), 'r', 'linewidth', 3);
            
            legend ('LE', 'ODA', 'HIND', 'Location', 'Southoutside', 'Orientation', 'Horizontal')
        %     title([cfg_ts.var, ', ', num2str(sta_lonlat{stai}(1)),'E, ', num2str(sta_lonlat{stai}(3)), 'N'])
            title([cfg_ts.var, ', ', num2str(sta_lonlat{stai}(1)),'E, ', ...
                    num2str(sta_lonlat{stai}(2)),'E, ', ...
                    num2str(sta_lonlat{stai}(3)), 'N, ', ...
                    num2str(sta_lonlat{stai}(4)), 'N'])
        
            grid minor
            xlim([1960 2025])
            set(gca, 'fontsize', 20)    
        
            dirs.figdir= [dirs.figroot, filesep, 'ens_all', filesep, cfg_ts.var, '_time_series_indv', filesep, 'l2_5'];
            if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
            cfg.figname=[dirs.figdir, filesep, 'ts_all_l2_5', '_', num2str(sta_lonlat{stai}(1)), 'E_', ...
                num2str(sta_lonlat{stai}(2)), 'E_', num2str(sta_lonlat{stai}(3)), 'N_', num2str(sta_lonlat{stai}(4)), 'N_',cfg_ts.var, '.tif'];
            print(fig_h, cfg.figname, '-dpng');
            RemoveWhiteSpace([], 'file', cfg.figname);
            close all;
        elseif sum(months>12)==0 % 1~3, 2~4, ... 10~12
            for ly=1:5
                fig_h = figure('name','ts','visible','off');
                fig_h.Position= [0,0,1000,500];
        
                tmp.lyear_str=num2str(ly);
                cfg.iyears_4ym=movmean(cfg.iyears, 4, 'Endpoints', 'discard');
                hold on
                tmp.reshp=reshape(m_data_lens2.(cfg_ts.var), [len_mem_lens2, 12, len_m/12]); % get seasonal mean
                m_data_lens2.([cfg_ts.var, '_sm'])=squeeze(mean(tmp.reshp(:,months,:),2));
                for mi=1:50
                    fig_ts.plots=plot(cfg.iyears, m_data_lens2.([cfg_ts.var, '_sm'])(mi,:), 'g');
                    set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                end
                plot(cfg.iyears, mean(m_data_lens2.([cfg_ts.var, '_sm']),1), 'g', 'linewidth', 3);
                
                tmp.reshp=reshape(m_data_hcst.(cfg_ts.var).(['ly',tmp.lyear_str]), [len_mem_oda, 12, len_m/12]); % get seasonal mean
                m_data_hcst.([cfg_ts.var, '_sm']).(['ly',tmp.lyear_str])=squeeze(mean(tmp.reshp(:,months,:),2));
                tmp.reshp=reshape(m_data_assm.(cfg_ts.var), [len_mem_oda, 12, len_m/12]); % get seasonal mean
                m_data_assm.([cfg_ts.var, '_sm'])=squeeze(mean(tmp.reshp(:,months,:),2));
                for mi=1:20
                    fig_ts.plots=plot(cfg.iyears+ly-1, m_data_hcst.([cfg_ts.var, '_sm']).(['ly',tmp.lyear_str])(mi,:), 'b');
                    set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                    fig_ts.plots=plot(cfg.iyears, m_data_assm.([cfg_ts.var, '_sm'])(mi,:), 'r');
                    set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');        
                end
                plot(cfg.iyears+ly-1, mean(m_data_hcst.([cfg_ts.var, '_sm']).(['ly',tmp.lyear_str]),1), 'b', 'linewidth', 3);
                
                plot(cfg.iyears, mean(m_data_assm.([cfg_ts.var, '_sm']),1), 'r', 'linewidth', 3);
                
                legend ('LE', 'ODA', 'HIND', 'Location', 'Southoutside', 'Orientation', 'Horizontal')
                title([cfg_ts.var, ', sm, ', num2str(sta_lonlat{stai}(1)),'E, ', ...
                    num2str(sta_lonlat{stai}(2)),'E, ', ...
                    num2str(sta_lonlat{stai}(3)), 'N, ', ...
                    num2str(sta_lonlat{stai}(4)), 'N'])
            
                grid minor
                xlim([1960 2025])
                set(gca, 'fontsize', 20)    
            
                dirs.figdir= [dirs.figroot, filesep, 'ens_all', filesep, cfg_ts.var, '_time_series_indv', filesep, 'l',tmp.lyear_str, '_', str_prepos];
                if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
                cfg.figname=[dirs.figdir, filesep, str_prepos, '_', num2str(min(months)), '_', num2str(max(months)), '_ts_all_l',tmp.lyear_str, '_', num2str(sta_lonlat{stai}(1)), 'E_', ...
                    num2str(sta_lonlat{stai}(2)), 'E_', num2str(sta_lonlat{stai}(3)), 'N_',num2str(sta_lonlat{stai}(4)), 'N_', cfg_ts.var, '.tif'];
                print(fig_h, cfg.figname, '-dpng'); 
                RemoveWhiteSpace([], 'file', cfg.figname);
                close all;
            end
        
            fig_h = figure('name','ts','visible','off');
            fig_h.Position= [0,0,1000,500];
        
            tmp.lyear_str=num2str(ly);
            hold on
            cfg.iyears_4ym=movmean(cfg.iyears, 4, 'Endpoints', 'discard');
            for mi=1:50
                fig_ts.plots=plot(cfg.iyears_4ym, movmean(m_data_lens2.([cfg_ts.var, '_sm'])(mi,:),4,'Endpoints', 'discard'), 'g');
                set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
            plot(cfg.iyears_4ym, movmean(mean(m_data_lens2.([cfg_ts.var, '_sm']),1),4,'Endpoints', 'discard'), 'g', 'linewidth', 3);
            for mi=1:20
        %         fig_ts.plots=plot(cfg.iyears+ly-1, m_data_hcst.([cfg_ts.var, '_ym']).(['ly',tmp.lyear_str])(mi,:), 'b');
        %         set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                
                tmp.data(mi,:)= m_data_hcst.([cfg_ts.var, '_sm']).ly2(mi,:) + ...
                    m_data_hcst.([cfg_ts.var, '_sm']).ly3(mi,:) + ...
                    m_data_hcst.([cfg_ts.var, '_sm']).ly4(mi,:) + ...
                    m_data_hcst.([cfg_ts.var, '_sm']).ly5(mi,:);
                tmp.data(mi,:)=tmp.data(mi,:)/4;
        
                fig_ts.plots=plot(cfg.iyears+2.5, tmp.data(mi,:), 'b');
                set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                
                fig_ts.plots=plot(cfg.iyears_4ym, movmean(m_data_assm.([cfg_ts.var, '_sm'])(mi,:),4,'Endpoints', 'discard'), 'r');
                set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');        
            end
            plot(cfg.iyears+2.5, mean(tmp.data,1), 'b', 'linewidth', 3);
            plot(cfg.iyears_4ym, movmean(mean(m_data_assm.([cfg_ts.var, '_sm']),1),4,'Endpoints', 'discard'), 'r', 'linewidth', 3);
            
            legend ('LE', 'ODA', 'HIND', 'Location', 'Southoutside', 'Orientation', 'Horizontal')
        %     title([cfg_ts.var, ', ', num2str(sta_lonlat{stai}(1)),'E, ', num2str(sta_lonlat{stai}(3)), 'N'])
            title([cfg_ts.var, ', sm, ', num2str(sta_lonlat{stai}(1)),'E, ', ...
                    num2str(sta_lonlat{stai}(2)),'E, ', ...
                    num2str(sta_lonlat{stai}(3)), 'N, ', ...
                    num2str(sta_lonlat{stai}(4)), 'N'])
        
            grid minor
            xlim([1960 2025])
            set(gca, 'fontsize', 20)    
        
            dirs.figdir= [dirs.figroot, filesep, 'ens_all', filesep, cfg_ts.var, '_time_series_indv', filesep, 'l2_5','_', str_prepos];
            if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
            cfg.figname=[dirs.figdir, filesep, str_prepos, '_', num2str(min(months)), '_', num2str(max(months)), '_ts_all_l2_5', '_', num2str(sta_lonlat{stai}(1)), 'E_', ...
                num2str(sta_lonlat{stai}(2)), 'E_', num2str(sta_lonlat{stai}(3)), 'N_', num2str(sta_lonlat{stai}(4)), 'N_',cfg_ts.var, '.tif'];
            print(fig_h, cfg.figname, '-dpng');
            RemoveWhiteSpace([], 'file', cfg.figname);
            close all;
    
        elseif sum(months>12)>0  % DJF(12~14), ...
            for ly=1:5
                fig_h = figure('name','ts','visible','off');
                fig_h.Position= [0,0,1000,500];
        
                tmp.lyear_str=num2str(ly);
                tmp.lyear_str2=num2str(ly+1);
                cfg.iyears_4ym=movmean(cfg.iyears, 4, 'Endpoints', 'discard');
                hold on
                %% LE
                tmp.reshp=reshape(m_data_lens2.(cfg_ts.var), [len_mem_lens2, 12, len_m/12]); % get seasonal mean
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                m_data_lens2.([cfg_ts.var, '_sm'])=(tmp.sm_m1+tmp.sm_m2)/length(months);
    
                for mi=1:50
                    fig_ts.plots=plot(cfg.iyears, m_data_lens2.([cfg_ts.var, '_sm'])(mi,:), 'g');
                    set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                end
                plot(cfg.iyears, mean(m_data_lens2.([cfg_ts.var, '_sm']),1), 'g', 'linewidth', 3);
                
                %% HCST
                tmp.reshp=reshape(m_data_hcst.(cfg_ts.var).(['ly',tmp.lyear_str]), [len_mem_oda, 12, len_m/12]); % get seasonal mean
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                if ly<=4
                    tmp.reshp=reshape(m_data_hcst.(cfg_ts.var).(['ly',tmp.lyear_str2]), [len_mem_oda, 12, len_m/12]); % get seasonal mean
                    tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,:),2));
                else
                    tmp.sm_m2=NaN(size(tmp.sm_m1));
                end
                m_data_hcst.([cfg_ts.var, '_sm']).(['ly',tmp.lyear_str])=(tmp.sm_m1+tmp.sm_m2)/length(months);
    
                %% ODA
                tmp.reshp=reshape(m_data_assm.(cfg_ts.var), [len_mem_oda, 12, len_m/12]); % get seasonal mean
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                m_data_assm.([cfg_ts.var, '_sm'])=(tmp.sm_m1+tmp.sm_m2)/length(months);
                
                for mi=1:20
                    fig_ts.plots=plot(cfg.iyears+ly-1, m_data_hcst.([cfg_ts.var, '_sm']).(['ly',tmp.lyear_str])(mi,:), 'b');
                    set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                    fig_ts.plots=plot(cfg.iyears, m_data_assm.([cfg_ts.var, '_sm'])(mi,:), 'r');
                    set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');        
                end
                plot(cfg.iyears+ly-1, mean(m_data_hcst.([cfg_ts.var, '_sm']).(['ly',tmp.lyear_str]),1), 'b', 'linewidth', 3);
                
                plot(cfg.iyears, mean(m_data_assm.([cfg_ts.var, '_sm']),1), 'r', 'linewidth', 3);
                
                legend ('LE', 'ODA', 'HIND', 'Location', 'Southoutside', 'Orientation', 'Horizontal')
                title([cfg_ts.var, ', sm, ', num2str(sta_lonlat{stai}(1)),'E, ', ...
                    num2str(sta_lonlat{stai}(2)),'E, ', ...
                    num2str(sta_lonlat{stai}(3)), 'N, ', ...
                    num2str(sta_lonlat{stai}(4)), 'N'])
            
                grid minor
                xlim([1960 2025])
                set(gca, 'fontsize', 20)    
            
                dirs.figdir= [dirs.figroot, filesep, 'ens_all', filesep, cfg_ts.var, '_time_series_indv', filesep, 'l',tmp.lyear_str, '_', str_prepos];
                if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
                cfg.figname=[dirs.figdir, filesep, str_prepos, '_', num2str(min(months)), '_', num2str(max(months)), '_ts_all_l',tmp.lyear_str, '_', num2str(sta_lonlat{stai}(1)), 'E_', ...
                    num2str(sta_lonlat{stai}(2)), 'E_', num2str(sta_lonlat{stai}(3)), 'N_',num2str(sta_lonlat{stai}(4)), 'N_', cfg_ts.var, '.tif'];
                print(fig_h, cfg.figname, '-dpng'); 
                RemoveWhiteSpace([], 'file', cfg.figname);
                close all;
            end
        
            fig_h = figure('name','ts','visible','off');
            fig_h.Position= [0,0,1000,500];
        
            tmp.lyear_str=num2str(ly);
            hold on
            cfg.iyears_4ym=movmean(cfg.iyears, 4, 'Endpoints', 'discard');
            for mi=1:50
                fig_ts.plots=plot(cfg.iyears_4ym, movmean(m_data_lens2.([cfg_ts.var, '_sm'])(mi,:),4,'Endpoints', 'discard'), 'g');
                set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
            plot(cfg.iyears_4ym, movmean(mean(m_data_lens2.([cfg_ts.var, '_sm']),1),4,'Endpoints', 'discard'), 'g', 'linewidth', 3);
            for mi=1:20
        %         fig_ts.plots=plot(cfg.iyears+ly-1, m_data_hcst.([cfg_ts.var, '_ym']).(['ly',tmp.lyear_str])(mi,:), 'b');
        %         set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                
                tmp.data(mi,:)= m_data_hcst.([cfg_ts.var, '_sm']).ly2(mi,:) + ...
                    m_data_hcst.([cfg_ts.var, '_sm']).ly3(mi,:) + ...
                    m_data_hcst.([cfg_ts.var, '_sm']).ly4(mi,:) + ...
                    m_data_hcst.([cfg_ts.var, '_sm']).ly5(mi,:);
                tmp.data(mi,:)=tmp.data(mi,:)/4;
        
                fig_ts.plots=plot(cfg.iyears+2.5, tmp.data(mi,:), 'b');
                set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                
                fig_ts.plots=plot(cfg.iyears_4ym, movmean(m_data_assm.([cfg_ts.var, '_sm'])(mi,:),4,'Endpoints', 'discard'), 'r');
                set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');        
            end
            plot(cfg.iyears+2.5, mean(tmp.data,1), 'b', 'linewidth', 3);
            plot(cfg.iyears_4ym, movmean(mean(m_data_assm.([cfg_ts.var, '_sm']),1),4,'Endpoints', 'discard'), 'r', 'linewidth', 3);
            
            legend ('LE', 'ODA', 'HIND', 'Location', 'Southoutside', 'Orientation', 'Horizontal')
        %     title([cfg_ts.var, ', ', num2str(sta_lonlat{stai}(1)),'E, ', num2str(sta_lonlat{stai}(3)), 'N'])
            title([cfg_ts.var, ', sm, ', num2str(sta_lonlat{stai}(1)),'E, ', ...
                    num2str(sta_lonlat{stai}(2)),'E, ', ...
                    num2str(sta_lonlat{stai}(3)), 'N, ', ...
                    num2str(sta_lonlat{stai}(4)), 'N'])
        
            grid minor
            xlim([1960 2025])
            set(gca, 'fontsize', 20)    
        
            dirs.figdir= [dirs.figroot, filesep, 'ens_all', filesep, cfg_ts.var, '_time_series_indv', filesep, 'l2_5','_', str_prepos];
            if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
            cfg.figname=[dirs.figdir, filesep, str_prepos, '_', num2str(min(months)), '_', num2str(max(months)), '_ts_all_l2_5', '_', num2str(sta_lonlat{stai}(1)), 'E_', ...
                num2str(sta_lonlat{stai}(2)), 'E_', num2str(sta_lonlat{stai}(3)), 'N_', num2str(sta_lonlat{stai}(4)), 'N_',cfg_ts.var, '.tif'];
            print(fig_h, cfg.figname, '-dpng');
            RemoveWhiteSpace([], 'file', cfg.figname);
            close all;
    
            
        end
    
    
        %% with obs
        if length(months)==12
            %% annual graph
            for ly=1:5
                fig_h = figure('name','ts','visible','off');
                fig_h.Position= [0,0,1000,500];
        
                tmp.lyear_str=num2str(ly);
                hold on
                for mi=1:50
                    fig_ts.plots=plot(cfg.iyears, m_data_lens2.([cfg_ts.var, '_ym'])(mi,:), 'g');
                    set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                end
                plot(cfg.iyears, mean(m_data_lens2.([cfg_ts.var, '_ym']),1), 'g', 'linewidth', 3);
                
                for mi=1:20
                    fig_ts.plots=plot(cfg.iyears+ly-1, m_data_hcst.([cfg_ts.var, '_ym']).(['ly',tmp.lyear_str])(mi,:), 'b');
                    set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                    
                    fig_ts.plots=plot(cfg.iyears, m_data_assm.([cfg_ts.var, '_ym'])(mi,:), 'r');
                    set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');        
                end
                plot(cfg.iyears+ly-1, mean(m_data_hcst.([cfg_ts.var, '_ym']).(['ly',tmp.lyear_str]),1), 'b', 'linewidth', 3);
                
                plot(cfg.iyears, mean(m_data_assm.([cfg_ts.var, '_ym']),1), 'r', 'linewidth', 3);
                
                mm_assm=mean(mean(m_data_assm.([cfg_ts.var, '_ym']),1),'omitnan');
                std_assm=std(mean(m_data_assm.([cfg_ts.var, '_ym']),1),'omitnan');
                mm_obs=mean(m_data_obs.([cfg_ts.var, '_ym']),'omitnan');
                if ~(isfinite(mm_obs) && mm_obs <= mm_assm+std_assm && mm_obs >= mm_assm-std_assm)
                    yyaxis right;
                end

                fig_ts.plots=plot(cfg.iyears, m_data_obs.([cfg_ts.var, '_ym'])(:), 'k', 'linewidth', 3);
    
                legend ('LE', 'ODA', 'HIND', 'OBS', 'Location', 'Southoutside', 'Orientation', 'Horizontal')
                title([cfg_ts.var, ', ', num2str(sta_lonlat{stai}(1)),'E, ', ...
                    num2str(sta_lonlat{stai}(2)),'E, ', ...
                    num2str(sta_lonlat{stai}(3)), 'N, ', ...
                    num2str(sta_lonlat{stai}(4)), 'N'])
            
                grid minor
                xlim([1960 2025])
                set(gca, 'fontsize', 20)    
            
                dirs.figdir= [dirs.figroot, filesep, 'ens_all', filesep, cfg_ts.var, '_time_series_indv_withobs', filesep, 'l',tmp.lyear_str];
                if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
                cfg.figname=[dirs.figdir, filesep, 'ts_all_l',tmp.lyear_str, '_', num2str(sta_lonlat{stai}(1)), 'E_', ...
                    num2str(sta_lonlat{stai}(2)), 'E_', num2str(sta_lonlat{stai}(3)), 'N_',num2str(sta_lonlat{stai}(4)), 'N_', cfg_ts.var, '.tif'];
                print(fig_h, cfg.figname, '-dpng');
                RemoveWhiteSpace([], 'file', cfg.figname);
                close all;
            end
        
            fig_h = figure('name','ts','visible','off');
            fig_h.Position= [0,0,1000,500];
        
            tmp.lyear_str=num2str(ly);
            hold on
            cfg.iyears_4ym=movmean(cfg.iyears, 4, 'Endpoints', 'discard');
            for mi=1:50
                fig_ts.plots=plot(cfg.iyears_4ym, movmean(m_data_lens2.([cfg_ts.var, '_ym'])(mi,:),4,'Endpoints', 'discard'), 'g');
                set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
            plot(cfg.iyears_4ym, movmean(mean(m_data_lens2.([cfg_ts.var, '_ym']),1),4,'Endpoints', 'discard'), 'g', 'linewidth', 3);
            for mi=1:20
        %         fig_ts.plots=plot(cfg.iyears+ly-1, m_data_hcst.([cfg_ts.var, '_ym']).(['ly',tmp.lyear_str])(mi,:), 'b');
        %         set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                
                tmp.data(mi,:)= m_data_hcst.([cfg_ts.var, '_ym']).ly2(mi,:) + ...
                    m_data_hcst.([cfg_ts.var, '_ym']).ly3(mi,:) + ...
                    m_data_hcst.([cfg_ts.var, '_ym']).ly4(mi,:) + ...
                    m_data_hcst.([cfg_ts.var, '_ym']).ly5(mi,:);
                tmp.data(mi,:)=tmp.data(mi,:)/4;
        
                fig_ts.plots=plot(cfg.iyears+2.5, tmp.data(mi,:), 'b');
                set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                
                fig_ts.plots=plot(cfg.iyears_4ym, movmean(m_data_assm.([cfg_ts.var, '_ym'])(mi,:),4,'Endpoints', 'discard'), 'r');
                set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');        
            end
        %     plot(cfg.iyears+ly-1, mean(m_data_hcst.([cfg_ts.var, '_ym']).(['ly',tmp.lyear_str]),1), 'b', 'linewidth', 3);
            plot(cfg.iyears+2.5, mean(tmp.data,1), 'b', 'linewidth', 3);
            plot(cfg.iyears_4ym, movmean(mean(m_data_assm.([cfg_ts.var, '_ym']),1),4,'Endpoints', 'discard'), 'r', 'linewidth', 3);
            
            abc=movmean(mean(m_data_assm.([cfg_ts.var, '_ym']),1),4,'Endpoints', 'discard');
            mm_assm=mean(abc,'omitnan');
            std_assm=std(abc,'omitnan');
            abcd=movmean(m_data_obs.([cfg_ts.var, '_ym'])(:),4,'Endpoints', 'discard');
            mm_obs=mean(abcd,'omitnan');
            if ~(isfinite(mm_obs) && mm_obs <= mm_assm+std_assm && mm_obs >= mm_assm-std_assm)
                yyaxis right;
            end

            plot(cfg.iyears_4ym, movmean(m_data_obs.([cfg_ts.var, '_ym'])(:),4,'Endpoints', 'discard'), 'k', 'linewidth', 3);
    
    
            legend ('LE', 'ODA', 'HIND', 'OBS', 'Location', 'Southoutside', 'Orientation', 'Horizontal')
        %     title([cfg_ts.var, ', ', num2str(sta_lonlat{stai}(1)),'E, ', num2str(sta_lonlat{stai}(3)), 'N'])
            title([cfg_ts.var, ', ', num2str(sta_lonlat{stai}(1)),'E, ', ...
                    num2str(sta_lonlat{stai}(2)),'E, ', ...
                    num2str(sta_lonlat{stai}(3)), 'N, ', ...
                    num2str(sta_lonlat{stai}(4)), 'N'])
        
            grid minor
            xlim([1960 2025])
            set(gca, 'fontsize', 20)    
        
            dirs.figdir= [dirs.figroot, filesep, 'ens_all', filesep, cfg_ts.var, '_time_series_indv_withobs', filesep, 'l2_5'];
            if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
            cfg.figname=[dirs.figdir, filesep, 'ts_all_l2_5', '_', num2str(sta_lonlat{stai}(1)), 'E_', ...
                num2str(sta_lonlat{stai}(2)), 'E_', num2str(sta_lonlat{stai}(3)), 'N_', num2str(sta_lonlat{stai}(4)), 'N_',cfg_ts.var, '.tif'];
            print(fig_h, cfg.figname, '-dpng');
            RemoveWhiteSpace([], 'file', cfg.figname);
            close all;
        elseif sum(months>12)==0 % 1~3, 2~4, ... 10~12
            for ly=1:5
                fig_h = figure('name','ts','visible','off');
                fig_h.Position= [0,0,1000,500];
        
                tmp.lyear_str=num2str(ly);
                cfg.iyears_4ym=movmean(cfg.iyears, 4, 'Endpoints', 'discard');
                hold on
                tmp.reshp=reshape(m_data_lens2.(cfg_ts.var), [len_mem_lens2, 12, len_m/12]); % get seasonal mean
                m_data_lens2.([cfg_ts.var, '_sm'])=squeeze(mean(tmp.reshp(:,months,:),2));
                for mi=1:50
                    fig_ts.plots=plot(cfg.iyears, m_data_lens2.([cfg_ts.var, '_sm'])(mi,:), 'g');
                    set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                end
                plot(cfg.iyears, mean(m_data_lens2.([cfg_ts.var, '_sm']),1), 'g', 'linewidth', 3);
                
                tmp.reshp=reshape(m_data_hcst.(cfg_ts.var).(['ly',tmp.lyear_str]), [len_mem_oda, 12, len_m/12]); % get seasonal mean
                m_data_hcst.([cfg_ts.var, '_sm']).(['ly',tmp.lyear_str])=squeeze(mean(tmp.reshp(:,months,:),2));
                tmp.reshp=reshape(m_data_assm.(cfg_ts.var), [len_mem_oda, 12, len_m/12]); % get seasonal mean
                m_data_assm.([cfg_ts.var, '_sm'])=squeeze(mean(tmp.reshp(:,months,:),2));
                for mi=1:20
                    fig_ts.plots=plot(cfg.iyears+ly-1, m_data_hcst.([cfg_ts.var, '_sm']).(['ly',tmp.lyear_str])(mi,:), 'b');
                    set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                    fig_ts.plots=plot(cfg.iyears, m_data_assm.([cfg_ts.var, '_sm'])(mi,:), 'r');
                    set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');        
                end
                plot(cfg.iyears+ly-1, mean(m_data_hcst.([cfg_ts.var, '_sm']).(['ly',tmp.lyear_str]),1), 'b', 'linewidth', 3);
                
                plot(cfg.iyears, mean(m_data_assm.([cfg_ts.var, '_sm']),1), 'r', 'linewidth', 3);
                mm_assm=mean(mean(m_data_assm.([cfg_ts.var, '_sm']),1),'omitnan');
                std_assm=std(mean(m_data_assm.([cfg_ts.var, '_sm']),1),'omitnan');
                mm_obs=mean(m_data_obs.([cfg_ts.var, '_ym']),'omitnan');
                if ~(isfinite(mm_obs) && mm_obs <= mm_assm+std_assm && mm_obs >= mm_assm-std_assm)
                    yyaxis right;
                end
                fig_ts.plots=plot(cfg.iyears, m_data_obs.([cfg_ts.var, '_ym'])(:), 'k', 'linewidth', 3);
    
                legend ('LE', 'ODA', 'HIND', 'OBS', 'Location', 'Southoutside', 'Orientation', 'Horizontal')
                title([cfg_ts.var, ', sm, ', num2str(sta_lonlat{stai}(1)),'E, ', ...
                    num2str(sta_lonlat{stai}(2)),'E, ', ...
                    num2str(sta_lonlat{stai}(3)), 'N, ', ...
                    num2str(sta_lonlat{stai}(4)), 'N'])
            
                grid minor
                xlim([1960 2025])
                set(gca, 'fontsize', 20)    
            
                dirs.figdir= [dirs.figroot, filesep, 'ens_all', filesep, cfg_ts.var, '_time_series_indv_withobs', filesep, 'l',tmp.lyear_str, '_', str_prepos];
                if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
                cfg.figname=[dirs.figdir, filesep, str_prepos, '_', num2str(min(months)), '_', num2str(max(months)), '_ts_all_l',tmp.lyear_str, '_', num2str(sta_lonlat{stai}(1)), 'E_', ...
                    num2str(sta_lonlat{stai}(2)), 'E_', num2str(sta_lonlat{stai}(3)), 'N_',num2str(sta_lonlat{stai}(4)), 'N_', cfg_ts.var, '.tif'];
                print(fig_h, cfg.figname, '-dpng'); 
                RemoveWhiteSpace([], 'file', cfg.figname);
                close all;
            end
        
            fig_h = figure('name','ts','visible','off');
            fig_h.Position= [0,0,1000,500];
        
            tmp.lyear_str=num2str(ly);
            hold on
            cfg.iyears_4ym=movmean(cfg.iyears, 4, 'Endpoints', 'discard');
            for mi=1:50
                fig_ts.plots=plot(cfg.iyears_4ym, movmean(m_data_lens2.([cfg_ts.var, '_sm'])(mi,:),4,'Endpoints', 'discard'), 'g');
                set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
            plot(cfg.iyears_4ym, movmean(mean(m_data_lens2.([cfg_ts.var, '_sm']),1),4,'Endpoints', 'discard'), 'g', 'linewidth', 3);
            for mi=1:20
        %         fig_ts.plots=plot(cfg.iyears+ly-1, m_data_hcst.([cfg_ts.var, '_ym']).(['ly',tmp.lyear_str])(mi,:), 'b');
        %         set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                
                tmp.data(mi,:)= m_data_hcst.([cfg_ts.var, '_sm']).ly2(mi,:) + ...
                    m_data_hcst.([cfg_ts.var, '_sm']).ly3(mi,:) + ...
                    m_data_hcst.([cfg_ts.var, '_sm']).ly4(mi,:) + ...
                    m_data_hcst.([cfg_ts.var, '_sm']).ly5(mi,:);
                tmp.data(mi,:)=tmp.data(mi,:)/4;
        
                fig_ts.plots=plot(cfg.iyears+2.5, tmp.data(mi,:), 'b');
                set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                
                fig_ts.plots=plot(cfg.iyears_4ym, movmean(m_data_assm.([cfg_ts.var, '_sm'])(mi,:),4,'Endpoints', 'discard'), 'r');
                set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');        
            end
            plot(cfg.iyears+2.5, mean(tmp.data,1), 'b', 'linewidth', 3);
            plot(cfg.iyears_4ym, movmean(mean(m_data_assm.([cfg_ts.var, '_sm']),1),4,'Endpoints', 'discard'), 'r', 'linewidth', 3);
            
            abc=movmean(mean(m_data_assm.([cfg_ts.var, '_sm']),1),4,'Endpoints', 'discard');
            mm_assm=mean(abc,'omitnan');
            std_assm=std(abc,'omitnan');
            abcd=movmean(m_data_obs.([cfg_ts.var, '_ym'])(:),4,'Endpoints', 'discard');
            mm_obs=mean(abcd,'omitnan');
            if ~(isfinite(mm_obs) && mm_obs <= mm_assm+std_assm && mm_obs >= mm_assm-std_assm)
                yyaxis right;
            end

            plot(cfg.iyears_4ym, movmean(m_data_obs.([cfg_ts.var, '_ym'])(:),4,'Endpoints', 'discard'), 'k', 'linewidth', 3);
    
            legend ('LE', 'ODA', 'HIND', 'OBS', 'Location', 'Southoutside', 'Orientation', 'Horizontal')
        %     title([cfg_ts.var, ', ', num2str(sta_lonlat{stai}(1)),'E, ', num2str(sta_lonlat{stai}(3)), 'N'])
            title([cfg_ts.var, ', sm, ', num2str(sta_lonlat{stai}(1)),'E, ', ...
                    num2str(sta_lonlat{stai}(2)),'E, ', ...
                    num2str(sta_lonlat{stai}(3)), 'N, ', ...
                    num2str(sta_lonlat{stai}(4)), 'N'])
        
            grid minor
            xlim([1960 2025])
            set(gca, 'fontsize', 20)    
        
            dirs.figdir= [dirs.figroot, filesep, 'ens_all', filesep, cfg_ts.var, '_time_series_indv_withobs', filesep, 'l2_5','_', str_prepos];
            if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
            cfg.figname=[dirs.figdir, filesep, str_prepos, '_', num2str(min(months)), '_', num2str(max(months)), '_ts_all_l2_5', '_', num2str(sta_lonlat{stai}(1)), 'E_', ...
                num2str(sta_lonlat{stai}(2)), 'E_', num2str(sta_lonlat{stai}(3)), 'N_', num2str(sta_lonlat{stai}(4)), 'N_',cfg_ts.var, '.tif'];
            print(fig_h, cfg.figname, '-dpng');
            RemoveWhiteSpace([], 'file', cfg.figname);
            close all;
    
        elseif sum(months>12)>0  % DJF(12~14), ...
            for ly=1:5
                fig_h = figure('name','ts','visible','off');
                fig_h.Position= [0,0,1000,500];
        
                tmp.lyear_str=num2str(ly);
                tmp.lyear_str2=num2str(ly+1);
                cfg.iyears_4ym=movmean(cfg.iyears, 4, 'Endpoints', 'discard');
                hold on
                %% LE
                tmp.reshp=reshape(m_data_lens2.(cfg_ts.var), [len_mem_lens2, 12, len_m/12]); % get seasonal mean
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                m_data_lens2.([cfg_ts.var, '_sm'])=(tmp.sm_m1+tmp.sm_m2)/length(months);
    
                for mi=1:50
                    fig_ts.plots=plot(cfg.iyears, m_data_lens2.([cfg_ts.var, '_sm'])(mi,:), 'g');
                    set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                end
                plot(cfg.iyears, mean(m_data_lens2.([cfg_ts.var, '_sm']),1), 'g', 'linewidth', 3);
                
                %% HCST
                tmp.reshp=reshape(m_data_hcst.(cfg_ts.var).(['ly',tmp.lyear_str]), [len_mem_oda, 12, len_m/12]); % get seasonal mean
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                if ly<=4
                    tmp.reshp=reshape(m_data_hcst.(cfg_ts.var).(['ly',tmp.lyear_str2]), [len_mem_oda, 12, len_m/12]); % get seasonal mean
                    tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,:),2));
                else
                    tmp.sm_m2=NaN(size(tmp.sm_m1));
                end
                m_data_hcst.([cfg_ts.var, '_sm']).(['ly',tmp.lyear_str])=(tmp.sm_m1+tmp.sm_m2)/length(months);
    
                %% ODA
                tmp.reshp=reshape(m_data_assm.(cfg_ts.var), [len_mem_oda, 12, len_m/12]); % get seasonal mean
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                m_data_assm.([cfg_ts.var, '_sm'])=(tmp.sm_m1+tmp.sm_m2)/length(months);
                
                for mi=1:20
                    fig_ts.plots=plot(cfg.iyears+ly-1, m_data_hcst.([cfg_ts.var, '_sm']).(['ly',tmp.lyear_str])(mi,:), 'b');
                    set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                    fig_ts.plots=plot(cfg.iyears, m_data_assm.([cfg_ts.var, '_sm'])(mi,:), 'r');
                    set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');        
                end
                plot(cfg.iyears+ly-1, mean(m_data_hcst.([cfg_ts.var, '_sm']).(['ly',tmp.lyear_str]),1), 'b', 'linewidth', 3);
                
                plot(cfg.iyears, mean(m_data_assm.([cfg_ts.var, '_sm']),1), 'r', 'linewidth', 3);
                mm_assm=mean(mean(m_data_assm.([cfg_ts.var, '_sm']),1),'omitnan');
                std_assm=std(mean(m_data_assm.([cfg_ts.var, '_sm']),1),'omitnan');
                mm_obs=mean(m_data_obs.([cfg_ts.var, '_ym']),'omitnan');
                if ~(isfinite(mm_obs) && mm_obs <= mm_assm+std_assm && mm_obs >= mm_assm-std_assm)
                    yyaxis right;
                end
                fig_ts.plots=plot(cfg.iyears, m_data_obs.([cfg_ts.var, '_ym'])(:), 'k', 'linewidth', 3);
    
    
                legend ('LE', 'ODA', 'HIND', 'OBS', 'Location', 'Southoutside', 'Orientation', 'Horizontal')
                title([cfg_ts.var, ', sm, ', num2str(sta_lonlat{stai}(1)),'E, ', ...
                    num2str(sta_lonlat{stai}(2)),'E, ', ...
                    num2str(sta_lonlat{stai}(3)), 'N, ', ...
                    num2str(sta_lonlat{stai}(4)), 'N'])
            
                grid minor
                xlim([1960 2025])
                set(gca, 'fontsize', 20)    
            
                dirs.figdir= [dirs.figroot, filesep, 'ens_all', filesep, cfg_ts.var, '_time_series_indv_withobs', filesep, 'l',tmp.lyear_str, '_', str_prepos];
                if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
                cfg.figname=[dirs.figdir, filesep, str_prepos, '_', num2str(min(months)), '_', num2str(max(months)), '_ts_all_l',tmp.lyear_str, '_', num2str(sta_lonlat{stai}(1)), 'E_', ...
                    num2str(sta_lonlat{stai}(2)), 'E_', num2str(sta_lonlat{stai}(3)), 'N_',num2str(sta_lonlat{stai}(4)), 'N_', cfg_ts.var, '.tif'];
                print(fig_h, cfg.figname, '-dpng'); 
                RemoveWhiteSpace([], 'file', cfg.figname);
                close all;
            end
        
            fig_h = figure('name','ts','visible','off');
            fig_h.Position= [0,0,1000,500];
        
            tmp.lyear_str=num2str(ly);
            hold on
            cfg.iyears_4ym=movmean(cfg.iyears, 4, 'Endpoints', 'discard');
            for mi=1:50
                fig_ts.plots=plot(cfg.iyears_4ym, movmean(m_data_lens2.([cfg_ts.var, '_sm'])(mi,:),4,'Endpoints', 'discard'), 'g');
                set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
            plot(cfg.iyears_4ym, movmean(mean(m_data_lens2.([cfg_ts.var, '_sm']),1),4,'Endpoints', 'discard'), 'g', 'linewidth', 3);
            for mi=1:20
        %         fig_ts.plots=plot(cfg.iyears+ly-1, m_data_hcst.([cfg_ts.var, '_ym']).(['ly',tmp.lyear_str])(mi,:), 'b');
        %         set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                
                tmp.data(mi,:)= m_data_hcst.([cfg_ts.var, '_sm']).ly2(mi,:) + ...
                    m_data_hcst.([cfg_ts.var, '_sm']).ly3(mi,:) + ...
                    m_data_hcst.([cfg_ts.var, '_sm']).ly4(mi,:) + ...
                    m_data_hcst.([cfg_ts.var, '_sm']).ly5(mi,:);
                tmp.data(mi,:)=tmp.data(mi,:)/4;
        
                fig_ts.plots=plot(cfg.iyears+2.5, tmp.data(mi,:), 'b');
                set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                
                fig_ts.plots=plot(cfg.iyears_4ym, movmean(m_data_assm.([cfg_ts.var, '_sm'])(mi,:),4,'Endpoints', 'discard'), 'r');
                set(get(get(fig_ts.plots,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');        
            end
            plot(cfg.iyears+2.5, mean(tmp.data,1), 'b', 'linewidth', 3);
            plot(cfg.iyears_4ym, movmean(mean(m_data_assm.([cfg_ts.var, '_sm']),1),4,'Endpoints', 'discard'), 'r', 'linewidth', 3);
            abc=movmean(mean(m_data_assm.([cfg_ts.var, '_sm']),1),4,'Endpoints', 'discard');
            mm_assm=mean(abc,'omitnan');
            std_assm=std(abc,'omitnan');
            abcd=movmean(m_data_obs.([cfg_ts.var, '_ym'])(:),4,'Endpoints', 'discard');
            mm_obs=mean(abcd,'omitnan');
            if ~(isfinite(mm_obs) && mm_obs <= mm_assm+std_assm && mm_obs >= mm_assm-std_assm)
                yyaxis right;
            end
            plot(cfg.iyears_4ym, movmean(m_data_obs.([cfg_ts.var, '_ym'])(:),4,'Endpoints', 'discard'), 'k', 'linewidth', 3);
    
            legend ('LE', 'ODA', 'HIND', 'OBS', 'Location', 'Southoutside', 'Orientation', 'Horizontal')
        %     title([cfg_ts.var, ', ', num2str(sta_lonlat{stai}(1)),'E, ', num2str(sta_lonlat{stai}(3)), 'N'])
            title([cfg_ts.var, ', sm, ', num2str(sta_lonlat{stai}(1)),'E, ', ...
                    num2str(sta_lonlat{stai}(2)),'E, ', ...
                    num2str(sta_lonlat{stai}(3)), 'N, ', ...
                    num2str(sta_lonlat{stai}(4)), 'N'])
        
            grid minor
            xlim([1960 2025])
            set(gca, 'fontsize', 20)    
        
            dirs.figdir= [dirs.figroot, filesep, 'ens_all', filesep, cfg_ts.var, '_time_series_indv_withobs', filesep, 'l2_5','_', str_prepos];
            if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
            cfg.figname=[dirs.figdir, filesep, str_prepos, '_', num2str(min(months)), '_', num2str(max(months)), '_ts_all_l2_5', '_', num2str(sta_lonlat{stai}(1)), 'E_', ...
                num2str(sta_lonlat{stai}(2)), 'E_', num2str(sta_lonlat{stai}(3)), 'N_', num2str(sta_lonlat{stai}(4)), 'N_',cfg_ts.var, '.tif'];
            print(fig_h, cfg.figname, '-dpng');
            RemoveWhiteSpace([], 'file', cfg.figname);
            close all;
    
            
        end
    
    
    
        
    %     plot(m_data_obs.([cfg_ts.var, '_ym']), 'k', 'linewidth', 2);
    %     hold off
        
        [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_20', tmp.dropboxpath);
        
        % abc(:,:,1)=rand(3,3);
        % abc(:,:,2)=rand(3,3);
        % abc(:,:,3)=rand(3,3);
        % abc(:,:,4)=rand(3,3);
        % plot(squeeze(mean(mean(abc,1),2)))
        % tmmean(1,1,1:4)=1;
        % abcd=abc-tmmean;
        % plot(squeeze(mean(mean(abcd,1),2)))    
        
        % ly=4;
        % close all;
        % hold on
        % for mi=1:50
        %     plot(1960:2020, m_data_lens2.([cfg_ts.var, '_ym'])(mi,:), 'g')
        % end
        % plot(1960:2020, mean( m_data_lens2.([cfg_ts.var, '_ym']),1), 'g', 'linewidth', 3)
        % for mi=1:20
        %     plot(1960+ly-1:2020+ly-1, m_data_hcst.([cfg_ts.var, '_ym']).(['ly',num2str(ly)])(mi,:), 'b')
        %     plot(1960:2020, m_data_assm.([cfg_ts.var, '_ym'])(mi,:), 'r')
        % end
        % plot(1960+ly-1:2020+ly-1, mean(m_data_hcst.([cfg_ts.var, '_ym']).(['ly',num2str(ly)]),1), 'b', 'linewidth', 3)
        % plot(1960:2020, mean(m_data_assm.([cfg_ts.var, '_ym']),1), 'r', 'linewidth', 3)
        % 
        % plot(m_data_obs.([cfg_ts.var, '_ym']), 'k', 'linewidth', 2)
        % hold off
        
        
        if once==1
            %% Indices read
            %% AMO
            load ('/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/clim_indices/clim_indices_SST_all_AMO_obs_ERSST.mat');
            %% ENSO
            load ('/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/clim_indices/clim_indices_SST_all_ENSO_obs_ERSST.mat');
            %% IOD
            load ('/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/clim_indices/clim_indices_SST_all_IOD_obs_ERSST.mat');
            %% IPO
            load ('/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/clim_indices/clim_indices_SST_all_IPO_obs_ERSST.mat');
            %% PDO
            load ('/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/clim_indices/clim_indices_SST_all_PDO_obs_ERSST.mat');
            %% PDO_l
            load ('/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/clim_indices/clim_indices_SST_all_PDO_l_obs_ERSST.mat');
            %% NPGO
            data_NPGO=load ('/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/clim_indices/clim_indices_SSH_all_PDO_obs_CMEMS.mat');
            %% SAM
            load ('/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/clim_indices/clim_indices_PSL_all_SAM_obs_ERA5.mat');
            %% SAM_d
            load ('/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/clim_indices/clim_indices_PSL_all_SAM_d_obs_ERA5.mat');
            %% ATL3
            load ('/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/clim_indices/clim_indices_SST_all_ATL3_obs_ERSST.mat');
            %% ATL3_d
            load ('/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/clim_indices/clim_indices_SST_all_ATL3_d_obs_ERSST.mat');
            %% TNA
            load ('/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/clim_indices/clim_indices_SST_all_TNA_obs_ERSST.mat');
            %% TSA
            load ('/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/clim_indices/clim_indices_SST_all_TSA_obs_ERSST.mat');
            %% KEI_O
            load ('/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/clim_indices/clim_indices_SSH_all_KEI_O_obs_CMEMS.mat');
            %% KEI_LE
            load ('/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/clim_indices/clim_indices_SSH_all_KEI_LE_obs_CMEMS.mat');
            %% NAO_STA
            load ('/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/clim_indices/clim_indices_PSL_all_NAO_STA_obs_ERA5.mat');
            %% NAO_PC
            load ('/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/clim_indices/clim_indices_PSL_all_NAO_PC_obs_ERA5.mat');
            
            %% yearly mean indices
            %% AMO
            len_m=size(data_AMO.obs_dseason,1);
            len_mem_oda=size(data_AMO.assm_dseason,1);
            len_mem_lens2=size(data_AMO.lens2_dseason,1);
            
            if sum(months>12)==0
                tmp.reshp=reshape(data_AMO.obs_dseason, [12, len_m/12]);
                data_AMO.obs_ym=squeeze(mean(tmp.reshp(months,:),1));
                data_AMO.obs_4ym=movmean(data_AMO.obs_ym, 4, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_AMO.assm_dseason, [len_mem_oda, 12, len_m/12]);
                data_AMO.assm_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                data_AMO.assm_4ym=movmean(data_AMO.assm_ym, 4, 2, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_AMO.lens2_dseason, [len_mem_lens2, 12, len_m/12]);
                data_AMO.lens2_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                data_AMO.lens2_4ym=movmean(data_AMO.lens2_ym, 4, 2, 'Endpoints', 'discard');
    
                tmp.reshp=reshape(data_AMO.hcst_dseason, [5, len_mem_oda, 12, len_m/12]);
                data_AMO.hcst_ym=squeeze(mean(tmp.reshp(:,:,months,:),3));
                data_AMO.hcst_4ym=squeeze(mean(data_AMO.hcst_ym(2:5,:,:),1));
    
    %         figure;
    %         for mi=1:len_mem_lens2
    %             hold on
    %             plot(data_AMO.lens2_ym(mi,:));
    %         end
    
            else
                tmp.reshp=reshape(data_AMO.obs_dseason, [12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(month1,:),1));
                tmp.sm_m2=squeeze(sum(tmp.reshp(month2-12,2:end),1));  tmp.sm_m2(end+1)=NaN;
                data_AMO.obs_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_AMO.obs_4ym=movmean(data_AMO.obs_ym, 4, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_AMO.assm_dseason, [len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_AMO.assm_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_AMO.assm_4ym=movmean(data_AMO.assm_ym, 4, 2, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_AMO.lens2_dseason, [len_mem_lens2, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_AMO.lens2_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_AMO.lens2_4ym=movmean(data_AMO.lens2_ym, 4, 2, 'Endpoints', 'discard');
    
                tmp.reshp=reshape(data_AMO.hcst_dseason, [5, len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,:,month1,:),3));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,:,month2-12,:),3));  
                data_AMO.hcst_ym=(tmp.sm_m1(1:4,:,:)+tmp.sm_m2(2:5,:,:))/length(months);
                data_AMO.hcst_ym(end+1,:,:)=NaN;
                data_AMO.hcst_4ym=squeeze(mean(data_AMO.hcst_ym(2:5,:,:),1));
            end
            
            %% ENSO
            if sum(months>12)==0
                tmp.reshp=reshape(data_ENSO.obs, [12, len_m/12]);
                data_ENSO.obs_ym=squeeze(mean(tmp.reshp(months,:),1));
                data_ENSO.obs_4ym=movmean(data_ENSO.obs_ym, 4, 'Endpoints', 'discard');
        
                tmp.reshp=reshape(data_ENSO.assm, [len_mem_oda, 12, len_m/12]);
                data_ENSO.assm_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                data_ENSO.assm_4ym=movmean(data_ENSO.assm_ym, 4, 2, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_ENSO.lens2, [len_mem_lens2, 12, len_m/12]);
                data_ENSO.lens2_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                data_ENSO.lens2_4ym=movmean(data_ENSO.lens2_ym, 4, 2, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_ENSO.hcst, [5, len_mem_oda, 12, len_m/12]);
                data_ENSO.hcst_ym=squeeze(mean(tmp.reshp(:,:,months,:),3));
                data_ENSO.hcst_4ym=squeeze(mean(data_ENSO.hcst_ym(2:5,:,:),1)); 
            else
                tmp.reshp=reshape(data_ENSO.obs, [12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(month1,:),1));
                tmp.sm_m2=squeeze(sum(tmp.reshp(month2-12,2:end),1));  tmp.sm_m2(end+1)=NaN;
                data_ENSO.obs_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_ENSO.obs_4ym=movmean(data_ENSO.obs_ym, 4, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_ENSO.assm, [len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_ENSO.assm_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_ENSO.assm_4ym=movmean(data_ENSO.assm_ym, 4, 2, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_ENSO.lens2, [len_mem_lens2, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_ENSO.lens2_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_ENSO.lens2_4ym=movmean(data_ENSO.lens2_ym, 4, 2, 'Endpoints', 'discard');
    
                tmp.reshp=reshape(data_ENSO.hcst, [5, len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,:,month1,:),3));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,:,month2-12,:),3));  
                data_ENSO.hcst_ym=(tmp.sm_m1(1:4,:,:)+tmp.sm_m2(2:5,:,:))/length(months);
                data_ENSO.hcst_ym(end+1,:,:)=NaN;
                data_ENSO.hcst_4ym=squeeze(mean(data_ENSO.hcst_ym(2:5,:,:),1));
            end
    
            %% IOD
            if sum(months>12)==0
                tmp.reshp=reshape(data_IOD.obs, [12, len_m/12]);
                data_IOD.obs_ym=squeeze(mean(tmp.reshp(months,:),1));
                
                tmp.reshp=reshape(data_IOD.assm, [len_mem_oda, 12, len_m/12]);
                data_IOD.assm_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                
                tmp.reshp=reshape(data_IOD.lens2, [len_mem_lens2, 12, len_m/12]);
                data_IOD.lens2_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                
                tmp.reshp=reshape(data_IOD.hcst, [5, len_mem_oda, 12, len_m/12]);
                data_IOD.hcst_ym=squeeze(mean(tmp.reshp(:,:,months,:),3));
        
                data_IOD.obs_4ym=movmean(data_IOD.obs_ym, 4, 'Endpoints', 'discard');
                data_IOD.assm_4ym=movmean(data_IOD.assm_ym, 4, 2, 'Endpoints', 'discard');
                data_IOD.lens2_4ym=movmean(data_IOD.lens2_ym, 4, 2, 'Endpoints', 'discard');
                data_IOD.hcst_4ym=squeeze(mean(data_IOD.hcst_ym(2:5,:,:),1));
            else
                tmp.reshp=reshape(data_IOD.obs, [12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(month1,:),1));
                tmp.sm_m2=squeeze(sum(tmp.reshp(month2-12,2:end),1));  tmp.sm_m2(end+1)=NaN;
                data_IOD.obs_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_IOD.obs_4ym=movmean(data_IOD.obs_ym, 4, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_IOD.assm, [len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_IOD.assm_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_IOD.assm_4ym=movmean(data_IOD.assm_ym, 4, 2, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_IOD.lens2, [len_mem_lens2, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_IOD.lens2_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_IOD.lens2_4ym=movmean(data_IOD.lens2_ym, 4, 2, 'Endpoints', 'discard');
    
                tmp.reshp=reshape(data_IOD.hcst, [5, len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,:,month1,:),3));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,:,month2-12,:),3));  
                data_IOD.hcst_ym=(tmp.sm_m1(1:4,:,:)+tmp.sm_m2(2:5,:,:))/length(months);
                data_IOD.hcst_ym(end+1,:,:)=NaN;
                data_IOD.hcst_4ym=squeeze(mean(data_IOD.hcst_ym(2:5,:,:),1));
            end

            %% IPO
            if sum(months>12)==0
                tmp.reshp=reshape(data_IPO.obs, [12, len_m/12]);
                data_IPO.obs_ym=squeeze(mean(tmp.reshp(months,:),1));
                
                tmp.reshp=reshape(data_IPO.assm, [len_mem_oda, 12, len_m/12]);
                data_IPO.assm_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                
                tmp.reshp=reshape(data_IPO.lens2, [len_mem_lens2, 12, len_m/12]);
                data_IPO.lens2_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                
                tmp.reshp=reshape(data_IPO.hcst, [5, len_mem_oda, 12, len_m/12]);
                data_IPO.hcst_ym=squeeze(mean(tmp.reshp(:,:,months,:),3));
        
                data_IPO.obs_4ym=movmean(data_IPO.obs_ym, 4, 'Endpoints', 'discard');
                data_IPO.assm_4ym=movmean(data_IPO.assm_ym, 4, 2, 'Endpoints', 'discard');
                data_IPO.lens2_4ym=movmean(data_IPO.lens2_ym, 4, 2, 'Endpoints', 'discard');
                data_IPO.hcst_4ym=squeeze(mean(data_IPO.hcst_ym(2:5,:,:),1));
            else
                tmp.reshp=reshape(data_IPO.obs, [12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(month1,:),1));
                tmp.sm_m2=squeeze(sum(tmp.reshp(month2-12,2:end),1));  tmp.sm_m2(end+1)=NaN;
                data_IPO.obs_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_IPO.obs_4ym=movmean(data_IPO.obs_ym, 4, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_IPO.assm, [len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_IPO.assm_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_IPO.assm_4ym=movmean(data_IPO.assm_ym, 4, 2, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_IPO.lens2, [len_mem_lens2, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_IPO.lens2_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_IPO.lens2_4ym=movmean(data_IPO.lens2_ym, 4, 2, 'Endpoints', 'discard');
    
                tmp.reshp=reshape(data_IPO.hcst, [5, len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,:,month1,:),3));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,:,month2-12,:),3));  
                data_IPO.hcst_ym=(tmp.sm_m1(1:4,:,:)+tmp.sm_m2(2:5,:,:))/length(months);
                data_IPO.hcst_ym(end+1,:,:)=NaN;
                data_IPO.hcst_4ym=squeeze(mean(data_IPO.hcst_ym(2:5,:,:),1));
            end    
    
            %% PDO
            lmode=1;
            if sum(months>12)==0
                tmp.reshp=reshape(data_PDO.pct_obs(:,lmode), [12, len_m/12]);
                data_PDO.obs_ym=squeeze(mean(tmp.reshp(months,:),1));
                
                tmp.reshp=reshape(data_PDO.pct_assm(:,:,lmode), [len_mem_oda, 12, len_m/12]);
                data_PDO.assm_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                
                tmp.reshp=reshape(data_PDO.pct_lens2(:,:,lmode), [len_mem_lens2, 12, len_m/12]);
                data_PDO.lens2_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                
                tmp.reshp=reshape(data_PDO.pct_hcst(:,:,:,lmode), [5, len_mem_oda, 12, len_m/12]);
                data_PDO.hcst_ym=squeeze(mean(tmp.reshp(:,:,months,:),3)); 
            else
                tmp.reshp=reshape(data_PDO.pct_obs(:,lmode), [12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(month1,:),1));
                tmp.sm_m2=squeeze(sum(tmp.reshp(month2-12,2:end),1));  tmp.sm_m2(end+1)=NaN;
                data_PDO.obs_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_PDO.obs_4ym=movmean(data_PDO.obs_ym, 4, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_PDO.pct_assm(:,:,lmode), [len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_PDO.assm_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_PDO.assm_4ym=movmean(data_PDO.assm_ym, 4, 2, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_PDO.pct_lens2(:,:,lmode), [len_mem_lens2, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_PDO.lens2_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_PDO.lens2_4ym=movmean(data_PDO.lens2_ym, 4, 2, 'Endpoints', 'discard');
    
                tmp.reshp=reshape(data_PDO.pct_hcst(:,:,:,lmode), [5, len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,:,month1,:),3));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,:,month2-12,:),3));  
                data_PDO.hcst_ym=(tmp.sm_m1(1:4,:,:)+tmp.sm_m2(2:5,:,:))/length(months);
                data_PDO.hcst_ym(end+1,:,:)=NaN;
                data_PDO.hcst_4ym=squeeze(mean(data_PDO.hcst_ym(2:5,:,:),1));
            end
            
            
            %% x: 50~70 y=10:25 sign + -> negative change
            %% PDO sign correction
            xrange=50:70;
            yrange=10:25;
            %% obs
            tmp.sign=data_PDO.lv_obs(xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
            if tmp.sign>0
                data_PDO.obs_ym=-data_PDO.obs_ym;
            end
            %% assm
            for mi=1:len_mem_oda
                tmp.sign=data_PDO.lv_assm(mi,xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
                if tmp.sign>0
        %             disp('0')
                    data_PDO.assm_ym(mi,:)=-data_PDO.assm_ym(mi,:);
                end
            end
            %% lens2
            for mi=1:len_mem_lens2
                tmp.sign=data_PDO.lv_lens2(mi,xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
                if tmp.sign>0
                    data_PDO.lens2_ym(mi,:)=-data_PDO.lens2_ym(mi,:);
                end
            end
            %% hcst
            for ly=1:5
                for mi=1:len_mem_oda
                    tmp.sign=data_PDO.lv_hcst(ly,mi,xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
                    if tmp.sign>0
                        data_PDO.hcst_ym(ly,mi,:)=-data_PDO.hcst_ym(ly,mi,:);
                    end
                end
            end
    
            data_PDO.obs_4ym=movmean(data_PDO.obs_ym, 4, 'Endpoints', 'discard');
            data_PDO.assm_4ym=movmean(data_PDO.assm_ym, 4, 2, 'Endpoints', 'discard');
            data_PDO.lens2_4ym=movmean(data_PDO.lens2_ym, 4, 2, 'Endpoints', 'discard');
            data_PDO.hcst_4ym=squeeze(mean(data_PDO.hcst_ym(2:5,:,:),1));
    
    
            %% PDO_l
            lmode=1;
            if sum(months>12)==0
            tmp.reshp=reshape(data_PDO_l.pct_obs(:,lmode), [12, len_m/12]);
            data_PDO_l.obs_ym=squeeze(mean(tmp.reshp(months,:),1));
            
            tmp.reshp=reshape(data_PDO_l.pct_assm(:,:,lmode), [len_mem_oda, 12, len_m/12]);
            data_PDO_l.assm_ym=squeeze(mean(tmp.reshp(:,months,:),2));
            
            tmp.reshp=reshape(data_PDO_l.pct_lens2(:,:,lmode), [len_mem_lens2, 12, len_m/12]);
            data_PDO_l.lens2_ym=squeeze(mean(tmp.reshp(:,months,:),2));
            
            tmp.reshp=reshape(data_PDO_l.pct_hcst(:,:,:,lmode), [5, len_mem_oda, 12, len_m/12]);
            data_PDO_l.hcst_ym=squeeze(mean(tmp.reshp(:,:,months,:),3)); 
            else
                tmp.reshp=reshape(data_PDO_l.pct_obs(:,lmode), [12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(month1,:),1));
                tmp.sm_m2=squeeze(sum(tmp.reshp(month2-12,2:end),1));  tmp.sm_m2(end+1)=NaN;
                data_PDO_l.obs_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_PDO_l.obs_4ym=movmean(data_PDO_l.obs_ym, 4, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_PDO_l.pct_assm(:,:,lmode), [len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_PDO_l.assm_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_PDO_l.assm_4ym=movmean(data_PDO_l.assm_ym, 4, 2, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_PDO_l.pct_lens2(:,:,lmode), [len_mem_lens2, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_PDO_l.lens2_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_PDO_l.lens2_4ym=movmean(data_PDO_l.lens2_ym, 4, 2, 'Endpoints', 'discard');
    
                tmp.reshp=reshape(data_PDO_l.pct_hcst(:,:,:,lmode), [5, len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,:,month1,:),3));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,:,month2-12,:),3));  
                data_PDO_l.hcst_ym=(tmp.sm_m1(1:4,:,:)+tmp.sm_m2(2:5,:,:))/length(months);
                data_PDO_l.hcst_ym(end+1,:,:)=NaN;
                data_PDO_l.hcst_4ym=squeeze(mean(data_PDO_l.hcst_ym(2:5,:,:),1));
            end
            
            
            %% x: 50~70 y=10:25 sign + -> negative change
            %% PDO_l sign correction
            xrange=50:70;
            yrange=10:25;
            %% obs
            tmp.sign=data_PDO_l.lv_obs(xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
            if tmp.sign>0
                data_PDO_l.obs_ym=-data_PDO_l.obs_ym;
            end
            %% assm
            for mi=1:len_mem_oda
                tmp.sign=data_PDO_l.lv_assm(mi,xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
                if tmp.sign>0
        %             disp('0')
                    data_PDO_l.assm_ym(mi,:)=-data_PDO_l.assm_ym(mi,:);
                end
            end
    %         figure;
    %         for mi=1:len_mem_oda
    %             hold on
    %             plot(data_PDO_l.assm_ym(mi,:));
    %         end
            %% lens2
            for mi=1:len_mem_lens2
                tmp.sign=data_PDO_l.lv_lens2(mi,xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
                if tmp.sign>0
                    data_PDO_l.lens2_ym(mi,:)=-data_PDO_l.lens2_ym(mi,:);
                end
            end
            %% hcst
            for ly=1:5
                for mi=1:len_mem_oda
                    tmp.sign=data_PDO_l.lv_hcst(ly,mi,xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
                    if tmp.sign>0
                        data_PDO_l.hcst_ym(ly,mi,:)=-data_PDO_l.hcst_ym(ly,mi,:);
                    end
                end
            end
    
            data_PDO_l.obs_4ym=movmean(data_PDO_l.obs_ym, 4, 'Endpoints', 'discard');
            data_PDO_l.assm_4ym=movmean(data_PDO_l.assm_ym, 4, 2, 'Endpoints', 'discard');
            data_PDO_l.lens2_4ym=movmean(data_PDO_l.lens2_ym, 4, 2, 'Endpoints', 'discard');
            data_PDO_l.hcst_4ym=squeeze(mean(data_PDO_l.hcst_ym(2:5,:,:),1));
    
    
        
            %% NPGO_SST
            lmode=2;
            if sum(months>12)==0
                tmp.reshp=reshape(data_PDO.pct_obs(:,lmode), [12, len_m/12]);
                data_NPGO_SST.obs_ym=squeeze(mean(tmp.reshp(months,:),1));
                
                tmp.reshp=reshape(data_PDO.pct_assm(:,:,lmode), [len_mem_oda, 12, len_m/12]);
                data_NPGO_SST.assm_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                
                tmp.reshp=reshape(data_PDO.pct_lens2(:,:,lmode), [len_mem_lens2, 12, len_m/12]);
                data_NPGO_SST.lens2_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                
                tmp.reshp=reshape(data_PDO.pct_hcst(:,:,:,lmode), [5, len_mem_oda, 12, len_m/12]);
                data_NPGO_SST.hcst_ym=squeeze(mean(tmp.reshp(:,:,months,:),3));
            else
                tmp.reshp=reshape(data_PDO.pct_obs(:,lmode), [12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(month1,:),1));
                tmp.sm_m2=squeeze(sum(tmp.reshp(month2-12,2:end),1));  tmp.sm_m2(end+1)=NaN;
                data_NPGO_SST.obs_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_NPGO_SST.obs_4ym=movmean(data_NPGO_SST.obs_ym, 4, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_PDO.pct_assm(:,:,lmode), [len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_NPGO_SST.assm_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_NPGO_SST.assm_4ym=movmean(data_NPGO_SST.assm_ym, 4, 2, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_PDO.pct_lens2(:,:,lmode), [len_mem_lens2, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_NPGO_SST.lens2_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_NPGO_SST.lens2_4ym=movmean(data_NPGO_SST.lens2_ym, 4, 2, 'Endpoints', 'discard');
    
                tmp.reshp=reshape(data_PDO.pct_hcst(:,:,:,lmode), [5, len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,:,month1,:),3));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,:,month2-12,:),3));  
                data_NPGO_SST.hcst_ym=(tmp.sm_m1(1:4,:,:)+tmp.sm_m2(2:5,:,:))/length(months);
                data_NPGO_SST.hcst_ym(end+1,:,:)=NaN;
                data_NPGO_SST.hcst_4ym=squeeze(mean(data_NPGO_SST.hcst_ym(2:5,:,:),1));
            end
    
            
            %% x: 70~90 y=20:30 sign + -> negative change
            xrange=40:60;
            yrange=1:15;
            %% obs
            tmp.sign=data_PDO.lv_obs(xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
            if tmp.sign<0
                data_NPGO_SST.obs_ym=-data_NPGO_SST.obs_ym;
            end
            %% assm
            for mi=1:len_mem_oda
                tmp.sign=data_PDO.lv_assm(mi,xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
                if tmp.sign<0
                    data_NPGO_SST.assm_ym(mi,:)=-data_NPGO_SST.assm_ym(mi,:);
                end
            end
            %% lens2
            for mi=1:len_mem_lens2
                tmp.sign=data_PDO.lv_lens2(mi,xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
                if tmp.sign<0
                    data_NPGO_SST.lens2_ym(mi,:)=-data_NPGO_SST.lens2_ym(mi,:);
                end
            end
            %% hcst
            for ly=1:5
                for mi=1:len_mem_oda
                    tmp.sign=data_PDO.lv_hcst(ly,mi,xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
                    if tmp.sign<0
                        data_NPGO_SST.hcst_ym(ly,mi,:)=-data_NPGO_SST.hcst_ym(ly,mi,:);
                    end
                end
            end
    
            data_NPGO_SST.obs_4ym=movmean(data_NPGO_SST.obs_ym, 4, 'Endpoints', 'discard');
            data_NPGO_SST.assm_4ym=movmean(data_NPGO_SST.assm_ym, 4, 2, 'Endpoints', 'discard');
            data_NPGO_SST.lens2_4ym=movmean(data_NPGO_SST.lens2_ym, 4, 2, 'Endpoints', 'discard');
            data_NPGO_SST.hcst_4ym=squeeze(mean(data_NPGO_SST.hcst_ym(2:5,:,:),1));
    
    
            %% NPGO_l_SST
            lmode=2;
            if sum(months>12)==0
                tmp.reshp=reshape(data_PDO_l.pct_obs(:,lmode), [12, len_m/12]);
                data_NPGO_l_SST.obs_ym=squeeze(mean(tmp.reshp(months,:),1));
                
                tmp.reshp=reshape(data_PDO_l.pct_assm(:,:,lmode), [len_mem_oda, 12, len_m/12]);
                data_NPGO_l_SST.assm_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                
                tmp.reshp=reshape(data_PDO_l.pct_lens2(:,:,lmode), [len_mem_lens2, 12, len_m/12]);
                data_NPGO_l_SST.lens2_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                
                tmp.reshp=reshape(data_PDO_l.pct_hcst(:,:,:,lmode), [5, len_mem_oda, 12, len_m/12]);
                data_NPGO_l_SST.hcst_ym=squeeze(mean(tmp.reshp(:,:,months,:),3));
            else
                tmp.reshp=reshape(data_PDO_l.pct_obs(:,lmode), [12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(month1,:),1));
                tmp.sm_m2=squeeze(sum(tmp.reshp(month2-12,2:end),1));  tmp.sm_m2(end+1)=NaN;
                data_NPGO_l_SST.obs_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_NPGO_l_SST.obs_4ym=movmean(data_NPGO_l_SST.obs_ym, 4, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_PDO_l.pct_assm(:,:,lmode), [len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_NPGO_l_SST.assm_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_NPGO_l_SST.assm_4ym=movmean(data_NPGO_l_SST.assm_ym, 4, 2, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_PDO_l.pct_lens2(:,:,lmode), [len_mem_lens2, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_NPGO_l_SST.lens2_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_NPGO_l_SST.lens2_4ym=movmean(data_NPGO_l_SST.lens2_ym, 4, 2, 'Endpoints', 'discard');
    
                tmp.reshp=reshape(data_PDO_l.pct_hcst(:,:,:,lmode), [5, len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,:,month1,:),3));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,:,month2-12,:),3));  
                data_NPGO_l_SST.hcst_ym=(tmp.sm_m1(1:4,:,:)+tmp.sm_m2(2:5,:,:))/length(months);
                data_NPGO_l_SST.hcst_ym(end+1,:,:)=NaN;
                data_NPGO_l_SST.hcst_4ym=squeeze(mean(data_NPGO_l_SST.hcst_ym(2:5,:,:),1));
            end
            
            %% x: 70~90 y=20:30 sign + -> negative change
    %         xrange=70:90;
    %         yrange=20:30;
            xrange=40:60;
            yrange=1:15;
            %% obs
            tmp.sign=data_PDO_l.lv_obs(xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
            if tmp.sign<0
                data_NPGO_l_SST.obs_ym=-data_NPGO_l_SST.obs_ym;
            end
            %% assm
            for mi=1:len_mem_oda
                tmp.sign=data_PDO_l.lv_assm(mi,xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
                if tmp.sign<0
                    data_NPGO_l_SST.assm_ym(mi,:)=-data_NPGO_l_SST.assm_ym(mi,:);
                end
            end
    %         pcolor(squeeze(data_PDO_l.lv_assm(1,:,:,2))'); shading flat; colorbar;
    %         figure;
    %         for mi=1:len_mem_oda
    %             hold on
    %             plot(data_NPGO_l_SST.assm_ym(mi,:));
    %         end
    
            %% lens2
            for mi=1:len_mem_lens2
                tmp.sign=data_PDO_l.lv_lens2(mi,xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
                if tmp.sign<0
                    data_NPGO_l_SST.lens2_ym(mi,:)=-data_NPGO_l_SST.lens2_ym(mi,:);
                end
            end
    
    %         figure;
    %         for mi=1:len_mem_lens2
    %             hold on
    %             plot(data_NPGO_l_SST.lens2_ym(mi,:));
    %         end
    
            %% hcst
            for ly=1:5
                for mi=1:len_mem_oda
                    tmp.sign=data_PDO_l.lv_hcst(ly,mi,xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
                    if tmp.sign<0
                        data_NPGO_l_SST.hcst_ym(ly,mi,:)=-data_NPGO_l_SST.hcst_ym(ly,mi,:);
                    end
                end
            end
    
            data_NPGO_l_SST.obs_4ym=movmean(data_NPGO_l_SST.obs_ym, 4, 'Endpoints', 'discard');
            data_NPGO_l_SST.assm_4ym=movmean(data_NPGO_l_SST.assm_ym, 4, 2, 'Endpoints', 'discard');
            data_NPGO_l_SST.lens2_4ym=movmean(data_NPGO_l_SST.lens2_ym, 4, 2, 'Endpoints', 'discard');
            data_NPGO_l_SST.hcst_4ym=squeeze(mean(data_NPGO_l_SST.hcst_ym(2:5,:,:),1));
    
    
        
            %% NPGO_SSH
            lmode=2;
            tmp.pct_obs=NaN(732,3);
            tmp.pct_obs(397:732,:)=data_NPGO.data_PDO.pct_obs(:,:);
            if sum(months>12)==0
                tmp.reshp=reshape(tmp.pct_obs(:,lmode), [12, len_m/12]);
                data_NPGO_SSH.obs_ym=squeeze(mean(tmp.reshp(months,:),1));
                
                tmp.reshp=reshape(data_NPGO.data_PDO.pct_assm(:,:,lmode), [len_mem_oda, 12, len_m/12]);
                data_NPGO_SSH.assm_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                
                tmp.reshp=reshape(data_NPGO.data_PDO.pct_lens2(:,:,lmode), [len_mem_lens2, 12, len_m/12]);
                data_NPGO_SSH.lens2_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                
                tmp.reshp=reshape(data_NPGO.data_PDO.pct_hcst(:,:,:,lmode), [5, len_mem_oda, 12, len_m/12]);
                data_NPGO_SSH.hcst_ym=squeeze(mean(tmp.reshp(:,:,months,:),3));
            else
                tmp.reshp=reshape(tmp.pct_obs(:,lmode), [12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(month1,:),1));
                tmp.sm_m2=squeeze(sum(tmp.reshp(month2-12,2:end),1));  tmp.sm_m2(end+1)=NaN;
                data_NPGO_SSH.obs_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_NPGO_SSH.obs_4ym=movmean(data_NPGO_SSH.obs_ym, 4, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_NPGO.data_PDO.pct_assm(:,:,lmode), [len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_NPGO_SSH.assm_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_NPGO_SSH.assm_4ym=movmean(data_NPGO_SSH.assm_ym, 4, 2, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_NPGO.data_PDO.pct_lens2(:,:,lmode), [len_mem_lens2, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_NPGO_SSH.lens2_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_NPGO_SSH.lens2_4ym=movmean(data_NPGO_SSH.lens2_ym, 4, 2, 'Endpoints', 'discard');
    
                tmp.reshp=reshape(data_NPGO.data_PDO.pct_hcst(:,:,:,lmode), [5, len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,:,month1,:),3));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,:,month2-12,:),3));  
                data_NPGO_SSH.hcst_ym=(tmp.sm_m1(1:4,:,:)+tmp.sm_m2(2:5,:,:))/length(months);
                data_NPGO_SSH.hcst_ym(end+1,:,:)=NaN;
                data_NPGO_SSH.hcst_4ym=squeeze(mean(data_NPGO_SSH.hcst_ym(2:5,:,:),1));
            end
        
            %% x: 70~90 y=20:30 sign - -> positive change
            xrange=20:40;
            yrange=1:20;
            %% obs
            tmp.sign=data_NPGO.data_PDO.lv_obs(xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
            if tmp.sign<0
                data_NPGO_SSH.obs_ym=-data_NPGO_SSH.obs_ym;
            end
            %% assm
            for mi=1:len_mem_oda
                tmp.sign=data_NPGO.data_PDO.lv_assm(mi,xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
                if tmp.sign<0
                    data_NPGO_SSH.assm_ym(mi,:)=-data_NPGO_SSH.assm_ym(mi,:);
                end
            end
            
    %         figure;
    %         for mi=1:len_mem_oda
    %             hold on
    %             plot(data_NPGO_SSH.assm_ym(mi,:));
    %         end
    
            %% lens2
            for mi=1:len_mem_lens2
                tmp.sign=data_NPGO.data_PDO.lv_lens2(mi,xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
                if tmp.sign<0
                    data_NPGO_SSH.lens2_ym(mi,:)=-data_NPGO_SSH.lens2_ym(mi,:);
                end
            end
    
    %         figure;
    %         for mi=1:len_mem_lens2
    %             hold on
    %             plot(data_NPGO_SSH.lens2_ym(mi,:));
    %         end
    
            %% hcst
            for ly=1:5
                for mi=1:len_mem_oda
                    tmp.sign=data_NPGO.data_PDO.lv_hcst(ly,mi,xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
                    if tmp.sign<0
                        data_NPGO_SSH.hcst_ym(ly,mi,:)=-data_NPGO_SSH.hcst_ym(ly,mi,:);
                    end
                end
            end
    
            data_NPGO_SSH.obs_4ym=movmean(data_NPGO_SSH.obs_ym, 4, 'Endpoints', 'discard');
            data_NPGO_SSH.assm_4ym=movmean(data_NPGO_SSH.assm_ym, 4, 2, 'Endpoints', 'discard');
            data_NPGO_SSH.lens2_4ym=movmean(data_NPGO_SSH.lens2_ym, 4, 2, 'Endpoints', 'discard');
            data_NPGO_SSH.hcst_4ym=squeeze(mean(data_NPGO_SSH.hcst_ym(2:5,:,:),1));
            
            %% SAM
            if sum(months>12)==0
                tmp.reshp=reshape(data_SAM.obs, [12, len_m/12]);
                data_SAM.obs_ym=squeeze(mean(tmp.reshp(months,:),1));
                
                tmp.reshp=reshape(data_SAM.assm, [len_mem_oda, 12, len_m/12]);
                data_SAM.assm_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                
                tmp.reshp=reshape(data_SAM.lens2, [len_mem_lens2, 12, len_m/12]);
                data_SAM.lens2_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                
                tmp.reshp=reshape(data_SAM.hcst, [5, len_mem_oda, 12, len_m/12]);
                data_SAM.hcst_ym=squeeze(mean(tmp.reshp(:,:,months,:),3));
        
                data_SAM.obs_4ym=movmean(data_SAM.obs_ym, 4, 'Endpoints', 'discard');
                data_SAM.assm_4ym=movmean(data_SAM.assm_ym, 4, 2, 'Endpoints', 'discard');
                data_SAM.lens2_4ym=movmean(data_SAM.lens2_ym, 4, 2, 'Endpoints', 'discard');
                data_SAM.hcst_4ym=squeeze(mean(data_SAM.hcst_ym(2:5,:,:),1));
            else
                tmp.reshp=reshape(data_SAM.obs, [12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(month1,:),1));
                tmp.sm_m2=squeeze(sum(tmp.reshp(month2-12,2:end),1));  tmp.sm_m2(end+1)=NaN;
                data_SAM.obs_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_SAM.obs_4ym=movmean(data_SAM.obs_ym, 4, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_SAM.assm, [len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_SAM.assm_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_SAM.assm_4ym=movmean(data_SAM.assm_ym, 4, 2, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_SAM.lens2, [len_mem_lens2, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_SAM.lens2_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_SAM.lens2_4ym=movmean(data_SAM.lens2_ym, 4, 2, 'Endpoints', 'discard');
    
                tmp.reshp=reshape(data_SAM.hcst, [5, len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,:,month1,:),3));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,:,month2-12,:),3));  
                data_SAM.hcst_ym=(tmp.sm_m1(1:4,:,:)+tmp.sm_m2(2:5,:,:))/length(months);
                data_SAM.hcst_ym(end+1,:,:)=NaN;
                data_SAM.hcst_4ym=squeeze(mean(data_SAM.hcst_ym(2:5,:,:),1));
            end
    
            %% SAM_d
            if sum(months>12)==0
                tmp.reshp=reshape(data_SAM_d.obs, [12, len_m/12]);
                data_SAM_d.obs_ym=squeeze(mean(tmp.reshp(months,:),1));
                
                tmp.reshp=reshape(data_SAM_d.assm, [len_mem_oda, 12, len_m/12]);
                data_SAM_d.assm_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                
                tmp.reshp=reshape(data_SAM_d.lens2, [len_mem_lens2, 12, len_m/12]);
                data_SAM_d.lens2_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                
                tmp.reshp=reshape(data_SAM_d.hcst, [5, len_mem_oda, 12, len_m/12]);
                data_SAM_d.hcst_ym=squeeze(mean(tmp.reshp(:,:,months,:),3));
        
                data_SAM_d.obs_4ym=movmean(data_SAM_d.obs_ym, 4, 'Endpoints', 'discard');
                data_SAM_d.assm_4ym=movmean(data_SAM_d.assm_ym, 4, 2, 'Endpoints', 'discard');
                data_SAM_d.lens2_4ym=movmean(data_SAM_d.lens2_ym, 4, 2, 'Endpoints', 'discard');
                data_SAM_d.hcst_4ym=squeeze(mean(data_SAM_d.hcst_ym(2:5,:,:),1));
            else
                tmp.reshp=reshape(data_SAM_d.obs, [12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(month1,:),1));
                tmp.sm_m2=squeeze(sum(tmp.reshp(month2-12,2:end),1));  tmp.sm_m2(end+1)=NaN;
                data_SAM_d.obs_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_SAM_d.obs_4ym=movmean(data_SAM_d.obs_ym, 4, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_SAM_d.assm, [len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_SAM_d.assm_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_SAM_d.assm_4ym=movmean(data_SAM_d.assm_ym, 4, 2, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_SAM_d.lens2, [len_mem_lens2, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_SAM_d.lens2_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_SAM_d.lens2_4ym=movmean(data_SAM_d.lens2_ym, 4, 2, 'Endpoints', 'discard');
    
                tmp.reshp=reshape(data_SAM_d.hcst, [5, len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,:,month1,:),3));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,:,month2-12,:),3));  
                data_SAM_d.hcst_ym=(tmp.sm_m1(1:4,:,:)+tmp.sm_m2(2:5,:,:))/length(months);
                data_SAM_d.hcst_ym(end+1,:,:)=NaN;
                data_SAM_d.hcst_4ym=squeeze(mean(data_SAM_d.hcst_ym(2:5,:,:),1));
            end
        
            %% ATL3
            if sum(months>12)==0
                tmp.reshp=reshape(data_ATL3.obs, [12, len_m/12]);
                data_ATL3.obs_ym=squeeze(mean(tmp.reshp(months,:),1));
                
                tmp.reshp=reshape(data_ATL3.assm, [len_mem_oda, 12, len_m/12]);
                data_ATL3.assm_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                
                tmp.reshp=reshape(data_ATL3.lens2, [len_mem_lens2, 12, len_m/12]);
                data_ATL3.lens2_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                
                tmp.reshp=reshape(data_ATL3.hcst, [5, len_mem_oda, 12, len_m/12]);
                data_ATL3.hcst_ym=squeeze(mean(tmp.reshp(:,:,months,:),3));
        
                data_ATL3.obs_4ym=movmean(data_ATL3.obs_ym, 4, 'Endpoints', 'discard');
                data_ATL3.assm_4ym=movmean(data_ATL3.assm_ym, 4, 2, 'Endpoints', 'discard');
                data_ATL3.lens2_4ym=movmean(data_ATL3.lens2_ym, 4, 2, 'Endpoints', 'discard');
                data_ATL3.hcst_4ym=squeeze(mean(data_ATL3.hcst_ym(2:5,:,:),1));
            else
                tmp.reshp=reshape(data_ATL3.obs, [12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(month1,:),1));
                tmp.sm_m2=squeeze(sum(tmp.reshp(month2-12,2:end),1));  tmp.sm_m2(end+1)=NaN;
                data_ATL3.obs_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_ATL3.obs_4ym=movmean(data_ATL3.obs_ym, 4, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_ATL3.assm, [len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_ATL3.assm_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_ATL3.assm_4ym=movmean(data_ATL3.assm_ym, 4, 2, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_ATL3.lens2, [len_mem_lens2, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_ATL3.lens2_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_ATL3.lens2_4ym=movmean(data_ATL3.lens2_ym, 4, 2, 'Endpoints', 'discard');
    
                tmp.reshp=reshape(data_ATL3.hcst, [5, len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,:,month1,:),3));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,:,month2-12,:),3));  
                data_ATL3.hcst_ym=(tmp.sm_m1(1:4,:,:)+tmp.sm_m2(2:5,:,:))/length(months);
                data_ATL3.hcst_ym(end+1,:,:)=NaN;
                data_ATL3.hcst_4ym=squeeze(mean(data_ATL3.hcst_ym(2:5,:,:),1));
            end
    
            %% ATL3_d
            if sum(months>12)==0
                tmp.reshp=reshape(data_ATL3_d.obs, [12, len_m/12]);
                data_ATL3_d.obs_ym=squeeze(mean(tmp.reshp(months,:),1));
                
                tmp.reshp=reshape(data_ATL3_d.assm, [len_mem_oda, 12, len_m/12]);
                data_ATL3_d.assm_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                
                tmp.reshp=reshape(data_ATL3_d.lens2, [len_mem_lens2, 12, len_m/12]);
                data_ATL3_d.lens2_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                
                tmp.reshp=reshape(data_ATL3_d.hcst, [5, len_mem_oda, 12, len_m/12]);
                data_ATL3_d.hcst_ym=squeeze(mean(tmp.reshp(:,:,months,:),3));
        
                data_ATL3_d.obs_4ym=movmean(data_ATL3_d.obs_ym, 4, 'Endpoints', 'discard');
                data_ATL3_d.assm_4ym=movmean(data_ATL3_d.assm_ym, 4, 2, 'Endpoints', 'discard');
                data_ATL3_d.lens2_4ym=movmean(data_ATL3_d.lens2_ym, 4, 2, 'Endpoints', 'discard');
                data_ATL3_d.hcst_4ym=squeeze(mean(data_ATL3_d.hcst_ym(2:5,:,:),1));
            else
                tmp.reshp=reshape(data_ATL3_d.obs, [12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(month1,:),1));
                tmp.sm_m2=squeeze(sum(tmp.reshp(month2-12,2:end),1));  tmp.sm_m2(end+1)=NaN;
                data_ATL3_d.obs_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_ATL3_d.obs_4ym=movmean(data_ATL3_d.obs_ym, 4, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_ATL3_d.assm, [len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_ATL3_d.assm_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_ATL3_d.assm_4ym=movmean(data_ATL3_d.assm_ym, 4, 2, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_ATL3_d.lens2, [len_mem_lens2, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_ATL3_d.lens2_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_ATL3_d.lens2_4ym=movmean(data_ATL3_d.lens2_ym, 4, 2, 'Endpoints', 'discard');
    
                tmp.reshp=reshape(data_ATL3_d.hcst, [5, len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,:,month1,:),3));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,:,month2-12,:),3));  
                data_ATL3_d.hcst_ym=(tmp.sm_m1(1:4,:,:)+tmp.sm_m2(2:5,:,:))/length(months);
                data_ATL3_d.hcst_ym(end+1,:,:)=NaN;
                data_ATL3_d.hcst_4ym=squeeze(mean(data_ATL3_d.hcst_ym(2:5,:,:),1));
            end
        
        
            %% TNA
            if sum(months>12)==0
                tmp.reshp=reshape(data_TNA.obs, [12, len_m/12]);
                data_TNA.obs_ym=squeeze(mean(tmp.reshp(months,:),1));
                
                tmp.reshp=reshape(data_TNA.assm, [len_mem_oda, 12, len_m/12]);
                data_TNA.assm_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                
                tmp.reshp=reshape(data_TNA.lens2, [len_mem_lens2, 12, len_m/12]);
                data_TNA.lens2_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                
                tmp.reshp=reshape(data_TNA.hcst, [5, len_mem_oda, 12, len_m/12]);
                data_TNA.hcst_ym=squeeze(mean(tmp.reshp(:,:,months,:),3));
        
                data_TNA.obs_4ym=movmean(data_TNA.obs_ym, 4, 'Endpoints', 'discard');
                data_TNA.assm_4ym=movmean(data_TNA.assm_ym, 4, 2, 'Endpoints', 'discard');
                data_TNA.lens2_4ym=movmean(data_TNA.lens2_ym, 4, 2, 'Endpoints', 'discard');
                data_TNA.hcst_4ym=squeeze(mean(data_TNA.hcst_ym(2:5,:,:),1));
            else
                tmp.reshp=reshape(data_TNA.obs, [12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(month1,:),1));
                tmp.sm_m2=squeeze(sum(tmp.reshp(month2-12,2:end),1));  tmp.sm_m2(end+1)=NaN;
                data_TNA.obs_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_TNA.obs_4ym=movmean(data_TNA.obs_ym, 4, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_TNA.assm, [len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_TNA.assm_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_TNA.assm_4ym=movmean(data_TNA.assm_ym, 4, 2, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_TNA.lens2, [len_mem_lens2, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_TNA.lens2_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_TNA.lens2_4ym=movmean(data_TNA.lens2_ym, 4, 2, 'Endpoints', 'discard');
    
                tmp.reshp=reshape(data_TNA.hcst, [5, len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,:,month1,:),3));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,:,month2-12,:),3));  
                data_TNA.hcst_ym=(tmp.sm_m1(1:4,:,:)+tmp.sm_m2(2:5,:,:))/length(months);
                data_TNA.hcst_ym(end+1,:,:)=NaN;
                data_TNA.hcst_4ym=squeeze(mean(data_TNA.hcst_ym(2:5,:,:),1));
            end
        
        
            %% TSA
            if sum(months>12)==0
                tmp.reshp=reshape(data_TSA.obs, [12, len_m/12]);
                data_TSA.obs_ym=squeeze(mean(tmp.reshp(months,:),1));
                
                tmp.reshp=reshape(data_TSA.assm, [len_mem_oda, 12, len_m/12]);
                data_TSA.assm_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                
                tmp.reshp=reshape(data_TSA.lens2, [len_mem_lens2, 12, len_m/12]);
                data_TSA.lens2_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                
                tmp.reshp=reshape(data_TSA.hcst, [5, len_mem_oda, 12, len_m/12]);
                data_TSA.hcst_ym=squeeze(mean(tmp.reshp(:,:,months,:),3));
        
                data_TSA.obs_4ym=movmean(data_TSA.obs_ym, 4, 'Endpoints', 'discard');
                data_TSA.assm_4ym=movmean(data_TSA.assm_ym, 4, 2, 'Endpoints', 'discard');
                data_TSA.lens2_4ym=movmean(data_TSA.lens2_ym, 4, 2, 'Endpoints', 'discard');
                data_TSA.hcst_4ym=squeeze(mean(data_TSA.hcst_ym(2:5,:,:),1));
            else
                tmp.reshp=reshape(data_TSA.obs, [12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(month1,:),1));
                tmp.sm_m2=squeeze(sum(tmp.reshp(month2-12,2:end),1));  tmp.sm_m2(end+1)=NaN;
                data_TSA.obs_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_TSA.obs_4ym=movmean(data_TSA.obs_ym, 4, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_TSA.assm, [len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_TSA.assm_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_TSA.assm_4ym=movmean(data_TSA.assm_ym, 4, 2, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_TSA.lens2, [len_mem_lens2, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_TSA.lens2_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_TSA.lens2_4ym=movmean(data_TSA.lens2_ym, 4, 2, 'Endpoints', 'discard');
    
                tmp.reshp=reshape(data_TSA.hcst, [5, len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,:,month1,:),3));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,:,month2-12,:),3));  
                data_TSA.hcst_ym=(tmp.sm_m1(1:4,:,:)+tmp.sm_m2(2:5,:,:))/length(months);
                data_TSA.hcst_ym(end+1,:,:)=NaN;
                data_TSA.hcst_4ym=squeeze(mean(data_TSA.hcst_ym(2:5,:,:),1));
            end
    
            %% KEI_O
            if sum(months>12)==0
                tmp.reshp=reshape(data_KEI_O.obs, [12, len_m/12]);
                data_KEI_O.obs_ym=squeeze(mean(tmp.reshp(months,:),1));
                
                tmp.reshp=reshape(data_KEI_O.assm, [len_mem_oda, 12, len_m/12]);
                data_KEI_O.assm_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                
                tmp.reshp=reshape(data_KEI_O.lens2, [len_mem_lens2, 12, len_m/12]);
                data_KEI_O.lens2_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                
                tmp.reshp=reshape(data_KEI_O.hcst, [5, len_mem_oda, 12, len_m/12]);
                data_KEI_O.hcst_ym=squeeze(mean(tmp.reshp(:,:,months,:),3));
        
                data_KEI_O.obs_4ym=movmean(data_KEI_O.obs_ym, 4, 'Endpoints', 'discard');
                data_KEI_O.assm_4ym=movmean(data_KEI_O.assm_ym, 4, 2, 'Endpoints', 'discard');
                data_KEI_O.lens2_4ym=movmean(data_KEI_O.lens2_ym, 4, 2, 'Endpoints', 'discard');
                data_KEI_O.hcst_4ym=squeeze(mean(data_KEI_O.hcst_ym(2:5,:,:),1));
            else
                tmp.reshp=reshape(data_KEI_O.obs, [12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(month1,:),1));
                tmp.sm_m2=squeeze(sum(tmp.reshp(month2-12,2:end),1));  tmp.sm_m2(end+1)=NaN;
                data_KEI_O.obs_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_KEI_O.obs_4ym=movmean(data_KEI_O.obs_ym, 4, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_KEI_O.assm, [len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_KEI_O.assm_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_KEI_O.assm_4ym=movmean(data_KEI_O.assm_ym, 4, 2, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_KEI_O.lens2, [len_mem_lens2, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_KEI_O.lens2_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_KEI_O.lens2_4ym=movmean(data_KEI_O.lens2_ym, 4, 2, 'Endpoints', 'discard');
    
                tmp.reshp=reshape(data_KEI_O.hcst, [5, len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,:,month1,:),3));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,:,month2-12,:),3));  
                data_KEI_O.hcst_ym=(tmp.sm_m1(1:4,:,:)+tmp.sm_m2(2:5,:,:))/length(months);
                data_KEI_O.hcst_ym(end+1,:,:)=NaN;
                data_KEI_O.hcst_4ym=squeeze(mean(data_KEI_O.hcst_ym(2:5,:,:),1));
            end
    
            %% KEI_LE
            if sum(months>12)==0
                tmp.reshp=reshape(data_KEI_LE.obs, [12, len_m/12]);
                data_KEI_LE.obs_ym=squeeze(mean(tmp.reshp(months,:),1));
                
                tmp.reshp=reshape(data_KEI_LE.assm, [len_mem_oda, 12, len_m/12]);
                data_KEI_LE.assm_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                
                tmp.reshp=reshape(data_KEI_LE.lens2, [len_mem_lens2, 12, len_m/12]);
                data_KEI_LE.lens2_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                
        %         figure;
        %         for mi=1:len_mem_lens2
        %             hold on
        %             plot(data_KEI_LE.lens2_ym(mi,:));
        %         end
        
                tmp.reshp=reshape(data_KEI_LE.hcst, [5, len_mem_oda, 12, len_m/12]);
                data_KEI_LE.hcst_ym=squeeze(mean(tmp.reshp(:,:,months,:),3));
        
                data_KEI_LE.obs_4ym=movmean(data_KEI_LE.obs_ym, 4, 'Endpoints', 'discard');
                data_KEI_LE.assm_4ym=movmean(data_KEI_LE.assm_ym, 4, 2, 'Endpoints', 'discard');
                data_KEI_LE.lens2_4ym=movmean(data_KEI_LE.lens2_ym, 4, 2, 'Endpoints', 'discard');
                data_KEI_LE.hcst_4ym=squeeze(mean(data_KEI_LE.hcst_ym(2:5,:,:),1));
            else
                tmp.reshp=reshape(data_KEI_LE.obs, [12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(month1,:),1));
                tmp.sm_m2=squeeze(sum(tmp.reshp(month2-12,2:end),1));  tmp.sm_m2(end+1)=NaN;
                data_KEI_LE.obs_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_KEI_LE.obs_4ym=movmean(data_KEI_LE.obs_ym, 4, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_KEI_LE.assm, [len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_KEI_LE.assm_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_KEI_LE.assm_4ym=movmean(data_KEI_LE.assm_ym, 4, 2, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_KEI_LE.lens2, [len_mem_lens2, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_KEI_LE.lens2_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_KEI_LE.lens2_4ym=movmean(data_KEI_LE.lens2_ym, 4, 2, 'Endpoints', 'discard');
    
                tmp.reshp=reshape(data_KEI_LE.hcst, [5, len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,:,month1,:),3));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,:,month2-12,:),3));  
                data_KEI_LE.hcst_ym=(tmp.sm_m1(1:4,:,:)+tmp.sm_m2(2:5,:,:))/length(months);
                data_KEI_LE.hcst_ym(end+1,:,:)=NaN;
                data_KEI_LE.hcst_4ym=squeeze(mean(data_KEI_LE.hcst_ym(2:5,:,:),1));
            end
    
            %% NAO_STA
            if sum(months>12)==0
                tmp.reshp=reshape(data_NAO_STA.obs, [12, len_m/12]);
                data_NAO_STA.obs_ym=squeeze(mean(tmp.reshp(months,:),1));
                
                tmp.reshp=reshape(data_NAO_STA.assm, [len_mem_oda, 12, len_m/12]);
                data_NAO_STA.assm_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                
                tmp.reshp=reshape(data_NAO_STA.lens2, [len_mem_lens2, 12, len_m/12]);
                data_NAO_STA.lens2_ym=squeeze(mean(tmp.reshp(:,months,:),2));
                
        %         figure;
        %         for mi=1:len_mem_lens2
        %             hold on
        %             plot(data_NAO_STA.lens2_ym(mi,:));
        %         end
        
                tmp.reshp=reshape(data_NAO_STA.hcst, [5, len_mem_oda, 12, len_m/12]);
                data_NAO_STA.hcst_ym=squeeze(mean(tmp.reshp(:,:,months,:),3));
        
                data_NAO_STA.obs_4ym=movmean(data_NAO_STA.obs_ym, 4, 'Endpoints', 'discard');
                data_NAO_STA.assm_4ym=movmean(data_NAO_STA.assm_ym, 4, 2, 'Endpoints', 'discard');
                data_NAO_STA.lens2_4ym=movmean(data_NAO_STA.lens2_ym, 4, 2, 'Endpoints', 'discard');
                data_NAO_STA.hcst_4ym=squeeze(mean(data_NAO_STA.hcst_ym(2:5,:,:),1));
            else
                tmp.reshp=reshape(data_NAO_STA.obs, [12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(month1,:),1));
                tmp.sm_m2=squeeze(sum(tmp.reshp(month2-12,2:end),1));  tmp.sm_m2(end+1)=NaN;
                data_NAO_STA.obs_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_NAO_STA.obs_4ym=movmean(data_NAO_STA.obs_ym, 4, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_NAO_STA.assm, [len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_NAO_STA.assm_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_NAO_STA.assm_4ym=movmean(data_NAO_STA.assm_ym, 4, 2, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_NAO_STA.lens2, [len_mem_lens2, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_NAO_STA.lens2_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_NAO_STA.lens2_4ym=movmean(data_NAO_STA.lens2_ym, 4, 2, 'Endpoints', 'discard');
    
                tmp.reshp=reshape(data_NAO_STA.hcst, [5, len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,:,month1,:),3));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,:,month2-12,:),3));  
                data_NAO_STA.hcst_ym=(tmp.sm_m1(1:4,:,:)+tmp.sm_m2(2:5,:,:))/length(months);
                data_NAO_STA.hcst_ym(end+1,:,:)=NaN;
                data_NAO_STA.hcst_4ym=squeeze(mean(data_NAO_STA.hcst_ym(2:5,:,:),1));
            end
    
    
    
            %% NAO_PC
            lmode=1;
            if sum(months>12)==0
            tmp.reshp=reshape(data_NAO_PC.pct_obs(:,lmode), [12, len_m/12]);
            data_NAO_PC.obs_ym=squeeze(mean(tmp.reshp(months,:),1));
            
            tmp.reshp=reshape(data_NAO_PC.pct_assm(:,:,lmode), [len_mem_oda, 12, len_m/12]);
            data_NAO_PC.assm_ym=squeeze(mean(tmp.reshp(:,months,:),2));
            
            tmp.reshp=reshape(data_NAO_PC.pct_lens2(:,:,lmode), [len_mem_lens2, 12, len_m/12]);
            data_NAO_PC.lens2_ym=squeeze(mean(tmp.reshp(:,months,:),2));
            
            tmp.reshp=reshape(data_NAO_PC.pct_hcst(:,:,:,lmode), [5, len_mem_oda, 12, len_m/12]);
            data_NAO_PC.hcst_ym=squeeze(mean(tmp.reshp(:,:,months,:),3)); 
            else
                tmp.reshp=reshape(data_NAO_PC.pct_obs(:,lmode), [12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(month1,:),1));
                tmp.sm_m2=squeeze(sum(tmp.reshp(month2-12,2:end),1));  tmp.sm_m2(end+1)=NaN;
                data_NAO_PC.obs_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_NAO_PC.obs_4ym=movmean(data_NAO_PC.obs_ym, 4, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_NAO_PC.pct_assm(:,:,lmode), [len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_NAO_PC.assm_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_NAO_PC.assm_4ym=movmean(data_NAO_PC.assm_ym, 4, 2, 'Endpoints', 'discard');
                
                tmp.reshp=reshape(data_NAO_PC.pct_lens2(:,:,lmode), [len_mem_lens2, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,month1,:),2));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,month2-12,2:end),2));  tmp.sm_m2(:,end+1)=NaN;
                data_NAO_PC.lens2_ym=(tmp.sm_m1+tmp.sm_m2)/length(months);
                data_NAO_PC.lens2_4ym=movmean(data_NAO_PC.lens2_ym, 4, 2, 'Endpoints', 'discard');
    
                tmp.reshp=reshape(data_NAO_PC.pct_hcst(:,:,:,lmode), [5, len_mem_oda, 12, len_m/12]);
                tmp.sm_m1=squeeze(sum(tmp.reshp(:,:,month1,:),3));
                tmp.sm_m2=squeeze(sum(tmp.reshp(:,:,month2-12,:),3));  
                data_NAO_PC.hcst_ym=(tmp.sm_m1(1:4,:,:)+tmp.sm_m2(2:5,:,:))/length(months);
                data_NAO_PC.hcst_ym(end+1,:,:)=NaN;
                data_NAO_PC.hcst_4ym=squeeze(mean(data_NAO_PC.hcst_ym(2:5,:,:),1));
            end
            
            
            %% x: 50~70 y=10:25 sign + -> negative change
            %% NAO_PC sign correction
            xrange=35:66;
            yrange=15:30;
            %% obs
            tmp.sign=data_NAO_PC.lv_obs(xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
            if tmp.sign<0
                data_NAO_PC.obs_ym=-data_NAO_PC.obs_ym;
            end
            %% assm
            for mi=1:len_mem_oda
                tmp.sign=data_NAO_PC.lv_assm(mi,xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
                if tmp.sign<0
        %             disp('0')
                    data_NAO_PC.assm_ym(mi,:)=-data_NAO_PC.assm_ym(mi,:);
                end
            end
    %         figure;
    %         for mi=1:len_mem_oda
    %             hold on
    %             plot(data_NAO_PC.assm_ym(mi,:));
    %         end
            %% lens2
            for mi=1:len_mem_lens2
                tmp.sign=data_NAO_PC.lv_lens2(mi,xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
                if tmp.sign<0
                    data_NAO_PC.lens2_ym(mi,:)=-data_NAO_PC.lens2_ym(mi,:);
                end
            end
            %% hcst
            for ly=1:5
                for mi=1:len_mem_oda
                    tmp.sign=data_NAO_PC.lv_hcst(ly,mi,xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
                    if tmp.sign<0
                        data_NAO_PC.hcst_ym(ly,mi,:)=-data_NAO_PC.hcst_ym(ly,mi,:);
                    end
                end
            end
    
            data_NAO_PC.obs_4ym=movmean(data_NAO_PC.obs_ym, 4, 'Endpoints', 'discard');
            data_NAO_PC.assm_4ym=movmean(data_NAO_PC.assm_ym, 4, 2, 'Endpoints', 'discard');
            data_NAO_PC.lens2_4ym=movmean(data_NAO_PC.lens2_ym, 4, 2, 'Endpoints', 'discard');
            data_NAO_PC.hcst_4ym=squeeze(mean(data_NAO_PC.hcst_ym(2:5,:,:),1));
    
            once=0;
        end
    
    
        
        %% corr between data and indices
        
    %     tmp.corr=corrcoef(squeeze(data_obs.([cfg_ts.var,'_ym'])), squeeze(data_AMO.obs_ym), 'Rows', 'complete');
        
    %     hold on
    %     for mi=1:20
    %         plot(data_PDO.assm_ym(mi,:));
    %     end
    
    % clim_indices={'AMO', 'ENSO', 'PDO', 'PDO_l', 'NPGO_SST', 'NPGO_l_SST', 'NPGO_SSH', 'SAM', 'ATL3', 'TNA', 'TSA', 'KEI_O', 'KEI_LE'};
    % clim_indices={'AMO', 'ENSO', 'PDO_l', 'NPGO_l_SST', 'NPGO_SSH', 'KEI_O', 'KEI_LE'};
    % clim_indices={'AMO', 'ENSO', 'IOD', 'PDO_l', 'NPGO_l_SST', 'NPGO_SSH', 'ATL3', 'TNA', 'TSA', 'NAO_PC', 'NAO_STA', 'SAM'};
    clim_indices={'IOD', 'ENSO', 'ATL3_d', 'NAO_STA', 'SAM_d', 'PDO_l', 'AMO', 'IPO'};
        
        %% corr_obs
        for ci=1:length(clim_indices)
            clim_indice=clim_indices{ci};
            eval_line= ['tmp.corr=corrcoef(squeeze(m_data_obs.([cfg_ts.var,', ...
            '''','_ym','''','])(:)), squeeze(data_', clim_indice, '.obs_ym(:)), ',...
            '''','Rows','''',', ','''','complete','''',');'];
            eval(eval_line);
            eval_line= ['corr_all.obs_ym_', clim_indice, '=tmp.corr(1,2);'];
            eval(eval_line);
        end
    
        
        %% 4y movmean
        tmp.data=squeeze(m_data_obs.([cfg_ts.var,'_ym'])(:));
        tmp.data_4ym=movmean(tmp.data,4,'Endpoints', 'discard');
            
        for ci=1:length(clim_indices)
            clim_indice=clim_indices{ci};
            eval_line= ['tmp.corr=corrcoef(tmp.data_4ym, squeeze(data_', clim_indice, '.obs_4ym(:)), ',...
            '''','Rows','''',', ','''','complete','''',');'];
            eval(eval_line);
            eval_line= ['corr_all.obs_4ym_', clim_indice, '=tmp.corr(1,2);'];
            eval(eval_line);
        end
    
    
    
        %% corr_assm
        for mi=1:len_mem_oda
    
            for ci=1:length(clim_indices)
                clim_indice=clim_indices{ci};
                eval_line= ['tmp.corr=corrcoef(squeeze(m_data_assm.([cfg_ts.var,', ...
                '''','_ym','''','])(mi,:)), squeeze(data_', clim_indice, '.assm_ym(mi,:)), ',...
                '''','Rows','''',', ','''','complete','''',');'];
                eval(eval_line);
                eval_line= ['corr_all.assm_ym_', clim_indice, '(mi)=tmp.corr(1,2);'];
                eval(eval_line);
            end
    
            
            %% 4y movmean
            tmp.data=squeeze(m_data_assm.([cfg_ts.var,'_ym'])(mi,:));
            tmp.data_4ym=movmean(tmp.data,4,'Endpoints', 'discard');
                
            for ci=1:length(clim_indices)
                clim_indice=clim_indices{ci};
                eval_line= ['tmp.corr=corrcoef(tmp.data_4ym, squeeze(data_', clim_indice, '.assm_4ym(mi,:)), ',...
                '''','Rows','''',', ','''','complete','''',');'];
                eval(eval_line);
                eval_line= ['corr_all.assm_4ym_', clim_indice, '(mi)=tmp.corr(1,2);'];
                eval(eval_line);
            end
    
        end
        
        %% corr_assm_em
        for ci=1:length(clim_indices)
            clim_indice=clim_indices{ci};
            eval_line= ['tmp.corr=corrcoef(squeeze(mean(m_data_assm.([cfg_ts.var,', ...
            '''','_ym','''',']),1)), squeeze(mean(data_', clim_indice, '.assm_ym,1)), ',...
            '''','Rows','''',', ','''','complete','''',');'];
            eval(eval_line);
            eval_line= ['corr_all_em.assm_ym_', clim_indice, '=tmp.corr(1,2);'];
            eval(eval_line);
        end
    
        
        %% corr_lens2
        for mi=1:len_mem_lens2
    
            for ci=1:length(clim_indices)
                clim_indice=clim_indices{ci};
                eval_line= ['tmp.corr=corrcoef(squeeze(m_data_lens2.([cfg_ts.var,', ...
                '''','_ym','''','])(mi,:)), squeeze(data_', clim_indice, '.lens2_ym(mi,:)), ',...
                '''','Rows','''',', ','''','complete','''',');'];
                eval(eval_line);
                eval_line= ['corr_all.lens2_ym_', clim_indice, '(mi)=tmp.corr(1,2);'];
                eval(eval_line);
            end
    
            %% 4y movmean
            tmp.data=squeeze(m_data_lens2.([cfg_ts.var,'_ym'])(mi,:));
            tmp.data_4ym=movmean(tmp.data,4,'Endpoints', 'discard');
    
            for ci=1:length(clim_indices)
                clim_indice=clim_indices{ci};
                eval_line= ['tmp.corr=corrcoef(tmp.data_4ym, squeeze(data_', clim_indice, '.lens2_4ym(mi,:)), ',...
                '''','Rows','''',', ','''','complete','''',');'];
                eval(eval_line);
                eval_line= ['corr_all.lens2_4ym_', clim_indice, '(mi)=tmp.corr(1,2);'];
                eval(eval_line);
            end
        end
        
        %% corr_lens2_em
        for ci=1:length(clim_indices)
            clim_indice=clim_indices{ci};
            eval_line= ['tmp.corr=corrcoef(squeeze(mean(m_data_lens2.([cfg_ts.var,', ...
            '''','_ym','''',']),1)), squeeze(mean(data_', clim_indice, '.lens2_ym,1)), ',...
            '''','Rows','''',', ','''','complete','''',');'];
            eval(eval_line);
            eval_line= ['corr_all_em.lens2_ym_', clim_indice, '=tmp.corr(1,2);'];
            eval(eval_line);
        end
        
        %% corr_hcst
        for ly=1:5
            ly_str=['ly',num2str(ly)];
            for mi=1:len_mem_oda
    
                for ci=1:length(clim_indices)
                    clim_indice=clim_indices{ci};
                    eval_line= ['tmp.corr=corrcoef(squeeze(m_data_hcst.([cfg_ts.var,', ...
                    '''','_ym','''',']).(ly_str)(mi,:)), squeeze(data_', clim_indice, '.hcst_ym(ly,mi,:)), ',...
                    '''','Rows','''',', ','''','complete','''',');'];
                    eval(eval_line);
                    eval_line= ['corr_all.hcst_ym_', clim_indice, '(ly,mi)=tmp.corr(1,2);'];
                    eval(eval_line);
                end
            end
            
            %% corr_hcst_em
            for ci=1:length(clim_indices)
                clim_indice=clim_indices{ci};
                eval_line= ['tmp.corr=corrcoef(squeeze(mean(m_data_hcst.([cfg_ts.var,', ...
                '''','_ym','''',']).(ly_str),1)), squeeze(mean(data_', clim_indice, '.hcst_ym(ly,:,:),2)), ',...
                '''','Rows','''',', ','''','complete','''',');'];
                eval(eval_line);
                eval_line= ['corr_all_em.hcst_ym_', clim_indice, '(ly)=tmp.corr(1,2);'];
                eval(eval_line);
            end
        end
        
        %% 4y movmean
    %     tmp=rmfield(tmp, 'data_4ym');
        for mi=1:len_mem_oda
            tmp.data_4ym_hc(mi,:)=squeeze(m_data_hcst.([cfg_ts.var,'_ym']).ly2(mi,:)) + ...
                squeeze(m_data_hcst.([cfg_ts.var,'_ym']).ly3(mi,:)) + ...
                squeeze(m_data_hcst.([cfg_ts.var,'_ym']).ly4(mi,:)) + ...
                squeeze(m_data_hcst.([cfg_ts.var,'_ym']).ly5(mi,:));
            tmp.data_4ym_hc(mi,:)=tmp.data_4ym_hc(mi,:)/4;
    
    %         tmp.data_4ym_hc(mi,:)=movmean(tmp.data(mi,:),4,2,'Endpoints', 'discard');
            
            for ci=1:length(clim_indices)
                clim_indice=clim_indices{ci};
                eval_line= ['tmp.corr=corrcoef(tmp.data_4ym_hc(mi,:), squeeze(data_', clim_indice, '.hcst_4ym(mi,:)), ',...
                '''','Rows','''',', ','''','complete','''',');'];
                eval(eval_line);
                eval_line= ['corr_all.hcst_4ym_', clim_indice, '(mi)=tmp.corr(1,2);'];
                eval(eval_line);
            end
        end
        
    %% raw, corrplot
        tmp.rows=1;
        for ci=1:length(clim_indices)
            clim_indice=clim_indices{ci};
            %% obs
            eval(['tmp.fig_mat(tmp.rows,1)=corr_all.obs_ym_', clim_indice, '; tmp.fig_mat(tmp.rows,2:52)=NaN;']);
            tmp.rows=tmp.rows+1;
            %% assm
            eval(['tmp.fig_mat(tmp.rows,1)=mean(corr_all.assm_ym_', clim_indice, '); tmp.fig_mat(tmp.rows,2)=NaN;']);
            eval(['tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_ym_', clim_indice, ';'])
            tmp.fig_mat(tmp.rows,23:52)=NaN;
            tmp.rows=tmp.rows+1;
            %% hcst
            for ly=1:5
                eval(['tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_', clim_indice, '(', num2str(ly), ',:)); tmp.fig_mat(tmp.rows,2)=NaN;']);
                eval(['tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_', clim_indice, '(', num2str(ly), ',:);']);
                tmp.fig_mat(tmp.rows,23:52)=NaN;
                tmp.rows=tmp.rows+1;
            end
            %% lens2
            eval(['tmp.fig_mat(tmp.rows,1)=mean(corr_all.lens2_ym_', clim_indice, '); tmp.fig_mat(tmp.rows,2)=NaN;']);
            eval(['tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_ym_', clim_indice, ';']);
            tmp.rows=tmp.rows+1;
            tmp.fig_mat(tmp.rows,1:52)=NaN;
            tmp.rows=tmp.rows+1;
        end
        tmp.fig_mat(tmp.rows,1:52)=NaN;  
        lc=length(clim_indices);
        fig_h = figure('name','ts','visible','off');
    %     fig_h.Position= [0,0,1000,800];
        fig_h.Position= [0,0,1000,lc*100];
    
%         pcolor(tmp.fig_mat);
        hsc_d=imagesc(tmp.fig_mat);
        set(hsc_d, 'AlphaData', ~isnan(tmp.fig_mat))
        set(gca,'TickDir','out');
%         shading flat;
        colorbar; caxis([-1 1]); colormap(fig_cfg.c_map);
        
        clear vec_yticks
        vec_yticks(1:8)=[1,2,3,4,5,6,7,8];
        for ci=1:length(clim_indices)-1
            vec_yticks(end+1:end+8)=vec_yticks(end-7:end)+9;
        end
        yticks(vec_yticks);
        
        clear vec_yticklabels
        for ci=1:length(clim_indices)
            clim_indice=clim_indices{ci};
            clim_indice_rep=strrep(clim_indice,'_','-');
            vec_yticklabels{(ci-1)*8+1}=['OBS-',clim_indice_rep];
            vec_yticklabels{(ci-1)*8+2}=['ODA-',clim_indice_rep];
            vec_yticklabels{(ci-1)*8+3}=['(LY1) HIND-',clim_indice_rep];
            vec_yticklabels{(ci-1)*8+4}=['(LY2) HIND-',clim_indice_rep];
            vec_yticklabels{(ci-1)*8+5}=['(LY3) HIND-',clim_indice_rep];
            vec_yticklabels{(ci-1)*8+6}=['(LY4) HIND-',clim_indice_rep];
            vec_yticklabels{(ci-1)*8+7}=['(LY5) HIND-',clim_indice_rep];
            vec_yticklabels{(ci-1)*8+8}=['LE-',clim_indice_rep];
        end
        yticklabels(vec_yticklabels);
    
        xticks([1, 7:5:52]);
        xticklabels({'mean', '5', '10', '15', '20', '25', '30', '35', '40', '45', '50'})
    
        xlabel('members');
        if length(months)==12
             title(['r with indices', ', ', cfg_ts.var, ', ', ...
            num2str(sta_lonlat{stai}(1)),'E, ', num2str(sta_lonlat{stai}(2)),'E, ', ...
            num2str(sta_lonlat{stai}(3)), 'N, ', num2str(sta_lonlat{stai}(4)), 'N']);
        elseif sum(months>12)==0
             title(['r with indices', ', sm, ', cfg_ts.var, ', ', ...
            num2str(sta_lonlat{stai}(1)),'E, ', num2str(sta_lonlat{stai}(2)),'E, ', ...
            num2str(sta_lonlat{stai}(3)), 'N, ', num2str(sta_lonlat{stai}(4)), 'N'])
        else
            title(['r with indices', ', sm, ', cfg_ts.var, ', ', ...
            num2str(sta_lonlat{stai}(1)),'E, ', num2str(sta_lonlat{stai}(2)),'E, ', ...
            num2str(sta_lonlat{stai}(3)), 'N, ', num2str(sta_lonlat{stai}(4)), 'N'])
        end
        if length(months)==12
            dirs.figdir= [dirs.figroot, filesep, 'ens_all', filesep, cfg_ts.var, '_corr_indices'];
        else
            dirs.figdir= [dirs.figroot, filesep, 'ens_all', filesep, cfg_ts.var, '_corr_indices','_', str_prepos];
        end
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        if length(months)==12
            cfg.figname=[dirs.figdir, filesep, 'r_indices_', num2str(sta_lonlat{stai}(1)), 'E_', num2str(sta_lonlat{stai}(2)), 'E_', ...
                num2str(sta_lonlat{stai}(3)), 'N_', num2str(sta_lonlat{stai}(4)), 'N_', cfg_ts.var, '.tif'];
        else
            cfg.figname=[dirs.figdir, filesep, 'r_indices_', str_prepos, '_', num2str(min(months)), '_', num2str(max(months)), '_', num2str(sta_lonlat{stai}(1)), 'E_', num2str(sta_lonlat{stai}(2)), 'E_', ...
                num2str(sta_lonlat{stai}(3)), 'N_', num2str(sta_lonlat{stai}(4)), 'N_', cfg_ts.var, '.tif'];
        end
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;


    %% 2-5y, corrplot
        tmp=rmfield(tmp, 'fig_mat');
        tmp.rows=1;
        for ci=1:length(clim_indices)
            clim_indice=clim_indices{ci};
            %% obs
            eval(['tmp.fig_mat(tmp.rows,1)=corr_all.obs_4ym_', clim_indice, '; tmp.fig_mat(tmp.rows,2:52)=NaN;']);
            tmp.rows=tmp.rows+1;
            %% assm
            eval(['tmp.fig_mat(tmp.rows,1)=mean(corr_all.assm_4ym_', clim_indice, '); tmp.fig_mat(tmp.rows,2)=NaN;']);
            eval(['tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_4ym_', clim_indice, ';'])
            tmp.fig_mat(tmp.rows,23:52)=NaN;
            tmp.rows=tmp.rows+1;
            %% hcst
            eval(['tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_4ym_', clim_indice, '); tmp.fig_mat(tmp.rows,2)=NaN;']);
            eval(['tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_4ym_', clim_indice, ';']);
            tmp.fig_mat(tmp.rows,23:52)=NaN;
            tmp.rows=tmp.rows+1;
            %% lens2
            eval(['tmp.fig_mat(tmp.rows,1)=mean(corr_all.lens2_4ym_', clim_indice, '); tmp.fig_mat(tmp.rows,2)=NaN;']);
            eval(['tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_4ym_', clim_indice, ';']);
            tmp.rows=tmp.rows+1;
            tmp.fig_mat(tmp.rows,1:52)=NaN;
            tmp.rows=tmp.rows+1;
        end
        tmp.fig_mat(tmp.rows,1:52)=NaN;  
        lc=length(clim_indices);
        fig_h = figure('name','ts','visible','off');
    %     fig_h.Position= [0,0,1000,800];
        fig_h.Position= [0,0,1000,lc*100];
    
%         pcolor(tmp.fig_mat);
        hsc_d=imagesc(tmp.fig_mat);
        set(hsc_d, 'AlphaData', ~isnan(tmp.fig_mat))
        set(gca,'TickDir','out');
%         shading flat;
        colorbar; caxis([-1 1]); colormap(fig_cfg.c_map);
        
        clear vec_yticks
        vec_yticks(1:4)=[1,2,3,4];
        for ci=1:length(clim_indices)-1
            vec_yticks(end+1:end+4)=vec_yticks(end-3:end)+5;
        end
        yticks(vec_yticks);
        
        clear vec_yticklabels
        for ci=1:length(clim_indices)
            clim_indice=clim_indices{ci};
            clim_indice_rep=strrep(clim_indice,'_','-');
            vec_yticklabels{(ci-1)*4+1}=['OBS-',clim_indice_rep];
            vec_yticklabels{(ci-1)*4+2}=['ODA-',clim_indice_rep];
            vec_yticklabels{(ci-1)*4+3}=['(LY2-5) HIND-',clim_indice_rep];
            vec_yticklabels{(ci-1)*4+4}=['LE-',clim_indice_rep];
        end
        yticklabels(vec_yticklabels);
    
        xticks([1, 7:5:52]);
        xticklabels({'mean', '5', '10', '15', '20', '25', '30', '35', '40', '45', '50'})
    
        xlabel('members');
        if length(months)==12
             title(['r with indices', ', ', cfg_ts.var, ', ', ...
            num2str(sta_lonlat{stai}(1)),'E, ', num2str(sta_lonlat{stai}(2)),'E, ', ...
            num2str(sta_lonlat{stai}(3)), 'N, ', num2str(sta_lonlat{stai}(4)), 'N']);
        elseif sum(months>12)==0
             title(['r with indices', ', sm, ', cfg_ts.var, ', ', ...
            num2str(sta_lonlat{stai}(1)),'E, ', num2str(sta_lonlat{stai}(2)),'E, ', ...
            num2str(sta_lonlat{stai}(3)), 'N, ', num2str(sta_lonlat{stai}(4)), 'N'])
        else
            title(['r with indices', ', sm, ', cfg_ts.var, ', ', ...
            num2str(sta_lonlat{stai}(1)),'E, ', num2str(sta_lonlat{stai}(2)),'E, ', ...
            num2str(sta_lonlat{stai}(3)), 'N, ', num2str(sta_lonlat{stai}(4)), 'N'])
        end
        if length(months)==12
            dirs.figdir= [dirs.figroot, filesep, 'ens_all', filesep, cfg_ts.var, '_corr_indices'];
        else
            dirs.figdir= [dirs.figroot, filesep, 'ens_all', filesep, cfg_ts.var, '_corr_indices','_', str_prepos];
        end
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        if length(months)==12
            cfg.figname=[dirs.figdir, filesep, 'r_4ym_indices_', num2str(sta_lonlat{stai}(1)), 'E_', num2str(sta_lonlat{stai}(2)), 'E_', ...
                num2str(sta_lonlat{stai}(3)), 'N_', num2str(sta_lonlat{stai}(4)), 'N_', cfg_ts.var, '.tif'];
        else
            cfg.figname=[dirs.figdir, filesep, 'r_4ym_indices_', str_prepos, '_', num2str(min(months)), '_', num2str(max(months)), '_', num2str(sta_lonlat{stai}(1)), 'E_', num2str(sta_lonlat{stai}(2)), 'E_', ...
                num2str(sta_lonlat{stai}(3)), 'N_', num2str(sta_lonlat{stai}(4)), 'N_', cfg_ts.var, '.tif'];
        end
        print(fig_h, cfg.figname, '-dpng');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    
    for iiii=1:1  %% folded
    % % % %% 4y movmean
    % % %     tmp=rmfield(tmp, 'fig_mat');
    % % % %% ASSM_AMO (1~2 row)
    % % %     tmp.rows=1;
    % % %     tmp.fig_mat(tmp.rows,1)=mean(corr_all.assm_4ym_AMO); tmp.fig_mat(tmp.rows,2)=NaN;
    % % %     tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_4ym_AMO;
    % % %     tmp.fig_mat(tmp.rows,23:52)=NaN;
    % % %     
    % % %     %% HCST_AMO (3~12 row)
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_4ym_AMO(:)); tmp.fig_mat(tmp.rows,2)=NaN;
    % % %     tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_4ym_AMO(:);
    % % %     tmp.fig_mat(tmp.rows,23:52)=NaN;
    % % %         
    % % %     %% LENS2_AMO (13~14 row)
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1)=mean(corr_all.lens2_4ym_AMO); tmp.fig_mat(tmp.rows,2)=NaN;    
    % % %     tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_4ym_AMO;
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1:52)=NaN;
    % % %     
    % % %     %% ASSM_ENSO 
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1)=mean(corr_all.assm_4ym_ENSO); tmp.fig_mat(tmp.rows,2)=NaN;    
    % % %     tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_4ym_ENSO;
    % % %     tmp.fig_mat(tmp.rows,23:52)=NaN;
    % % %     
    % % %     %% HCST_ENSO 
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_4ym_ENSO(:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    % % %     tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_4ym_ENSO(:);
    % % %     tmp.fig_mat(tmp.rows,23:52)=NaN;
    % % %     
    % % %     %% LENS2_ENSO 
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1)=mean(corr_all.lens2_4ym_ENSO); tmp.fig_mat(tmp.rows,2)=NaN;    
    % % %     tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_4ym_ENSO;
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1:52)=NaN;
    % % %     
    % % %     %% ASSM_PDO 
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1)=mean(corr_all.assm_4ym_PDO); tmp.fig_mat(tmp.rows,2)=NaN;    
    % % %     tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_4ym_PDO;
    % % %     tmp.fig_mat(tmp.rows,23:52)=NaN;
    % % %     
    % % %     %% HCST_PDO 
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_4ym_PDO(:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    % % %     tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_4ym_PDO(:);
    % % %     tmp.fig_mat(tmp.rows,23:52)=NaN;
    % % %     
    % % %     %% LENS2_PDO 
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1)=mean(corr_all.lens2_4ym_PDO); tmp.fig_mat(tmp.rows,2)=NaN;    
    % % %     tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_4ym_PDO;
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1:52)=NaN;
    % % % 
    % % % 
    % % %     %% ASSM_PDO_l 
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1)=mean(corr_all.assm_4ym_PDO_l); tmp.fig_mat(tmp.rows,2)=NaN;    
    % % %     tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_4ym_PDO_l;
    % % %     tmp.fig_mat(tmp.rows,23:52)=NaN;
    % % %     
    % % %     %% HCST_PDO_l 
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_4ym_PDO_l(:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    % % %     tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_4ym_PDO_l(:);
    % % %     tmp.fig_mat(tmp.rows,23:52)=NaN;
    % % %     
    % % %     %% LENS2_PDO_l 
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1)=mean(corr_all.lens2_4ym_PDO_l); tmp.fig_mat(tmp.rows,2)=NaN;    
    % % %     tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_4ym_PDO_l;
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1:52)=NaN;
    % % %     
    % % %     
    % % %     %% ASSM_NPGO_SST 
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1)=mean(corr_all.assm_4ym_NPGO_SST); tmp.fig_mat(tmp.rows,2)=NaN;    
    % % %     tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_4ym_NPGO_SST;
    % % %     tmp.fig_mat(tmp.rows,23:52)=NaN;
    % % %     
    % % %     %% HCST_NPGO_SST 
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_4ym_NPGO_SST(:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    % % %     tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_4ym_NPGO_SST(:);
    % % %     tmp.fig_mat(tmp.rows,23:52)=NaN;
    % % %    
    % % %     
    % % %     %% LENS2_NPGO_SST 
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1)=mean(corr_all.lens2_4ym_NPGO_SST); tmp.fig_mat(tmp.rows,2)=NaN;    
    % % %     tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_4ym_NPGO_SST;
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1:52)=NaN;
    % % %     
    % % %     
    % % %     %% ASSM_NPGO_SSH 
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1)=mean(corr_all.assm_4ym_NPGO_SSH); tmp.fig_mat(tmp.rows,2)=NaN;    
    % % %     tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_4ym_NPGO_SSH;
    % % %     tmp.fig_mat(tmp.rows,23:52)=NaN;
    % % % 
    % % % 
    % % %     %% HCST_NPGO_SSH 
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_4ym_NPGO_SSH(:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    % % %     tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_4ym_NPGO_SSH(:);
    % % %     tmp.fig_mat(tmp.rows,23:52)=NaN;
    % % %    
    % % %     
    % % %     %% LENS2_NPGO_SSH 
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1)=mean(corr_all.lens2_4ym_NPGO_SSH); tmp.fig_mat(tmp.rows,2)=NaN;    
    % % %     tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_4ym_NPGO_SSH;
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1:52)=NaN;
    % % %     
    % % %     
    % % %     %% ASSM_SAM 
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1)=mean(corr_all.assm_4ym_SAM); tmp.fig_mat(tmp.rows,2)=NaN;    
    % % %     tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_4ym_SAM;
    % % %     tmp.fig_mat(tmp.rows,23:52)=NaN;
    % % %     
    % % %     %% HCST_SAM 
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_4ym_SAM(:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    % % %     tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_4ym_SAM(:);
    % % %     tmp.fig_mat(tmp.rows,23:52)=NaN;
    % % %     
    % % %     
    % % %     %% LENS2_SAM 
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1)=mean(corr_all.lens2_4ym_SAM); tmp.fig_mat(tmp.rows,2)=NaN;    
    % % %     tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_4ym_SAM;
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1:52)=NaN;
    % % % 
    % % % 
    % % %     %% ASSM_ATL3 
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1)=mean(corr_all.assm_4ym_ATL3); tmp.fig_mat(tmp.rows,2)=NaN;    
    % % %     tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_4ym_ATL3;
    % % %     tmp.fig_mat(tmp.rows,23:52)=NaN;
    % % %     
    % % %     %% HCST_ATL3 
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_4ym_ATL3(:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    % % %     tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_4ym_ATL3(:);
    % % %     tmp.fig_mat(tmp.rows,23:52)=NaN;
    % % %     
    % % %     
    % % %     %% LENS2_ATL3 
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1)=mean(corr_all.lens2_4ym_ATL3); tmp.fig_mat(tmp.rows,2)=NaN;    
    % % %     tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_4ym_ATL3;
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1:52)=NaN;
    % % % 
    % % % 
    % % % 
    % % %     %% ASSM_TNA 
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1)=mean(corr_all.assm_4ym_TNA); tmp.fig_mat(tmp.rows,2)=NaN;    
    % % %     tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_4ym_TNA;
    % % %     tmp.fig_mat(tmp.rows,23:52)=NaN;
    % % %     
    % % %     %% HCST_TNA 
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_4ym_TNA(:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    % % %     tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_4ym_TNA(:);
    % % %     tmp.fig_mat(tmp.rows,23:52)=NaN;
    % % %     
    % % %     
    % % %     %% LENS2_TNA 
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1)=mean(corr_all.lens2_4ym_TNA); tmp.fig_mat(tmp.rows,2)=NaN;    
    % % %     tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_4ym_TNA;
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1:52)=NaN;
    % % % 
    % % % 
    % % %     %% ASSM_TSA 
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1)=mean(corr_all.assm_4ym_TSA); tmp.fig_mat(tmp.rows,2)=NaN;    
    % % %     tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_4ym_TSA;
    % % %     tmp.fig_mat(tmp.rows,23:52)=NaN;
    % % %     
    % % %     %% HCST_TSA 
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_4ym_TSA(:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    % % %     tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_4ym_TSA(:);
    % % %     tmp.fig_mat(tmp.rows,23:52)=NaN;
    % % %     
    % % %     
    % % %     %% LENS2_TSA 
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1)=mean(corr_all.lens2_4ym_TSA); tmp.fig_mat(tmp.rows,2)=NaN;    
    % % %     tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_4ym_TSA;
    % % %     tmp.rows=tmp.rows+1;
    % % %     tmp.fig_mat(tmp.rows,1:52)=NaN;
    % % %     
    % % % 
    % % %     fig_h = figure('name','ts','visible','off');
    % % %         fig_h.Position= [0,0,1000,500];
    % % % 
    % % %     pcolor(tmp.fig_mat); shading flat; colorbar; caxis([-1 1]); colormap(fig_cfg.c_map);
    % % %     yticks([1,2,3, 5,6,7, 9,10,11, 13,14,15, 17,18,19, 21,22,23, 25,26,27, 29,30,31, 33,34,35])
    % % %     yticklabels({'ODA-AMO', 'HIND-AMO','LE-AMO', ...
    % % %         'ODA-ENSO', 'HIND-ENSO','LE-ENSO', ...
    % % %         'ODA-PDO', 'HIND-PDO','LE-PDO', ...
    % % %         'ODA-PDO_l', 'HIND-PDO_l','LE-PDO_l', ...
    % % %         'ODA-NPGO-SST', 'HIND-NPGO-SST','LE-NPGO-SST', ...
    % % %         'ODA-NPGO-SSH', 'HIND-NPGO-SSH','LE-NPGO-SSH', ...
    % % %         'ODA-SAM', 'HIND-SAM','LE-SAM', ...
    % % %         'ODA-ATL3', 'HIND-ATL3','LE-ATL3', ...
    % % %         'ODA-TNA', 'HIND-TNA','LE-TNA', ...
    % % %         'ODA-TSA', 'HIND-TSA','LE-TSA'})
    % % %     xticks([1, 7:5:52]);
    % % %     xticklabels({'mean', '5', '10', '15', '20', '25', '30', '35', '40', '45', '50'})
    % % % 
    % % % 
    % % %     xlabel('members');
    % % %     
    % % %     title(['r with indices', ', ', cfg_ts.var, ', ', num2str(sta_lonlat{stai}(1)),'E, ', num2str(sta_lonlat{stai}(3)), 'N']);
    % % %     
    % % %     dirs.figdir= [dirs.figroot, filesep, 'ens_all', filesep, cfg_ts.var, '_corr_indices'];
    % % %     if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    % % %     cfg.figname=[dirs.figdir, filesep, 'r_4ym_indices_', num2str(sta_lonlat{stai}(1)), 'E_', num2str(sta_lonlat{stai}(3)), 'N_', cfg_ts.var, '.tif'];
    % % %     print(fig_h, cfg.figname, '-dpng');
    % % %     RemoveWhiteSpace([], 'file', cfg.figname);
    % % %     close all;
    % % % 
    end
    
        end
    end
    
    
    
    %% 1y (assm)
    fig_h = figure('name','ts','visible','off');
    for ci=1:length(clim_indices)
        clim_indice=clim_indices{ci};
        clim_indice_rep=strrep(clim_indice,'_','-');
        sa{ci}=subplot(length(clim_indices),1,ci);
        eval(['shade_anomaly(cfg.iyears, normalize( mean(data_', clim_indice, '.assm_ym,1) - mean(mean(data_', clim_indice, '.assm_ym,1),', '''','omitnan','''',') )', ...
            '''',', ','''', 'r', '''', ', ', '''', 'b', '''', ', ', '0.3, sa{', num2str(ci), '});']);
        sa{ci}.YLabel.String=clim_indice_rep;
    end
    
    sa9.XLabel.String='Years';
    set(gcf, 'Position', [0 0 1000 ci*100+200]);
    if length(months)==12
        str_prepos='annual';
        print(gcf, ...
        '/Volumes/kyy_raid/kimyy/Research/Postdoc/03_IBS/2022_predictability_assimilation_run/paper/Figureset_raw/Fin_Figure/sub/A_climate_indices_normalized.png', ...
         '-dpng');
    else
        str_prepos='seasonal';
        print(gcf, ...
        ['/Volumes/kyy_raid/kimyy/Research/Postdoc/03_IBS/2022_predictability_assimilation_run/', ...
        'paper/Figureset_raw/Fin_Figure/sub/A_climate_indices_normalized_', str_prepos, num2str(min(months)), '_', num2str(max(months)), '.png'], ...
         '-dpng');
    end
    close all;
    
    
    
    
    
    %% 1y (obs)
    fig_h = figure('name','ts','visible','off');
    for ci=1:length(clim_indices)
        clim_indice=clim_indices{ci};
        clim_indice_rep=strrep(clim_indice,'_','-');
        sa{ci}=subplot(length(clim_indices),1,ci);
        eval(['shade_anomaly(cfg.iyears, normalize( mean(data_', clim_indice, '.obs_ym,1) - mean(mean(data_', clim_indice, '.obs_ym,1),', '''','omitnan','''',') )', ...
            '''',', ','''', 'r', '''', ', ', '''', 'b', '''', ', ', '0.3, sa{', num2str(ci), '});']);
        sa{ci}.YLabel.String=clim_indice_rep;
    end
    
    sa9.XLabel.String='Years';
    set(gcf, 'Position', [0 0 1000 ci*100+200]);
    if length(months)==12
        str_prepos='annual';
        print(gcf, ...
        '/Volumes/kyy_raid/kimyy/Research/Postdoc/03_IBS/2022_predictability_assimilation_run/paper/Figureset_raw/Fin_Figure/sub/A_climate_indices_normalized_obs.png', ...
         '-dpng');
    else
        str_prepos='seasonal';
        print(gcf, ...
        ['/Volumes/kyy_raid/kimyy/Research/Postdoc/03_IBS/2022_predictability_assimilation_run/', ...
        'paper/Figureset_raw/Fin_Figure/sub/A_climate_indices_normalized_obs_', str_prepos, num2str(min(months)), '_', num2str(max(months)), '.png'], ...
         '-dpng');
    end
    close all;
    
    
    
    
    % % %% 1y (obs)
    % % sa1=subplot(9,1,1);
    % % shade_anomaly(cfg.iyears, normalize( data_AMO.obs_ym - mean(data_AMO.obs_ym) )', ...
    % %     'r', 'b', 0.3, sa1);
    % % sa1.YLabel.String='AMO';
    % % 
    % % sa2=subplot(9,1,2);
    % % shade_anomaly(cfg.iyears, normalize( data_PDO.obs_ym - mean(data_PDO.obs_ym) )', ...
    % %     'r', 'b', 0.3, sa2);
    % % sa2.YLabel.String='PDO';
    % % 
    % % sa3=subplot(9,1,3);
    % % shade_anomaly(cfg.iyears, normalize( data_ENSO.obs_ym - mean(data_ENSO.obs_ym) )', ...
    % %     'r', 'b', 0.3, sa3);
    % % sa3.YLabel.String='ENSO';
    % % 
    % % sa4=subplot(9,1,4);
    % % shade_anomaly(cfg.iyears, normalize( data_NPGO_SST.obs_ym - mean(data_NPGO_SST.obs_ym) )', ...
    % %     'r', 'b', 0.3, sa4);
    % % sa4.YLabel.String='NPGO-SST';
    % % 
    % % sa5=subplot(9,1,5);
    % % shade_anomaly(cfg.iyears, normalize( data_NPGO_SSH.obs_ym - mean(data_NPGO_SSH.obs_ym) )', ...
    % %     'r', 'b', 0.3, sa5);
    % % sa5.YLabel.String='NPGO-SSH';
    % % 
    % % sa6=subplot(9,1,6);
    % % shade_anomaly(cfg.iyears, normalize( data_SAM.obs_ym - mean(data_SAM.obs_ym) )', ...
    % %     'r', 'b', 0.3, sa6);
    % % sa6.YLabel.String='SAM';
    % % 
    % % sa7=subplot(9,1,7);
    % % shade_anomaly(cfg.iyears, normalize( data_ATL3.obs_ym - mean(data_ATL3.obs_ym) )', ...
    % %     'r', 'b', 0.3, sa7);
    % % sa7.YLabel.String='ATL3';
    % % sa7.XLabel.String='Years';
    % % 
    % % sa8=subplot(9,1,8);
    % % shade_anomaly(cfg.iyears, normalize( data_TNA.obs_ym - mean(data_TNA.obs_ym) )', ...
    % %     'r', 'b', 0.3, sa8);
    % % sa8.YLabel.String='TNA';
    % % sa8.XLabel.String='Years';
    % % 
    % % sa9=subplot(9,1,9);
    % % shade_anomaly(cfg.iyears, normalize( data_TSA.obs_ym - mean(data_TSA.obs_ym) )', ...
    % %     'r', 'b', 0.3, sa9);
    % % sa9.YLabel.String='TSA';
    % % sa9.XLabel.String='Years';
    % % set(gcf, 'Position', [0 0 1000 800]);
    % % print(gcf, ...
    % %     '/kyy_raid/kimyy/Research/Postdoc/03_IBS/2022_predictability_assimilation_run/paper/Figureset_raw/Fin_Figure/sub/A_climate_indices_normalized_obs.png', ...
    % %      '-dpng');
    % % 
    % % close all;
    
    
% % % %     cfg.iyears_4ym=movmean(cfg.iyears, 4, 'Endpoints', 'discard');    
% % % %     %% 4ym (assm)
% % % %     sa1=subplot(9,1,1);
% % % %     shade_anomaly(cfg.iyears_4ym, normalize( mean(data_AMO.assm_4ym,1) - mean(mean(data_AMO.assm_4ym,1)) )', ...
% % % %         'r', 'b', 0.3, sa1);
% % % %     sa1.YLabel.String='AMO';
% % % %     
% % % %     sa2=subplot(9,1,2);
% % % %     shade_anomaly(cfg.iyears_4ym, normalize( mean(data_PDO.assm_4ym,1) - mean(mean(data_PDO.assm_4ym,1)) )', ...
% % % %         'r', 'b', 0.3, sa2);
% % % %     sa2.YLabel.String='PDO';
% % % %     
% % % %     sa3=subplot(9,1,3);
% % % %     shade_anomaly(cfg.iyears_4ym, normalize( mean(data_ENSO.assm_4ym,1) - mean(mean(data_ENSO.assm_4ym,1)) )', ...
% % % %         'r', 'b', 0.3, sa3);
% % % %     sa3.YLabel.String='ENSO';
% % % %     
% % % %     sa4=subplot(9,1,4);
% % % %     shade_anomaly(cfg.iyears_4ym, normalize( mean(data_NPGO_SST.assm_4ym,1) - mean(mean(data_NPGO_SST.assm_4ym,1)) )', ...
% % % %         'r', 'b', 0.3, sa4);
% % % %     sa4.YLabel.String='NPGO-SST';
% % % %     
% % % %     sa5=subplot(9,1,5);
% % % %     shade_anomaly(cfg.iyears_4ym, normalize( mean(data_NPGO_SSH.assm_4ym,1) - mean(mean(data_NPGO_SSH.assm_4ym,1)) )', ...
% % % %         'r', 'b', 0.3, sa5);
% % % %     sa5.YLabel.String='NPGO-SSH';
% % % %     
% % % %     sa6=subplot(9,1,6);
% % % %     shade_anomaly(cfg.iyears_4ym, normalize( mean(data_SAM.assm_4ym,1) - mean(mean(data_SAM.assm_4ym,1)) )', ...
% % % %         'r', 'b', 0.3, sa6);
% % % %     sa6.YLabel.String='SAM';
% % % %     
% % % %     sa7=subplot(9,1,7);
% % % %     shade_anomaly(cfg.iyears_4ym, normalize( mean(data_ATL3.assm_4ym,1) - mean(mean(data_ATL3.assm_4ym,1)) )', ...
% % % %         'r', 'b', 0.3, sa7);
% % % %     sa7.YLabel.String='ATL3';
% % % %     sa7.XLabel.String='Years';
% % % %     
% % % %     sa8=subplot(9,1,8);
% % % %     shade_anomaly(cfg.iyears_4ym, normalize( mean(data_TNA.assm_4ym,1) - mean(mean(data_TNA.assm_4ym,1)) )', ...
% % % %         'r', 'b', 0.3, sa8);
% % % %     sa8.YLabel.String='TNA';
% % % %     sa8.XLabel.String='Years';
% % % %     
% % % %     sa9=subplot(9,1,9);
% % % %     shade_anomaly(cfg.iyears_4ym, normalize( mean(data_TSA.assm_4ym,1) - mean(mean(data_TSA.assm_4ym,1)) )', ...
% % % %         'r', 'b', 0.3, sa9);
% % % %     sa9.YLabel.String='TSA';
% % % %     sa9.XLabel.String='Years';
% % % %     set(gcf, 'Position', [0 0 1000 800]);
% % % %     print(gcf, ...
% % % %         '/kyy_raid/kimyy/Research/Postdoc/03_IBS/2022_predictability_assimilation_run/paper/Figureset_raw/Fin_Figure/sub/A_climate_indices_normalized_4ym.png', ...
% % % %          '-dpng');
% % % %     
% % % %     close all;
% % % %     
% % % %     
% % % %     %% 4ym (obs)
% % % %     sa1=subplot(9,1,1);
% % % %     shade_anomaly(cfg.iyears_4ym , normalize( data_AMO.obs_4ym - mean(data_AMO.obs_4ym) )', ...
% % % %         'r', 'b', 0.3, sa1);
% % % %     sa1.YLabel.String='AMO';
% % % %     
% % % %     sa2=subplot(9,1,2);
% % % %     shade_anomaly(cfg.iyears_4ym, normalize( data_PDO.obs_4ym - mean(data_PDO.obs_4ym) )', ...
% % % %         'r', 'b', 0.3, sa2);
% % % %     sa2.YLabel.String='PDO';
% % % %     
% % % %     sa3=subplot(9,1,3);
% % % %     shade_anomaly(cfg.iyears_4ym, normalize( data_ENSO.obs_4ym - mean(data_ENSO.obs_4ym) )', ...
% % % %         'r', 'b', 0.3, sa3);
% % % %     sa3.YLabel.String='ENSO';
% % % %     
% % % %     sa4=subplot(9,1,4);
% % % %     shade_anomaly(cfg.iyears_4ym, normalize( data_NPGO_SST.obs_4ym - mean(data_NPGO_SST.obs_4ym) )', ...
% % % %         'r', 'b', 0.3, sa4);
% % % %     sa4.YLabel.String='NPGO-SST';
% % % %     
% % % %     sa5=subplot(9,1,5);
% % % %     shade_anomaly(cfg.iyears_4ym, normalize( data_NPGO_SSH.obs_4ym - mean(data_NPGO_SSH.obs_4ym) )', ...
% % % %         'r', 'b', 0.3, sa5);
% % % %     sa5.YLabel.String='NPGO-SSH';
% % % %     
% % % %     sa6=subplot(9,1,6);
% % % %     shade_anomaly(cfg.iyears_4ym, normalize( data_SAM.obs_4ym - mean(data_SAM.obs_4ym) )', ...
% % % %         'r', 'b', 0.3, sa6);
% % % %     sa6.YLabel.String='SAM';
% % % %     
% % % %     sa7=subplot(9,1,7);
% % % %     shade_anomaly(cfg.iyears_4ym, normalize( data_ATL3.obs_4ym - mean(data_ATL3.obs_4ym) )', ...
% % % %         'r', 'b', 0.3, sa7);
% % % %     sa7.YLabel.String='ATL3';
% % % %     sa7.XLabel.String='Years';
% % % %     
% % % %     sa8=subplot(9,1,8);
% % % %     shade_anomaly(cfg.iyears_4ym, normalize( data_TNA.obs_4ym - mean(data_TNA.obs_4ym) )', ...
% % % %         'r', 'b', 0.3, sa8);
% % % %     sa8.YLabel.String='TNA';
% % % %     sa8.XLabel.String='Years';
% % % %     
% % % %     sa9=subplot(9,1,9);
% % % %     shade_anomaly(cfg.iyears_4ym, normalize( data_TSA.obs_4ym - mean(data_TSA.obs_4ym) )', ...
% % % %         'r', 'b', 0.3, sa9);
% % % %     sa9.YLabel.String='TSA';
% % % %     sa9.XLabel.String='Years';
% % % %     set(gcf, 'Position', [0 0 1000 800]);
% % % %     print(gcf, ...
% % % %         '/kyy_raid/kimyy/Research/Postdoc/03_IBS/2022_predictability_assimilation_run/paper/Figureset_raw/Fin_Figure/sub/A_climate_indices_normalized_obs_4ym.png', ...
% % % %          '-dpng');
% % % %     
% % % %     close all;


end


function obsname_simple = f_obs_name(varn)
    switch varn
        case 'SST'
            obsname_simple='ERSST';
        case 'PRECT'
            obsname_simple='GPCC';
        case 'RAIN'
            obsname_simple='GPCC';
        case 'PSL'
            obsname_simple='ERA5';
        case 'U'
            obsname_simple='ERA5';
        case 'V'
            obsname_simple='ERA5';
        case 'SOILWATER_10CM'
%             obsname_simple='CMEMS';
            obsname_simple='GLEAM';
        case 'TWS'
            obsname_simple='NOAA';
        case 'SSH'
            obsname_simple='CMEMS';
        case 'TS'
%             obsname_simple='HadCRUT5';
            obsname_simple='ERA5';
        case 'TREFHT'
%             obsname_simple='HadCRUT5';
            obsname_simple='ERA5';
        case 'sumChl'
            obsname_simple='OC_CCI';
        case 'TLAI'
            obsname_simple='NOAA'; % NOAA LAI
        case 'FAREA_BURNED'
            obsname_simple='MODIS'; % MODIS Fire_cci v5.1
            obsname_simple='AVHRR'; % AVHRR-LTDR
        case 'COL_FIRE_CLOSS'
            obsname_simple='GFED'; % GFED
        case 'photoC_TOT_zint'
%             obsname_simple='VGPM'; % VGPM
            obsname_simple='CMEMS'; %Globcolour            
        case 'photoC_TOT_zint_100m'
%             obsname_simple='VGPM'; % VGPM
            obsname_simple='CMEMS'; %Globcolour
        case 'GPP'
%             obsname_simple='ORNL_DAAC';
            obsname_simple='VODCA2GPP';
        case 'TEMP'
            obsname_simple='EN4';
%             obsname_simple='projdv7.3';
        case 'SALT'
            obsname_simple='EN4';
%             obsname_simple='projdv7.3';
        otherwise
            obsname_simple='nan';
    end
end
