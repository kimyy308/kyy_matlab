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
cfg_ts.vars={'SST'};
cfg_ts.vars={'DIC', 'TEMP', 'SALT', 'DIC_ALT_CO2'};

% cfg_ts.var='NO3';

% cfg_ts.var='SSH';

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

%% 145 30 should be rechecked


% cfg.vlayer=1; % surface, vertical slice
cfg.vlayer=24; % surface, vertical slice

% cfg.vlayer=1:10; % 10layer. don't put more than 15
cfg.vlayer_1st=min(cfg.vlayer);
cfg.vlayer_cnt=max(cfg.vlayer)-cfg.vlayer_1st+1;

vstr=['v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i')];

cfg.iyears=1960:2020;

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
        title([cfg_ts.var, ', ', num2str(sta_lonlat{stai}(1)),'E, ', num2str(sta_lonlat{stai}(3)), 'N'])
    
        grid minor
        xlim([1960 2025])
        set(gca, 'fontsize', 20)    
    
        dirs.figdir= [dirs.figroot, filesep, 'ens_all', filesep, cfg_ts.var, '_time_series_indv', filesep, 'l',tmp.lyear_str];
        if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
        cfg.figname=[dirs.figdir, filesep, 'ts_all_l',tmp.lyear_str, '_', num2str(sta_lonlat{stai}(1)), 'E_', num2str(sta_lonlat{stai}(3)), 'N_', cfg_ts.var, '.tif'];
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
    title([cfg_ts.var, ', ', num2str(sta_lonlat{stai}(1)),'E, ', num2str(sta_lonlat{stai}(3)), 'N'])

    grid minor
    xlim([1960 2025])
    set(gca, 'fontsize', 20)    

    dirs.figdir= [dirs.figroot, filesep, 'ens_all', filesep, cfg_ts.var, '_time_series_indv', filesep, 'l2_5'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'ts_all_l2_5', '_', num2str(sta_lonlat{stai}(1)), 'E_', num2str(sta_lonlat{stai}(3)), 'N_', cfg_ts.var, '.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;


    
    %% no observation
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
        %% PDO
        load ('/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/clim_indices/clim_indices_SST_all_PDO_obs_ERSST.mat');
        %% PDO_l
        load ('/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/clim_indices/clim_indices_SST_all_PDO_l_obs_ERSST.mat');
        %% NPGO
        data_NPGO=load ('/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/clim_indices/clim_indices_SSH_all_PDO_obs_CMEMS.mat');
        %% SAM
        load ('/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/clim_indices/clim_indices_PSL_all_SAM_obs_ERA5.mat');
        %% ATL3
        load ('/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/clim_indices/clim_indices_SST_all_ATL3_obs_ERSST.mat');
        %% TNA
        load ('/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/clim_indices/clim_indices_SST_all_TNA_obs_ERSST.mat');
        %% TSA
        load ('/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/clim_indices/clim_indices_SST_all_TSA_obs_ERSST.mat');
        
        %% yearly mean indices
        %% AMO
        len_m=size(data_AMO.obs_dseason,1);
        len_mem_oda=size(data_AMO.assm_dseason,1);
        len_mem_lens2=size(data_AMO.lens2_dseason,1);
        
        tmp.reshp=reshape(data_AMO.obs_dseason, [12, len_m/12]);
        data_AMO.obs_ym=squeeze(mean(tmp.reshp,1));
        data_AMO.obs_4ym=movmean(data_AMO.obs_ym, 4, 'Endpoints', 'discard');
        
        tmp.reshp=reshape(data_AMO.assm_dseason, [len_mem_oda, 12, len_m/12]);
        data_AMO.assm_ym=squeeze(mean(tmp.reshp,2));
        data_AMO.assm_4ym=movmean(data_AMO.assm_ym, 4, 2, 'Endpoints', 'discard');
        
        tmp.reshp=reshape(data_AMO.lens2_dseason, [len_mem_lens2, 12, len_m/12]);
        data_AMO.lens2_ym=squeeze(mean(tmp.reshp,2));
        data_AMO.lens2_4ym=movmean(data_AMO.lens2_ym, 4, 2, 'Endpoints', 'discard');
        
        tmp.reshp=reshape(data_AMO.hcst_dseason, [5, len_mem_oda, 12, len_m/12]);
        data_AMO.hcst_ym=squeeze(mean(tmp.reshp,3));
        data_AMO.hcst_4ym=squeeze(mean(data_AMO.hcst_ym(2:5,:,:),1));

        
        %% ENSO
        tmp.reshp=reshape(data_ENSO.obs, [12, len_m/12]);
        data_ENSO.obs_ym=squeeze(mean(tmp.reshp,1));
        data_ENSO.obs_4ym=movmean(data_ENSO.obs_ym, 4, 'Endpoints', 'discard');

        tmp.reshp=reshape(data_ENSO.assm, [len_mem_oda, 12, len_m/12]);
        data_ENSO.assm_ym=squeeze(mean(tmp.reshp,2));
        data_ENSO.assm_4ym=movmean(data_ENSO.assm_ym, 4, 2, 'Endpoints', 'discard');
        
        tmp.reshp=reshape(data_ENSO.lens2, [len_mem_lens2, 12, len_m/12]);
        data_ENSO.lens2_ym=squeeze(mean(tmp.reshp,2));
        data_ENSO.lens2_4ym=movmean(data_ENSO.lens2_ym, 4, 2, 'Endpoints', 'discard');
        
        tmp.reshp=reshape(data_ENSO.hcst, [5, len_mem_oda, 12, len_m/12]);
        data_ENSO.hcst_ym=squeeze(mean(tmp.reshp,3));
        data_ENSO.hcst_4ym=squeeze(mean(data_ENSO.hcst_ym(2:5,:,:),1));   

        %% PDO
        lmode=1;
        tmp.reshp=reshape(data_PDO.pct_obs(:,lmode), [12, len_m/12]);
        data_PDO.obs_ym=squeeze(mean(tmp.reshp,1));
        
        tmp.reshp=reshape(data_PDO.pct_assm(:,:,lmode), [len_mem_oda, 12, len_m/12]);
        data_PDO.assm_ym=squeeze(mean(tmp.reshp,2));
        
        tmp.reshp=reshape(data_PDO.pct_lens2(:,:,lmode), [len_mem_lens2, 12, len_m/12]);
        data_PDO.lens2_ym=squeeze(mean(tmp.reshp,2));
        
        tmp.reshp=reshape(data_PDO.pct_hcst(:,:,:,lmode), [5, len_mem_oda, 12, len_m/12]);
        data_PDO.hcst_ym=squeeze(mean(tmp.reshp,3)); 
        
        
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
        tmp.reshp=reshape(data_PDO_l.pct_obs(:,lmode), [12, len_m/12]);
        data_PDO_l.obs_ym=squeeze(mean(tmp.reshp,1));
        
        tmp.reshp=reshape(data_PDO_l.pct_assm(:,:,lmode), [len_mem_oda, 12, len_m/12]);
        data_PDO_l.assm_ym=squeeze(mean(tmp.reshp,2));
        
        tmp.reshp=reshape(data_PDO_l.pct_lens2(:,:,lmode), [len_mem_lens2, 12, len_m/12]);
        data_PDO_l.lens2_ym=squeeze(mean(tmp.reshp,2));
        
        tmp.reshp=reshape(data_PDO_l.pct_hcst(:,:,:,lmode), [5, len_mem_oda, 12, len_m/12]);
        data_PDO_l.hcst_ym=squeeze(mean(tmp.reshp,3)); 
        
        
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
        tmp.reshp=reshape(data_PDO.pct_obs(:,lmode), [12, len_m/12]);
        data_NPGO_SST.obs_ym=squeeze(mean(tmp.reshp,1));
        
        tmp.reshp=reshape(data_PDO.pct_assm(:,:,lmode), [len_mem_oda, 12, len_m/12]);
        data_NPGO_SST.assm_ym=squeeze(mean(tmp.reshp,2));
        
        tmp.reshp=reshape(data_PDO.pct_lens2(:,:,lmode), [len_mem_lens2, 12, len_m/12]);
        data_NPGO_SST.lens2_ym=squeeze(mean(tmp.reshp,2));
        
        tmp.reshp=reshape(data_PDO.pct_hcst(:,:,:,lmode), [5, len_mem_oda, 12, len_m/12]);
        data_NPGO_SST.hcst_ym=squeeze(mean(tmp.reshp,3));

        
        %% x: 70~90 y=20:30 sign + -> negative change
        xrange=70:90;
        yrange=20:30;
        %% obs
        tmp.sign=data_PDO.lv_obs(xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
        if tmp.sign>0
            data_NPGO_SST.obs_ym=-data_NPGO_SST.obs_ym;
        end
        %% assm
        for mi=1:len_mem_oda
            tmp.sign=data_PDO.lv_assm(mi,xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
            if tmp.sign>0
                data_NPGO_SST.assm_ym(mi,:)=-data_NPGO_SST.assm_ym(mi,:);
            end
        end
        %% lens2
        for mi=1:len_mem_lens2
            tmp.sign=data_PDO.lv_lens2(mi,xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
            if tmp.sign>0
                data_NPGO_SST.lens2_ym(mi,:)=-data_NPGO_SST.lens2_ym(mi,:);
            end
        end
        %% hcst
        for ly=1:5
            for mi=1:len_mem_oda
                tmp.sign=data_PDO.lv_hcst(ly,mi,xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
                if tmp.sign>0
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
        tmp.reshp=reshape(data_PDO_l.pct_obs(:,lmode), [12, len_m/12]);
        data_NPGO_l_SST.obs_ym=squeeze(mean(tmp.reshp,1));
        
        tmp.reshp=reshape(data_PDO_l.pct_assm(:,:,lmode), [len_mem_oda, 12, len_m/12]);
        data_NPGO_l_SST.assm_ym=squeeze(mean(tmp.reshp,2));
        
        tmp.reshp=reshape(data_PDO_l.pct_lens2(:,:,lmode), [len_mem_lens2, 12, len_m/12]);
        data_NPGO_l_SST.lens2_ym=squeeze(mean(tmp.reshp,2));
        
        tmp.reshp=reshape(data_PDO_l.pct_hcst(:,:,:,lmode), [5, len_mem_oda, 12, len_m/12]);
        data_NPGO_l_SST.hcst_ym=squeeze(mean(tmp.reshp,3));

        
        %% x: 70~90 y=20:30 sign + -> negative change
        xrange=70:90;
        yrange=20:30;
        %% obs
        tmp.sign=data_PDO_l.lv_obs(xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
        if tmp.sign>0
            data_NPGO_l_SST.obs_ym=-data_NPGO_l_SST.obs_ym;
        end
        %% assm
        for mi=1:len_mem_oda
            tmp.sign=data_PDO_l.lv_assm(mi,xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
            if tmp.sign>0
                data_NPGO_l_SST.assm_ym(mi,:)=-data_NPGO_l_SST.assm_ym(mi,:);
            end
        end
        %% lens2
        for mi=1:len_mem_lens2
            tmp.sign=data_PDO_l.lv_lens2(mi,xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
            if tmp.sign>0
                data_NPGO_l_SST.lens2_ym(mi,:)=-data_NPGO_l_SST.lens2_ym(mi,:);
            end
        end
        %% hcst
        for ly=1:5
            for mi=1:len_mem_oda
                tmp.sign=data_PDO_l.lv_hcst(ly,mi,xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
                if tmp.sign>0
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
        tmp.pct_obs=NaN(732,1);
        tmp.pct_obs(397:end)=data_NPGO.data_PDO.pct_obs(:,lmode);
        tmp.reshp=reshape(tmp.pct_obs, [12, len_m/12]);
        data_NPGO_SSH.obs_ym=squeeze(mean(tmp.reshp,1));
        
        tmp.reshp=reshape(data_NPGO.data_PDO.pct_assm(:,:,lmode), [len_mem_oda, 12, len_m/12]);
        data_NPGO_SSH.assm_ym=squeeze(mean(tmp.reshp,2));
        
        tmp.reshp=reshape(data_NPGO.data_PDO.pct_lens2(:,:,lmode), [len_mem_lens2, 12, len_m/12]);
        data_NPGO_SSH.lens2_ym=squeeze(mean(tmp.reshp,2));
        
        tmp.reshp=reshape(data_NPGO.data_PDO.pct_hcst(:,:,:,lmode), [5, len_mem_oda, 12, len_m/12]);
        data_NPGO_SSH.hcst_ym=squeeze(mean(tmp.reshp,3));

    
        %% x: 70~90 y=20:30 sign - -> positive change
        xrange=70:90;
        yrange=20:30;
        %% obs
        tmp.sign=data_NPGO.data_PDO.lv_obs(xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
        if tmp.sign>0
            data_NPGO_SSH.obs_ym=-data_NPGO_SSH.obs_ym;
        end
        %% assm
        for mi=1:len_mem_oda
            tmp.sign=data_NPGO.data_PDO.lv_assm(mi,xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
            if tmp.sign<0
                data_NPGO_SSH.assm_ym(mi,:)=-data_NPGO_SSH.assm_ym(mi,:);
            end
        end
        %% lens2
        for mi=1:len_mem_lens2
            tmp.sign=data_NPGO.data_PDO.lv_lens2(mi,xrange,yrange,lmode); tmp.sign=mean(tmp.sign(:), 'omitnan');
            if tmp.sign<0
                data_NPGO_SSH.lens2_ym(mi,:)=-data_NPGO_SSH.lens2_ym(mi,:);
            end
        end
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
        tmp.reshp=reshape(data_SAM.obs, [12, len_m/12]);
        data_SAM.obs_ym=squeeze(mean(tmp.reshp,1));
        
        tmp.reshp=reshape(data_SAM.assm, [len_mem_oda, 12, len_m/12]);
        data_SAM.assm_ym=squeeze(mean(tmp.reshp,2));
        
        tmp.reshp=reshape(data_SAM.lens2, [len_mem_lens2, 12, len_m/12]);
        data_SAM.lens2_ym=squeeze(mean(tmp.reshp,2));
        
        tmp.reshp=reshape(data_SAM.hcst, [5, len_mem_oda, 12, len_m/12]);
        data_SAM.hcst_ym=squeeze(mean(tmp.reshp,3));

        data_SAM.obs_4ym=movmean(data_SAM.obs_ym, 4, 'Endpoints', 'discard');
        data_SAM.assm_4ym=movmean(data_SAM.assm_ym, 4, 2, 'Endpoints', 'discard');
        data_SAM.lens2_4ym=movmean(data_SAM.lens2_ym, 4, 2, 'Endpoints', 'discard');
        data_SAM.hcst_4ym=squeeze(mean(data_SAM.hcst_ym(2:5,:,:),1));
    
        %% ATL3
        tmp.reshp=reshape(data_ATL3.obs, [12, len_m/12]);
        data_ATL3.obs_ym=squeeze(mean(tmp.reshp,1));
        
        tmp.reshp=reshape(data_ATL3.assm, [len_mem_oda, 12, len_m/12]);
        data_ATL3.assm_ym=squeeze(mean(tmp.reshp,2));
        
        tmp.reshp=reshape(data_ATL3.lens2, [len_mem_lens2, 12, len_m/12]);
        data_ATL3.lens2_ym=squeeze(mean(tmp.reshp,2));
        
        tmp.reshp=reshape(data_ATL3.hcst, [5, len_mem_oda, 12, len_m/12]);
        data_ATL3.hcst_ym=squeeze(mean(tmp.reshp,3));

        data_ATL3.obs_4ym=movmean(data_ATL3.obs_ym, 4, 'Endpoints', 'discard');
        data_ATL3.assm_4ym=movmean(data_ATL3.assm_ym, 4, 2, 'Endpoints', 'discard');
        data_ATL3.lens2_4ym=movmean(data_ATL3.lens2_ym, 4, 2, 'Endpoints', 'discard');
        data_ATL3.hcst_4ym=squeeze(mean(data_ATL3.hcst_ym(2:5,:,:),1));
    
    
        %% TNA
        tmp.reshp=reshape(data_TNA.obs, [12, len_m/12]);
        data_TNA.obs_ym=squeeze(mean(tmp.reshp,1));
        
        tmp.reshp=reshape(data_TNA.assm, [len_mem_oda, 12, len_m/12]);
        data_TNA.assm_ym=squeeze(mean(tmp.reshp,2));
        
        tmp.reshp=reshape(data_TNA.lens2, [len_mem_lens2, 12, len_m/12]);
        data_TNA.lens2_ym=squeeze(mean(tmp.reshp,2));
        
        tmp.reshp=reshape(data_TNA.hcst, [5, len_mem_oda, 12, len_m/12]);
        data_TNA.hcst_ym=squeeze(mean(tmp.reshp,3));

        data_TNA.obs_4ym=movmean(data_TNA.obs_ym, 4, 'Endpoints', 'discard');
        data_TNA.assm_4ym=movmean(data_TNA.assm_ym, 4, 2, 'Endpoints', 'discard');
        data_TNA.lens2_4ym=movmean(data_TNA.lens2_ym, 4, 2, 'Endpoints', 'discard');
        data_TNA.hcst_4ym=squeeze(mean(data_TNA.hcst_ym(2:5,:,:),1));
    
    
        %% TSA
        tmp.reshp=reshape(data_TSA.obs, [12, len_m/12]);
        data_TSA.obs_ym=squeeze(mean(tmp.reshp,1));
        
        tmp.reshp=reshape(data_TSA.assm, [len_mem_oda, 12, len_m/12]);
        data_TSA.assm_ym=squeeze(mean(tmp.reshp,2));
        
        tmp.reshp=reshape(data_TSA.lens2, [len_mem_lens2, 12, len_m/12]);
        data_TSA.lens2_ym=squeeze(mean(tmp.reshp,2));
        
        tmp.reshp=reshape(data_TSA.hcst, [5, len_mem_oda, 12, len_m/12]);
        data_TSA.hcst_ym=squeeze(mean(tmp.reshp,3));

        data_TSA.obs_4ym=movmean(data_TSA.obs_ym, 4, 'Endpoints', 'discard');
        data_TSA.assm_4ym=movmean(data_TSA.assm_ym, 4, 2, 'Endpoints', 'discard');
        data_TSA.lens2_4ym=movmean(data_TSA.lens2_ym, 4, 2, 'Endpoints', 'discard');
        data_TSA.hcst_4ym=squeeze(mean(data_TSA.hcst_ym(2:5,:,:),1));

        once=0;
    end



    
    %% corr between data and indices
    
%     tmp.corr=corrcoef(squeeze(data_obs.([cfg_ts.var,'_ym'])), squeeze(data_AMO.obs_ym), 'Rows', 'complete');
    
%     hold on
%     for mi=1:20
%         plot(data_PDO.assm_ym(mi,:));
%     end

clim_indices={'AMO', 'ENSO', 'PDO', 'PDO_l', 'NPGO_SST', 'NPGO_l_SST', 'NPGO_SSH', 'SAM', 'ATL3', 'TNA', 'TSA'};

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


%         tmp.corr=corrcoef(squeeze(m_data_assm.([cfg_ts.var,'_ym'])(mi,:)), squeeze(data_AMO.assm_ym(mi,:)), 'Rows', 'complete');
%         corr_all.assm_ym_AMO(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(squeeze(m_data_assm.([cfg_ts.var,'_ym'])(mi,:)), squeeze(data_ENSO.assm_ym(mi,:)), 'Rows', 'complete');
%         corr_all.assm_ym_ENSO(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(squeeze(m_data_assm.([cfg_ts.var,'_ym'])(mi,:)), squeeze(data_PDO.assm_ym(mi,:)), 'Rows', 'complete');
%         corr_all.assm_ym_PDO(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(squeeze(m_data_assm.([cfg_ts.var,'_ym'])(mi,:)), squeeze(data_PDO_l.assm_ym(mi,:)), 'Rows', 'complete');
%         corr_all.assm_ym_PDO_l(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(squeeze(m_data_assm.([cfg_ts.var,'_ym'])(mi,:)), squeeze(data_NPGO_SST.assm_ym(mi,:)), 'Rows', 'complete');
%         corr_all.assm_ym_NPGO_SST(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(squeeze(m_data_assm.([cfg_ts.var,'_ym'])(mi,:)), squeeze(data_NPGO_l_SST.assm_ym(mi,:)), 'Rows', 'complete');
%         corr_all.assm_ym_NPGO_l_SST(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(squeeze(m_data_assm.([cfg_ts.var,'_ym'])(mi,:)), squeeze(data_NPGO_SSH.assm_ym(mi,:)), 'Rows', 'complete');
%         corr_all.assm_ym_NPGO_SSH(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(squeeze(m_data_assm.([cfg_ts.var,'_ym'])(mi,:)), squeeze(data_SAM.assm_ym(mi,:)), 'Rows', 'complete');
%         corr_all.assm_ym_SAM(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(squeeze(m_data_assm.([cfg_ts.var,'_ym'])(mi,:)), squeeze(data_ATL3.assm_ym(mi,:)), 'Rows', 'complete');
%         corr_all.assm_ym_ATL3(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(squeeze(m_data_assm.([cfg_ts.var,'_ym'])(mi,:)), squeeze(data_TNA.assm_ym(mi,:)), 'Rows', 'complete');
%         corr_all.assm_ym_TNA(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(squeeze(m_data_assm.([cfg_ts.var,'_ym'])(mi,:)), squeeze(data_TSA.assm_ym(mi,:)), 'Rows', 'complete');
%         corr_all.assm_ym_TSA(mi)=tmp.corr(1,2);
        
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

%         tmp.corr=corrcoef(tmp.data_4ym, squeeze(data_AMO.assm_4ym(mi,:)), 'Rows', 'complete');
%         corr_all.assm_4ym_AMO(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(tmp.data_4ym, squeeze(data_ENSO.assm_4ym(mi,:)), 'Rows', 'complete');
%         corr_all.assm_4ym_ENSO(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(tmp.data_4ym, squeeze(data_PDO.assm_4ym(mi,:)), 'Rows', 'complete');
%         corr_all.assm_4ym_PDO(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(tmp.data_4ym, squeeze(data_PDO_l.assm_4ym(mi,:)), 'Rows', 'complete');
%         corr_all.assm_4ym_PDO_l(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(tmp.data_4ym, squeeze(data_NPGO_SST.assm_4ym(mi,:)), 'Rows', 'complete');
%         corr_all.assm_4ym_NPGO_SST(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(tmp.data_4ym, squeeze(data_NPGO_l_SST.assm_4ym(mi,:)), 'Rows', 'complete');
%         corr_all.assm_4ym_NPGO_l_SST(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(tmp.data_4ym, squeeze(data_NPGO_SSH.assm_4ym(mi,:)), 'Rows', 'complete');
%         corr_all.assm_4ym_NPGO_SSH(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(tmp.data_4ym, squeeze(data_SAM.assm_4ym(mi,:)), 'Rows', 'complete');
%         corr_all.assm_4ym_SAM(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(tmp.data_4ym, squeeze(data_ATL3.assm_4ym(mi,:)), 'Rows', 'complete');
%         corr_all.assm_4ym_ATL3(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(tmp.data_4ym, squeeze(data_TNA.assm_4ym(mi,:)), 'Rows', 'complete');
%         corr_all.assm_4ym_TNA(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(tmp.data_4ym, squeeze(data_TSA.assm_4ym(mi,:)), 'Rows', 'complete');
%         corr_all.assm_4ym_TSA(mi)=tmp.corr(1,2);
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

%     tmp.corr=corrcoef(squeeze(mean(m_data_assm.([cfg_ts.var,'_ym']),1)), squeeze(mean(data_AMO.assm_ym,1)), 'Rows', 'complete');
%     corr_all_em.assm_ym_AMO=tmp.corr(1,2);
%     tmp.corr=corrcoef(squeeze(mean(m_data_assm.([cfg_ts.var,'_ym']),1)), squeeze(mean(data_ENSO.assm_ym,1)), 'Rows', 'complete');
%     corr_all_em.assm_ym_ENSO=tmp.corr(1,2);
%     tmp.corr=corrcoef(squeeze(mean(m_data_assm.([cfg_ts.var,'_ym']),1)), squeeze(mean(data_PDO.assm_ym,1)), 'Rows', 'complete');
%     corr_all_em.assm_ym_PDO=tmp.corr(1,2);
%     tmp.corr=corrcoef(squeeze(mean(m_data_assm.([cfg_ts.var,'_ym']),1)), squeeze(mean(data_NPGO_SST.assm_ym,1)), 'Rows', 'complete');
%     corr_all_em.assm_ym_NPGO_SST=tmp.corr(1,2);
%     tmp.corr=corrcoef(squeeze(mean(m_data_assm.([cfg_ts.var,'_ym']),1)), squeeze(mean(data_NPGO_SSH.assm_ym,1)), 'Rows', 'complete');
%     corr_all_em.assm_ym_NPGO_SSH=tmp.corr(1,2);
%     tmp.corr=corrcoef(squeeze(mean(m_data_assm.([cfg_ts.var,'_ym']),1)), squeeze(mean(data_SAM.assm_ym,1)), 'Rows', 'complete');
%     corr_all_em.assm_ym_SAM=tmp.corr(1,2);
%     tmp.corr=corrcoef(squeeze(mean(m_data_assm.([cfg_ts.var,'_ym']),1)), squeeze(mean(data_ATL3.assm_ym,1)), 'Rows', 'complete');
%     corr_all_em.assm_ym_ATL3=tmp.corr(1,2);
%     tmp.corr=corrcoef(squeeze(mean(m_data_assm.([cfg_ts.var,'_ym']),1)), squeeze(mean(data_TNA.assm_ym,1)), 'Rows', 'complete');
%     corr_all_em.assm_ym_TNA=tmp.corr(1,2);
%     tmp.corr=corrcoef(squeeze(mean(m_data_assm.([cfg_ts.var,'_ym']),1)), squeeze(mean(data_TSA.assm_ym,1)), 'Rows', 'complete');
%     corr_all_em.assm_ym_TSA=tmp.corr(1,2);
    
    
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

%         tmp.corr=corrcoef(squeeze(m_data_lens2.([cfg_ts.var,'_ym'])(mi,:)), squeeze(data_AMO.lens2_ym(mi,:)), 'Rows', 'complete');
%         corr_all.lens2_ym_AMO(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(squeeze(m_data_lens2.([cfg_ts.var,'_ym'])(mi,:)), squeeze(data_ENSO.lens2_ym(mi,:)), 'Rows', 'complete');
%         corr_all.lens2_ym_ENSO(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(squeeze(m_data_lens2.([cfg_ts.var,'_ym'])(mi,:)), squeeze(data_PDO.lens2_ym(mi,:)), 'Rows', 'complete');
%         corr_all.lens2_ym_PDO(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(squeeze(m_data_lens2.([cfg_ts.var,'_ym'])(mi,:)), squeeze(data_NPGO_SST.lens2_ym(mi,:)), 'Rows', 'complete');
%         corr_all.lens2_ym_NPGO_SST(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(squeeze(m_data_lens2.([cfg_ts.var,'_ym'])(mi,:)), squeeze(data_NPGO_SSH.lens2_ym(mi,:)), 'Rows', 'complete');
%         corr_all.lens2_ym_NPGO_SSH(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(squeeze(m_data_lens2.([cfg_ts.var,'_ym'])(mi,:)), squeeze(data_SAM.lens2_ym(mi,:)), 'Rows', 'complete');
%         corr_all.lens2_ym_SAM(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(squeeze(m_data_lens2.([cfg_ts.var,'_ym'])(mi,:)), squeeze(data_ATL3.lens2_ym(mi,:)), 'Rows', 'complete');
%         corr_all.lens2_ym_ATL3(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(squeeze(m_data_lens2.([cfg_ts.var,'_ym'])(mi,:)), squeeze(data_TNA.lens2_ym(mi,:)), 'Rows', 'complete');
%         corr_all.lens2_ym_TNA(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(squeeze(m_data_lens2.([cfg_ts.var,'_ym'])(mi,:)), squeeze(data_TSA.lens2_ym(mi,:)), 'Rows', 'complete');
%         corr_all.lens2_ym_TSA(mi)=tmp.corr(1,2);

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

%         tmp.corr=corrcoef(tmp.data_4ym, squeeze(data_AMO.lens2_4ym(mi,:)), 'Rows', 'complete');
%         corr_all.lens2_4ym_AMO(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(tmp.data_4ym, squeeze(data_ENSO.lens2_4ym(mi,:)), 'Rows', 'complete');
%         corr_all.lens2_4ym_ENSO(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(tmp.data_4ym, squeeze(data_PDO.lens2_4ym(mi,:)), 'Rows', 'complete');
%         corr_all.lens2_4ym_PDO(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(tmp.data_4ym, squeeze(data_NPGO_SST.lens2_4ym(mi,:)), 'Rows', 'complete');
%         corr_all.lens2_4ym_NPGO_SST(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(tmp.data_4ym, squeeze(data_NPGO_SSH.lens2_4ym(mi,:)), 'Rows', 'complete');
%         corr_all.lens2_4ym_NPGO_SSH(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(tmp.data_4ym, squeeze(data_SAM.lens2_4ym(mi,:)), 'Rows', 'complete');
%         corr_all.lens2_4ym_SAM(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(tmp.data_4ym, squeeze(data_ATL3.lens2_4ym(mi,:)), 'Rows', 'complete');
%         corr_all.lens2_4ym_ATL3(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(tmp.data_4ym, squeeze(data_TNA.lens2_4ym(mi,:)), 'Rows', 'complete');
%         corr_all.lens2_4ym_TNA(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(tmp.data_4ym, squeeze(data_TSA.lens2_4ym(mi,:)), 'Rows', 'complete');
%         corr_all.lens2_4ym_TSA(mi)=tmp.corr(1,2);
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

%     tmp.corr=corrcoef(squeeze(mean(m_data_lens2.([cfg_ts.var,'_ym']),1)), squeeze(mean(data_AMO.lens2_ym,1)), 'Rows', 'complete');
%     corr_all_em.lens2_ym_AMO=tmp.corr(1,2);
%     tmp.corr=corrcoef(squeeze(mean(m_data_lens2.([cfg_ts.var,'_ym']),1)), squeeze(mean(data_ENSO.lens2_ym,1)), 'Rows', 'complete');
%     corr_all_em.lens2_ym_ENSO=tmp.corr(1,2);
%     tmp.corr=corrcoef(squeeze(mean(m_data_lens2.([cfg_ts.var,'_ym']),1)), squeeze(mean(data_PDO.lens2_ym,1)), 'Rows', 'complete');
%     corr_all_em.lens2_ym_PDO=tmp.corr(1,2);
%     tmp.corr=corrcoef(squeeze(mean(m_data_lens2.([cfg_ts.var,'_ym']),1)), squeeze(mean(data_NPGO_SST.lens2_ym,1)), 'Rows', 'complete');
%     corr_all_em.lens2_ym_NPGO_SST=tmp.corr(1,2);
%     tmp.corr=corrcoef(squeeze(mean(m_data_lens2.([cfg_ts.var,'_ym']),1)), squeeze(mean(data_NPGO_SSH.lens2_ym,1)), 'Rows', 'complete');
%     corr_all_em.lens2_ym_NPGO_SSH=tmp.corr(1,2);
%     tmp.corr=corrcoef(squeeze(mean(m_data_lens2.([cfg_ts.var,'_ym']),1)), squeeze(mean(data_SAM.lens2_ym,1)), 'Rows', 'complete');
%     corr_all_em.lens2_ym_SAM=tmp.corr(1,2);
%     tmp.corr=corrcoef(squeeze(mean(m_data_lens2.([cfg_ts.var,'_ym']),1)), squeeze(mean(data_ATL3.lens2_ym,1)), 'Rows', 'complete');
%     corr_all_em.lens2_ym_ATL3=tmp.corr(1,2);
%     tmp.corr=corrcoef(squeeze(mean(m_data_lens2.([cfg_ts.var,'_ym']),1)), squeeze(mean(data_TNA.lens2_ym,1)), 'Rows', 'complete');
%     corr_all_em.lens2_ym_TNA=tmp.corr(1,2);
%     tmp.corr=corrcoef(squeeze(mean(m_data_lens2.([cfg_ts.var,'_ym']),1)), squeeze(mean(data_TSA.lens2_ym,1)), 'Rows', 'complete');
%     corr_all_em.lens2_ym_TSA=tmp.corr(1,2);
    
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

%             tmp.corr=corrcoef(squeeze(m_data_hcst.([cfg_ts.var,'_ym']).(ly_str)(mi,:)), squeeze(data_AMO.hcst_ym(ly,mi,:)), 'Rows', 'complete');
%             corr_all.hcst_ym_AMO(ly,mi)=tmp.corr(1,2);
%             tmp.corr=corrcoef(squeeze(m_data_hcst.([cfg_ts.var,'_ym']).(ly_str)(mi,:)), squeeze(data_ENSO.hcst_ym(ly,mi,:)), 'Rows', 'complete');
%             corr_all.hcst_ym_ENSO(ly,mi)=tmp.corr(1,2);
%             tmp.corr=corrcoef(squeeze(m_data_hcst.([cfg_ts.var,'_ym']).(ly_str)(mi,:)), squeeze(data_PDO.hcst_ym(ly,mi,:)), 'Rows', 'complete');
%             corr_all.hcst_ym_PDO(ly,mi)=tmp.corr(1,2);
%             tmp.corr=corrcoef(squeeze(m_data_hcst.([cfg_ts.var,'_ym']).(ly_str)(mi,:)), squeeze(data_NPGO_SST.hcst_ym(ly,mi,:)), 'Rows', 'complete');
%             corr_all.hcst_ym_NPGO_SST(ly,mi)=tmp.corr(1,2);
%             tmp.corr=corrcoef(squeeze(m_data_hcst.([cfg_ts.var,'_ym']).(ly_str)(mi,:)), squeeze(data_NPGO_SSH.hcst_ym(ly,mi,:)), 'Rows', 'complete');
%             corr_all.hcst_ym_NPGO_SSH(ly,mi)=tmp.corr(1,2);
%             tmp.corr=corrcoef(squeeze(m_data_hcst.([cfg_ts.var,'_ym']).(ly_str)(mi,:)), squeeze(data_SAM.hcst_ym(ly,mi,:)), 'Rows', 'complete');
%             corr_all.hcst_ym_SAM(ly,mi)=tmp.corr(1,2);
%             tmp.corr=corrcoef(squeeze(m_data_hcst.([cfg_ts.var,'_ym']).(ly_str)(mi,:)), squeeze(data_ATL3.hcst_ym(ly,mi,:)), 'Rows', 'complete');
%             corr_all.hcst_ym_ATL3(ly,mi)=tmp.corr(1,2);
%             tmp.corr=corrcoef(squeeze(m_data_hcst.([cfg_ts.var,'_ym']).(ly_str)(mi,:)), squeeze(data_TNA.hcst_ym(ly,mi,:)), 'Rows', 'complete');
%             corr_all.hcst_ym_TNA(ly,mi)=tmp.corr(1,2);
%             tmp.corr=corrcoef(squeeze(m_data_hcst.([cfg_ts.var,'_ym']).(ly_str)(mi,:)), squeeze(data_TSA.hcst_ym(ly,mi,:)), 'Rows', 'complete');
%             corr_all.hcst_ym_TSA(ly,mi)=tmp.corr(1,2);

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

%         tmp.corr=corrcoef(squeeze(mean(m_data_hcst.([cfg_ts.var,'_ym']).(ly_str),1)), squeeze(mean(data_AMO.hcst_ym(ly,:,:),2)), 'Rows', 'complete');
%         corr_all_em.hcst_ym_AMO(ly)=tmp.corr(1,2);
%         tmp.corr=corrcoef(squeeze(mean(m_data_hcst.([cfg_ts.var,'_ym']).(ly_str),1)), squeeze(mean(data_ENSO.hcst_ym(ly,:,:),2)), 'Rows', 'complete');
%         corr_all_em.hcst_ym_ENSO(ly)=tmp.corr(1,2);
%         tmp.corr=corrcoef(squeeze(mean(m_data_hcst.([cfg_ts.var,'_ym']).(ly_str),1)), squeeze(mean(data_PDO.hcst_ym(ly,:,:),2)), 'Rows', 'complete');
%         corr_all_em.hcst_ym_PDO(ly)=tmp.corr(1,2);
%         tmp.corr=corrcoef(squeeze(mean(m_data_hcst.([cfg_ts.var,'_ym']).(ly_str),1)), squeeze(mean(data_NPGO_SST.hcst_ym(ly,:,:),2)), 'Rows', 'complete');
%         corr_all_em.hcst_ym_NPGO_SST(ly)=tmp.corr(1,2);
%         tmp.corr=corrcoef(squeeze(mean(m_data_hcst.([cfg_ts.var,'_ym']).(ly_str),1)), squeeze(mean(data_NPGO_SSH.hcst_ym(ly,:,:),2)), 'Rows', 'complete');
%         corr_all_em.hcst_ym_NPGO_SSH(ly)=tmp.corr(1,2);
%         tmp.corr=corrcoef(squeeze(mean(m_data_hcst.([cfg_ts.var,'_ym']).(ly_str),1)), squeeze(mean(data_SAM.hcst_ym(ly,:,:),2)), 'Rows', 'complete');
%         corr_all_em.hcst_ym_SAM(ly)=tmp.corr(1,2);
%         tmp.corr=corrcoef(squeeze(mean(m_data_hcst.([cfg_ts.var,'_ym']).(ly_str),1)), squeeze(mean(data_ATL3.hcst_ym(ly,:,:),2)), 'Rows', 'complete');
%         corr_all_em.hcst_ym_ATL3(ly)=tmp.corr(1,2);
%         tmp.corr=corrcoef(squeeze(mean(m_data_hcst.([cfg_ts.var,'_ym']).(ly_str),1)), squeeze(mean(data_TNA.hcst_ym(ly,:,:),2)), 'Rows', 'complete');
%         corr_all_em.hcst_ym_TNA(ly)=tmp.corr(1,2);
%         tmp.corr=corrcoef(squeeze(mean(m_data_hcst.([cfg_ts.var,'_ym']).(ly_str),1)), squeeze(mean(data_TSA.hcst_ym(ly,:,:),2)), 'Rows', 'complete');
%         corr_all_em.hcst_ym_TSA(ly)=tmp.corr(1,2);
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

%         tmp.corr=corrcoef(tmp.data_4ym_hc(mi,:), squeeze(data_AMO.hcst_4ym(mi,:)), 'Rows', 'complete');
%         corr_all.hcst_4ym_AMO(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(tmp.data_4ym_hc(mi,:), squeeze(data_ENSO.hcst_4ym(mi,:)), 'Rows', 'complete');
%         corr_all.hcst_4ym_ENSO(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(tmp.data_4ym_hc(mi,:), squeeze(data_PDO.hcst_4ym(mi,:)), 'Rows', 'complete');
%         corr_all.hcst_4ym_PDO(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(tmp.data_4ym_hc(mi,:), squeeze(data_NPGO_SST.hcst_4ym(mi,:)), 'Rows', 'complete');
%         corr_all.hcst_4ym_NPGO_SST(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(tmp.data_4ym_hc(mi,:), squeeze(data_NPGO_SSH.hcst_4ym(mi,:)), 'Rows', 'complete');
%         corr_all.hcst_4ym_NPGO_SSH(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(tmp.data_4ym_hc(mi,:), squeeze(data_SAM.hcst_4ym(mi,:)), 'Rows', 'complete');
%         corr_all.hcst_4ym_SAM(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(tmp.data_4ym_hc(mi,:), squeeze(data_ATL3.hcst_4ym(mi,:)), 'Rows', 'complete');
%         corr_all.hcst_4ym_ATL3(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(tmp.data_4ym_hc(mi,:), squeeze(data_TNA.hcst_4ym(mi,:)), 'Rows', 'complete');
%         corr_all.hcst_4ym_TNA(mi)=tmp.corr(1,2);
%         tmp.corr=corrcoef(tmp.data_4ym_hc(mi,:), squeeze(data_TSA.hcst_4ym(mi,:)), 'Rows', 'complete');
%         corr_all.hcst_4ym_TSA(mi)=tmp.corr(1,2);
    end
    
    % % % % % %% ASSM_AMO (1~2 row)
    % % % % % tmp.rows=1;
    % % % % % tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_ym_AMO;
    % % % % % tmp.fig_mat(tmp.rows,23:52)=NaN;
    % % % % % tmp.rows=tmp.rows+1;
    % % % % % tmp.fig_mat(tmp.rows,1)=corr_all_em.assm_ym_AMO;
    % % % % % tmp.fig_mat(tmp.rows,2:50)=NaN;
    % % % % % 
    % % % % % %% HCST_AMO (3~12 row)
    % % % % % tmp.rows=tmp.rows+1;
    % % % % % tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_AMO(1,:);
    % % % % % tmp.fig_mat(tmp.rows,23:52)=NaN;
    % % % % % tmp.rows=tmp.rows+1;
    % % % % % tmp.fig_mat(tmp.rows,1)=corr_all_em.hcst_ym_AMO(1);
    % % % % % tmp.fig_mat(tmp.rows,2:50)=NaN;
    % % % % % tmp.rows=tmp.rows+1;
    % % % % % tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_AMO(2,:);
    % % % % % tmp.fig_mat(tmp.rows,23:52)=NaN;
    % % % % % tmp.rows=tmp.rows+1;
    % % % % % tmp.fig_mat(tmp.rows,1)=corr_all_em.hcst_ym_AMO(2);
    % % % % % tmp.fig_mat(tmp.rows,2:50)=NaN;
    % % % % % tmp.rows=tmp.rows+1;
    % % % % % tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_AMO(3,:);
    % % % % % tmp.fig_mat(tmp.rows,23:52)=NaN;
    % % % % % tmp.rows=tmp.rows+1;
    % % % % % tmp.fig_mat(tmp.rows,1)=corr_all_em.hcst_ym_AMO(3);
    % % % % % tmp.fig_mat(tmp.rows,2:50)=NaN;
    % % % % % tmp.rows=tmp.rows+1;
    % % % % % tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_AMO(4,:);
    % % % % % tmp.fig_mat(tmp.rows,23:52)=NaN;
    % % % % % tmp.rows=tmp.rows+1;
    % % % % % tmp.fig_mat(tmp.rows,1)=corr_all_em.hcst_ym_AMO(4);
    % % % % % tmp.fig_mat(tmp.rows,2:50)=NaN;
    % % % % % tmp.rows=tmp.rows+1;
    % % % % % tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_AMO(5,:);
    % % % % % tmp.fig_mat(tmp.rows,23:52)=NaN;
    % % % % % tmp.rows=tmp.rows+1;
    % % % % % tmp.fig_mat(tmp.rows,1)=corr_all_em.hcst_ym_AMO(5);
    % % % % % tmp.fig_mat(tmp.rows,2:50)=NaN;
    % % % % % 
    % % % % % %% LENS2_AMO (13~14 row)
    % % % % % tmp.rows=tmp.rows+1;
    % % % % % tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_ym_AMO;
    % % % % % tmp.rows=tmp.rows+1;
    % % % % % tmp.fig_mat(tmp.rows,1)=corr_all_em.lens2_ym_AMO;
    % % % % % tmp.fig_mat(tmp.rows,2:50)=NaN;
    
%% raw
    %% ASSM_AMO (1~2 row)
    tmp.rows=1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.assm_ym_AMO); tmp.fig_mat(tmp.rows,2)=NaN;
    tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_ym_AMO;
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    %% HCST_AMO (3~12 row)
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_AMO(1,:)); tmp.fig_mat(tmp.rows,2)=NaN;
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_AMO(1,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_AMO(2,:)); tmp.fig_mat(tmp.rows,2)=NaN;
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_AMO(2,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_AMO(3,:)); tmp.fig_mat(tmp.rows,2)=NaN;
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_AMO(3,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_AMO(4,:)); tmp.fig_mat(tmp.rows,2)=NaN;
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_AMO(4,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_AMO(5,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_AMO(5,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    %% LENS2_AMO (13~14 row)
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.lens2_ym_AMO); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_ym_AMO;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1:52)=NaN;
    
    %% ASSM_ENSO 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.assm_ym_ENSO); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_ym_ENSO;
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    %% HCST_ENSO 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_ENSO(1,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_ENSO(1,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_ENSO(2,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_ENSO(2,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_ENSO(3,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_ENSO(3,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_ENSO(4,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_ENSO(4,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_ENSO(5,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_ENSO(5,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    %% LENS2_ENSO 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.lens2_ym_ENSO); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_ym_ENSO;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1:52)=NaN;
    
    
    %% ASSM_PDO 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.assm_ym_PDO); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_ym_PDO;
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    %% HCST_PDO 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_PDO(1,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_PDO(1,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_PDO(2,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_PDO(2,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_PDO(3,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_PDO(3,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_PDO(4,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_PDO(4,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_PDO(5,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_PDO(5,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    %% LENS2_PDO 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.lens2_ym_PDO); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_ym_PDO;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1:52)=NaN;
    
    
    %% ASSM_NPGO_SST 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.assm_ym_NPGO_SST); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_ym_NPGO_SST;
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    %% HCST_NPGO_SST 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_NPGO_SST(1,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_NPGO_SST(1,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_NPGO_SST(2,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_NPGO_SST(2,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_NPGO_SST(3,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_NPGO_SST(3,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_NPGO_SST(4,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_NPGO_SST(4,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_NPGO_SST(5,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_NPGO_SST(5,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    %% LENS2_NPGO_SST 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.lens2_ym_NPGO_SST); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_ym_NPGO_SST;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1:52)=NaN;
    
    
    %% ASSM_NPGO_SSH 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.assm_ym_NPGO_SSH); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_ym_NPGO_SSH;
    tmp.fig_mat(tmp.rows,23:52)=NaN;


    %% HCST_NPGO_SSH 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_NPGO_SSH(1,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_NPGO_SSH(1,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_NPGO_SSH(2,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_NPGO_SSH(2,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_NPGO_SSH(3,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_NPGO_SSH(3,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_NPGO_SSH(4,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_NPGO_SSH(4,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_NPGO_SSH(5,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_NPGO_SSH(5,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    %% LENS2_NPGO_SSH 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.lens2_ym_NPGO_SSH); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_ym_NPGO_SSH;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1:52)=NaN;
    
    
    %% ASSM_SAM 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.assm_ym_SAM); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_ym_SAM;
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    %% HCST_SAM 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_SAM(1,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_SAM(1,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_SAM(2,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_SAM(2,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_SAM(3,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_SAM(3,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_SAM(4,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_SAM(4,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_SAM(5,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_SAM(5,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    %% LENS2_SAM 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.lens2_ym_SAM); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_ym_SAM;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1:52)=NaN;


    %% ASSM_ATL3 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.assm_ym_ATL3); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_ym_ATL3;
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    %% HCST_ATL3 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_ATL3(1,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_ATL3(1,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_ATL3(2,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_ATL3(2,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_ATL3(3,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_ATL3(3,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_ATL3(4,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_ATL3(4,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_ATL3(5,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_ATL3(5,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    %% LENS2_ATL3 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.lens2_ym_ATL3); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_ym_ATL3;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1:52)=NaN;



    %% ASSM_TNA 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.assm_ym_TNA); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_ym_TNA;
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    %% HCST_TNA 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_TNA(1,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_TNA(1,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_TNA(2,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_TNA(2,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_TNA(3,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_TNA(3,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_TNA(4,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_TNA(4,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_TNA(5,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_TNA(5,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    %% LENS2_TNA 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.lens2_ym_TNA); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_ym_TNA;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1:52)=NaN;


    %% ASSM_TSA 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.assm_ym_TSA); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_ym_TSA;
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    %% HCST_TSA 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_TSA(1,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_TSA(1,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_TSA(2,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_TSA(2,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_TSA(3,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_TSA(3,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_TSA(4,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_TSA(4,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_ym_TSA(5,:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_ym_TSA(5,:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    %% LENS2_TSA 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.lens2_ym_TSA); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_ym_TSA;
    tmp.rows=tmp.rows+1;


    
    tmp.fig_mat(tmp.rows,1:52)=NaN;
    

    fig_h = figure('name','ts','visible','off');
        fig_h.Position= [0,0,1000,800];

    pcolor(tmp.fig_mat); shading flat; colorbar; caxis([-1 1]); colormap(fig_cfg.c_map);
    yticks([1,2,7, 9,10,15, 17,18,23, 25,26,31, 33,34,39, 41,42,47, 49,50,55, 57,58,63, 65,66,71])
    yticklabels({'ODA-AMO', 'HIND-AMO','LE-AMO', ...
        'ODA-ENSO', 'HIND-ENSO','LE-ENSO', ...
        'ODA-PDO', 'HIND-PDO','LE-PDO', ...
        'ODA-NPGO-SST', 'HIND-NPGO-SST','LE-NPGO-SST', ...
        'ODA-NPGO-SSH', 'HIND-NPGO-SSH','LE-NPGO-SSH', ...
        'ODA-SAM', 'HIND-SAM','LE-SAM', ...
        'ODA-ATL3', 'HIND-ATL3','LE-ATL3', ...
        'ODA-TNA', 'HIND-TNA','LE-TNA', ...
        'ODA-TSA', 'HIND-TSA','LE-TSA'})
    xticks([1, 7:5:52]);
    xticklabels({'mean', '5', '10', '15', '20', '25', '30', '35', '40', '45', '50'})


    xlabel(['members'])
%     title([cfg_ts.var, ', ', 'r with indices', ])
    
%     if length(sta_lonlat{stai})==2
        title(['r with indices', ', ', cfg_ts.var, ', ', num2str(sta_lonlat{stai}(1)),'E, ', num2str(sta_lonlat{stai}(3)), 'N'])
%     elseif length(sta_lonlat{stai})==4
%         title(['r with indices', ', ', cfg_ts.var, ', ', num2str(sta_lonlat{stai}(1)),'E~', num2str(sta_lonlat{stai}(2)),'E, ', ...
%             num2str(sta_lonlat{stai}(3)), 'N~',num2str(sta_lonlat{stai}(4)), 'N'])
%     end
    
    dirs.figdir= [dirs.figroot, filesep, 'ens_all', filesep, cfg_ts.var, '_corr_indices'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'r_indices_', num2str(sta_lonlat{stai}(1)), 'E_', num2str(sta_lonlat{stai}(3)), 'N_', cfg_ts.var, '.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;

%% 4y movmean
    tmp=rmfield(tmp, 'fig_mat');
%% ASSM_AMO (1~2 row)
    tmp.rows=1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.assm_4ym_AMO); tmp.fig_mat(tmp.rows,2)=NaN;
    tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_4ym_AMO;
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    %% HCST_AMO (3~12 row)
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_4ym_AMO(:)); tmp.fig_mat(tmp.rows,2)=NaN;
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_4ym_AMO(:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
        
    %% LENS2_AMO (13~14 row)
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.lens2_4ym_AMO); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_4ym_AMO;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1:52)=NaN;
    
    %% ASSM_ENSO 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.assm_4ym_ENSO); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_4ym_ENSO;
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    %% HCST_ENSO 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_4ym_ENSO(:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_4ym_ENSO(:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    %% LENS2_ENSO 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.lens2_4ym_ENSO); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_4ym_ENSO;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1:52)=NaN;
    
    %% ASSM_PDO 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.assm_4ym_PDO); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_4ym_PDO;
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    %% HCST_PDO 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_4ym_PDO(:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_4ym_PDO(:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    %% LENS2_PDO 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.lens2_4ym_PDO); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_4ym_PDO;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1:52)=NaN;


    %% ASSM_PDO_l 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.assm_4ym_PDO_l); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_4ym_PDO_l;
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    %% HCST_PDO_l 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_4ym_PDO_l(:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_4ym_PDO_l(:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    %% LENS2_PDO_l 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.lens2_4ym_PDO_l); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_4ym_PDO_l;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1:52)=NaN;
    
    
    %% ASSM_NPGO_SST 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.assm_4ym_NPGO_SST); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_4ym_NPGO_SST;
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    %% HCST_NPGO_SST 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_4ym_NPGO_SST(:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_4ym_NPGO_SST(:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
   
    
    %% LENS2_NPGO_SST 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.lens2_4ym_NPGO_SST); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_4ym_NPGO_SST;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1:52)=NaN;
    
    
    %% ASSM_NPGO_SSH 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.assm_4ym_NPGO_SSH); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_4ym_NPGO_SSH;
    tmp.fig_mat(tmp.rows,23:52)=NaN;


    %% HCST_NPGO_SSH 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_4ym_NPGO_SSH(:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_4ym_NPGO_SSH(:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
   
    
    %% LENS2_NPGO_SSH 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.lens2_4ym_NPGO_SSH); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_4ym_NPGO_SSH;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1:52)=NaN;
    
    
    %% ASSM_SAM 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.assm_4ym_SAM); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_4ym_SAM;
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    %% HCST_SAM 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_4ym_SAM(:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_4ym_SAM(:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    
    %% LENS2_SAM 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.lens2_4ym_SAM); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_4ym_SAM;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1:52)=NaN;


    %% ASSM_ATL3 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.assm_4ym_ATL3); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_4ym_ATL3;
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    %% HCST_ATL3 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_4ym_ATL3(:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_4ym_ATL3(:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    
    %% LENS2_ATL3 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.lens2_4ym_ATL3); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_4ym_ATL3;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1:52)=NaN;



    %% ASSM_TNA 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.assm_4ym_TNA); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_4ym_TNA;
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    %% HCST_TNA 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_4ym_TNA(:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_4ym_TNA(:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    
    %% LENS2_TNA 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.lens2_4ym_TNA); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_4ym_TNA;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1:52)=NaN;


    %% ASSM_TSA 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.assm_4ym_TSA); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.assm_4ym_TSA;
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    %% HCST_TSA 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.hcst_4ym_TSA(:)); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:22)=corr_all.hcst_4ym_TSA(:);
    tmp.fig_mat(tmp.rows,23:52)=NaN;
    
    
    %% LENS2_TSA 
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1)=mean(corr_all.lens2_4ym_TSA); tmp.fig_mat(tmp.rows,2)=NaN;    
    tmp.fig_mat(tmp.rows,3:52)=corr_all.lens2_4ym_TSA;
    tmp.rows=tmp.rows+1;
    tmp.fig_mat(tmp.rows,1:52)=NaN;
    

    fig_h = figure('name','ts','visible','off');
        fig_h.Position= [0,0,1000,500];

    pcolor(tmp.fig_mat); shading flat; colorbar; caxis([-1 1]); colormap(fig_cfg.c_map);
    yticks([1,2,3, 5,6,7, 9,10,11, 13,14,15, 17,18,19, 21,22,23, 25,26,27, 29,30,31, 33,34,35])
    yticklabels({'ODA-AMO', 'HIND-AMO','LE-AMO', ...
        'ODA-ENSO', 'HIND-ENSO','LE-ENSO', ...
        'ODA-PDO', 'HIND-PDO','LE-PDO', ...
        'ODA-PDO_l', 'HIND-PDO_l','LE-PDO_l', ...
        'ODA-NPGO-SST', 'HIND-NPGO-SST','LE-NPGO-SST', ...
        'ODA-NPGO-SSH', 'HIND-NPGO-SSH','LE-NPGO-SSH', ...
        'ODA-SAM', 'HIND-SAM','LE-SAM', ...
        'ODA-ATL3', 'HIND-ATL3','LE-ATL3', ...
        'ODA-TNA', 'HIND-TNA','LE-TNA', ...
        'ODA-TSA', 'HIND-TSA','LE-TSA'})
    xticks([1, 7:5:52]);
    xticklabels({'mean', '5', '10', '15', '20', '25', '30', '35', '40', '45', '50'})


    xlabel('members');
    
    title(['r with indices', ', ', cfg_ts.var, ', ', num2str(sta_lonlat{stai}(1)),'E, ', num2str(sta_lonlat{stai}(3)), 'N']);
    
    dirs.figdir= [dirs.figroot, filesep, 'ens_all', filesep, cfg_ts.var, '_corr_indices'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'r_4ym_indices_', num2str(sta_lonlat{stai}(1)), 'E_', num2str(sta_lonlat{stai}(3)), 'N_', cfg_ts.var, '.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;

end
end


% % % % % 
% % % % % plot(data_AMO.time, data_AMO.obs_dseason_lp, 'color', 'k', 'linewidth', 2)
% % % % % hold on
% % % % % % for mi=1:20
% % % % %     plot(squeeze(mean(data_AMO.assm_dseason_lp(:,:),1)), 'color','r', 'linewidth', 2)
% % % % %     plot(squeeze(mean(data_AMO.hcst_dseason_lp(1,:,:),2)), 'color','b', 'linewidth', 2)
% % % % %     plot(mean(data_AMO.lens2_dseason_lp(:,:),1), 'color','g', 'linewidth', 2)
% % % % % hold off
% % % % % % end
% % % % % 
% % % % % for mi=1:20
% % % % %     hold on
% % % % %     plot(data_AMO.time,squeeze(data_AMO.hcst_dseason_lp(1,mi,:)), 'color','b', 'linewidth', 2)
% % % % % end
% % % % % 
% % % % % 
% % % % % 
% % % % % corrcoef(squeeze(mean(data_AMO.assm_dseason_lp(:,:),1)), data_AMO.obs_dseason_lp)
% % % % % 
% % % % % corrcoef(squeeze(median(data_AMO.hcst_dseason_lp(1,:,:),2)), data_AMO.obs_dseason_lp)
% % % % % corrcoef(squeeze(median(data_AMO.hcst_dseason_lp(2,:,1:end-1),2)), data_AMO.obs_dseason_lp(2:end))
% % % % % corrcoef(squeeze(median(data_AMO.hcst_dseason_lp(3,:,1:end-2),2)), data_AMO.obs_dseason_lp(3:end))
% % % % % corrcoef(squeeze(median(data_AMO.hcst_dseason_lp(4,:,1:end-3),2)), data_AMO.obs_dseason_lp(4:end))
% % % % % corrcoef(squeeze(median(data_AMO.hcst_dseason_lp(5,:,1:end-4),2)), data_AMO.obs_dseason_lp(5:end))
% % % % % 


figure;
plot(cfg.iyears, normalize( mean(data_AMO.assm_ym,1) - mean(mean(data_AMO.assm_ym,1)) ), 'b', 'linewidth', 3)
hold on
plot(cfg.iyears, normalize( mean(data_PDO.assm_ym,1) - mean(mean(data_PDO.assm_ym,1)) ), 'r', 'linewidth', 3)
plot(cfg.iyears, normalize( mean(data_ENSO.assm_ym,1) - mean(mean(data_ENSO.assm_ym,1)) ), 'g', 'linewidth', 3)
plot(cfg.iyears, normalize( mean(data_NPGO_SST.assm_ym,1) - mean(mean(data_NPGO_SST.assm_ym,1)) ), 'color', [0.9 0.9 0.9], 'linewidth', 3)
plot(cfg.iyears, normalize( mean(data_NPGO_SSH.assm_ym,1) - mean(mean(data_NPGO_SSH.assm_ym,1)) ), 'magenta', 'linewidth', 3)
plot(cfg.iyears, normalize( mean(data_SAM.assm_ym,1) - mean(mean(data_SAM.assm_ym,1)) ), 'black', 'linewidth', 3)
legend({'AMO', 'PDO', 'ENSO', 'NPGO-SST', 'NPGO-SSH', 'SAM'});
hold off

% 
% stackedplot(normalize( mean(data_AMO.assm_ym,1) - mean(mean(data_AMO.assm_ym,1)) ), ...
%     normalize( mean(data_PDO.assm_ym,1) - mean(mean(data_PDO.assm_ym,1)) ), ...
%     normalize( mean(data_ENSO.assm_ym,1) - mean(mean(data_ENSO.assm_ym,1)) ), ...
%     normalize( mean(data_NPGO_SST.assm_ym,1) - mean(mean(data_NPGO_SST.assm_ym,1)) ), ...
%     normalize( mean(data_NPGO_SSH.assm_ym,1) - mean(mean(data_NPGO_SSH.assm_ym,1)) ), ...
%     normalize( mean(data_SAM.assm_ym,1) - mean(mean(data_SAM.assm_ym,1)) ));

tbl=table(cfg.iyears', ...
    normalize( mean(data_AMO.assm_ym,1) - mean(mean(data_AMO.assm_ym,1)) )', ...
    normalize( mean(data_PDO.assm_ym,1) - mean(mean(data_PDO.assm_ym,1)) )', ...
    normalize( mean(data_ENSO.assm_ym,1) - mean(mean(data_ENSO.assm_ym,1)) )', ...
    normalize( mean(data_NPGO_SST.assm_ym,1) - mean(mean(data_NPGO_SST.assm_ym,1)) )', ...
    normalize( mean(data_NPGO_SSH.assm_ym,1) - mean(mean(data_NPGO_SSH.assm_ym,1)) )', ...
    normalize( mean(data_SAM.assm_ym,1) - mean(mean(data_SAM.assm_ym,1)) )', ...
    'VariableNames', {'Years', 'AMO', 'PDO', 'ENSO', 'NPGO-SST', 'NPGO-SSH', 'SAM'});
stackedplot(tbl, 'Xvariable', 'Years')
grid on
set(gca, 'linewidth', 3)

%% 1y (assm)
sa1=subplot(9,1,1);
shade_anomaly(cfg.iyears, normalize( mean(data_AMO.assm_ym,1) - mean(mean(data_AMO.assm_ym,1)) )', ...
    'r', 'b', 0.3, sa1);
sa1.YLabel.String='AMO';

sa2=subplot(9,1,2);
shade_anomaly(cfg.iyears, normalize( mean(data_PDO.assm_ym,1) - mean(mean(data_PDO.assm_ym,1)) )', ...
    'r', 'b', 0.3, sa2);
sa2.YLabel.String='PDO';

sa3=subplot(9,1,3);
shade_anomaly(cfg.iyears, normalize( mean(data_ENSO.assm_ym,1) - mean(mean(data_ENSO.assm_ym,1)) )', ...
    'r', 'b', 0.3, sa3);
sa3.YLabel.String='ENSO';

sa4=subplot(9,1,4);
shade_anomaly(cfg.iyears, normalize( mean(data_NPGO_SST.assm_ym,1) - mean(mean(data_NPGO_SST.assm_ym,1)) )', ...
    'r', 'b', 0.3, sa4);
sa4.YLabel.String='NPGO-SST';

sa5=subplot(9,1,5);
shade_anomaly(cfg.iyears, normalize( mean(data_NPGO_SSH.assm_ym,1) - mean(mean(data_NPGO_SSH.assm_ym,1)) )', ...
    'r', 'b', 0.3, sa5);
sa5.YLabel.String='NPGO-SSH';

sa6=subplot(9,1,6);
shade_anomaly(cfg.iyears, normalize( mean(data_SAM.assm_ym,1) - mean(mean(data_SAM.assm_ym,1)) )', ...
    'r', 'b', 0.3, sa6);
sa6.YLabel.String='SAM';

sa7=subplot(9,1,7);
shade_anomaly(cfg.iyears, normalize( mean(data_ATL3.assm_ym,1) - mean(mean(data_ATL3.assm_ym,1)) )', ...
    'r', 'b', 0.3, sa7);
sa7.YLabel.String='ATL3';
sa7.XLabel.String='Years';

sa8=subplot(9,1,8);
shade_anomaly(cfg.iyears, normalize( mean(data_TNA.assm_ym,1) - mean(mean(data_TNA.assm_ym,1)) )', ...
    'r', 'b', 0.3, sa8);
sa8.YLabel.String='TNA';
sa8.XLabel.String='Years';

sa9=subplot(9,1,9);
shade_anomaly(cfg.iyears, normalize( mean(data_TSA.assm_ym,1) - mean(mean(data_TSA.assm_ym,1)) )', ...
    'r', 'b', 0.3, sa9);
sa9.YLabel.String='TSA';
sa9.XLabel.String='Years';
set(gcf, 'Position', [0 0 1000 800]);
print(gcf, ...
    '/Users/kimyy/Desktop/backup/Research/Postdoc/03_IBS/2022_predictability_assimilation_run/paper/Figureset_raw/Fin_Figure/sub/A_climate_indices_normalized.png', ...
     '-dpng');

close all;

%% 1y (obs)
sa1=subplot(9,1,1);
shade_anomaly(cfg.iyears, normalize( data_AMO.obs_ym - mean(data_AMO.obs_ym) )', ...
    'r', 'b', 0.3, sa1);
sa1.YLabel.String='AMO';

sa2=subplot(9,1,2);
shade_anomaly(cfg.iyears, normalize( data_PDO.obs_ym - mean(data_PDO.obs_ym) )', ...
    'r', 'b', 0.3, sa2);
sa2.YLabel.String='PDO';

sa3=subplot(9,1,3);
shade_anomaly(cfg.iyears, normalize( data_ENSO.obs_ym - mean(data_ENSO.obs_ym) )', ...
    'r', 'b', 0.3, sa3);
sa3.YLabel.String='ENSO';

sa4=subplot(9,1,4);
shade_anomaly(cfg.iyears, normalize( data_NPGO_SST.obs_ym - mean(data_NPGO_SST.obs_ym) )', ...
    'r', 'b', 0.3, sa4);
sa4.YLabel.String='NPGO-SST';

sa5=subplot(9,1,5);
shade_anomaly(cfg.iyears, normalize( data_NPGO_SSH.obs_ym - mean(data_NPGO_SSH.obs_ym) )', ...
    'r', 'b', 0.3, sa5);
sa5.YLabel.String='NPGO-SSH';

sa6=subplot(9,1,6);
shade_anomaly(cfg.iyears, normalize( data_SAM.obs_ym - mean(data_SAM.obs_ym) )', ...
    'r', 'b', 0.3, sa6);
sa6.YLabel.String='SAM';

sa7=subplot(9,1,7);
shade_anomaly(cfg.iyears, normalize( data_ATL3.obs_ym - mean(data_ATL3.obs_ym) )', ...
    'r', 'b', 0.3, sa7);
sa7.YLabel.String='ATL3';
sa7.XLabel.String='Years';

sa8=subplot(9,1,8);
shade_anomaly(cfg.iyears, normalize( data_TNA.obs_ym - mean(data_TNA.obs_ym) )', ...
    'r', 'b', 0.3, sa8);
sa8.YLabel.String='TNA';
sa8.XLabel.String='Years';

sa9=subplot(9,1,9);
shade_anomaly(cfg.iyears, normalize( data_TSA.obs_ym - mean(data_TSA.obs_ym) )', ...
    'r', 'b', 0.3, sa9);
sa9.YLabel.String='TSA';
sa9.XLabel.String='Years';
set(gcf, 'Position', [0 0 1000 800]);
print(gcf, ...
    '/Users/kimyy/Desktop/backup/Research/Postdoc/03_IBS/2022_predictability_assimilation_run/paper/Figureset_raw/Fin_Figure/sub/A_climate_indices_normalized_obs.png', ...
     '-dpng');

close all;



%% 4ym (assm)
sa1=subplot(9,1,1);
shade_anomaly(cfg.iyears_4ym, normalize( mean(data_AMO.assm_4ym,1) - mean(mean(data_AMO.assm_4ym,1)) )', ...
    'r', 'b', 0.3, sa1);
sa1.YLabel.String='AMO';

sa2=subplot(9,1,2);
shade_anomaly(cfg.iyears_4ym, normalize( mean(data_PDO.assm_4ym,1) - mean(mean(data_PDO.assm_4ym,1)) )', ...
    'r', 'b', 0.3, sa2);
sa2.YLabel.String='PDO';

sa3=subplot(9,1,3);
shade_anomaly(cfg.iyears_4ym, normalize( mean(data_ENSO.assm_4ym,1) - mean(mean(data_ENSO.assm_4ym,1)) )', ...
    'r', 'b', 0.3, sa3);
sa3.YLabel.String='ENSO';

sa4=subplot(9,1,4);
shade_anomaly(cfg.iyears_4ym, normalize( mean(data_NPGO_SST.assm_4ym,1) - mean(mean(data_NPGO_SST.assm_4ym,1)) )', ...
    'r', 'b', 0.3, sa4);
sa4.YLabel.String='NPGO-SST';

sa5=subplot(9,1,5);
shade_anomaly(cfg.iyears_4ym, normalize( mean(data_NPGO_SSH.assm_4ym,1) - mean(mean(data_NPGO_SSH.assm_4ym,1)) )', ...
    'r', 'b', 0.3, sa5);
sa5.YLabel.String='NPGO-SSH';

sa6=subplot(9,1,6);
shade_anomaly(cfg.iyears_4ym, normalize( mean(data_SAM.assm_4ym,1) - mean(mean(data_SAM.assm_4ym,1)) )', ...
    'r', 'b', 0.3, sa6);
sa6.YLabel.String='SAM';

sa7=subplot(9,1,7);
shade_anomaly(cfg.iyears_4ym, normalize( mean(data_ATL3.assm_4ym,1) - mean(mean(data_ATL3.assm_4ym,1)) )', ...
    'r', 'b', 0.3, sa7);
sa7.YLabel.String='ATL3';
sa7.XLabel.String='Years';

sa8=subplot(9,1,8);
shade_anomaly(cfg.iyears_4ym, normalize( mean(data_TNA.assm_4ym,1) - mean(mean(data_TNA.assm_4ym,1)) )', ...
    'r', 'b', 0.3, sa8);
sa8.YLabel.String='TNA';
sa8.XLabel.String='Years';

sa9=subplot(9,1,9);
shade_anomaly(cfg.iyears_4ym, normalize( mean(data_TSA.assm_4ym,1) - mean(mean(data_TSA.assm_4ym,1)) )', ...
    'r', 'b', 0.3, sa9);
sa9.YLabel.String='TSA';
sa9.XLabel.String='Years';
set(gcf, 'Position', [0 0 1000 800]);
print(gcf, ...
    '/Users/kimyy/Desktop/backup/Research/Postdoc/03_IBS/2022_predictability_assimilation_run/paper/Figureset_raw/Fin_Figure/sub/A_climate_indices_normalized_4ym.png', ...
     '-dpng');

close all;


%% 4ym (obs)
sa1=subplot(9,1,1);
shade_anomaly(cfg.iyears_4ym , normalize( data_AMO.obs_4ym - mean(data_AMO.obs_4ym) )', ...
    'r', 'b', 0.3, sa1);
sa1.YLabel.String='AMO';

sa2=subplot(9,1,2);
shade_anomaly(cfg.iyears_4ym, normalize( data_PDO.obs_4ym - mean(data_PDO.obs_4ym) )', ...
    'r', 'b', 0.3, sa2);
sa2.YLabel.String='PDO';

sa3=subplot(9,1,3);
shade_anomaly(cfg.iyears_4ym, normalize( data_ENSO.obs_4ym - mean(data_ENSO.obs_4ym) )', ...
    'r', 'b', 0.3, sa3);
sa3.YLabel.String='ENSO';

sa4=subplot(9,1,4);
shade_anomaly(cfg.iyears_4ym, normalize( data_NPGO_SST.obs_4ym - mean(data_NPGO_SST.obs_4ym) )', ...
    'r', 'b', 0.3, sa4);
sa4.YLabel.String='NPGO-SST';

sa5=subplot(9,1,5);
shade_anomaly(cfg.iyears_4ym, normalize( data_NPGO_SSH.obs_4ym - mean(data_NPGO_SSH.obs_4ym) )', ...
    'r', 'b', 0.3, sa5);
sa5.YLabel.String='NPGO-SSH';

sa6=subplot(9,1,6);
shade_anomaly(cfg.iyears_4ym, normalize( data_SAM.obs_4ym - mean(data_SAM.obs_4ym) )', ...
    'r', 'b', 0.3, sa6);
sa6.YLabel.String='SAM';

sa7=subplot(9,1,7);
shade_anomaly(cfg.iyears_4ym, normalize( data_ATL3.obs_4ym - mean(data_ATL3.obs_4ym) )', ...
    'r', 'b', 0.3, sa7);
sa7.YLabel.String='ATL3';
sa7.XLabel.String='Years';

sa8=subplot(9,1,8);
shade_anomaly(cfg.iyears_4ym, normalize( data_TNA.obs_4ym - mean(data_TNA.obs_4ym) )', ...
    'r', 'b', 0.3, sa8);
sa8.YLabel.String='TNA';
sa8.XLabel.String='Years';

sa9=subplot(9,1,9);
shade_anomaly(cfg.iyears_4ym, normalize( data_TSA.obs_4ym - mean(data_TSA.obs_4ym) )', ...
    'r', 'b', 0.3, sa9);
sa9.YLabel.String='TSA';
sa9.XLabel.String='Years';
set(gcf, 'Position', [0 0 1000 800]);
print(gcf, ...
    '/Users/kimyy/Desktop/backup/Research/Postdoc/03_IBS/2022_predictability_assimilation_run/paper/Figureset_raw/Fin_Figure/sub/A_climate_indices_normalized_obs_4ym.png', ...
     '-dpng');

close all;





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
