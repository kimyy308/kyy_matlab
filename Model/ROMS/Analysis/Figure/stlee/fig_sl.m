close all; clc; clear all;

warning off;

%% set path
[error_status, tmp.hostname] = system('hostname');
tmp.hostname=tmp.hostname(1:end-1);
switch tmp.hostname
    case 'Yong-Yubs-iMac-Pro.local'
        tmp.dropboxpath = '/Volumes/kyy_raid/kimyy/Dropbox';
    case {'da1', 'da2', 'da3', 'da4'}
        tmp.dropboxpath = '/mnt/lustre/proj/kimyy/Dropbox';
    otherwise
        tmp.dropboxpath = '/mnt/lustre/proj/kimyy/Dropbox';
end
tmp.fs=filesep;
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Common', tmp.fs, 'mca']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Common', tmp.fs, 'order']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);

% cfg.varnames={'zeta','zosto','zosto_thermo','zosto_halo','sustr','svstr','ubar','vbar','Uwind','Vwind'};
cfg.varnames={'zeta','zosto','zosto_thermo','zosto_halo','ubar','vbar','Uwind','Vwind', 'wcurl'};
cfg.varnames={'Uwind', 'Vwind', 'wcurl', 'sustr', 'svstr', 'wscurl'};

dir.figdir='/Users/kimyy/Desktop/backup/Research/Postdoc/03_IBS/2024_TS_KS_connectivity/Figure';



%% data load
load('/Volumes/kyy_raid/kimyy/Model/ROMS/drifter_ROMS/mat/stlee_drifter_model_depth_1982_2019_01_12.mat');
load('/Volumes/kyy_raid/kimyy/Model/ROMS/drifter_ROMS/mat/stlee_drifter_model_steric_ssh_1982_2019_01_12.mat');

xn=size(lon,1);
yn=size(lat,2);
tn=size(comb_data.zeta,3);

% h2=h;
% mask_100=h2; 
% mask_50(mask_50>50)=NaN; mask_50=mask_50./mask_5J0;
% mask_100(mask_100>100)=NaN; mask_100=mask_100./mask_100;
% mask_200(mask_200>200)=NaN; mask_200=mask_200./mask_200;

depths=[100, 200, 300];
% depths=[150];
depths=[300];

for depthi=1:length(depths)
    dd=depths(depthi);
    str_d=num2str(dd);
    masks.(['d', str_d])=h;
    masks.(['d', str_d])(masks.(['d', str_d])>dd)=NaN;
%     masks.(['d', str_d])(masks.(['d', str_d])<50)=NaN;
    masks.(['d', str_d])=masks.(['d', str_d])./masks.(['d', str_d]);
end

[comb_data.wcurl] = Func_0037_get_curl(lon,lat,comb_data.Uwind,comb_data.Vwind);
[comb_data.wscurl] = Func_0037_get_curl(lon,lat,comb_data.sustr,comb_data.svstr);


%% annual mean
for vi=1:length(cfg.varnames)
    tmp.varn=cfg.varnames{vi};
    tmp.var=reshape(comb_data.(tmp.varn), [xn yn 12 tn/12]);
    comb_data_an.(tmp.varn)=squeeze(mean(tmp.var,3));
end




% monrange=4:10;
monrange=3:8;
monrange=6:8;

%% seasonal mean (summer)
for vi=1:length(cfg.varnames)
    tmp.varn=cfg.varnames{vi};
    tmp.var=reshape(comb_data.(tmp.varn), [xn yn 12 tn/12]);
    comb_data_an.(tmp.varn)=squeeze(mean(tmp.var(:,:,monrange,:),3));
end

part=[1800, 1750, 2050, 2030, 2000, 2400, 2200, 2010, 2500, 1730, 2150, 1930, 2700, 2030, 1510, ...
    2230, 2400, 2500, 2600, 2700, 1980, 2700, 2150, 2430, 3150, ...
    1930, 2400, 2250, 2400, 2380, 2370, 2000, 1800, 1500, 1900, 1950, 1880, 1490];


regionname='stlee_drifter';
%% polygon set
%% stleedrifter
% stlee_drifter_polygon = ...
%     [121, 24;
%     122, 24.5;
%     124, 25;
%     126, 26;
%     129, 28;
%     131, 34;
%     122, 34;
%     122, 30;
%     120, 28;
%     120, 27];

% stlee_drifter_polygon = ...
%     [127, 28;
%     129, 28;
%     131, 32;
%     127, 32];

stlee_drifter_polygon = ...
    [122, 26;
    124, 26;
    126, 26.5;
    128, 28;
    128, 30;
    130, 32;
    130, 34;
    122, 34;
    122, 32;];  %% with zeta promising result

stlee_drifter_polygon = ...
    [122, 27;
    124, 27;
    126, 28;
    127, 29;
    128, 30;
    130, 32;
    130, 34;
    122, 34;
    122, 32;]; %% with zosto promising result

stlee_drifter_polygon = ...
    [123, 27;
    124, 27;
    126, 28;
    127, 29;
    128, 30;
    128, 32;
    126, 32;
    124, 30;
    123, 28];%% with zosto dominant result


stlee_drifter_polygon = ...
    [122, 27;
    124, 27;
    125, 28;
    125, 32;
    122, 32]; %with wind dominant result

stlee_drifter_polygon = ...
    [126, 28;
    130, 30;
    130, 32;
    126, 32]; %EKB result

stlee_drifter_polygon = ...
    [120, 25;
    121, 25;
    124, 26;
    126, 26.5;
    128, 28;
    128, 30;
    130, 32;
    130, 35;
    120, 35;
    120, 32;
    122, 32;
    122, 29;
    120, 28;];  %% with onshore

switch(regionname)
    case('stlee_drifter') %% for debugging
        refpolygon=stlee_drifter_polygon;
    otherwise
        ('?');
end

%% masking
switch(regionname)
    case('NWP') %% North western Pacific
        mask_model(1:size(lon,1),1:size(lon,2))=1;
    otherwise
        mask_model = double(inpolygon(lon,lat,refpolygon(:,1),refpolygon(:,2)));
        mask_model(mask_model==0)=NaN;
end

stlee_orig_polygon = ...
    [121, 24;
    122, 24.5;
    124, 25;
    126, 26;
    129, 28;
    131, 32;
    122, 32;
    122, 30;
    120, 28;
    120, 27];

lonlat(1)=min(stlee_orig_polygon(:,1));
lonlat(2)=max(stlee_orig_polygon(:,1));
lonlat(3)=min(stlee_orig_polygon(:,2));
lonlat(4)=max(stlee_orig_polygon(:,2));

% lonlat(1)=min(refpolygon(:,1));
% lonlat(2)=max(refpolygon(:,1));
% lonlat(3)=min(refpolygon(:,2));
% lonlat(4)=max(refpolygon(:,2));

lonlat(1)=min(min(lon));
lonlat(2)=max(max(lon));
lonlat(3)=min(min(lat));
lonlat(4)=max(max(lat));



% % std_data.zosto=std(comb_data.zosto,1,3);
% % pcolor(std_data.zosto(:,:)'); shading flat; colorbar;

cfg.years=1982:2019;
[fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_20', tmp.dropboxpath);


% % % % % tmp.svd_modes=4;
% % % % % for vi=1:length(cfg.varnames)
% % % % %     tmp.varn=cfg.varnames{vi};
% % % % %     [EOF_result.(tmp.varn).lv, ...
% % % % %         EOF_result.(tmp.varn).pcs, ...
% % % % %         EOF_result.(tmp.varn).var_exp] = ...
% % % % %         Func_0024_EOF_3d( comb_data_an.(tmp.varn).*mask_model, tmp.svd_modes);
% % % % % end
% % % % % 
% % % % % % tmp.varn='Vwind'
% % % % % for vi=1:length(cfg.varnames)
% % % % %     tmp.varn=cfg.varnames{vi};
% % % % %     
% % % % %     close all;
% % % % %     fig_cfg.fig_size=[0 0 14 16];
% % % % %     fig_h = figure('name', 'EOF','PaperUnits','inches', ...
% % % % %             'PaperPosition',(fig_cfg.fig_size),'position', fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','on');
% % % % %     amod=1;
% % % % %     subplot(4,2,1)
% % % % %     lv_l=EOF_result.(tmp.varn).lv(:,:,amod);
% % % % %     m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
% % % % %     m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
% % % % %     % m_grid('fontsize', m_grid_fontsize-5, 'tickdir', m_grid_tickdir_type, 'xticklabels', [], 'xtick',(120:20:160), 'parent', ax1_1);  
% % % % %     m_grid;  
% % % % %     m_gshhs_i('color',[1 1 1])  
% % % % %     m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
% % % % %     title([tmp.varn, ' lv, ', num2str(round(EOF_result.(tmp.varn).var_exp(amod),2)), '%']);
% % % % %     mcax=max(abs(lv_l(:)));
% % % % %     caxis([-mcax mcax])
% % % % %     colormap(fig_cfg.c_map)
% % % % %     
% % % % %     subplot(4,2,2);
% % % % %     plot(cfg.years,EOF_result.(tmp.varn).pcs(:,amod))
% % % % %     tmp.corr=corrcoef(part, EOF_result.(tmp.varn).pcs(:,amod));
% % % % %     title([tmp.varn, ' pcs, ', num2str(round(EOF_result.(tmp.varn).var_exp(amod),2)), '%, ', num2str(round(tmp.corr(1,2),2))]);
% % % % %     grid on
% % % % %     
% % % % %     
% % % % %     amod=2;
% % % % %     subplot(4,2,3)
% % % % %     lv_l=EOF_result.(tmp.varn).lv(:,:,amod);
% % % % %     m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
% % % % %     m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
% % % % %     m_grid;  
% % % % %     m_gshhs_i('color',[1 1 1])  
% % % % %     m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
% % % % %     title([tmp.varn, ' lv, ', num2str(round(EOF_result.(tmp.varn).var_exp(amod),2)), '%']);
% % % % %     mcax=max(abs(lv_l(:)));
% % % % %     caxis([-mcax mcax])
% % % % %     colormap(fig_cfg.c_map)
% % % % %     
% % % % %     subplot(4,2,4);
% % % % %     plot(cfg.years,EOF_result.(tmp.varn).pcs(:,amod))
% % % % %     tmp.corr=corrcoef(part, EOF_result.(tmp.varn).pcs(:,amod));
% % % % %     title([tmp.varn, ' pcs, ', num2str(round(EOF_result.(tmp.varn).var_exp(amod),2)), '%, ', num2str(round(tmp.corr(1,2),2))]);
% % % % %     grid on
% % % % %     
% % % % %     amod=3;
% % % % %     subplot(4,2,5)
% % % % %     lv_l=EOF_result.(tmp.varn).lv(:,:,amod);
% % % % %     m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
% % % % %     m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
% % % % %     m_grid;  
% % % % %     m_gshhs_i('color',[1 1 1])  
% % % % %     m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
% % % % %     title([tmp.varn, ' lv, ', num2str(round(EOF_result.(tmp.varn).var_exp(amod),2)), '%']);
% % % % %     mcax=max(abs(lv_l(:)));
% % % % %     caxis([-mcax mcax])
% % % % %     colormap(fig_cfg.c_map)
% % % % %     
% % % % %     subplot(4,2,6);
% % % % %     plot(cfg.years,EOF_result.(tmp.varn).pcs(:,amod))
% % % % %     tmp.corr=corrcoef(part, EOF_result.(tmp.varn).pcs(:,amod));
% % % % %     title([tmp.varn, ' pcs, ', num2str(round(EOF_result.(tmp.varn).var_exp(amod),2)), '%, ', num2str(round(tmp.corr(1,2),2))]);
% % % % %     grid on
% % % % %     
% % % % %     amod=4;
% % % % %     subplot(4,2,7)
% % % % %     lv_l=EOF_result.(tmp.varn).lv(:,:,amod);
% % % % %     m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
% % % % %     m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
% % % % %     m_grid;  
% % % % %     m_gshhs_i('color',[1 1 1])  
% % % % %     m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
% % % % %     title([tmp.varn, ' lv, ', num2str(round(EOF_result.(tmp.varn).var_exp(amod),2)), '%']);
% % % % %     mcax=max(abs(lv_l(:)));
% % % % %     caxis([-mcax mcax])
% % % % %     colormap(fig_cfg.c_map)
% % % % %     
% % % % %     subplot(4,2,8);
% % % % %     plot(cfg.years,EOF_result.(tmp.varn).pcs(:,amod))
% % % % %     tmp.corr=corrcoef(part, EOF_result.(tmp.varn).pcs(:,amod));
% % % % %     title([tmp.varn, ' pcs, ', num2str(round(EOF_result.(tmp.varn).var_exp(amod),2)), '%, ', num2str(round(tmp.corr(1,2),2))]);
% % % % %     grid on
% % % % %     
% % % % %     
% % % % %     dir.figtgdir=[dir.figdir, '/', 'EOF'];
% % % % %     mkdir(dir.figtgdir);
% % % % %     
% % % % %     cfg.figname=[dir.figtgdir, '/', 'EOF_',tmp.varn,'_',num2str(min(monrange)),'_',num2str(max(monrange)) '.tif'];
% % % % %     print(fig_h, cfg.figname, '-dpng');
% % % % % end
% % % % % 
% % % % % 
% % % % % %% timeseries
% % % % % for vi=1:length(cfg.varnames)
% % % % %     close all;
% % % % %     fig_cfg.fig_size=[0 0 10 6];
% % % % %     fig_h = figure('name', 'EOF','PaperUnits','inches', ...
% % % % %             'PaperPosition',(fig_cfg.fig_size),'position', fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','on');
% % % % % 
% % % % %     tmp.varn=cfg.varnames{vi};
% % % % %     loninfo=129;
% % % % %     latinfo=30;
% % % % %     [indw, inde, inds, indn]=Func_0012_findind_Y(1/10,[loninfo latinfo],lon,lat, 1);
% % % % %     plot(cfg.years,squeeze(comb_data_an.(tmp.varn)(indw,inds,:)))
% % % % %     tmp.corr=corrcoef(part, squeeze(comb_data_an.(tmp.varn)(indw,inds,:)));
% % % % %     title([tmp.varn, ' ts, lon:', num2str(loninfo), ' lat: ', num2str(latinfo), ' ', num2str(round(tmp.corr(1,2),2))]);
% % % % %     grid on;
% % % % %     
% % % % %     dir.figtgdir=[dir.figdir, '/', 'timeseries'];
% % % % %         mkdir(dir.figtgdir);
% % % % %     cfg.figname=[dir.figtgdir, '/', 'ts_',tmp.varn,'_', ...
% % % % %         num2str(loninfo),'E_', num2str(latinfo), 'N_', ...
% % % % %         num2str(min(monrange)),'_',num2str(max(monrange)) '.tif'];
% % % % %         print(fig_h, cfg.figname, '-dpng');
% % % % % end


%% depth masked
% tmp.varn='Vwind'

for depthi=1:length(depths)
    dd=depths(depthi);
    str_d=num2str(dd);
    masks.(['d', str_d]);

    tmp.svd_modes=4;
    for vi=1:length(cfg.varnames)
        tmp.varn=cfg.varnames{vi};
        
%         [EOF_result.(tmp.varn).lv, ...
%             EOF_result.(tmp.varn).pcs, ...
%             EOF_result.(tmp.varn).var_exp] = ...
%             Func_0024_EOF_3d( comb_data_an.(tmp.varn).*mask_model.*masks.(['d', str_d]), tmp.svd_modes);
% 

%             [EOF_result.(tmp.varn).lv, ...
%                 EOF_result.(tmp.varn).pcs, ...            
%                 EOF_result.(tmp.varn).var_exp] = ...
%             Func_0024_EOF_3d( comb_data_an.(tmp.varn).*masks.(['d', str_d]), tmp.svd_modes);

          [EOF_result.(tmp.varn).lv, ...
        EOF_result.(tmp.varn).pcs, ...
        EOF_result.(tmp.varn).var_exp] = ...
        Func_0024_EOF_3d( comb_data_an.(tmp.varn).*mask_model.*masks.(['d', str_d]), tmp.svd_modes, lat);


    end

%% sign check


for vi=1:length(cfg.varnames)
    tmp.varn=cfg.varnames{vi};
    for amodi=1:tmp.svd_modes
        if mean(EOF_result.(tmp.varn).pcs(:,amodi))<0
            EOF_result.(tmp.varn).lv(:,:,amodi)=-EOF_result.(tmp.varn).lv(:,:,amodi);
            EOF_result.(tmp.varn).pcs(:,amodi)=-EOF_result.(tmp.varn).pcs(:,amodi);
        end
    end
    close all;
    fig_cfg.fig_size=[0 0 14 16];
    fig_h = figure('name', 'EOF','PaperUnits','inches', ...
            'PaperPosition',(fig_cfg.fig_size),'position', fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','on');
    amod=1;
    subplot(4,2,1)
    lv_l=EOF_result.(tmp.varn).lv(:,:,amod);
    m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
    m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
    % m_grid('fontsize', m_grid_fontsize-5, 'tickdir', m_grid_tickdir_type, 'xticklabels', [], 'xtick',(120:20:160), 'parent', ax1_1);  
    m_grid;  
    m_gshhs_i('color',[1 1 1])  
    m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
    title([tmp.varn, ' lv, ', num2str(round(EOF_result.(tmp.varn).var_exp(amod),2)), '%']);
    mcax=max(abs(lv_l(:)));
    caxis([-mcax mcax])
    colormap(fig_cfg.c_map)
    
    subplot(4,2,2);
    plot(cfg.years,EOF_result.(tmp.varn).pcs(:,amod))
    tmp.corr=corrcoef(part, EOF_result.(tmp.varn).pcs(:,amod));
    title([tmp.varn, ' pcs, ', num2str(round(EOF_result.(tmp.varn).var_exp(amod),2)), '%, ', num2str(round(tmp.corr(1,2),2))]);
    grid on
    
    
    amod=2;
    subplot(4,2,3)
    lv_l=EOF_result.(tmp.varn).lv(:,:,amod);
    m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
    m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
    m_grid;  
    m_gshhs_i('color',[1 1 1])  
    m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
    title([tmp.varn, ' lv, ', num2str(round(EOF_result.(tmp.varn).var_exp(amod),2)), '%']);
    mcax=max(abs(lv_l(:)));
    caxis([-mcax mcax])
    colormap(fig_cfg.c_map)
    
    subplot(4,2,4);
    plot(cfg.years,EOF_result.(tmp.varn).pcs(:,amod))
    tmp.corr=corrcoef(part, EOF_result.(tmp.varn).pcs(:,amod));
    title([tmp.varn, ' pcs, ', num2str(round(EOF_result.(tmp.varn).var_exp(amod),2)), '%, ', num2str(round(tmp.corr(1,2),2))]);
    grid on
    
    amod=3;
    subplot(4,2,5)
    lv_l=EOF_result.(tmp.varn).lv(:,:,amod);
    m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
    m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
    m_grid;  
    m_gshhs_i('color',[1 1 1])  
    m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
    title([tmp.varn, ' lv, ', num2str(round(EOF_result.(tmp.varn).var_exp(amod),2)), '%']);
    mcax=max(abs(lv_l(:)));
    caxis([-mcax mcax])
    colormap(fig_cfg.c_map)
    
    subplot(4,2,6);
    plot(cfg.years,EOF_result.(tmp.varn).pcs(:,amod))
    tmp.corr=corrcoef(part, EOF_result.(tmp.varn).pcs(:,amod));
    title([tmp.varn, ' pcs, ', num2str(round(EOF_result.(tmp.varn).var_exp(amod),2)), '%, ', num2str(round(tmp.corr(1,2),2))]);
    grid on
    
    amod=4;
    subplot(4,2,7)
    lv_l=EOF_result.(tmp.varn).lv(:,:,amod);
    m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
    m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
    m_grid;  
    m_gshhs_i('color',[1 1 1])  
    m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
    title([tmp.varn, ' lv, ', num2str(round(EOF_result.(tmp.varn).var_exp(amod),2)), '%']);
    mcax=max(abs(lv_l(:)));
    caxis([-mcax mcax])
    colormap(fig_cfg.c_map)
    
    subplot(4,2,8);
    plot(cfg.years,EOF_result.(tmp.varn).pcs(:,amod))
    tmp.corr=corrcoef(part, EOF_result.(tmp.varn).pcs(:,amod));
    title([tmp.varn, ' pcs, ', num2str(round(EOF_result.(tmp.varn).var_exp(amod),2)), '%, ', num2str(round(tmp.corr(1,2),2))]);
    grid on
    
    
    dir.figtgdir=[dir.figdir, '/', 'EOF'];
    mkdir(dir.figtgdir);
    
    cfg.figname=[dir.figtgdir, '/', 'EOF_',tmp.varn,'_',str_d, '_', num2str(min(monrange)),'_',num2str(max(monrange)) '.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
end


%% timeseries
for vi=1:length(cfg.varnames)
    close all;
    fig_cfg.fig_size=[0 0 10 6];
    fig_h = figure('name', 'EOF','PaperUnits','inches', ...
            'PaperPosition',(fig_cfg.fig_size),'position', fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','on');

    tmp.varn=cfg.varnames{vi};
    loninfo=129;
    latinfo=30;
    [indw, inde, inds, indn]=Func_0012_findind_Y(1/10,[loninfo latinfo],lon,lat, 1);
    plot(cfg.years,squeeze(comb_data_an.(tmp.varn)(indw,inds,:)))
    tmp.corr=corrcoef(part, squeeze(comb_data_an.(tmp.varn)(indw,inds,:)));
    title([tmp.varn, ' ts, lon:', num2str(loninfo), ' lat: ', num2str(latinfo), ' ', num2str(round(tmp.corr(1,2),2))]);
    grid on;
    
    dir.figtgdir=[dir.figdir, '/', 'timeseries'];
        mkdir(dir.figtgdir);
    cfg.figname=[dir.figtgdir, '/', 'ts_',tmp.varn,'_', ...
        str_d, '_', num2str(loninfo),'E_', num2str(latinfo), 'N_', ...
        num2str(min(monrange)),'_',num2str(max(monrange)) '.tif'];
        print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
end


end




%% 95% significance low & high range
    conf = 0.90;
%     n=cfg.t_ini_max-cfg.t_ini_min+1;
    n=38;
    alpha = 1 - conf;
    pLo = alpha/2;
    pUp = 1 - alpha/2;
    crit = tinv([pLo pUp], n-1);
    % xbar = mean(r); % = 0
    xbar = 0; % = 0
    r_crit=sqrt((crit.^2)./(n-2+(crit).^2));













%% fig_zosto
fig_cfg.fig_size=[0 0 10 10];
fig_h = figure('name', 'EOF','PaperUnits','inches', ...
        'PaperPosition',(fig_cfg.fig_size),'position', fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','on');

for loni=1:size(comb_data_an.zeta,1)
    for lati=1:size(comb_data_an.zeta,2)
        tmp.corr=corrcoef(comb_data_an.zeta(loni,lati,:), comb_data_an.zosto(loni,lati,:));
        corr_zeta(loni,lati)=tmp.corr(1,2);
    end
end
m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
m_pcolor(lon',lat',corr_zeta'); shading flat; colorbar; caxis([-1 1]); colormap(fig_cfg.c_map);
m_grid;  
m_gshhs_i('color',[1 1 1])  
m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
dir.figtgdir=[dir.figdir, '/', 'check_steric'];
        mkdir(dir.figtgdir);
cfg.figname=[dir.figtgdir, '/', 'corr_zosto',  '.tif'];
    print(fig_h, cfg.figname, '-dpng');

%% fig_zosto_thermo
fig_cfg.fig_size=[0 0 10 10];
fig_h = figure('name', 'EOF','PaperUnits','inches', ...
        'PaperPosition',(fig_cfg.fig_size),'position', fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','on');

for loni=1:size(comb_data_an.zeta,1)
    for lati=1:size(comb_data_an.zeta,2)
        tmp.corr=corrcoef(comb_data_an.zeta(loni,lati,:), comb_data_an.zosto_thermo(loni,lati,:));
        corr_zeta(loni,lati)=tmp.corr(1,2);
    end
end
m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
m_pcolor(lon',lat',corr_zeta'); shading flat; colorbar; caxis([-1 1]); colormap(fig_cfg.c_map);
m_grid;  
m_gshhs_i('color',[1 1 1])  
m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
dir.figtgdir=[dir.figdir, '/', 'check_steric'];
        mkdir(dir.figtgdir);
cfg.figname=[dir.figtgdir, '/', 'corr_zosto_thermo',  '.tif'];
    print(fig_h, cfg.figname, '-dpng');
close all;

%% fig_zosto_halo
fig_cfg.fig_size=[0 0 10 10];
fig_h = figure('name', 'EOF','PaperUnits','inches', ...
        'PaperPosition',(fig_cfg.fig_size),'position', fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','on');

for loni=1:size(comb_data_an.zeta,1)
    for lati=1:size(comb_data_an.zeta,2)
        tmp.corr=corrcoef(comb_data_an.zeta(loni,lati,:), comb_data_an.zosto_halo(loni,lati,:));
        corr_zeta(loni,lati)=tmp.corr(1,2);
    end
end
m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
m_pcolor(lon',lat',corr_zeta'); shading flat; colorbar; caxis([-1 1]); colormap(fig_cfg.c_map);
m_grid;  
m_gshhs_i('color',[1 1 1])  
m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
dir.figtgdir=[dir.figdir, '/', 'check_steric'];
        mkdir(dir.figtgdir);
cfg.figname=[dir.figtgdir, '/', 'corr_zosto_halo',  '.tif'];
    print(fig_h, cfg.figname, '-dpng');
close all;





% % [lv, pc, var_exp] = Func_0024_EOF_3d(data,X);

var_l='zosto_thermo';
% var_l='Vwind';
var_r='zeta';
amod=4;
[mca_result.lv_left, ...
        mca_result.pcs_left, ...
        mca_result.lv_right, ...
        mca_result.pcs_right, ...
        mca_result.lambda, ...
        mca_result.scf] = ...
        mca( double(comb_data_an.(var_l).*mask_model), double(comb_data_an.(var_r).*mask_model), tmp.svd_modes);


close all;
subplot(2,2,1)
mca_result.scf(amod)
lv_l=mca_result.lv_left(:,:,amod);
m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
% m_grid('fontsize', m_grid_fontsize-5, 'tickdir', m_grid_tickdir_type, 'xticklabels', [], 'xtick',(120:20:160), 'parent', ax1_1);  
m_grid;  
m_gshhs_i('color',[1 1 1])  
m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
title([var_l, ' lv, ', num2str(round(mca_result.scf(amod).*100,2)), '%']);
mcax=max(abs(lv_l(:)));
caxis([-mcax mcax])
colormap(fig_cfg.c_map)

subplot(2,2,2);
plot(cfg.years,mca_result.pcs_left(amod,:))
title([var_l, ' pcs, ', num2str(round(mca_result.scf(amod).*100,2)), '%']);
grid on


subplot(2,2,3);
lv_r=mca_result.lv_right(:,:,amod);
m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
m_pcolor(lon',lat',lv_r'); shading flat; colorbar;
% m_grid('fontsize', m_grid_fontsize-5, 'tickdir', m_grid_tickdir_type, 'xticklabels', [], 'xtick',(120:20:160), 'parent', ax1_1);  
m_grid;  
m_gshhs_i('color',[1 1 1]);
m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
title([var_r, ' lv, ', num2str(round(mca_result.scf(amod).*100,2)), '%']);
mcax=max(abs(lv_r(:)));
caxis([-mcax mcax])
colormap(fig_cfg.c_map)

subplot(2,2,4);
plot(cfg.years,mca_result.pcs_right(amod,:))
title([var_r, ' pcs, ', num2str(round(mca_result.scf(amod).*100,2)), '%']);
grid on



% pcolor(lon',lat',mca_result.lv_right(:,:,1)'); shading flat; colorbar;
% pcolor(mask_model'); shading flat; colorbar;