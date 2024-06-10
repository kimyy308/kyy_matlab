% %  Created 07-Jul-2023 by Yong-Yub Kim
clc; clear all; close all;
warning off;

%% set path
[error_status, tmp.hostname] = system('hostname');
tmp.hostname=tmp.hostname(1:end-1);
switch tmp.hostname
    case 'Yong-Yubs-iMac-Pro.local'
        tmp.dropboxpath = '/Volumes/kyy_raid/kimyy/Dropbox';
        tmp.kimyypath = '/Volumes/kyy_raid/kimyy';
    case 'Yong-Yubui-MacBookPro.local'
        tmp.dropboxpath = '/Users/kimyy/Dropbox';
        tmp.kimyypath = '/Users/kimyy';
    case {'da1', 'da2', 'da3', 'da4'}
        tmp.dropboxpath = '/mnt/lustre/proj/kimyy/Dropbox';
end
tmp.fs=filesep;
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);

%% model configuration
% cfg.var='TS'; %SST PRECT PSL TS SSH sumChl
cfg.vars={'photoC_TOT_zint_100m'};
cfg.vars={'zooC'};

% cfg.vlayer=1:10; % 10layer. don't put more than 15
cfg.vlayer=1; % surface, vertical slice (UVEL55, TEMP145, ...)

cfg.vlayer_1st=min(cfg.vlayer);
cfg.vlayer_cnt=max(cfg.vlayer)-cfg.vlayer_1st+1;


for vari=1:length(cfg.vars)

cfg.var=cfg.vars{vari};
cfg.obs_name=f_obs_name(cfg.var);
cfg.obs_fname_mid=f_obs_name_mid(cfg.var);
cfg.obs_varname=f_obs_varname(cfg.var);
cfg.comp=Func_0025_CESM2_cmpname_var(cfg.var);
cfg.obs_fname_module=f_obs_fname_module(cfg.comp);
cfg.obs_iyears=1965:2020;

disp(cfg.var);
tic;

dirs.hcstmatroot=[tmp.kimyypath, '/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', cfg.var];
vstr=['v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i')];
dirs.figroot=[tmp.kimyypath,'/Figure/CESM2/ESP/HCST_EXP/archive_anomaly/', cfg.comp,'/', cfg.var, filesep, vstr];



cfg.iyears=cfg.obs_iyears;
cfg.gnm='f09_g17';
cfg.proj_year=5;

cfg.len_t_y = length(cfg.iyears);
cfg.casename_m = ['ens_all'];


% [tmp.error_status, tmp.value]=system(['ls ', dirs.datadir, '/*once*']);  % b.e21.BHISTsmbb.f09_g17.assm.oras4_ba-10p1.pop.h.once.nc
tmp.gridname = [tmp.kimyypath, '/Model/CESM2/ESP/HCST_EXP/archive_transfer/', cfg.comp, '/grid.nc'];
tmp.maskname = [tmp.kimyypath, '/Model/CESM2/ESP/HCST_EXP/archive_transfer/ocn/RECCAP2_region_masks_all_v20210412_POP2_grid.nc'];

switch cfg.comp
    case {'ocn', 'ice'}
        grid.tlong=ncread(tmp.gridname, 'TLONG');
        grid.tlat=ncread(tmp.gridname, 'TLAT');
        grid.mask_ocn=ncread(tmp.maskname, 'open_ocean');
        grid.mask_ocn(grid.mask_ocn<-10e10)=NaN;
        grid.mask_ocn=grid.mask_ocn./grid.mask_ocn;
        grid.tarea=ncread(tmp.gridname, 'TAREA')/1000000.0; %(m^2 -> km^2)
        grid.tarea_60=grid.tarea; grid.tarea_60(grid.tlat>60 | grid.tlat<-60)=NaN;
    case {'atm', 'lnd'}
        grid.lon=ncread(tmp.gridname, 'lon');
        grid.lat=ncread(tmp.gridname, 'lat');
        [grid.tlat grid.tlong]=meshgrid(grid.lat, grid.lon);
        grid.tarea=ncread(tmp.gridname, 'AREA');
        grid.tarea_60=grid.tarea; grid.tarea_60(grid.tlat>60 | grid.tlat<-60)=NaN;
end

grid.nlon=size(grid.tlong,1);
grid.nlat=size(grid.tlat,2);
% grid.ntime=cfg.proj_year.*12;


% % model filename example
% /mnt/lustre/proj/kimyy/Model/CESM2/ESP/HCST_EXP/archive/atm/SST/ens_all/ens_all_i2020
% SST_f09_g17.hcst.ens_all_i2020.cam.h0.2024-12.nc

S = shaperead('landareas.shp');

%% read & plot data
tmp.varname=cfg.var;

clear tmp.ydata tmp.ydata_lens2 tmp.ydata_obs tmp.ydata_assm
for lyear=0:cfg.proj_year-1
    tmp.lyear_str=num2str(lyear, '%02i');
    fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
            '_l', tmp.lyear_str, 'y.mat'];
    load(fig_cfg.mat_name, 'data', 'data2')
    data2_all.(['l',tmp.lyear_str])=data2;
end


% data2_all.l00.TEMPCLINE_corr_assm_int_l00(180,240)

pred_hor=NaN(size(data2_all.l00.([tmp.varname,'_corr_assm_int_ano_l00'])));
for lyear=0:cfg.proj_year-1
    tmp.lyear_str=num2str(lyear, '%02i');
    for loni=1:size(data2_all.l00.([tmp.varname,'_corr_assm_int_ano_l00']), 1)
        for lati=1:size(data2_all.l00.([tmp.varname,'_corr_assm_int_ano_l00']), 2)
            if (data2_all.(['l',tmp.lyear_str]).([tmp.varname,'_corr_obs_ano_l',tmp.lyear_str])(loni,lati)>=0 && ...
                    data2_all.(['l',tmp.lyear_str]).([tmp.varname,'_corr_obs_ano_p_l',tmp.lyear_str])(loni,lati) <=0.1)
                pred_hor(loni,lati)=lyear+1;
            end
        end
    end
end

pred_hor_int=NaN(size(data2_all.l00.([tmp.varname,'_corr_assm_int_ano_l00'])));
for lyear=0:cfg.proj_year-1
    tmp.lyear_str=num2str(lyear, '%02i');
    for loni=1:size(data2_all.l00.([tmp.varname,'_corr_assm_int_ano_l00']), 1)
        for lati=1:size(data2_all.l00.([tmp.varname,'_corr_assm_int_ano_l00']), 2)
            if (data2_all.(['l',tmp.lyear_str]).([tmp.varname,'_corr_obs_int_ano_l',tmp.lyear_str])(loni,lati)>=0 && ...
                    data2_all.(['l',tmp.lyear_str]).([tmp.varname,'_corr_obs_int_ano_p_l',tmp.lyear_str])(loni,lati) <=0.1)
                pred_hor_int(loni,lati)=lyear+1;
            end
        end
    end
end


% tmp.data=data2_all.(['l00']).([tmp.varname,'_corr_obs_ano_p_l00']);
% tmp.data(tmp.data>0.1)=NaN;
% tmp.data(data2_all.(['l00']).([tmp.varname,'_corr_obs_ano_l00'])<0)=NaN;
% pcolor(tmp.data'); shading flat; colorbar;


fig_cfg.name_rgn = 'Glob';
fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
fig_cfg.x_lim = [-180 180];
fig_cfg.y_lim = [-80 89];
fig_cfg.fig_size = [0,0,6,3.5];
fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
fig_cfg.title_pos = [0.5,0.93];
fig_cfg.p_lim =0.1;
fig_cfg.c_lim = [0.5 5.5];
[fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('wr_5', tmp.dropboxpath);

tmp.X=grid.tlong([end, 1:end],:);
tmp.Y=grid.tlat([end, 1:end],:);
tmp.C=pred_hor;
tmp.C=tmp.C([end, 1:end],:);

fig_cfg.fig_name=['pred_hor, ', tmp. varname];
fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
    'PaperPosition',fig_cfg.fig_size, ...
    'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0], ...
    'visible','off');
%% map setting
ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
    'fontname','freeserif'); 

axis off; 
hold on;
setm(ax_m,'origin',[0,180],'MapLatLimit',fig_cfg.y_lim);
set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
'fontsize',14,'fontname','freeserif','interpreter','none')

%% caxis & colorbar
caxis(ax_m, fig_cfg.c_lim); 
colormap(fig_cfg.c_map);
cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
set(cb, 'Ticks', 1:5);
title(cb,'Y','fontsize',12);

%% draw on ax_m
h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
shading flat;
geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);

%% frame and label setting
setm(ax_m,'frame','on','FLineWidth',1);

label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
mlabel; plabel;
label_y=plabel; label_x=mlabel;
for lxi=1:length(label_x)
    tmp.tmppos=label_x(lxi,1).Position;
    tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55;
    label_x(lxi,1).Position=tmp.tmppos;
    label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
end
for lyi=1:length(label_y)
    label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
end

%% save
dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_predictability_horizon_map', filesep, 'model'];
if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
cfg.figname=[dirs.figdir, filesep, 'pred_hor_obs_ano_map_', tmp.varname, '.tif'];

print(fig_h, cfg.figname, '-dpng');
RemoveWhiteSpace([], 'file', cfg.figname);
close all;


%% obs_int

fig_cfg.name_rgn = 'Glob';
fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
fig_cfg.x_lim = [-180 180];
fig_cfg.y_lim = [-80 89];
fig_cfg.fig_size = [0,0,6,3.5];
fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
fig_cfg.title_pos = [0.5,0.93];
fig_cfg.p_lim =0.1;
fig_cfg.c_lim = [0.5 5.5];
[fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('wr_5', tmp.dropboxpath);

tmp.X=grid.tlong([end, 1:end],:);
tmp.Y=grid.tlat([end, 1:end],:);
tmp.C=pred_hor_int;
tmp.C=tmp.C([end, 1:end],:);

fig_cfg.fig_name=['pred_hor, ', tmp. varname];
fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
    'PaperPosition',fig_cfg.fig_size, ...
    'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0], ...
    'visible','off');
%% map setting
ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
    'fontname','freeserif'); 

axis off; 
hold on;
setm(ax_m,'origin',[0,180],'MapLatLimit',fig_cfg.y_lim);
set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
'fontsize',14,'fontname','freeserif','interpreter','none')

%% caxis & colorbar
caxis(ax_m, fig_cfg.c_lim); 
colormap(fig_cfg.c_map);
cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
set(cb, 'Ticks', 1:5);
title(cb,'Y','fontsize',12);

%% draw on ax_m
h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
shading flat;
geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);

%% frame and label setting
setm(ax_m,'frame','on','FLineWidth',1);

label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
mlabel; plabel;
label_y=plabel; label_x=mlabel;
for lxi=1:length(label_x)
    tmp.tmppos=label_x(lxi,1).Position;
    tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55;
    label_x(lxi,1).Position=tmp.tmppos;
    label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
end
for lyi=1:length(label_y)
    label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
end

%% save
dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_predictability_horizon_map', filesep, 'model'];
if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
cfg.figname=[dirs.figdir, filesep, 'pred_hor_obs_int_ano_map_', tmp.varname, '.tif'];

print(fig_h, cfg.figname, '-dpng');
RemoveWhiteSpace([], 'file', cfg.figname);
close all;








tmp.rootpath='/Volumes/kyy_raid';
matname2=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/ocn/photoC_TOT_zint_100m/','all', '_NPP_l', '00', '.mat'];
tmp.phyto='sp';
tmp.lyear_str='00';
matname=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/ocn/photoC_TOT_zint_100m/',tmp.phyto, '_nut_NPP_l', tmp.lyear_str, '.mat'];
load(matname, 'data_mask_nut_dominant', 'data_min_freq')
% ind(1~4): Fe, N, SiO3, P
mask.nlim.(tmp.phyto)=data_min_freq.N_lim_phyto.(tmp.phyto).assm;
mask.nlim.(tmp.phyto)(mask.nlim.(tmp.phyto)<0.5)=NaN;
% pcolor(mask.nlim.(tmp.phyto)'); shading flat; colorbar;

tmp.phyto='diat';
tmp.lyear_str='00';
matname=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/ocn/photoC_TOT_zint_100m/',tmp.phyto, '_nut_NPP_l', tmp.lyear_str, '.mat'];
load(matname, 'data_mask_nut_dominant', 'data_min_freq')
% ind(1~4): Fe, N, SiO3, P
mask.nlim.(tmp.phyto)=data_min_freq.N_lim_phyto.(tmp.phyto).assm;
mask.nlim.(tmp.phyto)(mask.nlim.(tmp.phyto)<0.5)=NaN;
% pcolor(mask.nlim.(tmp.phyto)'); shading flat; colorbar;

matname=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/HCST_EXP/mat/ocn/photoC_diat_zint_100m/', 'hcst_corr_assm_photoC_diat_zint_100m_v01_v10_l00y.mat'];
load(matname);
tmp.phyto='diat';
ensmean_npp.(tmp.phyto)=mean(data.photoC_diat_zint_100m_assm,3);

tmp.phyto='sp';
matname=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/HCST_EXP/mat/ocn/photoC_sp_zint_100m/', 'hcst_corr_assm_photoC_sp_zint_100m_v01_v10_l00y.mat'];
load(matname);
ensmean_npp.(tmp.phyto)=mean(data.photoC_sp_zint_100m_assm,3);

for loni=1:size(ensmean_npp.(tmp.phyto),1)
    for lati=1:size(ensmean_npp.(tmp.phyto),2)
        if ensmean_npp.diat(loni,lati) > ensmean_npp.sp(loni,lati)
            mask.nlim.comb(loni,lati)=mask.nlim.diat(loni,lati);
        else
            mask.nlim.comb(loni,lati)=mask.nlim.sp(loni,lati);
        end
    end
end




%% n-limitation filtered corr
fig_cfg.name_rgn = 'Glob';
fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
fig_cfg.x_lim = [-180 180];
fig_cfg.y_lim = [-80 89];
fig_cfg.fig_size = [0,0,6,3.5];
fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
fig_cfg.title_pos = [0.5,0.93];
fig_cfg.p_lim =0.1;
fig_cfg.c_lim = [0.5 5.5];
[fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('wr_5', tmp.dropboxpath);

tmp.X=grid.tlong([end, 1:end],:);
tmp.Y=grid.tlat([end, 1:end],:);
tmp.C=pred_hor.*mask.nlim.comb;
tmp.C=tmp.C([end, 1:end],:);

fig_cfg.fig_name=['pred_hor, ', tmp. varname];
fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
    'PaperPosition',fig_cfg.fig_size, ...
    'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0], ...
    'visible','off');
%% map setting
ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
    'fontname','freeserif'); 

axis off; 
hold on;
setm(ax_m,'origin',[0,180],'MapLatLimit',fig_cfg.y_lim);
set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
'fontsize',14,'fontname','freeserif','interpreter','none')

%% caxis & colorbar
caxis(ax_m, fig_cfg.c_lim); 
colormap(fig_cfg.c_map);
cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
set(cb, 'Ticks', 1:5);
title(cb,'Y','fontsize',12);

%% draw on ax_m
h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
shading flat;
geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);

%% frame and label setting
setm(ax_m,'frame','on','FLineWidth',1);

label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
mlabel; plabel;
label_y=plabel; label_x=mlabel;
for lxi=1:length(label_x)
    tmp.tmppos=label_x(lxi,1).Position;
    tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55;
    label_x(lxi,1).Position=tmp.tmppos;
    label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
end
for lyi=1:length(label_y)
    label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
end

%% save
dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_predictability_horizon_map', filesep, 'model', filesep, 'mask_n'];
if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
cfg.figname=[dirs.figdir, filesep, 'pred_hor_obs_ano_map_', tmp.varname, '.tif'];

print(fig_h, cfg.figname, '-dpng');
RemoveWhiteSpace([], 'file', cfg.figname);
close all;



%% obs_int n-limit
%% n-limitation filtered corr
fig_cfg.name_rgn = 'Glob';
fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
fig_cfg.x_lim = [-180 180];
fig_cfg.y_lim = [-80 89];
fig_cfg.fig_size = [0,0,6,3.5];
fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
fig_cfg.title_pos = [0.5,0.93];
fig_cfg.p_lim =0.1;
fig_cfg.c_lim = [0.5 5.5];
[fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('wr_5', tmp.dropboxpath);

tmp.X=grid.tlong([end, 1:end],:);
tmp.Y=grid.tlat([end, 1:end],:);
tmp.C=pred_hor_int.*mask.nlim.comb;
tmp.C=tmp.C([end, 1:end],:);

fig_cfg.fig_name=['pred_hor, ', tmp. varname];
fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
    'PaperPosition',fig_cfg.fig_size, ...
    'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0], ...
    'visible','off');
%% map setting
ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
    'fontname','freeserif'); 

axis off; 
hold on;
setm(ax_m,'origin',[0,180],'MapLatLimit',fig_cfg.y_lim);
set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
'fontsize',14,'fontname','freeserif','interpreter','none')

%% caxis & colorbar
caxis(ax_m, fig_cfg.c_lim); 
colormap(fig_cfg.c_map);
cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
set(cb, 'Ticks', 1:5);
title(cb,'Y','fontsize',12);

%% draw on ax_m
h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
shading flat;
geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);

%% frame and label setting
setm(ax_m,'frame','on','FLineWidth',1);

label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
mlabel; plabel;
label_y=plabel; label_x=mlabel;
for lxi=1:length(label_x)
    tmp.tmppos=label_x(lxi,1).Position;
    tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55;
    label_x(lxi,1).Position=tmp.tmppos;
    label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
end
for lyi=1:length(label_y)
    label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
end

%% save
dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_predictability_horizon_map', filesep, 'model', filesep, 'mask_n'];
if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
cfg.figname=[dirs.figdir, filesep, 'pred_hor_obs_int_ano_map_', tmp.varname, '.tif'];

print(fig_h, cfg.figname, '-dpng');
RemoveWhiteSpace([], 'file', cfg.figname);
close all;




end

% 
% pcolor(pred_hor'); 
% shading flat; 
% hc=colorbar; 
% colormap(fig_cfg.c_map)
% caxis([0.5 5.5])
% set(hc, 'Ticks', 1:5);
% % set(hc, 'Ticklabel', {'1','2','3','4','5'});




disp('abc')


function varname_C = f_varname_C(varname)
    switch varname
        case 'temp'
            varname_C='TEMP';
        case 'salt'
            varname_C='SALT';
    end
end



function obsname_simple = f_obs_name(varn)
    switch varn
        case 'SST'
            obsname_simple='ERSST';
        case 'PRECT'
            obsname_simple='GPCP';
        case 'PSL'
            obsname_simple='ERA5';
        case 'SSH'
            obsname_simple='CMEMS';
        case 'TS'
            obsname_simple='HadCRUT5';
        case 'sumChl'
            obsname_simple='OC_CCI';
        otherwise
            obsname_simple='nan';
    end
end


function obsname_simple = f_obs_name_mid(varn)
    switch varn
        case 'SST'
            obsname_simple='ersst_reg_cesm2.v5.';
        case 'PRECT'
            obsname_simple='GPCP_reg_cesm2.v5.';
        case 'PSL'
            obsname_simple='ERA5_msl_reg_cesm2.';
        case 'SSH'
            obsname_simple='CMEMS_reg_cesm2.';
        case 'TS'
            obsname_simple='HadCRUT5_reg_cesm2.';
        case 'sumChl'
            obsname_simple='OC_CCI_reg_cesm2.';
        otherwise
            obsname_simple='nan';
    end
end

function obsname_simple = f_obs_varname(varn)
    switch varn
        case 'SST'
            obsname_simple='sst';
        case 'PRECT'
            obsname_simple='precip';
        case 'PSL'
            obsname_simple='msl';
        case 'SSH'
            obsname_simple='sla';
        case 'TS'
            obsname_simple='tas_mean';
        case 'sumChl'
            obsname_simple='chlor_a';
        otherwise
            obsname_simple='nan';
    end
end

function obsname_simple = f_obs_fname_module(comp)
    switch comp
        case 'ocn'
            obsname_simple='.pop.h.';
        case 'ice'
            obsname_simple='.cice.h.';
        case 'atm'
            obsname_simple='.cam.h0.';
        case 'lnd'
            obsname_simple='.clm2.h0.';
    end
end

function obsname_simple = f_obs_iyears(varn)
    switch varn
        case 'PRECT'
            obsname_simple=1979:2020;
        case 'SSH'
            obsname_simple=1993:2020;
        case 'sumChl'
            obsname_simple=1998:2020;
        otherwise
            obsname_simple=1970:2020;
    end
end