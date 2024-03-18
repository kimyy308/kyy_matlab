% %  Created 23-May-2023 by Yong-Yub Kim
clc; clear all; close all;
warning off;

%% set path
[error_status, tmp.hostname] = system('hostname');
tmp.hostname=tmp.hostname(1:end-1);
switch tmp.hostname
    case 'Yong-Yubs-iMac-Pro.local'
        tmp.dropboxpath = '/Volumes/kyy_raid/kimyy/Dropbox';
        tmp.rootpath = '/Volumes/kyy_raid/';
    case 'Yong-Yubui-MacBookPro.local'
        tmp.dropboxpath = '/Users/kimyy/Dropbox';
    case {'da1', 'da2', 'da3', 'da4'}
        tmp.dropboxpath = '/mnt/lustre/proj/kimyy/Dropbox';
end
tmp.fs=filesep;
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);

%% model configuration

cfg.vlayer=1:10; % 10layer. don't put more than 15
cfg.vlayer_1st=min(cfg.vlayer);
cfg.vlayer_cnt=max(cfg.vlayer)-cfg.vlayer_1st+1;

cfg.var='photoC_TOT_zint_100m';
cfg.obs_name=f_obs_name(cfg.var);
cfg.obs_fname_mid=f_obs_name_mid(cfg.var);
cfg.obs_varname=f_obs_varname(cfg.var);
cfg.comp=Func_0025_CESM2_cmpname_var(cfg.var);
cfg.obs_fname_module=f_obs_fname_module(cfg.comp);
cfg.obs_iyears=1960:2020;

disp(cfg.var);
tic;

dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', cfg.var];
dirs.figroot=[tmp.rootpath, '/kimyy/Figure/CESM2/ESP/HCST_EXP/archive/', cfg.comp,'/', cfg.var];

cfg.iyears=cfg.obs_iyears;
cfg.gnm='f09_g17';
cfg.proj_year=5;

cfg.len_t_y = length(cfg.iyears);
cfg.casename_m = ['ens_all'];

% [tmp.error_status, tmp.value]=system(['ls ', dirs.datadir, '/*once*']);  % b.e21.BHISTsmbb.f09_g17.assm.oras4_ba-10p1.pop.h.once.nc
tmp.gridname = [tmp.rootpath, 'kimyy/Model/CESM2/ESP/HCST_EXP/archive_transfer/', cfg.comp, '/grid.nc'];
tmp.maskname = [tmp.rootpath, 'kimyy/Model/CESM2/ESP/HCST_EXP/archive_transfer/ocn/RECCAP2_region_masks_all_v20210412_POP2_grid.nc'];

switch cfg.comp
    case 'ocn'
        grid.tlong=ncread(tmp.gridname, 'TLONG');
        grid.tlat=ncread(tmp.gridname, 'TLAT');
        grid.mask_ocn=ncread(tmp.maskname, 'open_ocean');
        grid.mask_ocn(grid.mask_ocn<-10e10)=NaN;
        grid.mask_ocn=grid.mask_ocn./grid.mask_ocn;
        grid.tarea=ncread(tmp.gridname, 'TAREA')/1000000.0; %(m^2 -> km^2)
        grid.tarea_60=grid.tarea; grid.tarea_60(grid.tlat>60 | grid.tlat<-60)=NaN;
    case 'atm'
        grid.lon=ncread(tmp.gridname, 'lon');
        grid.lat=ncread(tmp.gridname, 'lat');
        [grid.tlat grid.tlong]=meshgrid(grid.lat, grid.lon);
        grid.tarea=ncread(tmp.gridname, 'AREA');
        grid.tarea_60=grid.tarea; grid.tarea_60(grid.tlat>60 | grid.tlat<-60)=NaN;
end

grid.nlon=size(grid.tlong,1);
grid.nlat=size(grid.tlat,2);
% grid.ntime=cfg.proj_year.*12;

S = shaperead('landareas.shp');

%% read & plot data
tmp.varname=cfg.var;
% lyear = 3;

for lyear = 0:4

tmp.lyear_str=num2str(lyear, '%02i');
dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', cfg.var];
fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', cfg.var, ...
            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
            '_l', tmp.lyear_str, 'y.mat'];
load(fig_cfg.mat_name, 'data');
data_NPP=data;

% tmp.varname='Fe'; %NO3 SiO3 PO4
% dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', tmp.varname];
% fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
%             '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
%             '_l', tmp.lyear_str, 'y.mat'];
% load(fig_cfg.mat_name, 'data');
% data_Fe=data;
% 
% tmp.varname='NO3'; %NO3 SiO3 PO4
% dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', tmp.varname];
% fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
%             '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
%             '_l', tmp.lyear_str, 'y.mat'];
% load(fig_cfg.mat_name, 'data');
% data_NO3=data;
% 
% tmp.varname='SiO3'; %NO3 SiO3 PO4
% dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', tmp.varname];
% fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
%             '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
%             '_l', tmp.lyear_str, 'y.mat'];
% load(fig_cfg.mat_name, 'data');
% data_SiO3=data;
% 
% tmp.varname='PO4'; %NO3 SiO3 PO4
% dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', tmp.varname];
% fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
%             '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
%             '_l', tmp.lyear_str, 'y.mat'];
% load(fig_cfg.mat_name, 'data');
% data_PO4=data;


tmp.varname='diat_Fe_lim_Cweight_avg_100m'; %NO3 SiO3 PO4
dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', tmp.varname];
fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
            '_l', tmp.lyear_str, 'y.mat'];
load(fig_cfg.mat_name, 'data');
data_Fe=data;

tmp.varname='diat_N_lim_Cweight_avg_100m'; %NO3 SiO3 PO4
dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', tmp.varname];
fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
            '_l', tmp.lyear_str, 'y.mat'];
load(fig_cfg.mat_name, 'data');
data_NO3=data;

tmp.varname='diat_SiO3_lim_Cweight_avg_100m'; %NO3 SiO3 PO4
dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', tmp.varname];
fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
            '_l', tmp.lyear_str, 'y.mat'];
load(fig_cfg.mat_name, 'data');
data_SiO3=data;

tmp.varname='diat_P_lim_Cweight_avg_100m'; %NO3 SiO3 PO4
dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', tmp.varname];
fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
            '_l', tmp.lyear_str, 'y.mat'];
load(fig_cfg.mat_name, 'data');
data_PO4=data;

tmp.varname='diat_light_lim_Cweight_avg_100m'; %NO3 SiO3 PO4
dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', tmp.varname];
fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
            '_l', tmp.lyear_str, 'y.mat'];
load(fig_cfg.mat_name, 'data');
data_light=data;

tmp.varname='TEMP'; 
dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', tmp.varname];
fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
            '_l', tmp.lyear_str, 'y.mat'];
load(fig_cfg.mat_name, 'data');
data_TEMP=data;



%% Fe
for loni=1:size(grid.tlong,1)
    for lati=1:size(grid.tlat,2)
        if (data_Fe.(['diat_Fe_lim_Cweight_avg_100m_corr_assm_int_l', tmp.lyear_str])(loni,lati)>0 && ...
                data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_l', tmp.lyear_str])(loni,lati) >0 && ...
                data_Fe.(['diat_Fe_lim_Cweight_avg_100m_corr_assm_int_p_l', tmp.lyear_str])(loni,lati)<0.05 && ...
                data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_p_l', tmp.lyear_str])(loni,lati) <0.05)
            
            [tmp.corr, tmp.corr_p]= corrcoef(data_NPP.(['photoC_TOT_zint_100m_model_l',tmp.lyear_str])(loni,lati,:), ...
                data_Fe.(['diat_Fe_lim_Cweight_avg_100m_model_l', tmp.lyear_str])(loni,lati,:));
            nut_corrset(loni,lati,1) = tmp.corr(1,2);
%          [tmp.corr, tmp.corr_p]=corrcoef(tmp.data(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));
        else
            nut_corrset(loni,lati,1)=0;
        end
    end
end

%% N
for loni=1:size(grid.tlong,1)
    for lati=1:size(grid.tlat,2)
        if (data_NO3.(['diat_N_lim_Cweight_avg_100m_corr_assm_int_l', tmp.lyear_str])(loni,lati)>0 && ...
                data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_l', tmp.lyear_str])(loni,lati) >0 && ...
                data_NO3.(['diat_N_lim_Cweight_avg_100m_corr_assm_int_p_l', tmp.lyear_str])(loni,lati)<0.05 && ...
                data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_p_l', tmp.lyear_str])(loni,lati) <0.05)
            
            [tmp.corr, tmp.corr_p]= corrcoef(data_NPP.(['photoC_TOT_zint_100m_model_l',tmp.lyear_str])(loni,lati,:), ...
                data_NO3.(['diat_N_lim_Cweight_avg_100m_model_l', tmp.lyear_str])(loni,lati,:));
            nut_corrset(loni,lati,2) = tmp.corr(1,2);
%          [tmp.corr, tmp.corr_p]=corrcoef(tmp.data(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));
        else
            nut_corrset(loni,lati,2)=0;
        end
    end
end


%% SiO3
for loni=1:size(grid.tlong,1)
    for lati=1:size(grid.tlat,2)
        if (data_SiO3.(['diat_SiO3_lim_Cweight_avg_100m_corr_assm_int_l', tmp.lyear_str])(loni,lati)>0 && ...
                data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_l', tmp.lyear_str])(loni,lati) >0 && ...
                data_SiO3.(['diat_SiO3_lim_Cweight_avg_100m_corr_assm_int_p_l', tmp.lyear_str])(loni,lati)<0.05 && ...
                data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_p_l', tmp.lyear_str])(loni,lati) <0.05)
            
            [tmp.corr, tmp.corr_p]= corrcoef(data_NPP.(['photoC_TOT_zint_100m_model_l',tmp.lyear_str])(loni,lati,:), ...
                data_SiO3.(['diat_SiO3_lim_Cweight_avg_100m_model_l', tmp.lyear_str])(loni,lati,:));
            nut_corrset(loni,lati,3) = tmp.corr(1,2);
%          [tmp.corr, tmp.corr_p]=corrcoef(tmp.data(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));
        else
            nut_corrset(loni,lati,3)=0;
        end
    end
end

            
%% P
for loni=1:size(grid.tlong,1)
    for lati=1:size(grid.tlat,2)
        if (data_PO4.(['diat_P_lim_Cweight_avg_100m_corr_assm_int_l', tmp.lyear_str])(loni,lati)>0 && ...
                data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_l', tmp.lyear_str])(loni,lati) >0 && ...
                data_PO4.(['diat_P_lim_Cweight_avg_100m_corr_assm_int_p_l', tmp.lyear_str])(loni,lati)<0.05 && ...
                data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_p_l', tmp.lyear_str])(loni,lati) <0.05)
            
           [tmp.corr, tmp.corr_p]= corrcoef(data_NPP.(['photoC_TOT_zint_100m_model_l',tmp.lyear_str])(loni,lati,:), ...
                data_PO4.(['diat_P_lim_Cweight_avg_100m_model_l', tmp.lyear_str])(loni,lati,:));
            nut_corrset(loni,lati,4) = tmp.corr(1,2);
%          [tmp.corr, tmp.corr_p]=corrcoef(tmp.data(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));
        else
            nut_corrset(loni,lati,4)=0;
        end
    end
end

%% Light
for loni=1:size(grid.tlong,1)
    for lati=1:size(grid.tlat,2)
        if (data_light.(['diat_light_lim_Cweight_avg_100m_corr_assm_int_l', tmp.lyear_str])(loni,lati)>0 && ...
                data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_l', tmp.lyear_str])(loni,lati) >0 && ...
                data_light.(['diat_light_lim_Cweight_avg_100m_corr_assm_int_p_l', tmp.lyear_str])(loni,lati)<0.05 && ...
                data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_p_l', tmp.lyear_str])(loni,lati) <0.05)
            
           [tmp.corr, tmp.corr_p]= corrcoef(data_NPP.(['photoC_TOT_zint_100m_model_l',tmp.lyear_str])(loni,lati,:), ...
                data_light.(['diat_light_lim_Cweight_avg_100m_model_l', tmp.lyear_str])(loni,lati,:));
            nut_corrset(loni,lati,5) = tmp.corr(1,2);
%          [tmp.corr, tmp.corr_p]=corrcoef(tmp.data(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));
        else
            nut_corrset(loni,lati,5)=0;
        end
    end
end


%% TEMP
for loni=1:size(grid.tlong,1)
    for lati=1:size(grid.tlat,2)
        if (data_TEMP.(['TEMP_corr_assm_int_l', tmp.lyear_str])(loni,lati)>0 && ...
                data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_l', tmp.lyear_str])(loni,lati) >0 && ...
                data_TEMP.(['TEMP_corr_assm_int_p_l', tmp.lyear_str])(loni,lati)<0.05 && ...
                data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_p_l', tmp.lyear_str])(loni,lati) <0.05)
            
           [tmp.corr, tmp.corr_p]= corrcoef(data_NPP.(['photoC_TOT_zint_100m_model_l',tmp.lyear_str])(loni,lati,:), ...
                data_TEMP.(['TEMP_model_l', tmp.lyear_str])(loni,lati,:));
            nut_corrset(loni,lati,6) = tmp.corr(1,2);
%          [tmp.corr, tmp.corr_p]=corrcoef(tmp.data(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));
        else
            nut_corrset(loni,lati,6)=0;
        end
    end
end
            

         



% [maxval, maxind] = max(abs(nut_corrset), [], 3);
[maxval, maxind] = max(nut_corrset, [], 3);

maxind(isnan(data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_l', tmp.lyear_str])))=NaN;
maxind(maxval==0)=NaN;
maxval(maxval==0)=NaN;

corrval=NaN(size(maxval));
for loni=1:size(grid.tlong,1)
    for lati=1:size(grid.tlat,2)
        if isfinite(maxind(loni,lati))
            corrval(loni,lati)=nut_corrset(loni,lati,maxind(loni,lati));
        end
    end
end

% pcolor(maxval'); shading flat; colorbar;
% pcolor(maxind'); shading flat; colorbar;



%% maxind map --------------------------------------
    fig_cfg.name_rgn = 'Glob';
    fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
%     fig_cfg.map_proj = 'robinson';  % robinson, eqdcylin

    fig_cfg.x_lim = [-180 180];
    fig_cfg.y_lim = [-80 89];
    fig_cfg.fig_size = [0,0,6,3.5];
    fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
    fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
    fig_cfg.title_pos = [0.5,0.93];
    fig_cfg.p_lim =0.1;
    fig_cfg.c_lim = [-1 1];
    [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    tmp.C=maxind;
    tmp.C=tmp.C([end, 1:end],:);

    fig_cfg.fig_name=['ind'];
    fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
        'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
    %% map setting
    ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
        'fontname','freeserif'); 

    axis off; 
    hold on;
    setm(ax_m,'origin',[0,200],'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
    set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
    text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
    'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
    'fontsize',14,'fontname','freeserif','interpreter','none')

    %% caxis & colorbar
    caxis(ax_m, [1 6]); 
    colormap(jet(6));
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
%     title(cb,'R','fontsize',12);

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
        tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55; % y position correction
        label_x(lxi,1).Position=tmp.tmppos;
        label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
    end
    for lyi=1:length(label_y)
        label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
        tmp.tmppos=label_y(lyi,1).Position;
        tmp.tmppos(1)=-fig_cfg.ax_size(3)+2.6; % x position correction
        label_y(lyi,1).Position=tmp.tmppos;
    end

    %% save
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, cfg.var, '_dominant_N'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'dominant_nut_diat', '_l', tmp.lyear_str, 'y_v3.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;



%% corrval map --------------------------------------
    fig_cfg.name_rgn = 'Glob';
    fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
%     fig_cfg.map_proj = 'robinson';  % robinson, eqdcylin

    fig_cfg.x_lim = [-180 180];
    fig_cfg.y_lim = [-80 89];
    fig_cfg.fig_size = [0,0,6,3.5];
    fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
    fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
    fig_cfg.title_pos = [0.5,0.93];
    fig_cfg.p_lim =0.1;
    fig_cfg.c_lim = [-1 1];
    [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    tmp.C=corrval;
    tmp.C=tmp.C([end, 1:end],:);

    fig_cfg.fig_name=['ind'];
    fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
        'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
    %% map setting
    ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
        'fontname','freeserif'); 

    axis off; 
    hold on;
    setm(ax_m,'origin',[0,200],'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
    set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
    text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
    'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
    'fontsize',14,'fontname','freeserif','interpreter','none')

    %% caxis & colorbar
    caxis(ax_m, fig_cfg.c_lim); 
    colormap(fig_cfg.c_map);
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
%     title(cb,'R','fontsize',12);

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
        tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55; % y position correction
        label_x(lxi,1).Position=tmp.tmppos;
        label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
    end
    for lyi=1:length(label_y)
        label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
        tmp.tmppos=label_y(lyi,1).Position;
        tmp.tmppos(1)=-fig_cfg.ax_size(3)+2.6; % x position correction
        label_y(lyi,1).Position=tmp.tmppos;
    end

    %% save
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, cfg.var, '_dominant_N'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'dominant_nut_diat_corr', '_l', tmp.lyear_str, 'y_v3.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;



%% Small Phyto

tmp.varname='sp_Fe_lim_Cweight_avg_100m'; %NO3 SiO3 PO4
dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', tmp.varname];
fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
            '_l', tmp.lyear_str, 'y.mat'];
load(fig_cfg.mat_name, 'data');
data_Fe=data;

tmp.varname='sp_N_lim_Cweight_avg_100m'; %NO3 SiO3 PO4
dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', tmp.varname];
fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
            '_l', tmp.lyear_str, 'y.mat'];
load(fig_cfg.mat_name, 'data');
data_NO3=data;

% tmp.varname='sp_SiO3_lim_Cweight_avg_100m'; %NO3 SiO3 PO4
% dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', tmp.varname];
% fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
%             '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
%             '_l', tmp.lyear_str, 'y.mat'];
% load(fig_cfg.mat_name, 'data');
% data_SiO3=data;

tmp.varname='sp_P_lim_Cweight_avg_100m'; %NO3 SiO3 PO4
dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', tmp.varname];
fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
            '_l', tmp.lyear_str, 'y.mat'];
load(fig_cfg.mat_name, 'data');
data_PO4=data;

tmp.varname='sp_light_lim_Cweight_avg_100m'; %NO3 SiO3 PO4
dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', tmp.varname];
fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
            '_l', tmp.lyear_str, 'y.mat'];
load(fig_cfg.mat_name, 'data');
data_light=data;

tmp.varname='TEMP'; 
dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', tmp.varname];
fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
            '_l', tmp.lyear_str, 'y.mat'];
load(fig_cfg.mat_name, 'data');
data_TEMP=data;





%% Fe
for loni=1:size(grid.tlong,1)
    for lati=1:size(grid.tlat,2)
        if (data_Fe.(['sp_Fe_lim_Cweight_avg_100m_corr_assm_int_l', tmp.lyear_str])(loni,lati)>0 && ...
                data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_l', tmp.lyear_str])(loni,lati) >0 && ...
                data_Fe.(['sp_Fe_lim_Cweight_avg_100m_corr_assm_int_p_l', tmp.lyear_str])(loni,lati)<0.05 && ...
                data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_p_l', tmp.lyear_str])(loni,lati) <0.05)
            
            [tmp.corr, tmp.corr_p]= corrcoef(data_NPP.(['photoC_TOT_zint_100m_model_l',tmp.lyear_str])(loni,lati,:), ...
                data_Fe.(['sp_Fe_lim_Cweight_avg_100m_model_l', tmp.lyear_str])(loni,lati,:));
            nut_corrset(loni,lati,1) = tmp.corr(1,2);
%          [tmp.corr, tmp.corr_p]=corrcoef(tmp.data(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));
        else
            nut_corrset(loni,lati,1)=0;
        end
    end
end

%% N
for loni=1:size(grid.tlong,1)
    for lati=1:size(grid.tlat,2)
        if (data_NO3.(['sp_N_lim_Cweight_avg_100m_corr_assm_int_l', tmp.lyear_str])(loni,lati)>0 && ...
                data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_l', tmp.lyear_str])(loni,lati) >0 && ...
                data_NO3.(['sp_N_lim_Cweight_avg_100m_corr_assm_int_p_l', tmp.lyear_str])(loni,lati)<0.05 && ...
                data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_p_l', tmp.lyear_str])(loni,lati) <0.05)
            
            [tmp.corr, tmp.corr_p]= corrcoef(data_NPP.(['photoC_TOT_zint_100m_model_l',tmp.lyear_str])(loni,lati,:), ...
                data_NO3.(['sp_N_lim_Cweight_avg_100m_model_l', tmp.lyear_str])(loni,lati,:));
            nut_corrset(loni,lati,2) = tmp.corr(1,2);
%          [tmp.corr, tmp.corr_p]=corrcoef(tmp.data(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));
        else
            nut_corrset(loni,lati,2)=0;
        end
    end
end


%% SiO3
for loni=1:size(grid.tlong,1)
    for lati=1:size(grid.tlat,2)
        nut_corrset(loni,lati,3)=0;
    end
end

            
%% P
for loni=1:size(grid.tlong,1)
    for lati=1:size(grid.tlat,2)
        if (data_PO4.(['sp_P_lim_Cweight_avg_100m_corr_assm_int_l', tmp.lyear_str])(loni,lati)>0 && ...
                data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_l', tmp.lyear_str])(loni,lati) >0 && ...
                data_PO4.(['sp_P_lim_Cweight_avg_100m_corr_assm_int_p_l', tmp.lyear_str])(loni,lati)<0.05 && ...
                data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_p_l', tmp.lyear_str])(loni,lati) <0.05)
            
           [tmp.corr, tmp.corr_p]= corrcoef(data_NPP.(['photoC_TOT_zint_100m_model_l',tmp.lyear_str])(loni,lati,:), ...
                data_PO4.(['sp_P_lim_Cweight_avg_100m_model_l', tmp.lyear_str])(loni,lati,:));
            nut_corrset(loni,lati,4) = tmp.corr(1,2);
%          [tmp.corr, tmp.corr_p]=corrcoef(tmp.data(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));
        else
            nut_corrset(loni,lati,4)=0;
        end
    end
end

%% Light
for loni=1:size(grid.tlong,1)
    for lati=1:size(grid.tlat,2)
        if (data_light.(['sp_light_lim_Cweight_avg_100m_corr_assm_int_l', tmp.lyear_str])(loni,lati)>0 && ...
                data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_l', tmp.lyear_str])(loni,lati) >0 && ...
                data_light.(['sp_light_lim_Cweight_avg_100m_corr_assm_int_p_l', tmp.lyear_str])(loni,lati)<0.05 && ...
                data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_p_l', tmp.lyear_str])(loni,lati) <0.05)
            
           [tmp.corr, tmp.corr_p]= corrcoef(data_NPP.(['photoC_TOT_zint_100m_model_l',tmp.lyear_str])(loni,lati,:), ...
                data_light.(['sp_light_lim_Cweight_avg_100m_model_l', tmp.lyear_str])(loni,lati,:));
            nut_corrset(loni,lati,5) = tmp.corr(1,2);
%          [tmp.corr, tmp.corr_p]=corrcoef(tmp.data(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));
        else
            nut_corrset(loni,lati,5)=0;
        end
    end
end


%% TEMP
for loni=1:size(grid.tlong,1)
    for lati=1:size(grid.tlat,2)
        if (data_TEMP.(['TEMP_corr_assm_int_l', tmp.lyear_str])(loni,lati)>0 && ...
                data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_l', tmp.lyear_str])(loni,lati) >0 && ...
                data_TEMP.(['TEMP_corr_assm_int_p_l', tmp.lyear_str])(loni,lati)<0.05 && ...
                data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_p_l', tmp.lyear_str])(loni,lati) <0.05)
            
           [tmp.corr, tmp.corr_p]= corrcoef(data_NPP.(['photoC_TOT_zint_100m_model_l',tmp.lyear_str])(loni,lati,:), ...
                data_TEMP.(['TEMP_model_l', tmp.lyear_str])(loni,lati,:));
            nut_corrset(loni,lati,6) = tmp.corr(1,2);
%          [tmp.corr, tmp.corr_p]=corrcoef(tmp.data(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));
        else
            nut_corrset(loni,lati,6)=0;
        end
    end
end



% [maxval, maxind] = max(abs(nut_corrset), [], 3);
[maxval, maxind] = max(nut_corrset, [], 3);
maxind(isnan(data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_l', tmp.lyear_str])))=NaN;
maxind(maxval==0)=NaN;
maxval(maxval==0)=NaN;
% pcolor(maxval'); shading flat; colorbar;
% pcolor(maxind'); shading flat; colorbar;

corrval=NaN(size(maxval));
for loni=1:size(grid.tlong,1)
    for lati=1:size(grid.tlat,2)
        if isfinite(maxind(loni,lati))
            corrval(loni,lati)=nut_corrset(loni,lati,maxind(loni,lati));
        end
    end
end

    %% model & assm corr map --------------------------------------
    fig_cfg.name_rgn = 'Glob';
    fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
%     fig_cfg.map_proj = 'robinson';  % robinson, eqdcylin

    fig_cfg.x_lim = [-180 180];
    fig_cfg.y_lim = [-80 89];
    fig_cfg.fig_size = [0,0,6,3.5];
    fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
    fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
    fig_cfg.title_pos = [0.5,0.93];
    fig_cfg.p_lim =0.1;
    fig_cfg.c_lim = [-1 1];
    [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    tmp.C=maxind;
    tmp.C=tmp.C([end, 1:end],:);

    fig_cfg.fig_name=['ind'];
    fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
        'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
    %% map setting
    ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
        'fontname','freeserif'); 

    axis off; 
    hold on;
    setm(ax_m,'origin',[0,200],'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
    set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
    text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
    'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
    'fontsize',14,'fontname','freeserif','interpreter','none')

    %% caxis & colorbar
    caxis(ax_m, [1 6]); 
    colormap(jet(6));
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
%     title(cb,'R','fontsize',12);

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
        tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55; % y position correction
        label_x(lxi,1).Position=tmp.tmppos;
        label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
    end
    for lyi=1:length(label_y)
        label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
        tmp.tmppos=label_y(lyi,1).Position;
        tmp.tmppos(1)=-fig_cfg.ax_size(3)+2.6; % x position correction
        label_y(lyi,1).Position=tmp.tmppos;
    end

    %% save
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, cfg.var, '_dominant_N'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'dominant_nut_sp', '_l', tmp.lyear_str, 'y_v3.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;


%% corrval map --------------------------------------
    fig_cfg.name_rgn = 'Glob';
    fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
%     fig_cfg.map_proj = 'robinson';  % robinson, eqdcylin

    fig_cfg.x_lim = [-180 180];
    fig_cfg.y_lim = [-80 89];
    fig_cfg.fig_size = [0,0,6,3.5];
    fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
    fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
    fig_cfg.title_pos = [0.5,0.93];
    fig_cfg.p_lim =0.1;
    fig_cfg.c_lim = [-1 1];
    [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    tmp.C=corrval;
    tmp.C=tmp.C([end, 1:end],:);

    fig_cfg.fig_name=['ind'];
    fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
        'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
    %% map setting
    ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
        'fontname','freeserif'); 

    axis off; 
    hold on;
    setm(ax_m,'origin',[0,200],'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
    set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
    text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
    'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
    'fontsize',14,'fontname','freeserif','interpreter','none')

    %% caxis & colorbar
    caxis(ax_m, fig_cfg.c_lim); 
    colormap(fig_cfg.c_map);
    cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
%     title(cb,'R','fontsize',12);

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
        tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55; % y position correction
        label_x(lxi,1).Position=tmp.tmppos;
        label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
    end
    for lyi=1:length(label_y)
        label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
        tmp.tmppos=label_y(lyi,1).Position;
        tmp.tmppos(1)=-fig_cfg.ax_size(3)+2.6; % x position correction
        label_y(lyi,1).Position=tmp.tmppos;
    end

    %% save
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, cfg.var, '_dominant_N'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'dominant_nut_sp_corr', '_l', tmp.lyear_str, 'y_v3.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;

end


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
        case 'atm'
            obsname_simple='.cam.h0.';
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