% %  Created 30-Jan-2023 by Yong-Yub Kim
clc; clear all; close all;
warning off;

%% set path
[error_status, tmp.hostname] = system('hostname');
tmp.hostname=tmp.hostname(1:end-1);
switch tmp.hostname
    case 'Yong-Yubs-iMac-Pro.local'
        tmp.dropboxpath = '/Volumes/kyy_raid/kimyy/Dropbox';
    case {'da1', 'da2', 'da3', 'da4'}
        tmp.dropboxpath = '/mnt/lustre/proj/kimyy/Dropbox';
end
tmp.fs=filesep;
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);

%% model configuration
cfg.var='SST';

% dirs.root='/mnt/lustre/proj/earth.system.predictability/HCST_EXP';
% dirs.yoshi_root='/proj/yoshi/DATA/CESM2_ODA';
% dirs.archive=[dirs.root, filesep, 'archive'];
% dirs.saveroot='/mnt/lustre/proj/kimyy/Model/CESM2/ESP/HCST_EXP';
dirs.hcstroot=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/HCST_EXP/archive/atm/', cfg.var];
dirs.obsroot=['/Volumes/kyy_raid/kimyy/Observation/ERSST/monthly_reg_cam'];
dirs.figroot=['/Volumes/kyy_raid/kimyy/Figure/CESM2/ESP/HCST_EXP/archive/atm/', cfg.var];

cfg.iyears=1991:2020;
% cfg.months=1:12;
% cfg.scnm='HIST';
cfg.gnm='f09_g17';
% cfg.assm_factor='10';
% cfg.ens_member='1';
% cfg.proj_year=5;
cfg.season = {'MAM', 'JJA', 'SON', 'DJF'};
% cfg.season = {'AMJ', 'JAS', 'OND', 'JFM'};

% cfg.component='ocn';
% cfg.varnames={'temp', 'salt'};
cfg.len_t_y = length(cfg.iyears);
% cfg.len_t_m = length(cfg.months);
% cfg.len_t = cfg.len_t_y *cfg.len_t_m;

% %% grid set(mask from model)
% tmp.obsname=cfg.obsnames{1};
% iyear=min(cfg.iyears);
% cfg.casename_m=[cfg.gridname, '.hcst.', tmp.obsname, '-', cfg.assm_factor, 'p', cfg.ens_member];
% cfg.casename=[cfg.casename_m, '_i', num2str(iyear)];
% dirs.datadir= [dirs.hcstroot, filesep, cfg.casename_m, filesep, 'GMSV'];

% f09_g17.hcst.en4.2_ba-10p1_i2021.pop.h.once.nc

% [tmp.error_status, tmp.value]=system(['ls ', dirs.datadir, '/*once*']);  % b.e21.BHISTsmbb.f09_g17.assm.oras4_ba-10p1.pop.h.once.nc
tmp.gridname = [dirs.hcstroot, tmp.fs, 'grid.nc'];
% grid.region_mask=ncread(tmp.gridname, 'REGION_MASK'); 
% grid.ocean_mask=NaN(size(grid.region_mask));
% grid.ocean_mask(grid.region_mask>0)=1;
% grid.tarea = ncread(tmp.gridname, 'TAREA');

grid.lon=ncread(tmp.gridname, 'lon');
grid.lat=ncread(tmp.gridname, 'lat');
[grid.tlat grid.tlong]=meshgrid(grid.lat, grid.lon);

grid.nlon=size(grid.tlong,1);
grid.nlat=size(grid.tlat,2);
% grid.ntime=cfg.proj_year.*12;


% % model filename example
% /mnt/lustre/proj/kimyy/Model/CESM2/ESP/HCST_EXP/archive/atm/SST/ens_all/ens_all_i2020
% SST_f09_g17.hcst.ens_all_i2020.cam.h0.2024-12.nc

S = shaperead('landareas.shp');

%% read & plot data
tmp.varname=cfg.var;
for lss=1:length(cfg.season)
    tmp.season=cfg.season{lss};
    tmp.mons = f_season_mons(tmp.season);
    cfg.casename_m=['ens_all'];

    dirs.datadir= [dirs.hcstroot, filesep, cfg.casename_m, filesep];
    fprintf('%s_%s_%s  ',tmp.season,cfg.casename_m, tmp.varname); lap_time = tic;
    
    %% variables initialization
    data.([tmp.varname, '_model', '_', tmp.season])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
    data.([tmp.varname, '_bias', '_', tmp.season])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
    data.([tmp.varname, '_obs', '_', tmp.season])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
%     data.([tmp.varname, '_assm'])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
    
    %% read variables
    for iyear=min(cfg.iyears):max(cfg.iyears)
        tmp.iyear_str=num2str(iyear, '%04i');
        tmp.iind=iyear-min(cfg.iyears)+1;
        cfg.casename=[cfg.casename_m, '_i', tmp.iyear_str];
%         tmp.fy=iyear+lss;
%         tmp.fy_str=num2str(tmp.fy, '%04i');
        %% monthly filename
        for mon=1:length(tmp.mons)
            tmp.mon=tmp.mons(mon);
            if tmp.mon>12
                tmp.mon_str=num2str(tmp.mon-12, '%02i');
                tmp.fy_str=num2str(iyear+1, '%04i');
            else
                tmp.mon_str=num2str(tmp.mon, '%02i');
                tmp.fy_str=num2str(iyear, '%04i');
            end
            cfg.mod_fnm=[dirs.datadir, tmp.fs, cfg.casename, tmp.fs, ...
            tmp.varname, '_', cfg.gnm, '.hcst.', cfg.casename, '.cam.h0.', tmp.fy_str, '-', tmp.mon_str, '.nc'];
            tmp.sdata(:,:,mon)=ncread(cfg.mod_fnm, tmp.varname);
            
            cfg.obs_fnm=[dirs.obsroot, tmp.fs, 'ersst_reg_cesm2.v5.',tmp.iyear_str,tmp.mon_str, '.nc'];
            tmp.sdata_obs(:,:,mon)=ncread(cfg.obs_fnm, 'sst');
        end
        tmp.sdata(tmp.sdata==0)=NaN;
%         data.time=ncread(cfg.datafilename, 'time');
        tmp.smean= mean(tmp.sdata,3);
        data.([tmp.varname, '_model', '_', tmp.season])(:,:,tmp.iind)= tmp.smean;
        tmp.smean_obs= mean(tmp.sdata_obs,3);
        data.([tmp.varname, '_obs', '_', tmp.season])(:,:,tmp.iind)= tmp.smean_obs;
%         tmp.ymean= mean(ncread(cfg.datafilename, ['assm_', tmp.varname], [1, 1, 1], [inf, inf, 12]), 3);
%         data.([tmp.varname, '_assm'])(:,:,iyear-min(cfg.iyears)+1)= tmp.ymean;
    end
    
    %% get correlation coefficient
    for loni=1:grid.nlon
        for lati=1:grid.nlat
             if (isnan(data.([tmp.varname, '_obs', '_', tmp.season])(loni,lati,1))~=1 & nansum(data.([tmp.varname, '_model', '_', tmp.season])(loni,lati,:))~=0)
                 tmp.data = squeeze(data.([tmp.varname, '_model', '_', tmp.season])(loni,lati,1:end));
                 tmp.data_obs = data.([tmp.varname, '_obs', '_', tmp.season])(loni,lati,1:end);
                 [tmp.corr, tmp.corr_p]=corrcoef(tmp.data, tmp.data_obs);
                 
                 tmp.data_det = Func_0028_detrend_linear_1d(squeeze(data.([tmp.varname, '_model', '_', tmp.season])(loni,lati,1:end)));
                 tmp.data_obs_det = Func_0028_detrend_linear_1d(data.([tmp.varname, '_obs', '_', tmp.season])(loni,lati,1:end));
                 [tmp.corr_det, tmp.corr_det_p]=corrcoef(tmp.data_det, tmp.data_obs_det);
%                  [tmp.corr, tmp.corr_p]=corrcoef(squeeze(data.([tmp.varname, '_model', '_', tmp.season])(loni,lati,1+lss:end-cfg.proj_year+1)), ...
%                      squeeze(data.([tmp.varname, '_obs', '_', tmp.season])(loni,lati,1+lss:end)));
                 data.([tmp.varname, '_corr_obs', '_', tmp.season])(loni,lati)=tmp.corr(1,2);
                 data.([tmp.varname, '_corr_obs_p', '_', tmp.season])(loni,lati)=tmp.corr_p(1,2);
                 data.([tmp.varname, '_corr_obs_det', '_', tmp.season])(loni,lati)=tmp.corr_det(1,2);
                 data.([tmp.varname, '_corr_obs_det_p', '_', tmp.season])(loni,lati)=tmp.corr_det_p(1,2);
%                  [tmp.corr, tmp.corr_p]=corrcoef(squeeze(data.([tmp.varname, '_model', '_', tmp.season])(loni,lati,1+lss:end-cfg.proj_year+1)), ...
%                      squeeze(data.([tmp.varname, '_assm'])(loni,lati,1+lss:end)));
%                  data.([tmp.varname, '_corr_assm', '_', tmp.season])(loni,lati)=tmp.corr(1,2);
%                  data.([tmp.varname, '_corr_assm_p', '_', tmp.season])(loni,lati)=tmp.corr_p(1,2);
             else
                 data.([tmp.varname, '_corr_obs', '_', tmp.season])(loni,lati)=NaN;
                 data.([tmp.varname, '_corr_obs_p', '_', tmp.season])(loni,lati)=NaN;
                 data.([tmp.varname, '_corr_obs_det', '_', tmp.season])(loni,lati)=NaN;
                 data.([tmp.varname, '_corr_obs_det_p', '_', tmp.season])(loni,lati)=NaN;
%                  data.([tmp.varname, '_corr_assm', '_', tmp.season])(loni,lati)=NaN;
%                  data.([tmp.varname, '_corr_assm_p', '_', tmp.season])(loni,lati)=NaN;   
             end
        end
    end
    disp('abc')
end

for lss_ind=1:length(cfg.season)
%     tmp.lss_str=tmp.yearset{lss_ind};
    tmp.season=cfg.season{lss_ind};
    
    %% model & obs corr map --------------------------------------
    fig_cfg.name_rgn = 'Glob';
    fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
    fig_cfg.x_lim = [-180 180];
    fig_cfg.y_lim = [-80 89];
    fig_cfg.fig_size = [0,0,6,3.5];
    fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
    fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
    fig_cfg.title_pos = [0.5,0.93];
    fig_cfg.p_lim =0.1;
    fig_cfg.c_lim = [-1 1];
    [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('byr', tmp.dropboxpath);

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    tmp.C=data.([tmp.varname, '_corr_obs', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%             tmp.C_H=data.([tmp.varname, '_corr_obs_p', '_', tmp.season]);
%             tmp.C_H=tmp.C_H([end, 1:end],:);
%             tmp.C_2=tmp.C;
%             tmp.C_2(tmp.C_H<fig_cfg.p_lim)=NaN;
    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data.([tmp.varname, '_corr_obs', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', corr_ob_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
    fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
        'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
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
    title(cb,'R','fontsize',12);

    %% draw on ax_m
    h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
    shading flat;
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);

%             pp2 = pcolorm(tmp.Y,tmp.X,tmp.C_2, 'parent', ax_m);
%             set(pp2,'linestyle','none','Tag','HatchingRegion');
%             hp = findobj(pp2,'Tag','HatchingRegion');
%             hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','k','HatchLineWidth',0.5);

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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_obs_map'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'corr_obs_map_', tmp.varname, '_', tmp.season, '.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    fprintf('%7.1f sec\n', toc(lap_time) );
    close all;


%% model & obs corr map (det) --------------------------------------
    fig_cfg.name_rgn = 'Glob';
    fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
    fig_cfg.x_lim = [-180 180];
    fig_cfg.y_lim = [-80 89];
    fig_cfg.fig_size = [0,0,6,3.5];
    fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
    fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
    fig_cfg.title_pos = [0.5,0.93];
    fig_cfg.p_lim =0.1;
    fig_cfg.c_lim = [-1 1];
    [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('byr', tmp.dropboxpath);

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    tmp.C=data.([tmp.varname, '_corr_obs_det', '_', tmp.season]);
    tmp.C=tmp.C([end, 1:end],:);
%             tmp.C_H=data.([tmp.varname, '_corr_obs_p', '_', tmp.season]);
%             tmp.C_H=tmp.C_H([end, 1:end],:);
%             tmp.C_2=tmp.C;
%             tmp.C_2(tmp.C_H<fig_cfg.p_lim)=NaN;
    [tmp.mean_corr, tmp.err] = ...
        Func_0011_get_area_weighted_mean(data.([tmp.varname, '_corr_obs_det', '_', tmp.season]), grid.tlong, grid.tlat);
    fig_cfg.fig_name=[tmp.season, ', corr_ob_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
    fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
        'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
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
    title(cb,'R','fontsize',12);

    %% draw on ax_m
    h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
    shading flat;
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);

%             pp2 = pcolorm(tmp.Y,tmp.X,tmp.C_2, 'parent', ax_m);
%             set(pp2,'linestyle','none','Tag','HatchingRegion');
%             hp = findobj(pp2,'Tag','HatchingRegion');
%             hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','k','HatchLineWidth',0.5);

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
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_obs_map'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=[dirs.figdir, filesep, 'corr_obs_det_map_', tmp.varname, '_', tmp.season, '.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
    fprintf('%7.1f sec\n', toc(lap_time) );
    close all;

end

function obsname_simple = f_obs_simple(obsname)
    switch obsname
        case 'en4.2_ba'
            obsname_simple='en4';
        case 'projdv7.3'
            obsname_simple='projd';
    end
end

function gmsst = f_gm_var(var_2d, area)
    mask_var=NaN(size(var_2d));
    mask_var(isfinite(var_2d))=1;
    area=area.*mask_var;
    var_2d=squeeze(var_2d);
    var_2d_sum=var_2d.*area;
    var_2d_sum=sum(var_2d_sum(:), 'omitnan');
    gmsst=var_2d_sum ./ sum(area(:), 'omitnan');
end

function varname_C = f_varname_C(varname)
    switch varname
        case 'temp'
            varname_C='TEMP';
        case 'salt'
            varname_C='SALT';
    end
end

function varname_unit = f_varname_unit(varname)
    switch varname
        case 'temp'
            varname_unit='\circC';
        case 'salt'
            varname_unit='g/kg';
    end
end

function mons = f_season_mons(season)
    switch season
        case 'MAM'
            mons = [3,4,5];
        case 'AMJ'
            mons = [4,5,6];
        case 'JJA'
            mons = [6,7,8];
        case 'JAS'
            mons = [7,8,9];
        case 'SON'
            mons = [9,10,11];
        case 'OND'
            mons = [10,11,12];
        case 'DJF'
            mons = [12,13,14];
        case 'JFM'
            mons = [13,14,15];
    end
end