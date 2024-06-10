% %  Created 29-Feb-2024 by Yong-Yub Kim
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
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Common', tmp.fs, 'mca']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Common', tmp.fs, 'order']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);

cfg.ly_s = 1;
cfg.ly_e = 5;





%% model configuration

cfg.vars = {'TS', 'PSL', 'PRECT', 'SST'};
% cfg.vars = {'PRECT', 'PSL'};
% cfg.vars = {'SST'};
% cfg.vars = {'NPP', 'GPP', 'RAIN', 'TOTVEGC', 'FIRE'};
% cfg.vars = {'TWS'};
% cfg.vars = {'SOILWATER_10CM'};
cfg.vars = {'TS', 'PRECT', 'SST', 'PSL'};
% cfg.vars = {'PRECT', 'SST', 'PSL'};
% cfg.vars = {'TS'};
cfg.vars={'SOILWATER_10CM', 'TLAI', 'FAREA_BURNED', 'COL_FIRE_CLOSS', 'TWS'};
% cfg.vars={'SSH', 'photoC_TOT_zint', 'photoC_TOT_zint_100m', 'IRON_FLUX', 'HMXL', 'HBLT'};
cfg.vars={'SSH'};
% cfg.vars={'photoC_TOT_zint_100m'};
cfg.vars={'TLAI'};
% tmp.dimids= [1, 2, 4];
cfg.vars={'HMXL', 'HBLT'};
cfg.vars={'NO3', 'SALT'};
cfg.vars={'NO3', 'SALT', 'mul_VVEL_NO3', 'mul_WVEL_NO3','mul_UVEL_NO3'};
cfg.vars={'mul_VVEL_NO3', 'mul_WVEL_NO3','mul_UVEL_NO3'};

% cfg.vars={'TS', 'SST', 'PRECT', 'TWS', 'GPP', 'FAREA_BURNED', 'photoC_TOT_zint_100m'};
cfg.vars={'TS', 'SST', 'PRECT', 'PSL', 'TWS', 'GPP', ...
    'photoC_TOT_zint_100m', 'COL_FIRE_CLOSS', 'FAREA_BURNED', 'SOILWATER_10CM', 'TLAI', 'AEROD_v', 'SSH'};
% cfg.vars={'SSH'};
% cfg.vars={'FAREA_BURNED'};
% cfg.vars={'TS'};
% cfg.vars={'PRECT'};
% cfg.vars={'PSL'};
% cfg.vars={'AEROD_v'};
% cfg.vars={'GPP', 'FAREA_BURNED'};
cfg.vars={'TREFHT'};
cfg.vars={'GPP'};
cfg.vars={'PRECT'};
cfg.vars={'PSL'};
cfg.vars={'SST'};
% cfg.vars={'photoC_TOT_zint_100m'};

% cfg.vars={'aice'};
% cfg.vars={'SALT'};

cfg.vlayer=1; % surf, vertical slice 

% cfg.vlayer=1:10; % 10layer. don't put more than 15

cfg.vlayer_1st=min(cfg.vlayer);
cfg.vlayer_cnt=max(cfg.vlayer)-cfg.vlayer_1st+1;


% for vari=1:length(cfg.vars)

    vari=1;
    cfg.var=cfg.vars{vari};
    cfg.obs_name=f_obs_name(cfg.var);
    cfg.obs_fname_mid=f_obs_name_mid(cfg.var);
    cfg.obs_varname=f_obs_varname(cfg.var);
    cfg.comp=Func_0025_CESM2_cmpname_var(cfg.var);
    cfg.obs_fname_module=f_obs_fname_module(cfg.comp);
    cfg.obs_iyears=1965:2020;
    cfg.obs_iyears2=f_obs_iyears(cfg.var);
    
    if strcmp(cfg.comp, 'ocn') || strcmp(cfg.comp, 'ice')
        origin_std=[0,205];
    elseif strcmp(cfg.comp, 'atm') || strcmp(cfg.comp, 'lnd')
        origin_std=[0,160];
    end

    if strcmp(cfg.var, 'SST')
        origin_std=[0,205];
    end
    % dirs.hcstroot=['/mnt/lustre/proj/earth.system.predictability/HCST_EXP/archive/', cfg.comp, '/', cfg.var];
    % dirs.obsroot=['/Volumes/kyy_raid/kimyy/Observation/', cfg.obs_name, '/monthly_reg_', cfg.obs_fname_module(2:4)];
    dirs.figroot=['/Volumes/kyy_raid/kimyy/Figure/CESM2/ESP/HCST_EXP/archive/', cfg.comp,'/', cfg.var, '/raw'];
    % dirs.lens2root=['/mnt/lustre/proj/earth.system.predictability/LENS2/archive_analysis/', cfg.comp, '/', cfg.var];
    % dirs.assmroot=['/mnt/lustre/proj/earth.system.predictability/ASSM_EXP/archive_analysis/', cfg.comp, '/', cfg.var];
    
    dirs.hcstroot=['/mnt/lustre/proj/earth.system.predictability/HCST_EXP/archive_yearly_transfer/', cfg.comp, '/', cfg.var];
    dirs.obsroot=['/mnt/lustre/proj/kimyy/Observation/', cfg.obs_name];
    dirs.matroot=['/mnt/lustre/proj/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp,'/', cfg.var];
    dirs.lens2root=['/mnt/lustre/proj/earth.system.predictability/LENS2/archive_yearly_transfer/', cfg.comp, '/', cfg.var];
    dirs.assmroot=['/mnt/lustre/proj/earth.system.predictability/ASSM_EXP/archive_yearly_transfer/', cfg.comp, '/', cfg.var];    
    
    dirs.corrroot=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/regrid_5deg/statistics/corr'];
   
%     load([dirs.corrroot, '/', 'corr_5deg_all_', cfg.var, '.mat'], 'cfg_assm', 'cfg_hcst', 'cfg_lens2', 'corrval', 'grid');


% % % % clear all;
corrval_assm=load(['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/corr_raw/corr_assm_', cfg.var, '_v1_v1.mat']);
corrval_hcst=load(['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/corr_raw/corr_hcst_', cfg.var, '_v1_v1.mat']);
corrval_lens2=load(['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/corr_raw/corr_lens2_', cfg.var, '_v1_v1.mat']);

% % % % assm check
% % % % pcolor(corrval_assm.data.([cfg.var,'_ym'])(:,:,1)'); shading flat; colorbar;
% % % pcolor(corrval_assm.data.([cfg.var,'_ym'])(:,:,56)'); shading flat; colorbar;
% % % pcolor(corrval_assm.data.spr(:,:,1)'); shading flat; colorbar;
% % % pcolor(corrval_assm.data.spr(:,:,56)'); shading flat; colorbar;
% % % pcolor(squeeze(corrval_assm.obs_assm.val(1,:,:))'); shading flat; colorbar;
% % % pcolor(squeeze(corrval_assm.obs_assm.val(20,:,:))'); shading flat; colorbar;
% % % pcolor(squeeze(corrval_assm.obs_assm.val_mean(:,:))'); shading flat; colorbar;
% % % pcolor(squeeze(corrval_assm.obs_assm_em.val(:,:))'); shading flat; colorbar;
% % % 
% % % % recasted_val_lens2=NaN(size(corrval_assm.obs_assm.val_median));
% % % % for reci=1:length(corrval_assm.grid.valid_ind)
% % % %     recasted_val_lens2(corrval_assm.grid.valid_ind_i(reci), corrval_assm.grid.valid_ind_j(reci))=corrval_lens2.obs_lens2.val_median(1,reci);
% % % % end
% % % % recasted_val_lens2(recasted_val_lens2==0)=NaN;
% % % pcolor(corrval_assm.obs_assm.val_median'-corrval_lens2.obs_lens2.val_median'); shading flat; colorbar;
% % % % pcolor(recasted_val_lens2'); shading flat; colorbar;
% % % 
% % % [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_20', tmp.dropboxpath);
% % % colormap(ax_m, fig_cfg.c_map); caxis([-1 1])
% % % 
% % % 
% % % % hcst check
% % % pcolor(corrval_hcst.data.ly1.([cfg.var,'_ym'])(:,:,1)'); shading flat; colorbar;
% % % pcolor(corrval_hcst.data.ly5.([cfg.var,'_ym'])(:,:,56)'); shading flat; colorbar;
% % % for reci=1:length(corrval_hcst.grid.valid_ind)
% % %     recasted_val(corrval_hcst.grid.valid_ind_i(reci), corrval_hcst.grid.valid_ind_j(reci))=corrval_hcst.assm_hcst.ly1.val_mean(reci);
% % % end
% % % pcolor(recasted_val'); shading flat; colorbar;
% % % 
% % % for reci=1:length(corrval_hcst.grid.valid_ind)
% % %     recasted_val(corrval_hcst.grid.valid_ind_i(reci), corrval_hcst.grid.valid_ind_j(reci))=corrval_hcst.assm_hcst.ly5.val_mean(reci);
% % % end
% % % pcolor(recasted_val'); shading flat; colorbar;
% % % % plot(squeeze(corrval_hcst.data.ly5.NO3_ym(91,95,:)))
% % % 
% % % for reci=1:length(corrval_hcst.grid.valid_ind)
% % %     recasted_val(corrval_hcst.grid.valid_ind_i(reci), corrval_hcst.grid.valid_ind_j(reci))=corrval_hcst.assm_hcst.ly5.val_median(reci);
% % % end
% % % pcolor(recasted_val'); shading flat; colorbar;
% % % 
% % % for reci=1:length(corrval_hcst.grid.valid_ind)
% % %     recasted_val(corrval_hcst.grid.valid_ind_i(reci), corrval_hcst.grid.valid_ind_j(reci))=corrval_hcst.assm_hcst.ly5.p(1,reci);
% % % end
% % % pcolor(recasted_val'); shading flat; colorbar;
% % % 
% % % pcolor(corrval_hcst.assm_hcst_em.ly1.val'); shading flat; colorbar;
% % % pcolor(corrval_hcst.assm_hcst_em.ly1.p'); shading flat; colorbar;
% % % 
% % % pcolor(corrval_hcst.assm_hcst_em.ly5.val'); shading flat; colorbar;
% % % pcolor(corrval_hcst.assm_hcst_em.ly5.p'); shading flat; colorbar;
% % % pcolor(corrval_hcst.assm_hcst_em_int.ly5.val'); shading flat; colorbar;
% % % pcolor(corrval_hcst.assm_hcst_em_int.ly5.p'); shading flat; colorbar;
% % % 
% % % %ly3
% % % for reci=1:length(corrval_assm.grid.valid_ind)
% % %     recasted_val_hcst(corrval_assm.grid.valid_ind_i(reci), corrval_assm.grid.valid_ind_j(reci))=corrval_hcst.assm_hcst.ly3.val_median(1,reci);
% % % end
% % % for reci=1:length(corrval_assm.grid.valid_ind)
% % %     recasted_val_lens2(corrval_assm.grid.valid_ind_i(reci), corrval_assm.grid.valid_ind_j(reci))=corrval_lens2.assm_lens2.val_median(1,reci);
% % % end
% % % recasted_val_hcst(recasted_val_hcst==0)=NaN;
% % % recasted_val_lens2(recasted_val_hcst==0)=NaN;
% % % pcolor(recasted_val_hcst'-recasted_val_lens2'); shading flat; colorbar;
% % % [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
% % % colormap(ax_m, fig_cfg.c_map); caxis([-0.5 0.5])
% % % 
% % % %ly2
% % % for reci=1:length(corrval_assm.grid.valid_ind)
% % %     recasted_val_hcst(corrval_assm.grid.valid_ind_i(reci), corrval_assm.grid.valid_ind_j(reci))=corrval_hcst.assm_hcst.ly2.val_median(1,reci);
% % % end
% % % for reci=1:length(corrval_assm.grid.valid_ind)
% % %     recasted_val_lens2(corrval_assm.grid.valid_ind_i(reci), corrval_assm.grid.valid_ind_j(reci))=corrval_lens2.assm_lens2.val_median(1,reci);
% % % end
% % % recasted_val_hcst(recasted_val_hcst==0)=NaN;
% % % recasted_val_lens2(recasted_val_lens2==0)=NaN;
% % % pcolor(recasted_val_hcst'-recasted_val_lens2'); shading flat; colorbar;
% % % [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
% % % colormap(ax_m, fig_cfg.c_map); caxis([-0.5 0.5])
% % % 
% % % %ly1
% % % for reci=1:length(corrval_assm.grid.valid_ind)
% % %     recasted_val_hcst(corrval_assm.grid.valid_ind_i(reci), corrval_assm.grid.valid_ind_j(reci))=corrval_hcst.assm_hcst.ly1.val_median(1,reci);
% % % end
% % % for reci=1:length(corrval_assm.grid.valid_ind)
% % %     recasted_val_lens2(corrval_assm.grid.valid_ind_i(reci), corrval_assm.grid.valid_ind_j(reci))=corrval_lens2.assm_lens2.val_median(1,reci);
% % % end
% % % recasted_val_hcst(recasted_val_hcst==0)=NaN;
% % % recasted_val_lens2(recasted_val_lens2==0)=NaN;
% % % pcolor(recasted_val_hcst'-recasted_val_lens2'); shading flat; colorbar;
% % % [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
% % % colormap(ax_m, fig_cfg.c_map); caxis([-0.5 0.5])
% % % 
% % % %% lens2 check
% % % for reci=1:length(corrval_assm.grid.valid_ind)
% % %     recasted_val_lens2(corrval_assm.grid.valid_ind_i(reci), corrval_assm.grid.valid_ind_j(reci))=corrval_lens2.assm_lens2.val_median(reci);
% % % end
% % % recasted_val_lens2(recasted_val_lens2==0)=NaN;
% % % pcolor(recasted_val_lens2'); shading flat; colorbar;
% % % 
% % % for reci=1:length(corrval_assm.grid.valid_ind)
% % %     recasted_val_lens2(corrval_assm.grid.valid_ind_i(reci), corrval_assm.grid.valid_ind_j(reci))=corrval_lens2.assm_lens2.p(1,reci);
% % % end
% % % recasted_val_lens2(recasted_val_lens2==0)=NaN;
% % % pcolor(recasted_val_lens2'); shading flat; colorbar;
% % % pcolor(corrval_lens2.obs_lens2.val_median'); shading flat; colorbar;
% % % pcolor(corrval_lens2.assm_lens2_em.val'); shading flat; colorbar;
% % % pcolor(corrval_lens2.assm_lens2_em.p'); shading flat; colorbar;

%0.05=r*sqrt(n-2) / (1-r^2)


%% statistical significance test
% r=0.2630;
% r=corrval_lens2.assm_lens2.val(500,3000)
% n=size(corrval_assm.data.([cfg.var,'_ym']),3)
% t=r*sqrt(n-2)./sqrt((1-r.^2))
% p=2*(1-tcdf(t,n-2)) % r=positive
% corrval_lens2.assm_lens2.p(500,3000)

% r=corrval_hcst.assm_hcst.ly3.val(200,19000)
% n=size(corrval_assm.data.([cfg.var,'_ym']),3)
% t=r*sqrt(n-2)./sqrt((1-r.^2))
% p=2*(tcdf(t,n-2)) % r=negative
% corrval_hcst.assm_hcst.ly3.p(200,19000)


% %% 95% significance low & high range
% conf = 0.95;
% alpha = 1 - conf;
% pLo = alpha/2;
% pUp = 1 - alpha/2;
% crit = tinv([pLo pUp], n-1);
% % xbar = mean(r); % = 0
% xbar = 0; % = 0
% r_crit=sqrt((crit.^2)./(n-2+(crit).^2));
% r_crit=0.2630; % 95% significance range if n=56
% r_crit=0.2220; % 90% significance range if n=56

% p=0.1
% n=size(corrval_assm.data.([cfg.var,'_ym']),3)
% t=tinv(p,n)
% r_standard_plus=(-(1/t)*sqrt(n-2)+sqrt((1/t*sqrt(n-2)).^2+4))/2
% r_standard_minus=(-(1/t)*sqrt(n-2)-sqrt((1/t*sqrt(n-2)).^2+4))/2
% r=r_standard_plus
% t=r*sqrt(n-2)./sqrt((1-r.^2))
% p=tcdf(t,n)

% 
% p=0.95
% n=size(corrval_assm.data.([cfg.var,'_ym']),3)
% r_standard_plus=(-(1/p)*sqrt(n-2)+sqrt((1/p*sqrt(n-2)).^2+4))/2;
% r_standard_minus=(-(1/p)*sqrt(n-2)-sqrt((1/p*sqrt(n-2)).^2+4))/2;
% r=r_standard_minus
% t=r*sqrt(n-2)./sqrt((1-r.^2))
% p=tcdf(t,n)

% 
% r=0.632
% n=10
% t=r*sqrt(n-2)./sqrt((1-r.^2))
% tcdf(t,n-2)
% p=tpdf(t,n-2)


    grid=corrval_hcst.grid;
    cfg.iyears=cfg.obs_iyears;
    cfg.gnm='f09_g17';
    cfg.proj_year=5;

    cfg.max_iy=max(cfg.iyears);
    cfg.min_iy=min(cfg.iyears);
    cfg.len_t_y = length(cfg.iyears);
       
    % [tmp.error_status, tmp.value]=system(['ls ', dirs.datadir, '/*once*']);  % b.e21.BHISTsmbb.f09_g17.assm.oras4_ba-10p1.pop.h.once.nc
    tmp.gridname = [dirs.hcstroot, tmp.fs, '../grid.nc'];
    tmp.maskname = '/mnt/lustre/proj/kimyy/Model/CESM2/ESP/grids/ocn/RECCAP2_region_masks_all_v20210412_POP2_grid.nc';
    
    % plot set, S.
    S = shaperead('landareas.shp');

    %% read & plot data
    tmp.varname=cfg.var;
    
    fig_cfg.name_rgn = 'Glob';
    fig_cfg.map_proj = 'robinson';  % robinson, eqdcylin

    fig_cfg.x_lim = [-180 180];
    fig_cfg.y_lim = [-80 89];
    fig_cfg.c_lim = [-1 1];
    fig_cfg.c_lim2 = [-0.5 0.5];
    [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_20', tmp.dropboxpath);
    [fig_cfg.c_map2, tmp.err_stat] = Func_0009_get_colormaps('bwg_10', tmp.dropboxpath);
%     fig_cfg.p_lim =0.05; %95% significance
    fig_cfg.p_lim =0.1; %90% significance
%         fig_cfg.fig_size = [0,0,6.5,3.5]; %% paper size (original)

    

%% SUBPLOT(3,3,1); corr, ASSM <-> HCST, ensmean (LY1)
for subi=1:1
    fig_cfg.fig_size = [0,0,18,14]; %% paper size (original)
    fig_cfg.ax_size = [1.0, 10, 5.4, 2.7];
    fig_cfg.cb_size = [3, 1, 12, 0.3];
    fig_cfg.title_pos = [0.5,1.02];

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    
    tmp.C=squeeze(corrval_hcst.obs_hcst_em.ly1.val);
    tmp.C=tmp.C([end, 1:end],:);
    
    % significance test
    for loni=1:size(tmp.C,1)-1
        for lati=1:size(tmp.C,2)
            sig_n=size(corrval_assm.data_obs.([cfg.var,'_ym'])(loni,lati,:),3);
            sig_t=tmp.C(loni,lati)*sqrt(sig_n-2)./sqrt((1-tmp.C(loni,lati).^2));
            sig_tcdf=tcdf(sig_t,sig_n-2);

            if tmp.C(loni,lati)>=0
                tmp.p(loni,lati)=2*(1-sig_tcdf); % r=positive
            elseif tmp.C(loni,lati)<0
                tmp.p(loni,lati)=2*(sig_tcdf); % r=negative
            end
        end
    end
    tmp.p=tmp.p([end, 1:end],:);
%     tmp.C(tmp.p>0.1)=NaN; % 90% significant

%     pcolor(tmp.p'); shading flat; colorbar;

    fig_cfg.fig_name='$$ (a) \hspace{1mm}  r_{O,E(I)}^{\tau=1} $$';
    fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
        'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','on');
    %% map setting
    subplot(3,3,1);

    ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
        'fontname','freeserif'); 
    set(ax_m, 'Parent', fig_h);
%             subplot(3,3,1,ax_m);

    axis off; 
    hold on;
    setm(ax_m,'origin',origin_std,'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
    set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
    text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize', 20,'fontname','freeserif','interpreter','latex')

    %% draw on ax_m
    h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
    shading flat;
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);

    % hatch -> insignificant area
    tmp.C2=tmp.C;
    tmp.C2(tmp.p<=fig_cfg.p_lim)=NaN;
    if sum(isfinite(tmp.C2(:)))>0
        pp2 = pcolorm(tmp.Y,tmp.X,tmp.C2, 'parent', ax_m);
        set(pp2,'linestyle','none','Tag','HatchingRegion');
        hp = findobj(pp2,'Tag','HatchingRegion');
        hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','w','HatchLineWidth',0.5);
    end

    %% frame and label setting
    setm(ax_m,'frame','off','FLineWidth',1);

    label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
    label_x=mlabel('MLabelParallel','south', 'MLineLocation',30, 'MLabelLocation',90, 'labelrotation','on');
    mlabel; plabel;
    label_y=plabel; label_x=mlabel;
    for lxi=1:length(label_x)
        tmp.tmppos=label_x(lxi,1).Position;
%         tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.67; % y position correction
        label_x(lxi,1).Position=tmp.tmppos;
        label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
        label_x(lxi,1).String{2} = ['$$ ', label_x(lxi,1).String{2}, ' $$']; % latex grammar
        set(label_x,'Interpreter','latex');
    end

    labelcorr=[3.5, 3.2, 2.9, 2.75, 2.7, 2.75, 2.9, 3.2, 3.5];
    for lyi=1:length(label_y)
        label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
        tmp.tmppos=label_y(lyi,1).Position;
        tmp.tmppos(1)=-fig_cfg.ax_size(3)+labelcorr(lyi); % x position correction 2.7=avg, for 0 deg
        label_y(lyi,1).Position=tmp.tmppos;
        label_y(lyi,1).String = ['$$ ', label_y(lyi,1).String, ' $$']; % latex grammar
        set(label_y,'Interpreter','latex');
    end

    %% color set
    caxis(ax_m, fig_cfg.c_lim);
    colormap(ax_m, fig_cfg.c_map);
end

%% SUBPLOT(3,3,2); corr, ASSM <-> HCST, median(individual) (LY1)
for subi=1:1
    fig_cfg.ax_size = [6.5, 10, 5.4, 2.7];

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    
%     recasted_val_hcst=NaN(size(corrval_hcst.assm_hcst_em.ly1.val));
%     for reci=1:length(corrval_assm.grid.valid_ind)
%         recasted_val_hcst(corrval_assm.grid.valid_ind_i(reci), corrval_assm.grid.valid_ind_j(reci))=corrval_hcst.assm_hcst.ly1.val_median(1,reci);
%     end
    
%     corrval_hcst.assm_hcst.ly1.val_median

    tmp.C=squeeze(corrval_hcst.obs_hcst.ly1.val_median);
    tmp.C=tmp.C([end, 1:end],:);
    
    tmp=rmfield(tmp,'p');
    % significance test
    for loni=1:size(tmp.C,1)-1
        for lati=1:size(tmp.C,2)
            sig_n=size(corrval_assm.data_obs.([cfg.var,'_ym'])(loni,lati,:),3);
            sig_t=tmp.C(loni,lati)*sqrt(sig_n-2)./sqrt((1-tmp.C(loni,lati).^2));
            sig_tcdf=tcdf(sig_t,sig_n-2);

            if tmp.C(loni,lati)>=0
                tmp.p(loni,lati)=2*(1-sig_tcdf); % r=positive
            elseif tmp.C(loni,lati)<0
                tmp.p(loni,lati)=2*(sig_tcdf); % r=negative
            end
        end
    end
    tmp.p=tmp.p([end, 1:end],:);


    fig_cfg.fig_name='$$ (b) \hspace{1mm}  M(r_{O,I}^{\tau=1}) $$';

    %% map setting
    subplot(3,3,2);
    ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
        'fontname','freeserif'); 
    set(ax_m, 'Parent', fig_h);
%     subplot(3,3,2,ax_m);
    axis off; 
    hold on;
    setm(ax_m,'origin',origin_std,'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
    set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
        text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize',20,'fontname','freeserif','interpreter','latex')

    %% draw on ax_m
    h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
    shading flat;
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);

    % hatch -> insignificant area
    tmp.C2=tmp.C;
    tmp.C2(tmp.p<=fig_cfg.p_lim)=NaN;
    pp2 = pcolorm(tmp.Y,tmp.X,tmp.C2, 'parent', ax_m);
    set(pp2,'linestyle','none','Tag','HatchingRegion');
    hp = findobj(pp2,'Tag','HatchingRegion');
    hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','w','HatchLineWidth',0.5);

    %% frame and label setting
    setm(ax_m,'frame','off','FLineWidth',1);

    label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
    label_x=mlabel('MLabelParallel','south', 'MLineLocation',30, 'MLabelLocation',90, 'labelrotation','on');
    mlabel; plabel;
    label_y=plabel; label_x=mlabel;
    for lxi=1:length(label_x)
        tmp.tmppos=label_x(lxi,1).Position;
%         tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.67; % y position correction
        label_x(lxi,1).Position=tmp.tmppos;
        label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
        label_x(lxi,1).String{2} = ['$$ ', label_x(lxi,1).String{2}, ' $$']; % latex grammar
        set(label_x,'Interpreter','latex');
    end

    labelcorr=[3.4, 3.1, 2.8, 2.7, 2.7, 2.7, 2.8, 3.1, 3.4];
    for lyi=1:length(label_y)
        label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
        tmp.tmppos=label_y(lyi,1).Position;
        tmp.tmppos(1)=-fig_cfg.ax_size(3)+labelcorr(lyi); % x position correction 2.7=avg, for 0 deg
        label_y(lyi,1).Position=tmp.tmppos;
        label_y(lyi,1).String = ['$$ ', label_y(lyi,1).String, ' $$']; % latex grammar
        set(label_y,'Interpreter','latex');
    end

    %% color set
    caxis(ax_m, fig_cfg.c_lim);
    colormap(ax_m, fig_cfg.c_map);
end

%% SUBPLOT(3,3,3); corr, ASSM <-> HCST, median(individual)  (LY2~5) !!! should be changed
for subi=1:1
    fig_cfg.ax_size = [12, 10, 5.4, 2.7];

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    
%     recasted_val_hcst=NaN(size(corrval_hcst.assm_hcst_em.ly1.val));
%     for reci=1:length(corrval_assm.grid.valid_ind)
%         recasted_val_hcst(corrval_assm.grid.valid_ind_i(reci), corrval_assm.grid.valid_ind_j(reci))=corrval_hcst.assm_hcst_4ym.val_median(1,reci);
%     end

    tmp.C=squeeze(corrval_hcst.obs_hcst_4ym.val_median);
    tmp.C=tmp.C([end, 1:end],:);

    tmp=rmfield(tmp,'p');
    % significance test
    for loni=1:size(tmp.C,1)-1
        for lati=1:size(tmp.C,2)
            sig_n=size(corrval_assm.data_obs.([cfg.var,'_4ym'])(loni,lati,:),3);
            sig_t=tmp.C(loni,lati)*sqrt(sig_n-2)./sqrt((1-tmp.C(loni,lati).^2));
            sig_tcdf=tcdf(sig_t,sig_n-2);

            if tmp.C(loni,lati)>=0
                tmp.p(loni,lati)=2*(1-sig_tcdf); % r=positive
            elseif tmp.C(loni,lati)<0
                tmp.p(loni,lati)=2*(sig_tcdf); % r=negative
            end
        end
    end
    tmp.p=tmp.p([end, 1:end],:);



    fig_cfg.fig_name='$$ (c) \hspace{1mm}  M(r_{O,I}^{\tau=2 \textendash 5}) $$';

    %% map setting
    subplot(3,3,3);
    ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
        'fontname','freeserif'); 
    set(ax_m, 'Parent', fig_h);
%     subplot(3,3,2,ax_m);
    axis off; 
    hold on;
    setm(ax_m,'origin',origin_std,'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
    set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
        text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize',20,'fontname','freeserif','interpreter','latex')

    %% draw on ax_m
    h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
    shading flat;
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);

    % hatch -> insignificant area
    tmp.C2=tmp.C;
    tmp.C2(tmp.p<=fig_cfg.p_lim)=NaN;
    pp2 = pcolorm(tmp.Y,tmp.X,tmp.C2, 'parent', ax_m);
    set(pp2,'linestyle','none','Tag','HatchingRegion');
    hp = findobj(pp2,'Tag','HatchingRegion');
    hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','w','HatchLineWidth',0.5);
    
    %% frame and label setting
    setm(ax_m,'frame','off','FLineWidth',1);

    label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
    label_x=mlabel('MLabelParallel','south', 'MLineLocation',30, 'MLabelLocation',90, 'labelrotation','on');
    mlabel; plabel;
    label_y=plabel; label_x=mlabel;
    for lxi=1:length(label_x)
        tmp.tmppos=label_x(lxi,1).Position;
%         tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.67; % y position correction
        label_x(lxi,1).Position=tmp.tmppos;
        label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
        label_x(lxi,1).String{2} = ['$$ ', label_x(lxi,1).String{2}, ' $$']; % latex grammar
        set(label_x,'Interpreter','latex');
    end

    labelcorr=[3.4, 3.1, 2.8, 2.7, 2.7, 2.7, 2.8, 3.1, 3.4];
    for lyi=1:length(label_y)
        label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
        tmp.tmppos=label_y(lyi,1).Position;
        tmp.tmppos(1)=-fig_cfg.ax_size(3)+labelcorr(lyi); % x position correction 2.7=avg, for 0 deg
        label_y(lyi,1).Position=tmp.tmppos;
        label_y(lyi,1).String = ['$$ ', label_y(lyi,1).String, ' $$']; % latex grammar
        set(label_y,'Interpreter','latex');
    end

    %% color set
    caxis(ax_m, fig_cfg.c_lim);
    colormap(ax_m, fig_cfg.c_map);
end

%% SUBPLOT(3,3,4); corr, ASSM <-> LENS2, ensmean (LY1)
for subi=1:1
    fig_cfg.ax_size = [1.0, 6.5, 5.4, 2.7];

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    
    tmp.C=squeeze(corrval_lens2.obs_lens2_em.val);
    tmp.C=tmp.C([end, 1:end],:);
    
    tmp=rmfield(tmp,'p');
    % significance test
    for loni=1:size(tmp.C,1)-1
        for lati=1:size(tmp.C,2)
            sig_n=size(corrval_assm.data_obs.([cfg.var,'_ym'])(loni,lati,:),3);
            sig_t=tmp.C(loni,lati)*sqrt(sig_n-2)./sqrt((1-tmp.C(loni,lati).^2));
            sig_tcdf=tcdf(sig_t,sig_n-2);

            if tmp.C(loni,lati)>=0
                tmp.p(loni,lati)=2*(1-sig_tcdf); % r=positive
            elseif tmp.C(loni,lati)<0
                tmp.p(loni,lati)=2*(sig_tcdf); % r=negative
            end
        end
    end
    tmp.p=tmp.p([end, 1:end],:);

    fig_cfg.fig_name='$$ (d) \hspace{1mm}  r_{O,E(U)}^{\tau=1} $$';
    %% map setting
    subplot(3,3,4);

    ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
        'fontname','freeserif'); 
    set(ax_m, 'Parent', fig_h);
%             subplot(3,3,1,ax_m);

    axis off; 
    hold on;
    setm(ax_m,'origin',origin_std,'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
    set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
    text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize', 20,'fontname','freeserif','interpreter','latex')

    %% draw on ax_m
    h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
    shading flat;
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);

    % hatch -> insignificant area
    tmp.C2=tmp.C;
    tmp.C2(tmp.p<=fig_cfg.p_lim)=NaN;
    pp2 = pcolorm(tmp.Y,tmp.X,tmp.C2, 'parent', ax_m);
    set(pp2,'linestyle','none','Tag','HatchingRegion');
    hp = findobj(pp2,'Tag','HatchingRegion');
    hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','w','HatchLineWidth',0.5);

    %% frame and label setting
    setm(ax_m,'frame','off','FLineWidth',1);

    label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
    label_x=mlabel('MLabelParallel','south', 'MLineLocation',30, 'MLabelLocation',90, 'labelrotation','on');
    mlabel; plabel;
    label_y=plabel; label_x=mlabel;
    for lxi=1:length(label_x)
        tmp.tmppos=label_x(lxi,1).Position;
%         tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.67; % y position correction
        label_x(lxi,1).Position=tmp.tmppos;
        label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
        label_x(lxi,1).String{2} = ['$$ ', label_x(lxi,1).String{2}, ' $$']; % latex grammar
        set(label_x,'Interpreter','latex');
    end

    labelcorr=[3.5, 3.2, 2.9, 2.75, 2.7, 2.75, 2.9, 3.2, 3.5];
    for lyi=1:length(label_y)
        label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
        tmp.tmppos=label_y(lyi,1).Position;
        tmp.tmppos(1)=-fig_cfg.ax_size(3)+labelcorr(lyi); % x position correction 2.7=avg, for 0 deg
        label_y(lyi,1).Position=tmp.tmppos;
        label_y(lyi,1).String = ['$$ ', label_y(lyi,1).String, ' $$']; % latex grammar
        set(label_y,'Interpreter','latex');
    end

    %% color set
    caxis(ax_m, fig_cfg.c_lim);
    colormap(ax_m, fig_cfg.c_map);
end

%% SUBPLOT(3,3,5); corr, ASSM <-> LENS2, median(individual) (LY1)
for subi=1:1
    fig_cfg.ax_size = [6.5, 6.5, 5.4, 2.7];

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    
%     recasted_val_lens2=NaN(size(corrval_lens2.assm_lens2_em.val));
%     for reci=1:length(corrval_assm.grid.valid_ind)
%         recasted_val_lens2(corrval_assm.grid.valid_ind_i(reci), corrval_assm.grid.valid_ind_j(reci))=corrval_lens2.assm_lens2.val_median(1,reci);
%     end
%     corrval_hcst.assm_hcst.ly1.val_median

    tmp.C=squeeze(corrval_lens2.obs_lens2.val_median);
    tmp.C=tmp.C([end, 1:end],:);
    
    tmp=rmfield(tmp,'p');
    % significance test
    for loni=1:size(tmp.C,1)-1
        for lati=1:size(tmp.C,2)
            sig_n=size(corrval_assm.data_obs.([cfg.var,'_ym'])(loni,lati,:),3);
            sig_t=tmp.C(loni,lati)*sqrt(sig_n-2)./sqrt((1-tmp.C(loni,lati).^2));
            sig_tcdf=tcdf(sig_t,sig_n-2);

            if tmp.C(loni,lati)>=0
                tmp.p(loni,lati)=2*(1-sig_tcdf); % r=positive
            elseif tmp.C(loni,lati)<0
                tmp.p(loni,lati)=2*(sig_tcdf); % r=negative
            end
        end
    end
    tmp.p=tmp.p([end, 1:end],:);
    
    fig_cfg.fig_name='$$ (e) \hspace{1mm}  M(r_{O,U}^{\tau=1}) $$';

    %% map setting
    subplot(3,3,5);
    ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
        'fontname','freeserif'); 
    set(ax_m, 'Parent', fig_h);
%     subplot(3,3,2,ax_m);
    axis off; 
    hold on;
    setm(ax_m,'origin',origin_std,'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
    set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
        text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize',20,'fontname','freeserif','interpreter','latex')

    %% draw on ax_m
    h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
    shading flat;
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
    
    % hatch -> insignificant area
    tmp.C2=tmp.C;
    tmp.C2(tmp.p<=fig_cfg.p_lim)=NaN;
    pp2 = pcolorm(tmp.Y,tmp.X,tmp.C2, 'parent', ax_m);
    set(pp2,'linestyle','none','Tag','HatchingRegion');
    hp = findobj(pp2,'Tag','HatchingRegion');
    hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','w','HatchLineWidth',0.5);

    %% frame and label setting
    setm(ax_m,'frame','off','FLineWidth',1);

    label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
    label_x=mlabel('MLabelParallel','south', 'MLineLocation',30, 'MLabelLocation',90, 'labelrotation','on');
    mlabel; plabel;
    label_y=plabel; label_x=mlabel;
    for lxi=1:length(label_x)
        tmp.tmppos=label_x(lxi,1).Position;
%         tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.67; % y position correction
        label_x(lxi,1).Position=tmp.tmppos;
        label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
        label_x(lxi,1).String{2} = ['$$ ', label_x(lxi,1).String{2}, ' $$']; % latex grammar
        set(label_x,'Interpreter','latex');
    end

    labelcorr=[3.4, 3.1, 2.8, 2.7, 2.7, 2.7, 2.8, 3.1, 3.4];
    for lyi=1:length(label_y)
        label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
        tmp.tmppos=label_y(lyi,1).Position;
        tmp.tmppos(1)=-fig_cfg.ax_size(3)+labelcorr(lyi); % x position correction 2.7=avg, for 0 deg
        label_y(lyi,1).Position=tmp.tmppos;
        label_y(lyi,1).String = ['$$ ', label_y(lyi,1).String, ' $$']; % latex grammar
        set(label_y,'Interpreter','latex');
    end

    %% color set
    caxis(ax_m, fig_cfg.c_lim);
    colormap(ax_m, fig_cfg.c_map);
end

%% SUBPLOT(3,3,6); corr, ASSM <-> LENS2, median(individual)  (LY2~5) !!! should be changed
for subi=1:1
    fig_cfg.ax_size = [12, 6.5, 5.4, 2.7];

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    
%     recasted_val_lens2=NaN(size(corrval_lens2.assm_lens2_em.val));
%     for reci=1:length(corrval_assm.grid.valid_ind)
%         recasted_val_lens2(corrval_assm.grid.valid_ind_i(reci), corrval_assm.grid.valid_ind_j(reci))=corrval_lens2.assm_lens2_4ym.val_median(1,reci);
%     end

    tmp.C=squeeze(corrval_lens2.obs_lens2_4ym.val_median);
    tmp.C=tmp.C([end, 1:end],:);
    
    tmp=rmfield(tmp,'p');
    % significance test
    for loni=1:size(tmp.C,1)-1
        for lati=1:size(tmp.C,2)
            sig_n=size(corrval_assm.data_obs.([cfg.var,'_4ym'])(loni,lati,:),3);
            sig_t=tmp.C(loni,lati)*sqrt(sig_n-2)./sqrt((1-tmp.C(loni,lati).^2));
            sig_tcdf=tcdf(sig_t,sig_n-2);

            if tmp.C(loni,lati)>=0
                tmp.p(loni,lati)=2*(1-sig_tcdf); % r=positive
            elseif tmp.C(loni,lati)<0
                tmp.p(loni,lati)=2*(sig_tcdf); % r=negative
            end
        end
    end
    tmp.p=tmp.p([end, 1:end],:);

    fig_cfg.fig_name='$$ (f) \hspace{1mm}  M(r_{O,U}^{\tau=2 \textendash 5}) $$';

    %% map setting
    subplot(3,3,6);
    ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
        'fontname','freeserif'); 
    set(ax_m, 'Parent', fig_h);
%     subplot(3,3,2,ax_m);
    axis off; 
    hold on;
    setm(ax_m,'origin',origin_std,'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
    set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
        text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize',20,'fontname','freeserif','interpreter','latex')

    %% draw on ax_m
    h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
    shading flat;
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);

    % hatch -> insignificant area
    tmp.C2=tmp.C;
    tmp.C2(tmp.p<=fig_cfg.p_lim)=NaN;
    pp2 = pcolorm(tmp.Y,tmp.X,tmp.C2, 'parent', ax_m);
    set(pp2,'linestyle','none','Tag','HatchingRegion');
    hp = findobj(pp2,'Tag','HatchingRegion');
    hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','w','HatchLineWidth',0.5);

    %% frame and label setting
    setm(ax_m,'frame','off','FLineWidth',1);

    label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
    label_x=mlabel('MLabelParallel','south', 'MLineLocation',30, 'MLabelLocation',90, 'labelrotation','on');
    mlabel; plabel;
    label_y=plabel; label_x=mlabel;
    for lxi=1:length(label_x)
        tmp.tmppos=label_x(lxi,1).Position;
%         tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.67; % y position correction
        label_x(lxi,1).Position=tmp.tmppos;
        label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
        label_x(lxi,1).String{2} = ['$$ ', label_x(lxi,1).String{2}, ' $$']; % latex grammar
        set(label_x,'Interpreter','latex');
    end

    labelcorr=[3.4, 3.1, 2.8, 2.7, 2.7, 2.7, 2.8, 3.1, 3.4];
    for lyi=1:length(label_y)
        label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
        tmp.tmppos=label_y(lyi,1).Position;
        tmp.tmppos(1)=-fig_cfg.ax_size(3)+labelcorr(lyi); % x position correction 2.7=avg, for 0 deg
        label_y(lyi,1).Position=tmp.tmppos;
        label_y(lyi,1).String = ['$$ ', label_y(lyi,1).String, ' $$']; % latex grammar
        set(label_y,'Interpreter','latex');
    end

    %% color set
    caxis(ax_m, fig_cfg.c_lim);
    colormap(ax_m, fig_cfg.c_map);
end

%% caxis & colorbar (ACC)
    caxis(ax_m, fig_cfg.c_lim); 
    colormap(ax_m, fig_cfg.c_map);
    cb = colorbar(ax_m,'units','inches','Location', 'northoutside', 'position',fig_cfg.cb_size + [0, 1, 0, 0]);
    set(cb,'fontsize',15,'fontname','freeserif','TickDir','both');
%     cb.TickLabelPosition = 'top';
%     cb_title=title(cb,'$$ r $$','fontsize', 22, 'Position', [880, 0, 0]); % hor, ver, ?
    cb_title=title(cb,'$$ r $$','fontsize', 22, 'Position', [880, 0, 0]); % hor, ver, ?
    set(cb_title, 'interpreter', 'latex');

%% SUBPLOT(3,3,7); corr, ASSM<->HCST - ASSM <-> LENS2, ensmean (LY1)
for subi=1:1
    fig_cfg.ax_size = [1.0, 3, 5.4, 2.7];

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    
%     tic;
%     for loni=1:size(corrval_lens2.data.TREFHT_ym,1)
%         for lati=1:size(corrval_lens2.data.TREFHT_ym,2)
%             tmp.btstrp=bootstrp(10,'corr',corrval_assm.data.TREFHT_ym(loni,lati,:),corrval_lens2.data.TREFHT_ym(loni,lati,:));
%             tmp.B_low(loni,lati)=quantile(tmp.btstrp,0.05);
%             tmp.B_hi(loni,lati)=quantile(tmp.btstrp,0.95);
%         end
%     end
%     toc;
    
%     tmp.tt=NaN(size(corrval_lens2.assm_lens2_em.val)); 
%     tmp.a=corrval_hcst.assm_hcst_em.ly1.val;
%     tmp.b=corrval_lens2.assm_lens2_em.val;
%     tmp.b_low=corrval_lens2.assm_lens2_em.val_050;
%     tmp.b_upper=corrval_lens2.assm_lens2_em.val_950;
%     tmp.ind_f=find(tmp.a>tmp.b_low & tmp.a<tmp.b_upper);
%     tmp.tt(tmp.ind_f)=tmp.a(tmp.ind_f)-tmp.b(tmp.ind_f);
    
%     %% double-sample t-test
%     tmp.tt=NaN(size(corrval_hcst.assm_hcst_em.ly1.val));
%     for loni=1:size(corrval_hcst.assm_hcst_em.ly1.val,1)
%         for lati=1:size(corrval_hcst.assm_hcst_em.ly1.val,2)
%             tmp.a=corrval_hcst.assm_hcst_em.val_btst(loni,lati,:);
%             tmp.b=corrval_lens2.assm_lens2_em.val_btst(loni,lati,:);
%             tmp.tt(loni, lati)=ttest2(tmp.a,tmp.b,'alpha', fig_cfg.p_lim);
%         end
%     end
    
%     corrval_hcst.assm_hcst_em.val_btst
    
    tmp.A=corrval_hcst.obs_hcst_em.ly1.val;
    tmp.B=corrval_lens2.obs_lens2_em.val;
    
     % normal distributed sample -> ttest
    for loni=1:size(corrval_hcst.obs_hcst_em.ly1.val,1)
        for lati=1:size(corrval_hcst.obs_hcst_em.ly1.val,2)
            [tmp.hc_corr,tmp.hc_p,tmp.hc_rl,tmp.hc_ru] = ...
                corrcoef(corrval_hcst.data.ly1.([cfg.var,'_ym'])(loni,lati,:), ...
                corrval_assm.data_obs.([cfg.var,'_ym'])(loni,lati,:),'alpha', 0.10);
            [tmp.le_corr,tmp.le_p,tmp.le_rl,tmp.le_ru] = ...
                corrcoef(corrval_assm.data_obs.([cfg.var,'_ym'])(loni,lati,:), ...
                corrval_lens2.data.([cfg.var,'_ym'])(loni,lati,:), 'alpha', 0.10);
            tmp.hc_sigma=(tmp.hc_ru-tmp.hc_rl)/2;
            tmp.le_sigma=(tmp.le_ru-tmp.le_rl)/2;

%             for i=1:100
%                 tmp.hc_rnd(i)=normrnd(tmp.hc_corr(1,2),tmp.hc_sigma(1,2));
%                 tmp.le_rnd(i)=normrnd(tmp.le_corr(1,2),tmp.le_sigma(1,2));
%             end
            tmp.hc_rnd=tmp.hc_corr(1,2)+(rand(1,20)-0.5).*2.*(tmp.hc_sigma(1,2));
            tmp.le_rnd=tmp.le_corr(1,2)+(rand(1,50)-0.5).*2.*(tmp.le_sigma(1,2));

            tmp.tt(loni, lati)=ttest2(tmp.hc_rnd,tmp.le_rnd,'alpha', fig_cfg.p_lim);
        end
    end


%     [tmp.h, tmp.p]=ttest(tmp.A(:), tmp.B(:), 'Alpha', 0.1);
%     [tmp.p1, tmp.p2, tmp.z, tmp.za, tmp.zb] = ...
%         Func_0036_corr_diff_ttest(tmp.A, tmp.B, size(corrval_hcst.data.ly1.([cfg.var,'_ym']),3), size(corrval_lens2.data.([cfg.var,'_ym']),3));
    
    tmp.C=squeeze(tmp.A-tmp.B);
    tmp.C2=tmp.C;
    tmp.C2(tmp.tt==1)=NaN;

    tmp.C=tmp.C([end, 1:end],:);
    tmp.C2=tmp.C2([end, 1:end],:);
%     tmp.C(tmp.p2>0.1)=NaN;
    

    fig_cfg.fig_name='$$ (g) \hspace{1mm}  r_{O,E(I)}^{\tau=1} - r_{O,E(U)}^{\tau=1} $$';
    %% map setting
    subplot(3,3,7);

    ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
        'fontname','freeserif'); 
    set(ax_m, 'Parent', fig_h);
%             subplot(3,3,1,ax_m);

    axis off; 
    hold on;
    setm(ax_m,'origin',origin_std,'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
    set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
    text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize', 20,'fontname','freeserif','interpreter','latex')

    %% draw on ax_m
    h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
    shading flat;
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);

    % hatch -> insignificant area
    pp2 = pcolorm(tmp.Y,tmp.X,tmp.C2, 'parent', ax_m);
    set(pp2,'linestyle','none','Tag','HatchingRegion');
    hp = findobj(pp2,'Tag','HatchingRegion');
    hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','w','HatchLineWidth',0.5);

    %% frame and label setting
    setm(ax_m,'frame','off','FLineWidth',1);

    label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
    label_x=mlabel('MLabelParallel','south', 'MLineLocation',30, 'MLabelLocation',90, 'labelrotation','on');
    mlabel; plabel;
    label_y=plabel; label_x=mlabel;
    for lxi=1:length(label_x)
        tmp.tmppos=label_x(lxi,1).Position;
%         tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.67; % y position correction
        label_x(lxi,1).Position=tmp.tmppos;
        label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
        label_x(lxi,1).String{2} = ['$$ ', label_x(lxi,1).String{2}, ' $$']; % latex grammar
        set(label_x,'Interpreter','latex');
    end

    labelcorr=[3.5, 3.2, 2.9, 2.75, 2.7, 2.75, 2.9, 3.2, 3.5];
    for lyi=1:length(label_y)
        label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
        tmp.tmppos=label_y(lyi,1).Position;
        tmp.tmppos(1)=-fig_cfg.ax_size(3)+labelcorr(lyi); % x position correction 2.7=avg, for 0 deg
        label_y(lyi,1).Position=tmp.tmppos;
        label_y(lyi,1).String = ['$$ ', label_y(lyi,1).String, ' $$']; % latex grammar
        set(label_y,'Interpreter','latex');
    end

    %% color set
    caxis(ax_m, fig_cfg.c_lim2);
    colormap(ax_m, fig_cfg.c_map2);
end

%% SUBPLOT(3,3,8); corr, ASSM<->HCST - ASSM <-> LENS2, median(individual) (LY1)
for subi=1:1
    fig_cfg.ax_size = [6.5, 3, 5.4, 2.7];

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    
%     tmp.A=NaN(size(corrval_lens2.obs_lens2_em.val));
%     for reci=1:length(corrval_assm.grid.valid_ind)
%         tmp.A(corrval_assm.grid.valid_ind_i(reci), corrval_assm.grid.valid_ind_j(reci))=corrval_hcst.obs_hcst.ly1.val_median(1,reci);
%     end
%     tmp.B=NaN(size(corrval_lens2.obs_lens2_em.val));
%     for reci=1:length(corrval_assm.grid.valid_ind)
%         tmp.B(corrval_assm.grid.valid_ind_i(reci), corrval_assm.grid.valid_ind_j(reci))=corrval_lens2.obs_lens2.val_median(1,reci);
%     end
    
    tmp.A=corrval_hcst.obs_hcst.ly1.val_median;
    tmp.B=corrval_lens2.obs_lens2.val_median;

    %% double-sample t-test
%     tmp.tt=NaN(size(corrval_lens2.assm_lens2_em.val));
%     for reci=1:length(corrval_assm.grid.valid_ind)
%         tmp.a=corrval_hcst.assm_hcst.ly1.val(:,reci);
%         tmp.b=corrval_lens2.assm_lens2.val(:,reci);
%         tmp.tt(corrval_assm.grid.valid_ind_i(reci), corrval_assm.grid.valid_ind_j(reci))=ttest2(tmp.a,tmp.b,'alpha', fig_cfg.p_lim);
%     end

    for loni=1:size(corrval_hcst.obs_hcst.ly1.val,2)
        for lati=1:size(corrval_hcst.obs_hcst.ly1.val,3)
            tmp.a=corrval_hcst.obs_hcst.ly1.val(:,loni,lati);
            tmp.b=corrval_lens2.obs_lens2.val(:,loni,lati);
            tmp.tt(loni,lati)=ttest2(tmp.a,tmp.b,'alpha', fig_cfg.p_lim);
        end
    end

%     %% check significant range of lens2
%     tmp.ttt=NaN(size(corrval_assm.grid.valid_ind));
%     tmp.a=corrval_hcst.assm_hcst.ly1.val_median;
%     tmp.b=corrval_lens2.assm_lens2.val_median;
%     tmp.b_low=quantile(corrval_lens2.assm_lens2.val,0.05,1);
%     tmp.b_upper=quantile(corrval_lens2.assm_lens2.val,0.95,1);
%     tmp.ind_f=find(tmp.a>tmp.b_low & tmp.a<tmp.b_upper);
%     tmp.ttt(tmp.ind_f)=tmp.a(tmp.ind_f)-tmp.b(tmp.ind_f);
%     tmp.tt=NaN(size(corrval_lens2.assm_lens2_em.val));
%     for reci=1:length(corrval_assm.grid.valid_ind)
%         tmp.tt(corrval_assm.grid.valid_ind_i(reci), corrval_assm.grid.valid_ind_j(reci))=tmp.ttt(reci);
%     end

    
    %% pairwise subtraction
%     aimax=size(corrval_hcst.assm_hcst.ly1.val,1);
%     bimax=size(corrval_lens2.assm_lens2.val,1);
%     for reci=1:length(corrval_assm.grid.valid_ind)
%         for ai=1:aimax
%             for bi=1:bimax
%                 tmp.c((ai-1)*bimax+bi)=corrval_hcst.assm_hcst.ly1.val(ai,reci)-corrval_lens2.assm_lens2.val(bi,reci);
%             end
%         end
%         tmp.C(corrval_assm.grid.valid_ind_i(reci), corrval_assm.grid.valid_ind_j(reci))=median(tmp.c);
%     end

    %% get correlation p based on normal distribution, DOF
%     [tmp.p1, tmp.p2, tmp.z, tmp.za, tmp.zb] = ...
%         Func_0036_corr_diff_ttest(tmp.A, tmp.B, size(corrval_hcst.data.ly1.([cfg.var,'_ym']),3), size(corrval_lens2.data.([cfg.var,'_ym']),3));

    tmp.C=squeeze(tmp.A-tmp.B);
    tmp.C2=tmp.C;
    tmp.C2(tmp.tt==1)=NaN;
%     tmp.C2=tmp.tt;
    tmp.C=tmp.C([end, 1:end],:);
    tmp.C2=tmp.C2([end, 1:end],:);
%     tmp.C(tmp.p2>0.1)=NaN;
    


    fig_cfg.fig_name='$$ (h) \hspace{1mm}  M(r_{O,I}^{\tau=1}) -  M(r_{O,U}^{\tau=1}) $$';

    %% map setting
    subplot(3,3,8);
    ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
        'fontname','freeserif'); 
    set(ax_m, 'Parent', fig_h);
%     subplot(3,3,2,ax_m);
    axis off; 
    hold on;
    setm(ax_m,'origin',origin_std,'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
    set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
        text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize',20,'fontname','freeserif','interpreter','latex')

    %% draw on ax_m
    h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
    shading flat;
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);

    % hatch -> insignificant area
    pp2 = pcolorm(tmp.Y,tmp.X,tmp.C2, 'parent', ax_m);
    set(pp2,'linestyle','none','Tag','HatchingRegion');
    hp = findobj(pp2,'Tag','HatchingRegion');
    hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','w','HatchLineWidth',0.5);

    %% frame and label setting
    setm(ax_m,'frame','off','FLineWidth',1);

    label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
    label_x=mlabel('MLabelParallel','south', 'MLineLocation',30, 'MLabelLocation',90, 'labelrotation','on');
    mlabel; plabel;
    label_y=plabel; label_x=mlabel;
    for lxi=1:length(label_x)
        tmp.tmppos=label_x(lxi,1).Position;
%         tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.67; % y position correction
        label_x(lxi,1).Position=tmp.tmppos;
        label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
        label_x(lxi,1).String{2} = ['$$ ', label_x(lxi,1).String{2}, ' $$']; % latex grammar
        set(label_x,'Interpreter','latex');
    end

    labelcorr=[3.4, 3.1, 2.8, 2.7, 2.7, 2.7, 2.8, 3.1, 3.4];
    for lyi=1:length(label_y)
        label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
        tmp.tmppos=label_y(lyi,1).Position;
        tmp.tmppos(1)=-fig_cfg.ax_size(3)+labelcorr(lyi); % x position correction 2.7=avg, for 0 deg
        label_y(lyi,1).Position=tmp.tmppos;
        label_y(lyi,1).String = ['$$ ', label_y(lyi,1).String, ' $$']; % latex grammar
        set(label_y,'Interpreter','latex');
    end

    %% color set
    caxis(ax_m, fig_cfg.c_lim2);
    colormap(ax_m, fig_cfg.c_map2);
end

%% SUBPLOT(3,3,9); corr, ASSM<->HCST - ASSM <-> LENS2, median(individual) (LY2~5) !!! should be changed
for subi=1:1
    fig_cfg.ax_size = [12, 3, 5.4, 2.7];

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    
%     tmp.A=NaN(size(corrval_lens2.assm_lens2_em.val));
%     for reci=1:length(corrval_assm.grid.valid_ind)
%         tmp.A(corrval_assm.grid.valid_ind_i(reci), corrval_assm.grid.valid_ind_j(reci))=corrval_hcst.assm_hcst_4ym.val_median(1,reci);
%     end
%     tmp.B=NaN(size(corrval_lens2.assm_lens2_em.val));
%     for reci=1:length(corrval_assm.grid.valid_ind)
%         tmp.B(corrval_assm.grid.valid_ind_i(reci), corrval_assm.grid.valid_ind_j(reci))=corrval_lens2.assm_lens2_4ym.val_median(1,reci);
%     end

    % double-sample ttest
%     tmp.tt=NaN(size(corrval_lens2.assm_lens2_em.val));
%     for reci=1:length(corrval_assm.grid.valid_ind)
%         tmp.a=corrval_hcst.assm_hcst_4ym.val(:,reci);
%         tmp.b=corrval_lens2.assm_lens2_4ym.val(:,reci);
%         tmp.tt(corrval_assm.grid.valid_ind_i(reci), corrval_assm.grid.valid_ind_j(reci))=ttest2(tmp.a,tmp.b,'alpha', fig_cfg.p_lim);
%     end
    
    tmp.A=corrval_hcst.obs_hcst_4ym.val_median;
    tmp.B=corrval_lens2.obs_lens2_4ym.val_median;

    for loni=1:size(corrval_hcst.obs_hcst.ly1.val,2)
        for lati=1:size(corrval_hcst.obs_hcst.ly1.val,3)
            tmp.a=corrval_hcst.obs_hcst_4ym.val(:,loni,lati);
            tmp.b=corrval_lens2.obs_lens2_4ym.val(:,loni,lati);
            tmp.tt(loni,lati)=ttest2(tmp.a,tmp.b,'alpha', fig_cfg.p_lim);
        end
    end

%     %% check significant range of lens2
%     tmp.ttt=NaN(size(corrval_assm.grid.valid_ind));
%     tmp.a=corrval_hcst.assm_hcst_4ym.val_median;
%     tmp.b=corrval_lens2.assm_lens2_4ym.val_median;
%     tmp.b_low=quantile(corrval_lens2.assm_lens2_4ym.val,0.05,1);
%     tmp.b_upper=quantile(corrval_lens2.assm_lens2_4ym.val,0.95,1);
%     tmp.ind_f=find(tmp.a>tmp.b_low & tmp.a<tmp.b_upper);
%     tmp.ttt(tmp.ind_f)=tmp.a(tmp.ind_f)-tmp.b(tmp.ind_f);
%     tmp.tt=NaN(size(corrval_lens2.assm_lens2_em.val));
%     for reci=1:length(corrval_assm.grid.valid_ind)
%         tmp.tt(corrval_assm.grid.valid_ind_i(reci), corrval_assm.grid.valid_ind_j(reci))=tmp.ttt(reci);
%     end

%     corrval_hcst.assm_hcst.ly1.val_median
    
%     [tmp.p1, tmp.p2, tmp.z, tmp.za, tmp.zb] = ...
%         Func_0036_corr_diff_ttest(tmp.A, tmp.B, size(corrval_hcst.data.ly1.([cfg.var,'_ym']),3), size(corrval_lens2.data.([cfg.var,'_ym']),3));

    tmp.C=squeeze(tmp.A-tmp.B);
    tmp.C2=tmp.C;
    tmp.C2(tmp.tt==1)=NaN;
%     tmp.C2=tmp.tt;
    tmp.C=tmp.C([end, 1:end],:);
    tmp.C2=tmp.C2([end, 1:end],:);

    fig_cfg.fig_name='$$ (i) \hspace{1mm}  M(r_{O,I}^{\tau=2 \textendash 5}) - M(r_{O,U}^{\tau=2 \textendash 5})$$';

    %% map setting
    subplot(3,3,9);
    ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
        'fontname','freeserif'); 
    set(ax_m, 'Parent', fig_h);
%     subplot(3,3,2,ax_m);
    axis off; 
    hold on;
    setm(ax_m,'origin',origin_std,'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
    set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
        text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize',20,'fontname','freeserif','interpreter','latex')

    %% draw on ax_m
    h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
    shading flat;
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);

    % hatch -> insignificant area
    pp2 = pcolorm(tmp.Y,tmp.X,tmp.C2, 'parent', ax_m);
    set(pp2,'linestyle','none','Tag','HatchingRegion');
    hp = findobj(pp2,'Tag','HatchingRegion');
    hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','w','HatchLineWidth',0.5);

    %% frame and label setting
    setm(ax_m,'frame','off','FLineWidth',1);

    label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
    label_x=mlabel('MLabelParallel','south', 'MLineLocation',30, 'MLabelLocation',90, 'labelrotation','on');
    mlabel; plabel;
    label_y=plabel; label_x=mlabel;
    for lxi=1:length(label_x)
        tmp.tmppos=label_x(lxi,1).Position;
%         tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.67; % y position correction
        label_x(lxi,1).Position=tmp.tmppos;
        label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
        label_x(lxi,1).String{2} = ['$$ ', label_x(lxi,1).String{2}, ' $$']; % latex grammar
        set(label_x,'Interpreter','latex');
    end

    labelcorr=[3.4, 3.1, 2.8, 2.7, 2.7, 2.7, 2.8, 3.1, 3.4];
    for lyi=1:length(label_y)
        label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
        tmp.tmppos=label_y(lyi,1).Position;
        tmp.tmppos(1)=-fig_cfg.ax_size(3)+labelcorr(lyi); % x position correction 2.7=avg, for 0 deg
        label_y(lyi,1).Position=tmp.tmppos;
        label_y(lyi,1).String = ['$$ ', label_y(lyi,1).String, ' $$']; % latex grammar
        set(label_y,'Interpreter','latex');
    end

    %% color set
    caxis(ax_m, fig_cfg.c_lim2);
    colormap(ax_m, fig_cfg.c_map2);
end

    %% caxis & colorbar (dACC)
    caxis(ax_m, fig_cfg.c_lim2); 
    colormap(ax_m, fig_cfg.c_map2);
    cb = colorbar(ax_m,'units','inches','Location', 'northoutside', 'position',fig_cfg.cb_size);
    set(cb,'fontsize',15,'fontname','freeserif','TickDir','both');
%     cb.TickLabelPosition = 'top';
    cb_title=title(cb,'$$ \Delta r $$','fontsize', 22, 'Position', [885, -4, 0]); % hor, ver, ?
    set(cb_title, 'interpreter', 'latex');    

    %% save
    dirs.figdir= [dirs.figroot, filesep, tmp.varname, '_corr_assm_map', filesep, 'lens2'];
    if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
    cfg.figname=['/Users/kimyy/Desktop/backup/Research/Postdoc/03_IBS/2022_predictability_assimilation_run/paper', ...
        filesep, 'Figureset_raw', filesep, 'fig1','_obs_',cfg.var, '.tif'];
    print(fig_h, cfg.figname, '-dpng');
%     RemoveWhiteSpace([], 'file', cfg.figname);
    close all;

% end




function scenname = f_scen(year)
    if year<=2014
        scenname='BHISTsmbb';
    else
        scenname='BSSP370smbb';
    end
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
        otherwise
            obsname_simple='nan';
    end
end


function obsname_simple = f_obs_name_mid(varn)
    switch varn
        case 'SST'
            obsname_simple='ersst_reg_5deg.v5.';
        case 'PRECT'
            obsname_simple='GPCC_reg_5deg.v5.';
        case 'RAIN'
            obsname_simple='GPCC_reg_5deg.v5.';
        case 'PSL'
            obsname_simple='ERA5_msl_reg_5deg.';
        case 'SOILWATER_10CM'
%             obsname_simple='SM_reg_5deg.';
            obsname_simple='GLEAM_reg_5deg.v5.';
        case 'TWS'
            obsname_simple='TSW_reg_5deg.';
        case 'SSH'
            obsname_simple='CMEMS_reg_5deg.';
        case 'TS'
%             obsname_simple='HadCRUT5_reg_5deg.';
            obsname_simple='ERA5_t2m_reg_5deg.';
        case 'sumChl'
            obsname_simple='OC_CCI_reg_5deg.';
        case 'TLAI'
            obsname_simple='LAI_reg_5deg.'; % NOAA LAI
        case 'FAREA_BURNED'
            obsname_simple='burned_area_reg_5deg.'; % MODIS Fire_cci v5.1
            obsname_simple='AVHRR-LTDR_reg_5deg.'; % AVHRR-LTDR
        case 'COL_FIRE_CLOSS'
            obsname_simple='FIRE_CLOSS_reg_5deg.'; % GFED
        case 'photoC_TOT_zint'
%             obsname_simple='s_vgpm_reg_5deg.'; % VGPM
%             obsname_simple='ensmean_reg_5deg.'; % VGPM   
            obsname_simple='CMEMS_reg_5deg.'; %Globcolour;
        case 'photoC_TOT_zint_100m'
%             obsname_simple='s_vgpm_reg_5deg.'; % VGPM
%             obsname_simple='ensmean_reg_5deg.'; % VGPM
            obsname_simple='CMEMS_reg_5deg.'; %Globcolour;
        case 'GPP'
%             obsname_simple='ORNL_DAAC_reg_5deg.'; % ORNL_DAAC
            obsname_simple='VODCA2GPP_reg_5deg.'; % VODCA2GPP
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
        case 'RAIN'
            obsname_simple='precip';  % GPCC
        case 'PSL'
            obsname_simple='msl';
        case 'SOILWATER_10CM'
%             obsname_simple='sm';
            obsname_simple='SMsurf'; %GLEAM
        case 'TWS'
            obsname_simple='w';
        case 'SSH'
            obsname_simple='sla';
        case 'TS'
%             obsname_simple='tas_mean'; %HadCRUT5
            obsname_simple='t2m';  %ERA5
        case 'sumChl'
            obsname_simple='chlor_a';
        case 'TLAI'
            obsname_simple='LAI'; % NOAA LAI
        case 'FAREA_BURNED'
            obsname_simple='burned_area'; % MODIS Fire_cci v5.1, AVHRR-LTDR
        case 'COL_FIRE_CLOSS'
            obsname_simple='C'; % GFED
        case 'photoC_TOT_zint'
%             obsname_simple='npp'; % VGPM
            obsname_simple='PP'; %Globcolour
        case 'photoC_TOT_zint_100m'
%             obsname_simple='npp'; % VGPM
            obsname_simple='PP'; %Globcolour
        case 'GPP'
            obsname_simple='GPP';
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
        case 'lnd'
            obsname_simple='.clm2.h0.';
        case 'ice'
            obsname_simple='.cice.h.';
    end
end

function obsname_simple = f_obs_iyears(varn)
    switch varn
        case 'PRECT'
            obsname_simple=1979:2020; % GPCP
            obsname_simple=1960:2019; % GPCC
        case 'RAIN'
            obsname_simple=1960:2019;
        case 'SSH'
            obsname_simple=1993:2020;
        case 'sumChl'
            obsname_simple=1998:2020;
        case 'SOILWATER_10CM'
%             obsname_simple=1979:2020; % C3S Surface Soil Moisture
            obsname_simple=1980:2020; % GLEAM SMsurf
        case 'TLAI'
            obsname_simple=1982:2018; % NOAA LAI
        case 'FAREA_BURNED'
%             obsname_simple=2001:2020; % MODIS Fire_cci v5.1
            obsname_simple=1982:2018; % AVHRR-LTDR
        case 'COL_FIRE_CLOSS'
            obsname_simple=1997:2020; % GFED
        case 'photoC_TOT_zint'
            obsname_simple=2003:2020; % VGPM (ens)
            obsname_simple=1998:2020; % GlobColour            
        case 'photoC_TOT_zint_100m'
%             obsname_simple=2003:2020; % VGPM (ens)
            obsname_simple=1998:2020; % GlobColour
        case 'GPP'
%             obsname_simple=1982:2016; % ORNL_DAAC
            obsname_simple=1989:2019; % VODCA2GPP
        otherwise
            obsname_simple=1960:2020;
    end
end
