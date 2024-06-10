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


hatchflag=1;


%% model configuration
cfg.vlayer=1; % surf, vertical slice 

cfg.vlayer_1st=min(cfg.vlayer);
cfg.vlayer_cnt=max(cfg.vlayer)-cfg.vlayer_1st+1;


cfg.var1='TREFHT';
cfg.var2='PSL';
cfg.var3='PRECT';
    


cfg.var=cfg.var1;
corrval_assm=load(['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/corr_raw/corr_assm_', cfg.var, '_v1_v1.mat']);
corrval_hcst=load(['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/corr_raw/corr_hcst_', cfg.var, '_v1_v1.mat']);
corrval_lens2=load(['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/corr_raw/corr_lens2_', cfg.var, '_v1_v1.mat']);


    grid=corrval_hcst.grid;
    grid.lmask=grid.lfrac;
    grid.lmask(grid.lmask==0)=NaN;
    grid.lmask(isfinite(grid.lmask))=1;
    cfg.gnm='f09_g17';
       
    % [tmp.error_status, tmp.value]=system(['ls ', dirs.datadir, '/*once*']);  % b.e21.BHISTsmbb.f09_g17.assm.oras4_ba-10p1.pop.h.once.nc
    cfg.comp=Func_0025_CESM2_cmpname_var(cfg.var);
    dirs.hcstroot=['/mnt/lustre/proj/earth.system.predictability/HCST_EXP/archive_yearly_transfer/', cfg.comp, '/', cfg.var];
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

    fig_cfg.p_lim =0.05; %95% significance
    fig_cfg.p_lim =0.1; %90% significance
%         fig_cfg.fig_size = [0,0,6.5,3.5]; %% paper size (original)

    

loc_column_first=1;
loc_row_first=10;

%% SUBPLOT(3,2,1); corr, ASSM <-> HCST, median(individual) (LY1)
for subi=1:1
    fig_cfg.fig_size = [0,0,13,14]; %% paper size (original)
    fig_cfg.ax_size = [loc_column_first, loc_row_first, 5.4, 2.7];
    fig_cfg.cb_size = [2, 1, 9, 0.3];
    fig_cfg.title_pos = [0.5,1.02];

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    
    tmp.C=squeeze(corrval_hcst.obs_hcst.ly1.val_median.*grid.lmask);
    tmp.C=tmp.C([end, 1:end],:);

    % significance test
%     sig_n=size(corrval_hcst.data.ly1.([cfg.var,'_ym']),3);
    tmp.D=sum(isfinite(corrval_assm.data_obs.([cfg.var,'_ym'])),3);
    sig_n=tmp.D([end, 1:end],:);
    sig_t=tmp.C.*sqrt(sig_n-2)./sqrt((1-tmp.C.^2));
    sig_tcdf=tcdf(sig_t,sig_n-2);
    for loni=1:size(tmp.C,1)
        for lati=1:size(tmp.C,2)
            if tmp.C(loni,lati)>=0
                tmp.p(loni,lati)=2*(1-sig_tcdf(loni,lati)); % r=positive
            elseif tmp.C(loni,lati)<0
                tmp.p(loni,lati)=2*(sig_tcdf(loni,lati)); % r=negative
            elseif isnan(tmp.C(loni,lati))
                tmp.p(loni,lati)=0;
            end
        end
    end
%     tmp.C(tmp.p>0.1)=NaN; % 90% significant


    fig_cfg.fig_name='$$ (a) \hspace{1mm}  M(r_{O,I}^{\tau=1}) $$';
    fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
        'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','on');

    %% map setting
    subplot(3,2,1);
    ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
        'fontname','freeserif'); 
    set(ax_m, 'Parent', fig_h);
%     subplot(3,3,2,ax_m);
    axis off; 
    hold on;
    setm(ax_m,'origin',[0,160],'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
    set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
        text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize',20,'fontname','freeserif','interpreter','latex')

    %% draw on ax_m
    h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
    shading flat;
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);

    % hatch -> insignificant area
    if hatchflag==1
        tmp.C2=tmp.C;
        tmp.C2(tmp.p<=fig_cfg.p_lim)=NaN;
        pp2 = pcolorm(tmp.Y,tmp.X,tmp.C2, 'parent', ax_m);
        set(pp2,'linestyle','none','Tag','HatchingRegion');
        hp = findobj(pp2,'Tag','HatchingRegion');
        hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','w','HatchLineWidth',0.5);
    end
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
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
    colormap(ax_m,fig_cfg.c_map);
end

%% SUBPLOT(3,2,2); corr, ASSM<->HCST - ASSM <-> LENS2, median(individual) (LY1)
for subi=1:1
    fig_cfg.ax_size = [loc_column_first+5.5, loc_row_first, 5.4, 2.7];

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    
    tmp.A=corrval_hcst.obs_hcst.ly1.val_median;
    tmp.B=corrval_lens2.obs_lens2.val_median;
    
    %% double-sample t-test
    tmp.tt=NaN(size(corrval_lens2.assm_lens2_em.val));
    for loni=1:size(corrval_hcst.obs_hcst.ly1.val,2)
        for lati=1:size(corrval_hcst.obs_hcst.ly1.val,3)
            tmp.a=corrval_hcst.obs_hcst.ly1.val(:,loni,lati);
            tmp.b=corrval_lens2.obs_lens2.val(:,loni,lati);
            tmp.tt(loni,lati)=ttest2(tmp.a,tmp.b,'alpha', fig_cfg.p_lim);
        end
    end

    tmp.tt(isnan(grid.lmask))=1;
    %% get correlation p based on normal distribution, DOF

    tmp.C=squeeze(tmp.A-tmp.B).*grid.lmask;
    tmp.C2=tmp.C;
    tmp.C2(tmp.tt==1)=NaN;
    tmp.C=tmp.C([end, 1:end],:);
    tmp.C2=tmp.C2([end, 1:end],:);
    


    fig_cfg.fig_name='$$ (b) \hspace{1mm}  M(r_{O,I}^{\tau=1}) -  M(r_{O,U}^{\tau=1}) $$';

    %% map setting
    subplot(3,2,2);
    ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
        'fontname','freeserif'); 
    set(ax_m, 'Parent', fig_h);
    axis off; 
    hold on;
    setm(ax_m,'origin',[0,160],'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
    set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
        text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize',20,'fontname','freeserif','interpreter','latex')

    %% draw on ax_m
    h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
    shading flat;
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);

    % hatch -> insignificant area
    if hatchflag ==1
        pp2 = pcolorm(tmp.Y,tmp.X,tmp.C2, 'parent', ax_m);
        set(pp2,'linestyle','none','Tag','HatchingRegion');
        hp = findobj(pp2,'Tag','HatchingRegion');
        hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','w','HatchLineWidth',0.5);
    end
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
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
    colormap(ax_m,fig_cfg.c_map2);
end


cfg.var=cfg.var2;

corrval_assm=load(['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/corr_raw/corr_assm_', cfg.var, '_v1_v1.mat']);
corrval_hcst=load(['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/corr_raw/corr_hcst_', cfg.var, '_v1_v1.mat']);
corrval_lens2=load(['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/corr_raw/corr_lens2_', cfg.var, '_v1_v1.mat']);


%% SUBPLOT(3,2,3); corr, ASSM <-> HCST, median(individual) (LY1)
for subi=1:1
    fig_cfg.ax_size = [loc_column_first, loc_row_first-3.5, 5.4, 2.7];

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);

    tmp.C=squeeze(corrval_hcst.obs_hcst.ly1.val_median);
    tmp.C=tmp.C([end, 1:end],:);

    % significance test
%     sig_n=size(corrval_hcst.data.ly1.([cfg.var,'_ym']),3);
    tmp.D=sum(isfinite(corrval_assm.data_obs.([cfg.var,'_ym'])),3);
    sig_n=tmp.D([end, 1:end],:);    sig_t=tmp.C.*sqrt(sig_n-2)./sqrt((1-tmp.C.^2));
    sig_tcdf=tcdf(sig_t,sig_n-2);
    for loni=1:size(tmp.C,1)
        for lati=1:size(tmp.C,2)
            if tmp.C(loni,lati)>=0
                tmp.p(loni,lati)=2*(1-sig_tcdf(loni,lati)); % r=positive
            elseif tmp.C(loni,lati)<0
                tmp.p(loni,lati)=2*(sig_tcdf(loni,lati)); % r=negative
            end
        end
    end

    fig_cfg.fig_name='$$ (c) \hspace{1mm}  M(r_{O,I}^{\tau=1}) $$';

    %% map setting
    subplot(3,2,3);
    ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
        'fontname','freeserif'); 
    set(ax_m, 'Parent', fig_h);
    axis off; 
    hold on;
    setm(ax_m,'origin',[0,160],'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
    set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
        text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize',20,'fontname','freeserif','interpreter','latex')

    %% draw on ax_m
    h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
    shading flat;
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);

    % hatch -> insignificant area
    if hatchflag==1
        tmp.C2=tmp.C;
        tmp.C2(tmp.p<=fig_cfg.p_lim)=NaN;
        pp2 = pcolorm(tmp.Y,tmp.X,tmp.C2, 'parent', ax_m);
        set(pp2,'linestyle','none','Tag','HatchingRegion');
        hp = findobj(pp2,'Tag','HatchingRegion');
        hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','w','HatchLineWidth',0.5);
    end
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
    %% frame and label setting
    setm(ax_m,'frame','off','FLineWidth',1);

    label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
    label_x=mlabel('MLabelParallel','south', 'MLineLocation',30, 'MLabelLocation',90, 'labelrotation','on');
    mlabel; plabel;
    label_y=plabel; label_x=mlabel;
    for lxi=1:length(label_x)
        tmp.tmppos=label_x(lxi,1).Position;
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
    colormap(ax_m,fig_cfg.c_map);
end

%% caxis & colorbar (ACC)
    caxis(ax_m, fig_cfg.c_lim); 
    colormap(ax_m,fig_cfg.c_map);
    cb = colorbar(ax_m,'units','inches','Location', 'southoutside', 'position',fig_cfg.cb_size + [0, 1, 0, 0]);
    set(cb,'fontsize',15,'fontname','freeserif','TickDir','both');
    cb_title=title(cb,'$$ r $$','fontsize', 22, 'Position', [660, 0, 0]); % hor, ver, ?
    set(cb_title, 'interpreter', 'latex');


%% SUBPLOT(3,2,4); corr, ASSM<->HCST - ASSM <-> LENS2, median(individual) (LY1)
for subi=1:1
    fig_cfg.ax_size = [loc_column_first+5.5, loc_row_first-3.5, 5.4, 2.7];

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);
    tmp.A=corrval_hcst.obs_hcst.ly1.val_median;
    tmp.B=corrval_lens2.obs_lens2.val_median;

    %% double-sample t-test
    tmp.tt=NaN(size(corrval_lens2.assm_lens2_em.val));
    for loni=1:size(corrval_hcst.obs_hcst.ly1.val,2)
        for lati=1:size(corrval_hcst.obs_hcst.ly1.val,3)
            tmp.a=corrval_hcst.obs_hcst.ly1.val(:,loni,lati);
            tmp.b=corrval_lens2.obs_lens2.val(:,loni,lati);
            tmp.tt(loni,lati)=ttest2(tmp.a,tmp.b,'alpha', fig_cfg.p_lim);
        end
    end

    %% get correlation p based on normal distribution, DOF
    tmp.C=squeeze(tmp.A-tmp.B);
    tmp.C2=tmp.C;
    tmp.C2(tmp.tt==1)=NaN;
    tmp.C=tmp.C([end, 1:end],:);
    tmp.C2=tmp.C2([end, 1:end],:);
    


    fig_cfg.fig_name='$$ (d) \hspace{1mm}  M(r_{O,I}^{\tau=1}) -  M(r_{O,U}^{\tau=1}) $$';

    %% map setting
    subplot(3,2,4);
    ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
        'fontname','freeserif'); 
    set(ax_m, 'Parent', fig_h);
%     subplot(3,3,2,ax_m);
    axis off; 
    hold on;
    setm(ax_m,'origin',[0,160],'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
    set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
        text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize',20,'fontname','freeserif','interpreter','latex')

    %% draw on ax_m
    h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
    shading flat;
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
    
    % hatch -> insignificant area
    if hatchflag==1
        pp2 = pcolorm(tmp.Y,tmp.X,tmp.C2, 'parent', ax_m);
        set(pp2,'linestyle','none','Tag','HatchingRegion');
        hp = findobj(pp2,'Tag','HatchingRegion');
        hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','w','HatchLineWidth',0.5);
    end
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
    %% frame and label setting
    setm(ax_m,'frame','off','FLineWidth',1);

    label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
    label_x=mlabel('MLabelParallel','south', 'MLineLocation',30, 'MLabelLocation',90, 'labelrotation','on');
    mlabel; plabel;
    label_y=plabel; label_x=mlabel;
    for lxi=1:length(label_x)
        tmp.tmppos=label_x(lxi,1).Position;
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
    colormap(ax_m,fig_cfg.c_map2);
end


cfg.var=cfg.var3;

corrval_assm=load(['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/corr_raw/corr_assm_', cfg.var, '_v1_v1.mat']);
corrval_hcst=load(['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/corr_raw/corr_hcst_', cfg.var, '_v1_v1.mat']);
corrval_lens2=load(['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/corr_raw/corr_lens2_', cfg.var, '_v1_v1.mat']);

%% SUBPLOT(3,2,5); corr, ASSM <-> HCST, median(individual) (LY1)
for subi=1:1
    fig_cfg.ax_size = [loc_column_first, loc_row_first-7, 5.4, 2.7];

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);

    tmp.C=squeeze(corrval_hcst.obs_hcst.ly1.val_median);
    tmp.C=tmp.C([end, 1:end],:);

    % significance test
%     sig_n=size(corrval_hcst.data.ly1.([cfg.var,'_ym']),3);
    tmp.D=sum(isfinite(corrval_assm.data_obs.([cfg.var,'_ym'])),3);
    sig_n=tmp.D([end, 1:end],:);    
    sig_t=tmp.C.*sqrt(sig_n-2)./sqrt((1-tmp.C.^2));
    sig_tcdf=tcdf(sig_t,sig_n-2);
    for loni=1:size(tmp.C,1)
        for lati=1:size(tmp.C,2)
            if tmp.C(loni,lati)>=0
                tmp.p(loni,lati)=2*(1-sig_tcdf(loni,lati)); % r=positive
            elseif tmp.C(loni,lati)<0
                tmp.p(loni,lati)=2*(sig_tcdf(loni,lati)); % r=negative
            end
        end
    end

    fig_cfg.fig_name='$$ (e) \hspace{1mm}  M(r_{O,I}^{\tau=1}) $$';

    %% map setting
    subplot(3,2,5);
    ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
        'fontname','freeserif'); 
    set(ax_m, 'Parent', fig_h);
%     subplot(3,3,2,ax_m);
    axis off; 
    hold on;
    setm(ax_m,'origin',[0,160],'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
    set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
        text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize',20,'fontname','freeserif','interpreter','latex')

    %% draw on ax_m
    h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
    shading flat;
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);

    % hatch -> insignificant area
    if hatchflag==1
        tmp.C2=tmp.C;
        tmp.C2(tmp.p<=fig_cfg.p_lim)=NaN;
        pp2 = pcolorm(tmp.Y,tmp.X,tmp.C2, 'parent', ax_m);
        set(pp2,'linestyle','none','Tag','HatchingRegion');
        hp = findobj(pp2,'Tag','HatchingRegion');
        hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','w','HatchLineWidth',0.5);
    end
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
    %% frame and label setting
    setm(ax_m,'frame','off','FLineWidth',1);

    label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
    label_x=mlabel('MLabelParallel','south', 'MLineLocation',30, 'MLabelLocation',90, 'labelrotation','on');
    mlabel; plabel;
    label_y=plabel; label_x=mlabel;
    for lxi=1:length(label_x)
        tmp.tmppos=label_x(lxi,1).Position;
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
    colormap(ax_m,fig_cfg.c_map);
end

%% SUBPLOT(3,2,6); corr, ASSM<->HCST - ASSM <-> LENS2, median(individual) (LY1)
for subi=1:1
    fig_cfg.ax_size = [loc_column_first+5.5, loc_row_first-7, 5.4, 2.7];

    tmp.X=grid.tlong([end, 1:end],:);
    tmp.Y=grid.tlat([end, 1:end],:);

    tmp.A=corrval_hcst.obs_hcst.ly1.val_median;
    tmp.B=corrval_lens2.obs_lens2.val_median;
    
    %% double-sample t-test
    tmp.tt=NaN(size(corrval_lens2.assm_lens2_em.val));

    for loni=1:size(corrval_hcst.obs_hcst.ly1.val,2)
        for lati=1:size(corrval_hcst.obs_hcst.ly1.val,3)
            tmp.a=corrval_hcst.obs_hcst.ly1.val(:,loni,lati);
            tmp.b=corrval_lens2.obs_lens2.val(:,loni,lati);
            tmp.tt(loni,lati)=ttest2(tmp.a,tmp.b,'alpha', fig_cfg.p_lim);
        end
    end
    

    %% get correlation p based on normal distribution, DOF

    tmp.C=squeeze(tmp.A-tmp.B);
    tmp.C2=tmp.C;
    tmp.C2(tmp.tt==1)=NaN;
    tmp.C=tmp.C([end, 1:end],:);
    tmp.C2=tmp.C2([end, 1:end],:);

    fig_cfg.fig_name='$$ (f) \hspace{1mm}  M(r_{O,I}^{\tau=1}) -  M(r_{O,U}^{\tau=1}) $$';

    %% map setting
    subplot(3,2,6);
    ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
        'fontname','freeserif'); 
    set(ax_m, 'Parent', fig_h);
    axis off; 
    hold on;
    setm(ax_m,'origin',[0,160],'MapLatLimit',fig_cfg.y_lim);  % lat origin(middle point), lon origin (middle point)
    set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
        text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
        'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
        'fontsize',20,'fontname','freeserif','interpreter','latex')

    %% draw on ax_m
    h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
    shading flat;
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);

    % hatch -> insignificant area
    if hatchflag==1
        pp2 = pcolorm(tmp.Y,tmp.X,tmp.C2, 'parent', ax_m);
        set(pp2,'linestyle','none','Tag','HatchingRegion');
        hp = findobj(pp2,'Tag','HatchingRegion');
        hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','w','HatchLineWidth',0.5);
    end
    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
    %% frame and label setting
    setm(ax_m,'frame','off','FLineWidth',1);

    label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
    label_x=mlabel('MLabelParallel','south', 'MLineLocation',30, 'MLabelLocation',90, 'labelrotation','on');
    mlabel; plabel;
    label_y=plabel; label_x=mlabel;
    for lxi=1:length(label_x)
        tmp.tmppos=label_x(lxi,1).Position;
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
    colormap(ax_m,fig_cfg.c_map2);
end

%% caxis & colorbar (dACC)
    caxis(ax_m, fig_cfg.c_lim2); 
    colormap(ax_m,fig_cfg.c_map2);
    cb = colorbar(ax_m,'units','inches','Location', 'southoutside', 'position',fig_cfg.cb_size + [0 0 0 0]);
    set(cb,'fontsize',15,'fontname','freeserif','TickDir','both');
    cb_title=title(cb,'$$ \Delta r $$','fontsize', 22, 'Position', [670, -4, 0]); % hor, ver, ?
    set(cb_title, 'interpreter', 'latex');    

    %% save
    cfg.figname=['/Users/kimyy/Desktop/backup/Research/Postdoc/03_IBS/2022_predictability_assimilation_run/paper', ...
        filesep, 'Figureset_raw', filesep, 'fig2','_ATM_obs', '.tif'];
    print(fig_h, cfg.figname, '-dpng');
    close all;



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
