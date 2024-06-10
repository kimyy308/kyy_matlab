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

% pcolor(mask.nlim.comb'); shading flat; colorbar;


%% pp figures


cfg.vars={'photoC_TOT_zint_100m', 'subt_tend_zint_100m_NO3_Jint_100m_NO3', 'tend_zint_100m_NO3', 'Jint_100m_NO3'};
cfg.vars={'NO3'};

cfg.vlayer=1:10; % surface, vertical slice
cfg.vlayer_1st=min(cfg.vlayer);
cfg.vlayer_cnt=max(cfg.vlayer)-cfg.vlayer_1st+1;

for vari=1:length(cfg.vars)

    cfg.var=cfg.vars{vari};
    cfg.obs_name=f_obs_name(cfg.var);
    cfg.obs_fname_mid=f_obs_name_mid(cfg.var);
    cfg.obs_varname=f_obs_varname(cfg.var);
    cfg.comp=Func_0025_CESM2_cmpname_var(cfg.var);
    cfg.obs_fname_module=f_obs_fname_module(cfg.comp);
    cfg.obs_iyears=1960:2020;
    
    disp(cfg.var);
    tic;
    
    dirs.hcstmatroot=[tmp.kimyypath, '/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', cfg.var];
    vstr=['v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i')];
    dirs.figroot=[tmp.kimyypath, '/Figure/CESM2/ESP/HCST_EXP/archive_anomaly/', cfg.comp,'/', cfg.var, filesep, vstr];
    mkdir(dirs.figroot)
    
    cfg.iyears=cfg.obs_iyears;
    cfg.gnm='f09_g17';
    cfg.proj_year=5;  %%%%%%%%%%%%%% LY 0,1 only
    
    cfg.len_t_y = length(cfg.iyears);
    cfg.casename_m = ['ens_all'];
    
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
    
    S = shaperead('landareas.shp');
    
    fig_flags(1:100)=1;
    
    tmp.varname=cfg.var;
    for lyear=0:cfg.proj_year-1
        tmp.lyear_str=num2str(lyear, '%02i');
        fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
                '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
                '_l', tmp.lyear_str, 'y.mat'];
        load(fig_cfg.mat_name, 'data', 'data2')
        %% clim time range
        switch cfg.obs_name
           case 'GPCC'
               cfg.clim_ys=1965;
               cfg.clim_ye=2019;
           otherwise
               cfg.clim_ys=1965;
               cfg.clim_ye=2020;
       end
       cfg.clim_tlen = (cfg.clim_ys-1959)-lyear:(cfg.clim_ye-2020)+cfg.len_t_y-lyear;
       cfg.clim_tlen2=length(cfg.clim_tlen);

%% model & assm corr map --------------------------------------
        if fig_flags(11)==1
            for fake=1:1
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
                [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
            
                tmp.X=grid.tlong([end, 1:end],:);
                tmp.Y=grid.tlat([end, 1:end],:);
                tmp.C=data2.([tmp.varname, '_corr_assm_ano', '_l', tmp.lyear_str]).* mask.nlim.comb;
                tmp.C=tmp.C([end, 1:end],:);
            
                [tmp.mean_corr, tmp.err] = ...
                    Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_corr_assm_ano', '_l', tmp.lyear_str]), grid.tlong, grid.tlat);
                fig_cfg.fig_name=['l',tmp.lyear_str, ', corr_assm_ano_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
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
                title(cb,'R','fontsize',12);
            
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
                dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_assm_map', filesep, 'model', filesep, 'mask_n'];
                if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
                cfg.figname=[dirs.figdir, filesep, 'corr_assm_ano_mask_n_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
                print(fig_h, cfg.figname, '-dpng');
                RemoveWhiteSpace([], 'file', cfg.figname);
                close all;
            end
        end



%% model-lens2 & assm corr, model-lens2 & assm-lens2 map --------------------------------------
        if fig_flags(13)==1
        
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
            [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
        
            tmp.X=grid.tlong([end, 1:end],:);
            tmp.Y=grid.tlat([end, 1:end],:);
            tmp.C=data2.([tmp.varname, '_corr_assm_int_ano', '_l', tmp.lyear_str]).* mask.nlim.comb;
            tmp.C=tmp.C([end, 1:end],:);

            [tmp.mean_corr, tmp.err] = ...
                Func_0011_get_area_weighted_mean(data2.([tmp.varname, '_corr_assm_int_ano', '_l', tmp.lyear_str]), grid.tlong, grid.tlat);
            fig_cfg.fig_name=['l',tmp.lyear_str, ', corr_assm_int_ano_, ', tmp. varname, ', ', num2str(round(tmp.mean_corr,2))];
            fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
                'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
            %% map setting
            ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
                'fontname','freeserif'); 
        
            axis off; 
            hold on;
            setm(ax_m,'origin',[0,200],'MapLatLimit',fig_cfg.y_lim);
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
                tmp.tmppos=label_y(lyi,1).Position;
                tmp.tmppos(1)=-fig_cfg.ax_size(3)+2.6; % x position correction
                label_y(lyi,1).Position=tmp.tmppos;
            end
        
            %% save
            dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, tmp.varname, '_corr_assm_int_map', filesep, 'model', filesep, 'mask_n'];
            if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
            cfg.figname=[dirs.figdir, filesep, 'corr_assm_int_ano_mask_n_map_', tmp.varname, '_l', tmp.lyear_str, 'y.tif'];
        
            print(fig_h, cfg.figname, '-dpng');
            RemoveWhiteSpace([], 'file', cfg.figname);
            close all;
        
        end



    end
end

tmp.varname='photoC_TOT_zint_100m';
dirs.hcstmatroot=[tmp.kimyypath, '/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', tmp.varname];
fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
                '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
                '_l', tmp.lyear_str, 'y.mat'];
load(fig_cfg.mat_name);
npp_assm=data.photoC_TOT_zint_100m_assm;

tmp.varname='Jint_100m_NO3';
dirs.hcstmatroot=[tmp.kimyypath, '/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', tmp.varname];
fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
                '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
                '_l', tmp.lyear_str, 'y.mat'];
load(fig_cfg.mat_name);

NO3sink_assm=data.Jint_100m_NO3_assm;

for loni=1:size(ensmean_npp.(tmp.phyto),1)
    for lati=1:size(ensmean_npp.(tmp.phyto),2)
        tmp.corr=corrcoef(npp_assm(loni,lati,:), NO3sink_assm(loni,lati,:), 'Rows', 'complete');
        corr_sink(loni,lati)=tmp.corr(1,2);
    end
end

pcolor(corr_sink'); shading flat; colorbar; caxis([-1 1]); colormap(fig_cfg.c_map)


tmp.varname='subt_tend_zint_100m_NO3_Jint_100m_NO3';
dirs.hcstmatroot=[tmp.kimyypath, '/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', tmp.varname];
fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
                '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
                '_l', tmp.lyear_str, 'y.mat'];
load(fig_cfg.mat_name);
subt_assm=data.subt_tend_zint_100m_NO3_Jint_100m_NO3_assm;

for loni=1:size(ensmean_npp.(tmp.phyto),1)
    for lati=1:size(ensmean_npp.(tmp.phyto),2)
        tmp.corr=corrcoef(subt_assm(loni,lati,:), NO3sink_assm(loni,lati,:), 'Rows', 'complete');
        corr_circ(loni,lati)=tmp.corr(1,2);
    end
end

pcolor(corr_circ'); shading flat; colorbar; caxis([-1 1]); colormap(fig_cfg.c_map)





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
%             obsname_simple='GPCP';
            obsname_simple='GPCC';
        case 'PSL'
            obsname_simple='ERA5';
        case 'SSH'
            obsname_simple='CMEMS';
        case 'TS'
%             obsname_simple='HadCRUT5';
            obsname_simple='ERA5';
        case 'sumChl'
            obsname_simple='OC_CCI';
        case 'RAIN'
            obsname_simple='GPCC';
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
        case 'lnd'
            obsname_simple='.clm2.h0.';
        case 'ice'
            obsname_simple='.cice.h.';
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