% %  Created 23-May-2023 by Yong-Yub Kim
% %  Modified 25-Nov-2023 by Yong-Yub Kim
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

cfg.vlayer=1:10; % 10layer. don't put more than 15ㅅ
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

dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', cfg.var];
dirs.figroot=[tmp.rootpath, '/kimyy/Figure/CESM2/ESP/HCST_EXP/archive/', cfg.comp,'/', 'NPP_rec'];

cfg.iyears=cfg.obs_iyears;
cfg.gnm='f09_g17';
cfg.proj_year=5;

cfg.len_t_y = length(cfg.iyears);
cfg.casename_m = ['ens_all'];

% [tmp.error_status, tmp.value]=system(['ls ', dirs.datadir, '/*once*']);  % b.e21.BHISTsmbb.f09_g17.assm.oras4_ba-10p1.pop.h.once.nc
tmp.gridname = [tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/archive_transfer/', cfg.comp, '/grid.nc'];
tmp.maskname = [tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/archive_transfer/ocn/RECCAP2_region_masks_all_v20210412_POP2_grid.nc'];

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

phytoset={'diat', 'diaz', 'sp'};
% phytoset={'sp'};

%% loop start
for pind=1:length(phytoset)
%     clear data_mask_nut_dominant_ind data_mask_nut_dominant_val data_min_lim_phyto data_min_lim_ind_pytho ...
%         data_NPP_rec data_NPP_rec_corr data_NPP_rec_Fe data_NPP_rec_Fe_corr data_NPP_rec_L data_NPP_rec_L_corr ...
%         data_NPP_rec_N data_NPP_rec_N_corr data_NPP_rec_P data_NPP_rec_P_corr ...
%         data_NPP_rec_SiO3 data_NPP_rec_SiO3_corr data_NPP_rec_T data_NPP_rec_T_corr ...
%         data_NPP_rec_V data_NPP_rec_V_corr
    tmp.phyto=phytoset{pind};

    %% mu_ref
    if pind==2
        mu_ref(pind) = 2.5;
    else
        mu_ref(pind) = 5;
    end
    
    %% year loop start
    for lyear = 0:0  % 0:4

        tmp.lyear_str=num2str(lyear, '%02i');
        
        clear data_mask_nut_dominant data_min_lim_phyto data_NPP_rec data_NPP_rec_corr data_lim_set_phyto data_min_freq

        matname=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/ocn/photoC_TOT_zint_100m/',tmp.phyto, '_nut_NPP_l', tmp.lyear_str, '.mat'];
        %% data read
        if ~exist(matname)
%         if lyear<99
            disp(matname)
            tic;
            dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', cfg.var];
            fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', cfg.var, ...
                        '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
                        '_l', tmp.lyear_str, 'y.mat'];
            load(fig_cfg.mat_name, 'data');
            data_NPP=data;
            
            tmp.varname='TEMP'; 
            dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', tmp.varname];
            fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
                        '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
                        '_l', tmp.lyear_str, 'y.mat'];
            load(fig_cfg.mat_name, 'data');
            data_TEMP=data;

            %% photoC_phyto
            tmp.varname=['photoC_',tmp.phyto,'_zint_100m'];
            dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', tmp.varname];
            fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
                        '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
                        '_l', tmp.lyear_str, 'y.mat'];
            load(fig_cfg.mat_name, 'data');
            data_NPP_phyto.(tmp.phyto)=data;
            
            %% biomass
            tmp.varname=[tmp.phyto,'C'];
            dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', tmp.varname];
            fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
                        '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
                        '_l', tmp.lyear_str, 'y.mat'];
            load(fig_cfg.mat_name, 'data');
            data_C_phyto.(tmp.phyto)=data;
        
            %% Fe
            tmp.varname=[tmp.phyto,'_Fe_lim_Cweight_avg_100m'];
            dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', tmp.varname];
            fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
                        '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
                        '_l', tmp.lyear_str, 'y.mat'];
            load(fig_cfg.mat_name, 'data');
            data_Fe_lim_phyto.(tmp.phyto)=data;
        
            %% N
            if strcmp(tmp.phyto, 'diaz')~=1
                tmp.varname=[tmp.phyto,'_N_lim_Cweight_avg_100m'];
                dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', tmp.varname];
                fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
                            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
                            '_l', tmp.lyear_str, 'y.mat'];
                load(fig_cfg.mat_name, 'data');
                data_N_lim_phyto.(tmp.phyto)=data;
            end
        
            %% SiO3
            if strcmp(tmp.phyto, 'diat')==1
                tmp.varname=[tmp.phyto,'_SiO3_lim_Cweight_avg_100m'];
                dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', tmp.varname];
                fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
                            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
                            '_l', tmp.lyear_str, 'y.mat'];
                load(fig_cfg.mat_name, 'data');
                data_SiO3_lim_phyto.(tmp.phyto)=data;
            end
            
            %% P
            tmp.varname=[tmp.phyto,'_P_lim_Cweight_avg_100m'];
            dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', tmp.varname];
            fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
                        '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
                        '_l', tmp.lyear_str, 'y.mat'];
            load(fig_cfg.mat_name, 'data');
            data_P_lim_phyto.(tmp.phyto)=data;
        
            %% light 
            tmp.varname=[tmp.phyto,'_light_lim_Cweight_avg_100m'];
            dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', tmp.varname];
            fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
                        '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
                        '_l', tmp.lyear_str, 'y.mat'];
            load(fig_cfg.mat_name, 'data');
            data_light_lim_phyto.(tmp.phyto)=data;
            
        
            %% initialization for minimum limitation set
            data_lim_set_phyto.(tmp.phyto).assm = ...
                NaN(size(grid.tlong,1), size(grid.tlong,2), length(cfg.iyears), 4);
            data_lim_set_phyto.(tmp.phyto).(['lens2_l',tmp.lyear_str]) = ...
                NaN(size(grid.tlong,1), size(grid.tlong,2), length(cfg.iyears), 4);
            data_lim_set_phyto.(tmp.phyto).(['model_l',tmp.lyear_str]) = ...
                NaN(size(grid.tlong,1), size(grid.tlong,2), length(cfg.iyears), 4);
            
            %% Fe
            data_lim_set_phyto.(tmp.phyto).assm(:,:,:,1)=...
                data_Fe_lim_phyto.(tmp.phyto).([tmp.phyto,'_Fe_lim_Cweight_avg_100m_assm']);
            data_lim_set_phyto.(tmp.phyto).(['lens2_l',tmp.lyear_str])(:,:,:,1)=...
                data_Fe_lim_phyto.(tmp.phyto).([tmp.phyto,'_Fe_lim_Cweight_avg_100m_lens2_l', tmp.lyear_str]);
            data_lim_set_phyto.(tmp.phyto).(['model_l',tmp.lyear_str])(:,:,:,1)=...
                data_Fe_lim_phyto.(tmp.phyto).([tmp.phyto,'_Fe_lim_Cweight_avg_100m_model_l', tmp.lyear_str]);
            
            %% N
            if strcmp(tmp.phyto, 'diaz')~=1
                data_lim_set_phyto.(tmp.phyto).assm(:,:,:,2)=...
                    data_N_lim_phyto.(tmp.phyto).([tmp.phyto,'_N_lim_Cweight_avg_100m_assm']);
                data_lim_set_phyto.(tmp.phyto).(['lens2_l',tmp.lyear_str])(:,:,:,2)=...
                    data_N_lim_phyto.(tmp.phyto).([tmp.phyto,'_N_lim_Cweight_avg_100m_lens2_l', tmp.lyear_str]);
                data_lim_set_phyto.(tmp.phyto).(['model_l',tmp.lyear_str])(:,:,:,2)=...
                    data_N_lim_phyto.(tmp.phyto).([tmp.phyto,'_N_lim_Cweight_avg_100m_model_l', tmp.lyear_str]);
            else
                data_lim_set_phyto.(tmp.phyto).assm(:,:,:,2) = NaN;
                data_lim_set_phyto.(tmp.phyto).(['lens2_l',tmp.lyear_str])(:,:,:,2) = NaN;
                data_lim_set_phyto.(tmp.phyto).(['model_l',tmp.lyear_str])(:,:,:,2) = NaN;
            end
            %% SiO3
            if strcmp(tmp.phyto, 'diat')==1
                data_lim_set_phyto.(tmp.phyto).assm(:,:,:,3)=...
                    data_SiO3_lim_phyto.(tmp.phyto).([tmp.phyto,'_SiO3_lim_Cweight_avg_100m_assm']);
                data_lim_set_phyto.(tmp.phyto).(['lens2_l',tmp.lyear_str])(:,:,:,3)=...
                    data_SiO3_lim_phyto.(tmp.phyto).([tmp.phyto,'_SiO3_lim_Cweight_avg_100m_lens2_l', tmp.lyear_str]);
                data_lim_set_phyto.(tmp.phyto).(['model_l',tmp.lyear_str])(:,:,:,3)=...
                    data_SiO3_lim_phyto.(tmp.phyto).([tmp.phyto,'_SiO3_lim_Cweight_avg_100m_model_l', tmp.lyear_str]);
            else
                data_lim_set_phyto.(tmp.phyto).assm(:,:,:,3) = NaN;
                data_lim_set_phyto.(tmp.phyto).(['lens2_l',tmp.lyear_str])(:,:,:,3) = NaN;
                data_lim_set_phyto.(tmp.phyto).(['model_l',tmp.lyear_str])(:,:,:,3) = NaN;
            end
            %% P
            data_lim_set_phyto.(tmp.phyto).assm(:,:,:,4)=...
                data_P_lim_phyto.(tmp.phyto).([tmp.phyto,'_P_lim_Cweight_avg_100m_assm']);
            data_lim_set_phyto.(tmp.phyto).(['lens2_l',tmp.lyear_str])(:,:,:,4)=...
                data_P_lim_phyto.(tmp.phyto).([tmp.phyto,'_P_lim_Cweight_avg_100m_lens2_l', tmp.lyear_str]);
            data_lim_set_phyto.(tmp.phyto).(['model_l',tmp.lyear_str])(:,:,:,4)=...
                data_P_lim_phyto.(tmp.phyto).([tmp.phyto,'_P_lim_Cweight_avg_100m_model_l', tmp.lyear_str]);
    
            
    
            %% get V(minimum nutirient limitation)
            [data_min_lim_phyto.val.(tmp.phyto).assm, data_min_lim_phyto.ind.(tmp.phyto).assm]= ...
                min(data_lim_set_phyto.(tmp.phyto).assm,[], 4);
            [data_min_lim_phyto.val.(tmp.phyto).(['lens2_l',tmp.lyear_str]), ...
                data_min_lim_phyto.ind.(tmp.phyto).(['lens2_l',tmp.lyear_str])]= ...
                min(data_lim_set_phyto.(tmp.phyto).(['lens2_l',tmp.lyear_str]),[], 4);
            [data_min_lim_phyto.val.(tmp.phyto).(['model_l',tmp.lyear_str]), ...
                data_min_lim_phyto.ind.(tmp.phyto).(['model_l',tmp.lyear_str])]= ...
                min(data_lim_set_phyto.(tmp.phyto).(['model_l',tmp.lyear_str]),[], 4);
            
    
            % % pcolor(lim_set_assm_ind(:,:,5)'); shading flat; colorbar; colormap(jet(4))
            % % pcolor(mean(lim_set_assm_ind,3)'.*grid.mask_ocn'); shading flat; colorbar; colormap(jet(4)); caxis([0.5 4.5]);
            % % pcolor(std(lim_set_assm_ind,[],3)'.*grid.mask_ocn'); shading flat; colorbar; colormap(jet(256));
            

            %% reconstructed NPP
            data_NPP_rec.ALL.(tmp.phyto).assm = ...
                data_light_lim_phyto.(tmp.phyto).([tmp.phyto, '_light_lim_Cweight_avg_100m_assm']) .* ...
                data_min_lim_phyto.val.(tmp.phyto).assm .* ...
                1.7.^((data_TEMP.TEMP_assm-30.0)./10.0) .* ...
                mu_ref(pind) .* ...
                data_C_phyto.(tmp.phyto).([tmp.phyto,'C_assm']) ...
                ./ 86400.0 .* 10000.0;
            data_NPP_rec.ALL.(tmp.phyto).(['lens2_l',tmp.lyear_str]) = ...
                data_light_lim_phyto.(tmp.phyto).([tmp.phyto, '_light_lim_Cweight_avg_100m_lens2_l', tmp.lyear_str]) .* ...
                data_min_lim_phyto.val.(tmp.phyto).(['lens2_l',tmp.lyear_str]) .* ...
                1.7.^((data_TEMP.(['TEMP_','lens2_l',tmp.lyear_str])-30.0)./10.0) .* ...
                mu_ref(pind) .* ...
                data_C_phyto.(tmp.phyto).([tmp.phyto,'C_lens2_l', tmp.lyear_str]) ...
                ./ 86400.0 .* 10000.0;
            data_NPP_rec.ALL.(tmp.phyto).(['model_l',tmp.lyear_str]) = ...
                data_light_lim_phyto.(tmp.phyto).([tmp.phyto, '_light_lim_Cweight_avg_100m_model_l', tmp.lyear_str]) .* ...
                data_min_lim_phyto.val.(tmp.phyto).(['model_l',tmp.lyear_str]) .* ...
                1.7.^((data_TEMP.(['TEMP_','model_l',tmp.lyear_str])-30.0)./10.0) .* ...
                mu_ref(pind) .* ...
                data_C_phyto.(tmp.phyto).([tmp.phyto,'C_model_l', tmp.lyear_str]) ...
                ./ 86400.0 .* 10000.0;
            
    
            %% reconstructed NPP: T driven
            data_NPP_rec.T.(tmp.phyto).assm = ...
                mean(data_light_lim_phyto.(tmp.phyto).([tmp.phyto, '_light_lim_Cweight_avg_100m_assm']),3,'omitnan') .* ...
                mean(data_min_lim_phyto.val.(tmp.phyto).assm,3,'omitnan') .* ...
                1.7.^((data_TEMP.TEMP_assm-30.0)./10.0) .* ...
                mu_ref(pind) .* ...
                mean(data_C_phyto.(tmp.phyto).([tmp.phyto,'C_assm']),3,'omitnan') ...
                ./ 86400.0 .* 10000.0;
            data_NPP_rec.T.(tmp.phyto).(['lens2_l',tmp.lyear_str]) = ...
                mean(data_light_lim_phyto.(tmp.phyto).([tmp.phyto, '_light_lim_Cweight_avg_100m_lens2_l', tmp.lyear_str]),3,'omitnan') .* ...
                mean(data_min_lim_phyto.val.(tmp.phyto).(['lens2_l',tmp.lyear_str]),3,'omitnan') .* ...
                1.7.^((data_TEMP.(['TEMP_','lens2_l',tmp.lyear_str])-30.0)./10.0) .* ...
                mu_ref(pind) .* ...
                mean(data_C_phyto.(tmp.phyto).([tmp.phyto,'C_lens2_l', tmp.lyear_str]),3,'omitnan') ...
                ./ 86400.0 .* 10000.0;
            data_NPP_rec.T.(tmp.phyto).(['model_l',tmp.lyear_str]) = ...
                mean(data_light_lim_phyto.(tmp.phyto).([tmp.phyto, '_light_lim_Cweight_avg_100m_model_l', tmp.lyear_str]),3,'omitnan') .* ...
                mean(data_min_lim_phyto.val.(tmp.phyto).(['model_l',tmp.lyear_str]),3,'omitnan') .* ...
                1.7.^((data_TEMP.(['TEMP_','model_l',tmp.lyear_str])-30.0)./10.0) .* ...
                mu_ref(pind) .* ...
                mean(data_C_phyto.(tmp.phyto).([tmp.phyto,'C_model_l', tmp.lyear_str]),3,'omitnan') ...
                ./ 86400.0 .* 10000.0;

            
    
            %% reconstructed NPP: L driven
            data_NPP_rec.L.(tmp.phyto).assm = ...
                data_light_lim_phyto.(tmp.phyto).([tmp.phyto, '_light_lim_Cweight_avg_100m_assm']) .* ...
                mean(data_min_lim_phyto.val.(tmp.phyto).assm,3,'omitnan') .* ...
                mean(1.7.^((data_TEMP.TEMP_assm-30.0)./10.0),3,'omitnan') .* ...
                mu_ref(pind) .* ...
                mean(data_C_phyto.(tmp.phyto).([tmp.phyto,'C_assm']),3,'omitnan') ...
                ./ 86400.0 .* 10000.0;
            data_NPP_rec.L.(tmp.phyto).(['lens2_l',tmp.lyear_str]) = ...
                data_light_lim_phyto.(tmp.phyto).([tmp.phyto, '_light_lim_Cweight_avg_100m_lens2_l', tmp.lyear_str]) .* ...
                mean(data_min_lim_phyto.val.(tmp.phyto).(['lens2_l',tmp.lyear_str]),3,'omitnan') .* ...
                mean(1.7.^((data_TEMP.(['TEMP_','lens2_l',tmp.lyear_str])-30.0)./10.0),3, 'omitnan') .* ...
                mu_ref(pind) .* ...
                mean(data_C_phyto.(tmp.phyto).([tmp.phyto,'C_lens2_l', tmp.lyear_str]),3,'omitnan') ...
                ./ 86400.0 .* 10000.0;
            data_NPP_rec.L.(tmp.phyto).(['model_l',tmp.lyear_str]) = ...
                data_light_lim_phyto.(tmp.phyto).([tmp.phyto, '_light_lim_Cweight_avg_100m_model_l', tmp.lyear_str]) .* ...
                mean(data_min_lim_phyto.val.(tmp.phyto).(['model_l',tmp.lyear_str]),3,'omitnan') .* ...
                mean(1.7.^((data_TEMP.(['TEMP_','model_l',tmp.lyear_str])-30.0)./10.0),3,'omitnan') .* ...
                mu_ref(pind) .* ...
                mean(data_C_phyto.(tmp.phyto).([tmp.phyto,'C_model_l', tmp.lyear_str]),3,'omitnan') ...
                ./ 86400.0 .* 10000.0;
            
            
            %% reconstructed NPP: V driven (nutrient limiting factor)
            data_NPP_rec.V.(tmp.phyto).assm = ...
                mean(data_light_lim_phyto.(tmp.phyto).([tmp.phyto, '_light_lim_Cweight_avg_100m_assm']),3,'omitnan') .* ...
                data_min_lim_phyto.val.(tmp.phyto).assm .* ...
                mean(1.7.^((data_TEMP.TEMP_assm-30.0)./10.0),3,'omitnan') .* ...
                mu_ref(pind) .* ...
                mean(data_C_phyto.(tmp.phyto).([tmp.phyto,'C_assm']),3,'omitnan') ...
                ./ 86400.0 .* 10000.0;
            data_NPP_rec.V.(tmp.phyto).(['lens2_l',tmp.lyear_str]) = ...
                mean(data_light_lim_phyto.(tmp.phyto).([tmp.phyto, '_light_lim_Cweight_avg_100m_lens2_l', tmp.lyear_str]),3,'omitnan') .* ...
                data_min_lim_phyto.val.(tmp.phyto).(['lens2_l',tmp.lyear_str]) .* ...
                mean(1.7.^((data_TEMP.(['TEMP_','lens2_l',tmp.lyear_str])-30.0)./10.0),3,'omitnan') .* ...
                mu_ref(pind) .* ...
                mean(data_C_phyto.(tmp.phyto).([tmp.phyto,'C_lens2_l', tmp.lyear_str]),3,'omitnan') ...
                ./ 86400.0 .* 10000.0;
            data_NPP_rec.V.(tmp.phyto).(['model_l',tmp.lyear_str]) = ...
                mean(data_light_lim_phyto.(tmp.phyto).([tmp.phyto, '_light_lim_Cweight_avg_100m_model_l', tmp.lyear_str]),3,'omitnan') .* ...
                data_min_lim_phyto.val.(tmp.phyto).(['model_l',tmp.lyear_str]) .* ...
                mean(1.7.^((data_TEMP.(['TEMP_','model_l',tmp.lyear_str])-30.0)./10.0),3,'omitnan') .* ...
                mu_ref(pind) .* ...
                mean(data_C_phyto.(tmp.phyto).([tmp.phyto,'C_model_l', tmp.lyear_str]),3,'omitnan') ...
                ./ 86400.0 .* 10000.0;
            
    
            %% reconstructed NPP: Fe driven (nutrient limiting factor)
            data_NPP_rec.Fe.(tmp.phyto).assm = ...
                mean(data_light_lim_phyto.(tmp.phyto).([tmp.phyto, '_light_lim_Cweight_avg_100m_assm']),3,'omitnan') .* ...
                data_lim_set_phyto.(tmp.phyto).assm(:,:,:,1) .* ...
                mean(1.7.^((data_TEMP.TEMP_assm-30.0)./10.0),3,'omitnan') .* ...
                mu_ref(pind) .* ...
                mean(data_C_phyto.(tmp.phyto).([tmp.phyto,'C_assm']),3,'omitnan') ...
                ./ 86400.0 .* 10000.0;
            data_NPP_rec.Fe.(tmp.phyto).(['lens2_l',tmp.lyear_str]) = ...
                mean(data_light_lim_phyto.(tmp.phyto).([tmp.phyto, '_light_lim_Cweight_avg_100m_lens2_l', tmp.lyear_str]),3,'omitnan') .* ...
                data_lim_set_phyto.(tmp.phyto).(['lens2_l',tmp.lyear_str])(:,:,:,1) .* ...
                mean(1.7.^((data_TEMP.(['TEMP_','lens2_l',tmp.lyear_str])-30.0)./10.0),3,'omitnan') .* ...
                mu_ref(pind) .* ...
                mean(data_C_phyto.(tmp.phyto).([tmp.phyto,'C_lens2_l', tmp.lyear_str]),3,'omitnan') ...
                ./ 86400.0 .* 10000.0;
            data_NPP_rec.Fe.(tmp.phyto).(['model_l',tmp.lyear_str]) = ...
                mean(data_light_lim_phyto.(tmp.phyto).([tmp.phyto, '_light_lim_Cweight_avg_100m_model_l', tmp.lyear_str]),3,'omitnan') .* ...
                data_lim_set_phyto.(tmp.phyto).(['model_l',tmp.lyear_str])(:,:,:,1) .* ...
                mean(1.7.^((data_TEMP.(['TEMP_','model_l',tmp.lyear_str])-30.0)./10.0),3,'omitnan') .* ...
                mu_ref(pind) .* ...
                mean(data_C_phyto.(tmp.phyto).([tmp.phyto,'C_model_l', tmp.lyear_str]),3,'omitnan') ...
                ./ 86400.0 .* 10000.0;
            
    
            %% reconstructed NPP: N driven (nutrient limiting factor)
            data_NPP_rec.N.(tmp.phyto).assm = ...
                mean(data_light_lim_phyto.(tmp.phyto).([tmp.phyto, '_light_lim_Cweight_avg_100m_assm']),3,'omitnan') .* ...
                data_lim_set_phyto.(tmp.phyto).assm(:,:,:,2) .* ...
                mean(1.7.^((data_TEMP.TEMP_assm-30.0)./10.0),3,'omitnan') .* ...
                mu_ref(pind) .* ...
                mean(data_C_phyto.(tmp.phyto).([tmp.phyto,'C_assm']),3,'omitnan') ...
                ./ 86400.0 .* 10000.0;
            data_NPP_rec.N.(tmp.phyto).(['lens2_l',tmp.lyear_str]) = ...
                mean(data_light_lim_phyto.(tmp.phyto).([tmp.phyto, '_light_lim_Cweight_avg_100m_lens2_l', tmp.lyear_str]),3,'omitnan') .* ...
                data_lim_set_phyto.(tmp.phyto).(['lens2_l',tmp.lyear_str])(:,:,:,2) .* ...
                mean(1.7.^((data_TEMP.(['TEMP_','lens2_l',tmp.lyear_str])-30.0)./10.0),3,'omitnan') .* ...
                mu_ref(pind) .* ...
                mean(data_C_phyto.(tmp.phyto).([tmp.phyto,'C_lens2_l', tmp.lyear_str]),3,'omitnan') ...
                ./ 86400.0 .* 10000.0;
            data_NPP_rec.N.(tmp.phyto).(['model_l',tmp.lyear_str]) = ...
                mean(data_light_lim_phyto.(tmp.phyto).([tmp.phyto, '_light_lim_Cweight_avg_100m_model_l', tmp.lyear_str]),3,'omitnan') .* ...
                data_lim_set_phyto.(tmp.phyto).(['model_l',tmp.lyear_str])(:,:,:,2) .* ...
                mean(1.7.^((data_TEMP.(['TEMP_','model_l',tmp.lyear_str])-30.0)./10.0),3,'omitnan') .* ...
                mu_ref(pind) .* ...
                mean(data_C_phyto.(tmp.phyto).([tmp.phyto,'C_model_l', tmp.lyear_str]),3,'omitnan') ...
                ./ 86400.0 .* 10000.0;
            
    
            %% reconstructed NPP: SiO3 driven (nutrient limiting factor)
            data_NPP_rec.SiO3.(tmp.phyto).assm = ...
                mean(data_light_lim_phyto.(tmp.phyto).([tmp.phyto, '_light_lim_Cweight_avg_100m_assm']),3,'omitnan') .* ...
                data_lim_set_phyto.(tmp.phyto).assm(:,:,:,3) .* ...
                mean(1.7.^((data_TEMP.TEMP_assm-30.0)./10.0),3,'omitnan') .* ...
                mu_ref(pind) .* ...
                mean(data_C_phyto.(tmp.phyto).([tmp.phyto,'C_assm']),3,'omitnan') ...
                ./ 86400.0 .* 10000.0;
            data_NPP_rec.SiO3.(tmp.phyto).(['lens2_l',tmp.lyear_str]) = ...
                mean(data_light_lim_phyto.(tmp.phyto).([tmp.phyto, '_light_lim_Cweight_avg_100m_lens2_l', tmp.lyear_str]),3,'omitnan') .* ...
                data_lim_set_phyto.(tmp.phyto).(['lens2_l',tmp.lyear_str])(:,:,:,3) .* ...
                mean(1.7.^((data_TEMP.(['TEMP_','lens2_l',tmp.lyear_str])-30.0)./10.0),3,'omitnan') .* ...
                mu_ref(pind) .* ...
                mean(data_C_phyto.(tmp.phyto).([tmp.phyto,'C_lens2_l', tmp.lyear_str]),3,'omitnan') ...
                ./ 86400.0 .* 10000.0;
            data_NPP_rec.SiO3.(tmp.phyto).(['model_l',tmp.lyear_str]) = ...
                mean(data_light_lim_phyto.(tmp.phyto).([tmp.phyto, '_light_lim_Cweight_avg_100m_model_l', tmp.lyear_str]),3,'omitnan') .* ...
                data_lim_set_phyto.(tmp.phyto).(['model_l',tmp.lyear_str])(:,:,:,3) .* ...
                mean(1.7.^((data_TEMP.(['TEMP_','model_l',tmp.lyear_str])-30.0)./10.0),3,'omitnan') .* ...
                mu_ref(pind) .* ...
                mean(data_C_phyto.(tmp.phyto).([tmp.phyto,'C_model_l', tmp.lyear_str]),3,'omitnan') ...
                ./ 86400.0 .* 10000.0;
            
    
            %% reconstructed NPP: P driven (nutrient limiting factor)
            data_NPP_rec.P.(tmp.phyto).assm = ...
                mean(data_light_lim_phyto.(tmp.phyto).([tmp.phyto, '_light_lim_Cweight_avg_100m_assm']),3,'omitnan') .* ...
                data_lim_set_phyto.(tmp.phyto).assm(:,:,:,4) .* ...
                mean(1.7.^((data_TEMP.TEMP_assm-30.0)./10.0),3,'omitnan') .* ...
                mu_ref(pind) .* ...
                mean(data_C_phyto.(tmp.phyto).([tmp.phyto,'C_assm']),3,'omitnan') ...
                ./ 86400.0 .* 10000.0;
            data_NPP_rec.P.(tmp.phyto).(['lens2_l',tmp.lyear_str]) = ...
                mean(data_light_lim_phyto.(tmp.phyto).([tmp.phyto, '_light_lim_Cweight_avg_100m_lens2_l', tmp.lyear_str]),3,'omitnan') .* ...
                data_lim_set_phyto.(tmp.phyto).(['lens2_l',tmp.lyear_str])(:,:,:,4) .* ...
                mean(1.7.^((data_TEMP.(['TEMP_','lens2_l',tmp.lyear_str])-30.0)./10.0),3,'omitnan') .* ...
                mu_ref(pind) .* ...
                mean(data_C_phyto.(tmp.phyto).([tmp.phyto,'C_lens2_l', tmp.lyear_str]),3,'omitnan') ...
                ./ 86400.0 .* 10000.0;
            data_NPP_rec.P.(tmp.phyto).(['model_l',tmp.lyear_str]) = ...
                mean(data_light_lim_phyto.(tmp.phyto).([tmp.phyto, '_light_lim_Cweight_avg_100m_model_l', tmp.lyear_str]),3,'omitnan') .* ...
                data_lim_set_phyto.(tmp.phyto).(['model_l',tmp.lyear_str])(:,:,:,4) .* ...
                mean(1.7.^((data_TEMP.(['TEMP_','model_l',tmp.lyear_str])-30.0)./10.0),3,'omitnan') .* ...
                mu_ref(pind) .* ...
                mean(data_C_phyto.(tmp.phyto).([tmp.phyto,'C_model_l', tmp.lyear_str]),3,'omitnan') ...
                ./ 86400.0 .* 10000.0;
            
    
            %% corr NPP_rec, NPP
    
            %% corr initialization
    
            strset={'assm', ['lens2_l',tmp.lyear_str], ['model_l',tmp.lyear_str], ...
                ['PP_model_l',tmp.lyear_str], ['PP_hcst_int_l',tmp.lyear_str], ['PP_lens2_l',tmp.lyear_str]};
            for stri=1:length(strset)
                tmpstr=strset{stri};
                data_NPP_rec_corr.ALL.(tmp.phyto).(tmpstr)=NaN(size(grid.tlong));
                data_NPP_rec_corr.T.(tmp.phyto).(tmpstr)=NaN(size(grid.tlong));
                data_NPP_rec_corr.L.(tmp.phyto).(tmpstr)=NaN(size(grid.tlong));
                data_NPP_rec_corr.V.(tmp.phyto).(tmpstr)=NaN(size(grid.tlong));
                data_NPP_rec_corr.Fe.(tmp.phyto).(tmpstr)=NaN(size(grid.tlong));
                data_NPP_rec_corr.N.(tmp.phyto).(tmpstr)=NaN(size(grid.tlong));
                data_NPP_rec_corr.SiO3.(tmp.phyto).(tmpstr)=NaN(size(grid.tlong));
                data_NPP_rec_corr.P.(tmp.phyto).(tmpstr)=NaN(size(grid.tlong));
            end
    
            for loni=1:size(grid.tlong,1)
                for lati=1:size(grid.tlat,2)
                    if(isfinite(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,1)))
                        %% rec validation (raw <-> reconstructed)
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.ALL.(tmp.phyto).assm(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.ALL.(tmp.phyto).assm(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_lens2_l', tmp.lyear_str])(loni,lati,:), ...
                            data_NPP_rec.ALL.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.ALL.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_model_l', tmp.lyear_str])(loni,lati,:), ...
                            data_NPP_rec.ALL.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.ALL.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
    
                        %% rec Potential Predictability (assm raw <-> model reconstructed)
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.ALL.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.ALL.(tmp.phyto).(['PP_model_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.ALL.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati,:) - ...
                            data_NPP_rec.ALL.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.ALL.(tmp.phyto).(['PP_hcst_int_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.ALL.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.ALL.(tmp.phyto).(['PP_lens2_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
    
                        %% rec_T validation (raw <-> reconstructed)
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.T.(tmp.phyto).assm(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.T.(tmp.phyto).assm(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_lens2_l', tmp.lyear_str])(loni,lati,:), ...
                            data_NPP_rec.T.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.T.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_model_l', tmp.lyear_str])(loni,lati,:), ...
                            data_NPP_rec.T.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.T.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
                        %% rec_T PP
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.T.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.T.(tmp.phyto).(['PP_model_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.T.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati,:) - ...
                            data_NPP_rec.T.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.T.(tmp.phyto).(['PP_hcst_int_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.T.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.T.(tmp.phyto).(['PP_lens2_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
    
                        %% rec_L validation (raw <-> reconstructed)
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.L.(tmp.phyto).assm(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.L.(tmp.phyto).assm(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_lens2_l', tmp.lyear_str])(loni,lati,:), ...
                            data_NPP_rec.L.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.L.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_model_l', tmp.lyear_str])(loni,lati,:), ...
                            data_NPP_rec.L.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.L.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
                        %% rec_L PP
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.L.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.L.(tmp.phyto).(['PP_model_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.L.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati,:) - ...
                            data_NPP_rec.L.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.L.(tmp.phyto).(['PP_hcst_int_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.L.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.L.(tmp.phyto).(['PP_lens2_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
    
                        %% rec_V validation (raw <-> reconstructed)
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.V.(tmp.phyto).assm(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.V.(tmp.phyto).assm(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_lens2_l', tmp.lyear_str])(loni,lati,:), ...
                            data_NPP_rec.V.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.V.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_model_l', tmp.lyear_str])(loni,lati,:), ...
                            data_NPP_rec.V.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.V.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
                        %% rec_V PP
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.V.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.V.(tmp.phyto).(['PP_model_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.V.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati,:) - ...
                            data_NPP_rec.V.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.V.(tmp.phyto).(['PP_hcst_int_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.V.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.V.(tmp.phyto).(['PP_lens2_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
    
                        %% rec_Fe validation (raw <-> reconstructed)
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.Fe.(tmp.phyto).assm(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.Fe.(tmp.phyto).assm(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_lens2_l', tmp.lyear_str])(loni,lati,:), ...
                            data_NPP_rec.Fe.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.Fe.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_model_l', tmp.lyear_str])(loni,lati,:), ...
                            data_NPP_rec.Fe.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.Fe.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
                        %% rec_Fe PP
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.Fe.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.Fe.(tmp.phyto).(['PP_model_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.Fe.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati,:) - ...
                            data_NPP_rec.Fe.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.Fe.(tmp.phyto).(['PP_hcst_int_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.Fe.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.Fe.(tmp.phyto).(['PP_lens2_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
    
                        %% rec_N validation (raw <-> reconstructed)
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.N.(tmp.phyto).assm(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.N.(tmp.phyto).assm(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_lens2_l', tmp.lyear_str])(loni,lati,:), ...
                            data_NPP_rec.N.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.N.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_model_l', tmp.lyear_str])(loni,lati,:), ...
                            data_NPP_rec.N.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.N.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
                        %% rec_N PP
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.N.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.N.(tmp.phyto).(['PP_model_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.N.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati,:) - ...
                            data_NPP_rec.N.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.N.(tmp.phyto).(['PP_hcst_int_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.N.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.N.(tmp.phyto).(['PP_lens2_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
    
                        %% rec_SiO3 validation (raw <-> reconstructed)
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.SiO3.(tmp.phyto).assm(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.SiO3.(tmp.phyto).assm(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_lens2_l', tmp.lyear_str])(loni,lati,:), ...
                            data_NPP_rec.SiO3.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.SiO3.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_model_l', tmp.lyear_str])(loni,lati,:), ...
                            data_NPP_rec.SiO3.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.SiO3.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
    
                        %% rec_SiO3 PP
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.SiO3.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.SiO3.(tmp.phyto).(['PP_model_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.SiO3.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati,:) - ...
                            data_NPP_rec.SiO3.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.SiO3.(tmp.phyto).(['PP_hcst_int_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.SiO3.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.SiO3.(tmp.phyto).(['PP_lens2_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
    
                        %% rec_P validation (raw <-> reconstructed)
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.P.(tmp.phyto).assm(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.P.(tmp.phyto).assm(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_lens2_l', tmp.lyear_str])(loni,lati,:), ...
                            data_NPP_rec.P.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.P.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_model_l', tmp.lyear_str])(loni,lati,:), ...
                            data_NPP_rec.P.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.P.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
                        %% rec_P PP 
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.P.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.P.(tmp.phyto).(['PP_model_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.P.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati,:) - ...
                            data_NPP_rec.P.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.P.(tmp.phyto).(['PP_hcst_int_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_NPP_phyto.(tmp.phyto).(['photoC_',tmp.phyto,'_zint_100m_assm'])(loni,lati,:), ...
                            data_NPP_rec.P.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati,:), 'rows','complete');
                        data_NPP_rec_corr.P.(tmp.phyto).(['PP_lens2_l',tmp.lyear_str])(loni,lati)=tmp_corr_r(1,2);
                    end
                end
            end
    
            %% contribution of each nutrient in V
            data_mask_nut_dominant.ind.(tmp.phyto).assm=NaN(size(grid.tlong,1), size(grid.tlat,2), 4);
            data_mask_nut_dominant.ind.(tmp.phyto).(['lens2_l',tmp.lyear_str])=NaN(size(grid.tlong,1), size(grid.tlat,2), 4);
            data_mask_nut_dominant.ind.(tmp.phyto).(['model_l',tmp.lyear_str])=NaN(size(grid.tlong,1), size(grid.tlat,2), 4);
            data_mask_nut_dominant.val.(tmp.phyto).assm=NaN(size(grid.tlong,1), size(grid.tlat,2), 4);
            data_mask_nut_dominant.val.(tmp.phyto).(['lens2_l',tmp.lyear_str])=NaN(size(grid.tlong,1), size(grid.tlat,2), 4);
            data_mask_nut_dominant.val.(tmp.phyto).(['model_l',tmp.lyear_str])=NaN(size(grid.tlong,1), size(grid.tlat,2), 4);

    %         [data_min_lim_phyto.val.(tmp.phyto).assm, data_min_lim_phyto.ind.(tmp.phyto).assm]= ...
    %             min(data_lim_set_phyto.(tmp.phyto).assm,[], 4);
    %         [data_min_lim_phyto.val.(tmp.phyto).(['lens2_l',tmp.lyear_str]), ...
    %             data_min_lim_phyto.ind.(tmp.phyto).(['lens2_l',tmp.lyear_str])]= ...
    %             min(data_lim_set_phyto.(tmp.phyto).(['lens2_l',tmp.lyear_str]),[], 4);
    %         [data_min_lim_phyto.val.(tmp.phyto).(['model_l',tmp.lyear_str]), ...
    %             data_min_lim_phyto.ind.(tmp.phyto).(['model_l',tmp.lyear_str])]= ...
    %             min(data_lim_set_phyto.(tmp.phyto).(['model_l',tmp.lyear_str]),[], 4);
    
            for loni=1:size(grid.tlong,1)
                for lati=1:size(grid.tlat,2)
                 %% Fe
                    data_min_freq.Fe_lim_phyto.(tmp.phyto).assm(loni,lati)= ...
                        length(find(data_min_lim_phyto.ind.(tmp.phyto).assm(loni,lati,:)==1))/length(cfg.iyears);
                    data_min_freq.Fe_lim_phyto.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati)= ...
                        length(find(data_min_lim_phyto.ind.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati,:)==1))/length(cfg.iyears);
                    data_min_freq.Fe_lim_phyto.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati)= ...
                        length(find(data_min_lim_phyto.ind.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati,:)==1))/length(cfg.iyears);
    
                %% N
                    data_min_freq.N_lim_phyto.(tmp.phyto).assm(loni,lati)= ...
                        length(find(data_min_lim_phyto.ind.(tmp.phyto).assm(loni,lati,:)==2))/length(cfg.iyears);
                    data_min_freq.N_lim_phyto.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati)= ...
                        length(find(data_min_lim_phyto.ind.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati,:)==2))/length(cfg.iyears);
                    data_min_freq.N_lim_phyto.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati)= ...
                        length(find(data_min_lim_phyto.ind.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati,:)==2))/length(cfg.iyears);
    
                %% SiO3
                    data_min_freq.SiO3_lim_phyto.(tmp.phyto).assm(loni,lati)= ...
                        length(find(data_min_lim_phyto.ind.(tmp.phyto).assm(loni,lati,:)==3))/length(cfg.iyears);
                    data_min_freq.SiO3_lim_phyto.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati)= ...
                        length(find(data_min_lim_phyto.ind.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati,:)==3))/length(cfg.iyears);
                    data_min_freq.SiO3_lim_phyto.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati)= ...
                        length(find(data_min_lim_phyto.ind.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati,:)==3))/length(cfg.iyears);
    
                %% P
                    data_min_freq.P_lim_phyto.(tmp.phyto).assm(loni,lati)= ...
                        length(find(data_min_lim_phyto.ind.(tmp.phyto).assm(loni,lati,:)==4))/length(cfg.iyears);
                    data_min_freq.P_lim_phyto.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati)= ...
                        length(find(data_min_lim_phyto.ind.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati,:)==4))/length(cfg.iyears);
                    data_min_freq.P_lim_phyto.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati)= ...
                        length(find(data_min_lim_phyto.ind.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati,:)==4))/length(cfg.iyears);
                    
                %% comparison
                    tmp.min_freq = [data_min_freq.Fe_lim_phyto.(tmp.phyto).assm(loni,lati);...
                        data_min_freq.N_lim_phyto.(tmp.phyto).assm(loni,lati); ...
                        data_min_freq.SiO3_lim_phyto.(tmp.phyto).assm(loni,lati); ...
                        data_min_freq.P_lim_phyto.(tmp.phyto).assm(loni,lati)];
                    [max_minfreq, maxind_minfreq]=max(tmp.min_freq);
                    data_mask_nut_dominant.val.(tmp.phyto).assm(loni,lati) = max_minfreq;
                    data_mask_nut_dominant.ind.(tmp.phyto).assm(loni,lati) = maxind_minfreq;

                    tmp.min_freq = [data_min_freq.Fe_lim_phyto.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati);...
                        data_min_freq.N_lim_phyto.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati); ...
                        data_min_freq.SiO3_lim_phyto.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati); ...
                        data_min_freq.P_lim_phyto.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati)];
                    [max_minfreq, maxind_minfreq]=max(tmp.min_freq);
                    data_mask_nut_dominant.val.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati) = max_minfreq;
                    data_mask_nut_dominant.ind.(tmp.phyto).(['lens2_l',tmp.lyear_str])(loni,lati) = maxind_minfreq;
    
                    tmp.min_freq = [data_min_freq.Fe_lim_phyto.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati);...
                        data_min_freq.N_lim_phyto.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati); ...
                        data_min_freq.SiO3_lim_phyto.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati); ...
                        data_min_freq.P_lim_phyto.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati)];
                    [max_minfreq, maxind_minfreq]=max(tmp.min_freq);
                    data_mask_nut_dominant.val.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati) = max_minfreq;
                    data_mask_nut_dominant.ind.(tmp.phyto).(['model_l',tmp.lyear_str])(loni,lati) = maxind_minfreq;
    
                end
            end
            save(matname, 'data_mask_nut_dominant', 'data_min_lim_phyto', 'data_NPP_rec', 'data_NPP_rec_corr', 'data_lim_set_phyto', 'data_min_freq')
            toc;
        else
            load(matname)
            disp(matname)
        end
    

% % % % % %% pictures - each phyto
% % % % % %%
% % % % % %%
% % % % % %%
% % % % % %%
        
% % % % %         cfg.vstrset={'assm', ['lens2_l',tmp.lyear_str], ['model_l',tmp.lyear_str]};
% % % % %         for stri=1:length(cfg.vstrset)
% % % % %             tmp.vstr=cfg.vstrset{stri};
% % % % % 
% % % % %             cfg.nuts={'Fe', 'N', 'SiO3', 'P'};
% % % % %             for ni=1:length(cfg.nuts)
% % % % %                 tmp.nut=cfg.nuts{ni};
% % % % %         
% % % % %         %% Each nutrient occupied rate of minimum map --------------------------------------
% % % % %                 fig_cfg.name_rgn = 'Glob';
% % % % %                 fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
% % % % %                 fig_cfg.x_lim = [-180 180];
% % % % %                 fig_cfg.y_lim = [-80 89];
% % % % %                 fig_cfg.fig_size = [0,0,6,3.5];
% % % % %                 fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
% % % % %                 fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
% % % % %                 fig_cfg.title_pos = [0.5,0.93];
% % % % %                 fig_cfg.p_lim =0.1;
% % % % %                 fig_cfg.c_lim = [0 1];
% % % % %                 [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('wr_10', tmp.dropboxpath);
% % % % %             
% % % % %                 tmp.X=grid.tlong([end, 1:end],:);
% % % % %                 tmp.Y=grid.tlat([end, 1:end],:);
% % % % %                 tmp.C=data_min_freq.([tmp.nut, '_lim_phyto']).(tmp.phyto).(tmp.vstr).*grid.mask_ocn;
% % % % %                 tmp.C=tmp.C([end, 1:end],:); 
% % % % %             
% % % % %                 fig_cfg.fig_name='data_min_freq';
% % % % %                 fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
% % % % %                     'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
% % % % %                 %% map setting
% % % % %                 ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
% % % % %                     'fontname','freeserif'); 
% % % % %             
% % % % %                 axis off; 
% % % % %                 hold on;
% % % % %                 setm(ax_m,'origin',[0,180],'MapLatLimit',fig_cfg.y_lim);
% % % % %                 set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
% % % % %                 text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
% % % % %                 'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
% % % % %                 'fontsize',14,'fontname','freeserif','interpreter','none')
% % % % %             
% % % % %                 %% caxis & colorbar
% % % % %                 caxis(ax_m, fig_cfg.c_lim); 
% % % % %                 colormap(fig_cfg.c_map);
% % % % %                 cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
% % % % %                 set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
% % % % %                 title(cb,'R','fontsize',12);
% % % % %             
% % % % %                 %% draw on ax_m
% % % % %                 h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
% % % % %                 shading flat;
% % % % %                 geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
% % % % %             
% % % % %             
% % % % %             %% frame and label setting
% % % % %                 setm(ax_m,'frame','on','FLineWidth',1);
% % % % %             
% % % % %                 label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
% % % % %                 label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
% % % % %                 mlabel; plabel;
% % % % %                 label_y=plabel; label_x=mlabel;
% % % % %                 for lxi=1:length(label_x)
% % % % %                     tmp.tmppos=label_x(lxi,1).Position;
% % % % %                     tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55;
% % % % %                     label_x(lxi,1).Position=tmp.tmppos;
% % % % %                     label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
% % % % %                 end
% % % % %                 for lyi=1:length(label_y)
% % % % %                     label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
% % % % %                 end
% % % % %             
% % % % %                 %% save
% % % % %                 dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, 'min_freq', filesep, tmp.nut, filesep, tmp.phyto];
% % % % %                 if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
% % % % %                 cfg.figname=[dirs.figdir, filesep, 'min_freq_map_', tmp.nut, '_', tmp.phyto, '_', tmp.vstr, '.tif'];
% % % % %                 print(fig_h, cfg.figname, '-dpng');
% % % % %                 disp(cfg.figname)
% % % % %                 RemoveWhiteSpace([], 'file', cfg.figname);
% % % % %                 close all;
% % % % %     
% % % % %     
% % % % %             end
% % % % % 
% % % % %             cfg.limv={'ALL', 'T', 'L', 'V', 'Fe', 'N', 'SiO3', 'P'};
% % % % %             for ni=1:length(cfg.limv)
% % % % %                 tmp.limv=cfg.limv{ni};
% % % % % 
% % % % %          %% correlation between NPP <-> NPP_rec in same dataset map --------------------------------------
% % % % %                 fig_cfg.c_lim = [-1 1];
% % % % %                 [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
% % % % %             
% % % % %                 tmp.X=grid.tlong([end, 1:end],:);
% % % % %                 tmp.Y=grid.tlat([end, 1:end],:);
% % % % % 
% % % % %                 %% data set
% % % % %                 tmp.C=data_NPP_rec_corr.(tmp.limv).(tmp.phyto).(tmp.vstr) .* grid.mask_ocn;
% % % % %                 switch tmp.limv
% % % % %                     case {'Fe', 'N', 'SiO3', 'P'}
% % % % %                     tmp.mask=data_min_freq.([tmp.limv, '_lim_phyto']).(tmp.phyto).(tmp.vstr).*grid.mask_ocn;
% % % % %                     tmp.mask(tmp.mask<0.8)=NaN;
% % % % %                     tmp.C=tmp.C.*tmp.mask;
% % % % %                 end
% % % % %                 tmp.C=tmp.C([end, 1:end],:); 
% % % % %             
% % % % %                 fig_cfg.fig_name='data_NPP_rec_corr_explained';
% % % % %                 fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
% % % % %                     'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
% % % % %                 %% map setting
% % % % %                 ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
% % % % %                     'fontname','freeserif'); 
% % % % %             
% % % % %                 axis off; 
% % % % %                 hold on;
% % % % %                 setm(ax_m,'origin',[0,180],'MapLatLimit',fig_cfg.y_lim);
% % % % %                 set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
% % % % %                 text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
% % % % %                 'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
% % % % %                 'fontsize',14,'fontname','freeserif','interpreter','none')
% % % % %             
% % % % %                 %% caxis & colorbar
% % % % %                 caxis(ax_m, fig_cfg.c_lim); 
% % % % %                 colormap(fig_cfg.c_map);
% % % % %                 cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
% % % % %                 set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
% % % % %                 title(cb,'R','fontsize',12);
% % % % %             
% % % % %                 %% draw on ax_m
% % % % %                 h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
% % % % %                 shading flat;
% % % % %                 geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
% % % % %             
% % % % %             
% % % % %                 %% frame and label setting
% % % % %                 setm(ax_m,'frame','on','FLineWidth',1);
% % % % %             
% % % % %                 label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
% % % % %                 label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
% % % % %                 mlabel; plabel;
% % % % %                 label_y=plabel; label_x=mlabel;
% % % % %                 for lxi=1:length(label_x)
% % % % %                     tmp.tmppos=label_x(lxi,1).Position;
% % % % %                     tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55;
% % % % %                     label_x(lxi,1).Position=tmp.tmppos;
% % % % %                     label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
% % % % %                 end
% % % % %                 for lyi=1:length(label_y)
% % % % %                     label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
% % % % %                 end
% % % % %             
% % % % %                 %% save
% % % % %                 dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, 'corr_explained', filesep, tmp.limv, filesep, tmp.phyto];
% % % % %                 if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
% % % % %                 cfg.figname=[dirs.figdir, filesep, 'corr_exp_map_', tmp.limv, '_', tmp.phyto, '_', tmp.vstr, '.tif'];
% % % % %                 print(fig_h, cfg.figname, '-dpng');
% % % % %                 disp(cfg.figname)
% % % % %                 RemoveWhiteSpace([], 'file', cfg.figname);
% % % % %                 close all;
% % % % % 
% % % % %             end
% % % % %         end
% % % % %         

% % % % % % data_NPP_rec.ALL.(tmp.phyto).assm

% % % % %         cfg.vstrset={['lens2_l',tmp.lyear_str], ['model_l',tmp.lyear_str], ['hcst_int_l',tmp.lyear_str]};
% % % % %         for stri=1:length(cfg.vstrset)
% % % % %             tmp.vstr=cfg.vstrset{stri};
% % % % % 
% % % % %             cfg.limv={'ALL', 'T', 'L', 'V', 'Fe', 'N', 'SiO3', 'P'};
% % % % %             for ni=1:length(cfg.limv)
% % % % %                 tmp.limv=cfg.limv{ni};
% % % % % 
% % % % %         %% Potential Predictability of NPP, but by NPP_rec map --------------------------------------
% % % % %                 fig_cfg.c_lim = [-1 1];
% % % % %                 [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
% % % % %             
% % % % %                 tmp.X=grid.tlong([end, 1:end],:);
% % % % %                 tmp.Y=grid.tlat([end, 1:end],:);
% % % % % 
% % % % %                 %% data set
% % % % %                 tmp.C=data_NPP_rec_corr.(tmp.limv).(tmp.phyto).(['PP_', tmp.vstr]) .* grid.mask_ocn;
% % % % %                 switch tmp.limv
% % % % %                     case {'Fe', 'N', 'SiO3', 'P'}
% % % % %                         if strcmp(tmp.vstr, ['hcst_int_l',tmp.lyear_str])==1
% % % % %                             tmp.mask=data_min_freq.([tmp.limv, '_lim_phyto']).(tmp.phyto).(['model_l',tmp.lyear_str]).*grid.mask_ocn;
% % % % %                         else
% % % % %                             tmp.mask=data_min_freq.([tmp.limv, '_lim_phyto']).(tmp.phyto).(tmp.vstr).*grid.mask_ocn;
% % % % %                         end
% % % % %                         tmp.mask(tmp.mask<0.8)=NaN;
% % % % %                         tmp.C=tmp.C.*tmp.mask;
% % % % %                 end
% % % % %                 tmp.C=tmp.C([end, 1:end],:); 
% % % % %             
% % % % %                 fig_cfg.fig_name='data_NPP_rec_PP';
% % % % %                 fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
% % % % %                     'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
% % % % %                 %% map setting
% % % % %                 ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
% % % % %                     'fontname','freeserif'); 
% % % % %             
% % % % %                 axis off; 
% % % % %                 hold on;
% % % % %                 setm(ax_m,'origin',[0,180],'MapLatLimit',fig_cfg.y_lim);
% % % % %                 set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
% % % % %                 text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
% % % % %                 'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
% % % % %                 'fontsize',14,'fontname','freeserif','interpreter','none')
% % % % %             
% % % % %                 %% caxis & colorbar
% % % % %                 caxis(ax_m, fig_cfg.c_lim); 
% % % % %                 colormap(fig_cfg.c_map);
% % % % %                 cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
% % % % %                 set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
% % % % %                 title(cb,'R','fontsize',12);
% % % % %             
% % % % %                 %% draw on ax_m
% % % % %                 h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
% % % % %                 shading flat;
% % % % %                 geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
% % % % %             
% % % % %             
% % % % %                 %% frame and label setting
% % % % %                 setm(ax_m,'frame','on','FLineWidth',1);
% % % % %             
% % % % %                 label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
% % % % %                 label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
% % % % %                 mlabel; plabel;
% % % % %                 label_y=plabel; label_x=mlabel;
% % % % %                 for lxi=1:length(label_x)
% % % % %                     tmp.tmppos=label_x(lxi,1).Position;
% % % % %                     tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55;
% % % % %                     label_x(lxi,1).Position=tmp.tmppos;
% % % % %                     label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
% % % % %                 end
% % % % %                 for lyi=1:length(label_y)
% % % % %                     label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
% % % % %                 end
% % % % %             
% % % % %                 %% save
% % % % %                 dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, 'corr_PP', filesep, tmp.limv, filesep, tmp.phyto];
% % % % %                 if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
% % % % %                 cfg.figname=[dirs.figdir, filesep, 'corr_PP_map_', tmp.limv, '_', tmp.phyto, '_', tmp.vstr, '.tif'];
% % % % %                 print(fig_h, cfg.figname, '-dpng');
% % % % %                 disp(cfg.figname)
% % % % %                 RemoveWhiteSpace([], 'file', cfg.figname);
% % % % %                 close all;
% % % % %             end
% % % % %         end

    end
end


% % % % % %% sum of phyto
% % % % % for lyear = 0:0  % 0:4
% % % % %     tic
% % % % %     tmp.lyear_str=num2str(lyear, '%02i');
% % % % %     cfg.vstrset={'assm', ['lens2_l',tmp.lyear_str], ['model_l',tmp.lyear_str]};
% % % % %     
% % % % %     for pind=1:length(phytoset)
% % % % %         tmp.phyto=phytoset{pind};
% % % % %         matname=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/ocn/photoC_TOT_zint_100m/',tmp.phyto, '_nut_NPP_l', tmp.lyear_str, '.mat'];
% % % % %         load(matname);
% % % % %         data_all.TOT.(tmp.phyto)=data_NPP_rec;
% % % % %         data_all_min_freq.TOT.(tmp.phyto)=data_min_freq;
% % % % % % 
% % % % % %         %% photoC_phyto
% % % % % %         tmp.varname=['photoC_',tmp.phyto,'_zint_100m'];
% % % % % %         dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', tmp.varname];
% % % % % %         fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
% % % % % %                     '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
% % % % % %                     '_l', tmp.lyear_str, 'y.mat'];
% % % % % %         load(fig_cfg.mat_name, 'data');
% % % % % %         data_all_NPP_phyto.(tmp.phyto)=data;
% % % % %     end
% % % % %     tmp.varname=['photoC_','TOT','_zint_100m'];
% % % % %     dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', tmp.varname];
% % % % %     fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
% % % % %                 '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
% % % % %                 '_l', tmp.lyear_str, 'y.mat'];
% % % % %     load(fig_cfg.mat_name, 'data');
% % % % %     data_all_NPP_phyto.TOT=data;
% % % % % 
% % % % %     for stri=1:length(cfg.vstrset)
% % % % %         tmp.vstr=cfg.vstrset{stri};
% % % % %         cfg.limv={'ALL', 'T', 'L', 'V', 'Fe', 'N', 'SiO3', 'P'};
% % % % %         for ni=1:length(cfg.limv)
% % % % %             tmp.limv=cfg.limv{ni};
% % % % %             switch tmp.limv
% % % % %                 case {'Fe', 'N', 'SiO3', 'P'}
% % % % %                     for pind=1:length(phytoset)
% % % % %                         tmp.phyto=phytoset{pind};
% % % % %                         tmp.mask=data_all_min_freq.TOT.(tmp.phyto).([tmp.limv, '_lim_phyto']).(tmp.phyto).(tmp.vstr).*grid.mask_ocn;
% % % % %                         tmp.mask(tmp.mask<0.8)=NaN;
% % % % %                         tmp.comb(:,:,:,pind)=data_all.TOT.(tmp.phyto).(tmp.limv).(tmp.phyto).(tmp.vstr).*tmp.mask;
% % % % %                     end
% % % % %                 otherwise
% % % % %                     tmp.comb(:,:,:,1)=data_all.TOT.diat.(tmp.limv).diat.(tmp.vstr);
% % % % %                     tmp.comb(:,:,:,2)=data_all.TOT.diaz.(tmp.limv).diaz.(tmp.vstr);
% % % % %                     tmp.comb(:,:,:,3)=data_all.TOT.sp.(tmp.limv).sp.(tmp.vstr);
% % % % %             end
% % % % %             
% % % % %             data_NPP_rec.(tmp.limv).TOT.(tmp.vstr)= sum(tmp.comb, 4, 'omitnan');
% % % % % %             data_NPP_rec.(tmp.limv).TOT.(tmp.vstr) = ...
% % % % % %                 data_all.TOT.diat.(tmp.limv).diat.(tmp.vstr)  ...
% % % % % %                 + data_all.TOT.diaz.(tmp.limv).diaz.(tmp.vstr) ...
% % % % % %                 + data_all.TOT.sp.(tmp.limv).sp.(tmp.vstr);
% % % % % 
% % % % %             data_NPP_rec_corr.(tmp.limv).TOT.(tmp.vstr)=NaN(size(grid.tlong));
% % % % %             data_NPP_rec_corr.(tmp.limv).TOT.(['PP_', tmp.vstr])=NaN(size(grid.tlong));
% % % % % 
% % % % %             for loni=1:size(grid.tlong,1)
% % % % %                 for lati=1:size(grid.tlat,2)
% % % % %                     if(isfinite(data_all_NPP_phyto.TOT.photoC_TOT_zint_100m_assm(loni,lati,1)))
% % % % %                         %% rec validation (raw <-> reconstructed)
% % % % %                         [tmp_corr_r,tmp_corr_p] = ...
% % % % %                             corrcoef(data_all_NPP_phyto.TOT.(['photoC_TOT_zint_100m_',tmp.vstr])(loni,lati,:), ...
% % % % %                             data_NPP_rec.(tmp.limv).TOT.(tmp.vstr)(loni,lati,:), 'rows','complete');
% % % % %                         data_NPP_rec_corr.(tmp.limv).TOT.(tmp.vstr)(loni,lati)=tmp_corr_r(1,2);
% % % % %     
% % % % %                         %% rec Potential Predictability (assm raw <-> model reconstructed)
% % % % %                         [tmp_corr_r,tmp_corr_p] = ...
% % % % %                             corrcoef(data_all_NPP_phyto.TOT.(['photoC_TOT_zint_100m_assm'])(loni,lati,:), ...
% % % % %                             data_NPP_rec.(tmp.limv).TOT.(tmp.vstr)(loni,lati,:), 'rows','complete');
% % % % %                         data_NPP_rec_corr.(tmp.limv).TOT.(['PP_', tmp.vstr])(loni,lati)=tmp_corr_r(1,2);
% % % % %                     end
% % % % %                 end
% % % % %             end
% % % % % 
% % % % %             
% % % % %             %% correlation between NPP <-> NPP_rec in same dataset map --------------------------------------
% % % % %                 fig_cfg.c_lim = [-1 1];
% % % % %                 fig_cfg.name_rgn = 'Glob';
% % % % %                 fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
% % % % %                 fig_cfg.x_lim = [-180 180];
% % % % %                 fig_cfg.y_lim = [-80 89];
% % % % %                 fig_cfg.fig_size = [0,0,6,3.5];
% % % % %                 fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
% % % % %                 fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
% % % % %                 fig_cfg.title_pos = [0.5,0.93];
% % % % %                 fig_cfg.p_lim =0.1;
% % % % %                 [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
% % % % %             
% % % % %                 tmp.X=grid.tlong([end, 1:end],:);
% % % % %                 tmp.Y=grid.tlat([end, 1:end],:);
% % % % % 
% % % % %                 %% data set
% % % % %                 tmp.C=data_NPP_rec_corr.(tmp.limv).TOT.(tmp.vstr) .* grid.mask_ocn;
% % % % % 
% % % % %                 tmp.C=tmp.C([end, 1:end],:); 
% % % % %             
% % % % %                 fig_cfg.fig_name='data_NPP_rec_corr_explained';
% % % % %                 fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
% % % % %                     'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
% % % % %                 %% map setting
% % % % %                 ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
% % % % %                     'fontname','freeserif'); 
% % % % %             
% % % % %                 axis off; 
% % % % %                 hold on;
% % % % %                 setm(ax_m,'origin',[0,180],'MapLatLimit',fig_cfg.y_lim);
% % % % %                 set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
% % % % %                 text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
% % % % %                 'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
% % % % %                 'fontsize',14,'fontname','freeserif','interpreter','none')
% % % % %             
% % % % %                 %% caxis & colorbar
% % % % %                 caxis(ax_m, fig_cfg.c_lim); 
% % % % %                 colormap(fig_cfg.c_map);
% % % % %                 cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
% % % % %                 set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
% % % % %                 title(cb,'R','fontsize',12);
% % % % %             
% % % % %                 %% draw on ax_m
% % % % %                 h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
% % % % %                 shading flat;
% % % % %                 geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
% % % % %             
% % % % %             
% % % % %                 %% frame and label setting
% % % % %                 setm(ax_m,'frame','on','FLineWidth',1);
% % % % %             
% % % % %                 label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
% % % % %                 label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
% % % % %                 mlabel; plabel;
% % % % %                 label_y=plabel; label_x=mlabel;
% % % % %                 for lxi=1:length(label_x)
% % % % %                     tmp.tmppos=label_x(lxi,1).Position;
% % % % %                     tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55;
% % % % %                     label_x(lxi,1).Position=tmp.tmppos;
% % % % %                     label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
% % % % %                 end
% % % % %                 for lyi=1:length(label_y)
% % % % %                     label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
% % % % %                 end
% % % % %             
% % % % %                 %% save
% % % % %                 dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, 'corr_explained', filesep, tmp.limv, filesep, 'TOT'];
% % % % %                 if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
% % % % %                 cfg.figname=[dirs.figdir, filesep, 'corr_exp_map_', tmp.limv, '_', 'TOT', '_', tmp.vstr, '.tif'];
% % % % %                 print(fig_h, cfg.figname, '-dpng');
% % % % %                 disp(cfg.figname)
% % % % %                 RemoveWhiteSpace([], 'file', cfg.figname);
% % % % %                 close all;
% % % % % 
% % % % % 
% % % % %         %% Potential Predictability of NPP, but by NPP_rec map --------------------------------------
% % % % %                 fig_cfg.c_lim = [-1 1];
% % % % %                 [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
% % % % %             
% % % % %                 tmp.X=grid.tlong([end, 1:end],:);
% % % % %                 tmp.Y=grid.tlat([end, 1:end],:);
% % % % % 
% % % % %                 %% data set
% % % % %                 tmp.C=data_NPP_rec_corr.(tmp.limv).TOT.(['PP_', tmp.vstr]) .* grid.mask_ocn;
% % % % % %                 switch tmp.limv
% % % % % %                     case {'Fe', 'N', 'SiO3', 'P'}
% % % % % %                         if strcmp(tmp.vstr, ['hcst_int_l',tmp.lyear_str])==1
% % % % % %                             tmp.mask=data_min_freq.([tmp.limv, '_lim_phyto']).TOT.(['model_l',tmp.lyear_str]).*grid.mask_ocn;
% % % % % %                         else
% % % % % %                             tmp.mask=data_min_freq.([tmp.limv, '_lim_phyto']).TOT.(tmp.vstr).*grid.mask_ocn;
% % % % % %                         end
% % % % % %                         tmp.mask(tmp.mask<0.8)=NaN;
% % % % % %                         tmp.C=tmp.C.*tmp.mask;
% % % % % %                 end
% % % % %                 tmp.C=tmp.C([end, 1:end],:); 
% % % % %             
% % % % %                 fig_cfg.fig_name='data_NPP_rec_PP';
% % % % %                 fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
% % % % %                     'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
% % % % %                 %% map setting
% % % % %                 ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
% % % % %                     'fontname','freeserif'); 
% % % % %             
% % % % %                 axis off; 
% % % % %                 hold on;
% % % % %                 setm(ax_m,'origin',[0,180],'MapLatLimit',fig_cfg.y_lim);
% % % % %                 set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
% % % % %                 text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
% % % % %                 'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
% % % % %                 'fontsize',14,'fontname','freeserif','interpreter','none')
% % % % %             
% % % % %                 %% caxis & colorbar
% % % % %                 caxis(ax_m, fig_cfg.c_lim); 
% % % % %                 colormap(fig_cfg.c_map);
% % % % %                 cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
% % % % %                 set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
% % % % %                 title(cb,'R','fontsize',12);
% % % % %             
% % % % %                 %% draw on ax_m
% % % % %                 h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
% % % % %                 shading flat;
% % % % %                 geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
% % % % %             
% % % % %             
% % % % %                 %% frame and label setting
% % % % %                 setm(ax_m,'frame','on','FLineWidth',1);
% % % % %             
% % % % %                 label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
% % % % %                 label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
% % % % %                 mlabel; plabel;
% % % % %                 label_y=plabel; label_x=mlabel;
% % % % %                 for lxi=1:length(label_x)
% % % % %                     tmp.tmppos=label_x(lxi,1).Position;
% % % % %                     tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55;
% % % % %                     label_x(lxi,1).Position=tmp.tmppos;
% % % % %                     label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
% % % % %                 end
% % % % %                 for lyi=1:length(label_y)
% % % % %                     label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
% % % % %                 end
% % % % %             
% % % % %                 %% save
% % % % %                 dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, 'corr_PP', filesep, tmp.limv, filesep, 'TOT'];
% % % % %                 if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
% % % % %                 cfg.figname=[dirs.figdir, filesep, 'corr_PP_map_', tmp.limv, '_', 'TOT', '_', tmp.vstr, '.tif'];
% % % % %                 print(fig_h, cfg.figname, '-dpng');
% % % % %                 disp(cfg.figname)
% % % % %                 RemoveWhiteSpace([], 'file', cfg.figname);
% % % % %                 close all;
% % % % % 
% % % % % 
% % % % %         end
% % % % %     end
% % % % %     toc;
% % % % % end



%% data_TOT
for lyear = [0] % [0,4]
    tic
    tmp.lyear_str=num2str(lyear, '%02i');
    cfg.vstrset={'assm', ['lens2_l',tmp.lyear_str], ['model_l',tmp.lyear_str]};
    
    for pind=1:length(phytoset)
        tmp.phyto=phytoset{pind};
        matname=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/ocn/photoC_TOT_zint_100m/',tmp.phyto, '_nut_NPP_l', tmp.lyear_str, '.mat'];
        load(matname);
        data_all.TOT.(tmp.phyto)=data_NPP_rec;
        data_all_min_freq.TOT.(tmp.phyto)=data_min_freq;
% 
%         %% photoC_phyto
%         tmp.varname=['photoC_',tmp.phyto,'_zint_100m'];
%         dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', tmp.varname];
%         fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
%                     '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
%                     '_l', tmp.lyear_str, 'y.mat'];
%         load(fig_cfg.mat_name, 'data');
%         data_all_NPP_phyto.(tmp.phyto)=data;

    end
    tmp.varname=['photoC_','TOT','_zint_100m'];
    dirs.hcstmatroot=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp, '/', tmp.varname];
    fig_cfg.mat_name=[dirs.hcstmatroot, filesep, 'hcst_corr_assm_', tmp.varname, ...
                '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
                '_l', tmp.lyear_str, 'y.mat'];
    load(fig_cfg.mat_name, 'data');
    data_all_NPP_phyto.TOT=data;

    for stri=1:length(cfg.vstrset)
        tmp.vstr=cfg.vstrset{stri};
        cfg.limv={'ALL', 'T', 'L', 'V', 'Fe', 'N', 'SiO3', 'P'};
        for ni=1:length(cfg.limv)
            tmp.limv=cfg.limv{ni};
            switch tmp.limv
                case {'Fe', 'N', 'SiO3', 'P'}
                    for pind=1:length(phytoset)
                        tmp.phyto=phytoset{pind};
                        tmp.mask=data_all_min_freq.TOT.(tmp.phyto).([tmp.limv, '_lim_phyto']).(tmp.phyto).(tmp.vstr).*grid.mask_ocn;                            
                        tmp.mask(tmp.mask<0.8)=NaN;
                        tmp.comb(:,:,:,pind)=data_all.TOT.(tmp.phyto).(tmp.limv).(tmp.phyto).(tmp.vstr).*tmp.mask;
                    end
                otherwise
                    tmp.comb(:,:,:,1)=data_all.TOT.diat.(tmp.limv).diat.(tmp.vstr);
                    tmp.comb(:,:,:,2)=data_all.TOT.diaz.(tmp.limv).diaz.(tmp.vstr);
                    tmp.comb(:,:,:,3)=data_all.TOT.sp.(tmp.limv).sp.(tmp.vstr);
            end
            
            data_all_NPP_rec.(tmp.limv).TOT.(tmp.vstr)= sum(tmp.comb, 4, 'omitnan');
%             data_NPP_rec.(tmp.limv).TOT.(tmp.vstr) = ...
%                 data_all.TOT.diat.(tmp.limv).diat.(tmp.vstr)  ...
%                 + data_all.TOT.diaz.(tmp.limv).diaz.(tmp.vstr) ...
%                 + data_all.TOT.sp.(tmp.limv).sp.(tmp.vstr);

            data_all_NPP_rec_corr.(tmp.limv).TOT.(tmp.vstr)=NaN(size(grid.tlong));
            data_all_NPP_rec_corr.(tmp.limv).TOT.(['PP_', tmp.vstr])=NaN(size(grid.tlong));
            
            data_all_NPP_phyto.TOT.trend.('photoC_TOT_zint_100m_assm')=NaN(size(grid.tlong));
            data_all_NPP_rec.(tmp.limv).TOT.trend.(tmp.vstr)=NaN(size(grid.tlong));
            data_all_NPP_rec.(tmp.limv).TOT.ratio.(tmp.vstr)=NaN(size(grid.tlong));
            for loni=1:size(grid.tlong,1)
                for lati=1:size(grid.tlat,2)
                    if(isfinite(data_all_NPP_phyto.TOT.photoC_TOT_zint_100m_assm(loni,lati,1)))
                        %% rec validation (raw <-> reconstructed)
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_all_NPP_phyto.TOT.(['photoC_TOT_zint_100m_',tmp.vstr])(loni,lati,:), ...
                            data_all_NPP_rec.(tmp.limv).TOT.(tmp.vstr)(loni,lati,:), 'rows','complete');
                        data_all_NPP_rec_corr.(tmp.limv).TOT.(tmp.vstr)(loni,lati)=tmp_corr_r(1,2);
                        data_all_NPP_rec_corr.(tmp.limv).TOT.(['p_', tmp.vstr])(loni,lati)=tmp_corr_p(1,2);
                        if tmp_corr_p(1,2) <= 0.05
                            data_all_NPP_rec_corr.(tmp.limv).TOT.(['sig_',tmp.vstr])(loni,lati)=tmp_corr_r(1,2);
                        else
                            data_all_NPP_rec_corr.(tmp.limv).TOT.(['sig_',tmp.vstr])(loni,lati)=NaN;
                        end

                        %% rec Potential Predictability (assm raw <-> model reconstructed)
                        [tmp_corr_r,tmp_corr_p] = ...
                            corrcoef(data_all_NPP_phyto.TOT.(['photoC_TOT_zint_100m_assm'])(loni,lati,:), ...
                            data_all_NPP_rec.(tmp.limv).TOT.(tmp.vstr)(loni,lati,:), 'rows','complete');
                        data_all_NPP_rec_corr.(tmp.limv).TOT.(['PP_', tmp.vstr])(loni,lati)=tmp_corr_r(1,2);                        
                        data_all_NPP_rec_corr.(tmp.limv).TOT.(['PP_p_', tmp.vstr])(loni,lati)=tmp_corr_p(1,2);

                        if tmp_corr_p(1,2) <= 0.05
                            data_all_NPP_rec_corr.(tmp.limv).TOT.(['sig_PP_', tmp.vstr])(loni,lati)=tmp_corr_r(1,2);                        
                        else
                            data_all_NPP_rec_corr.(tmp.limv).TOT.(['sig_PP_', tmp.vstr])(loni,lati)=NaN;                        
                        end

                        %% original and rec trend
                        [tmp.data_det, data_all_NPP_phyto.TOT.trend.('photoC_TOT_zint_100m_assm')(loni,lati)] = ...
                            Func_0028_detrend_linear_1d(data_all_NPP_phyto.TOT.('photoC_TOT_zint_100m_assm')(loni,lati,:));
                        [tmp.data_det, data_all_NPP_rec.(tmp.limv).TOT.trend.(tmp.vstr)(loni,lati)] = ...
                            Func_0028_detrend_linear_1d(data_all_NPP_rec.(tmp.limv).TOT.(tmp.vstr)(loni,lati,:));
                        
                        data_all_NPP_rec.(tmp.limv).TOT.ratio.(tmp.vstr)(loni,lati) = ...
                            data_all_NPP_phyto.TOT.trend.('photoC_TOT_zint_100m_assm')(loni,lati) ./ ...
                            data_all_NPP_rec.(tmp.limv).TOT.trend.(tmp.vstr)(loni,lati);



                    else
                        data_all_NPP_rec_corr.(tmp.limv).TOT.(tmp.vstr)(loni,lati)=NaN;
                        data_all_NPP_rec_corr.(tmp.limv).TOT.(['p_', tmp.vstr])(loni,lati)=NaN;
                        data_all_NPP_rec_corr.(tmp.limv).TOT.(['sig_',tmp.vstr])(loni,lati)=NaN;
                        data_all_NPP_rec_corr.(tmp.limv).TOT.(['PP_', tmp.vstr])(loni,lati)=NaN;
                        data_all_NPP_rec_corr.(tmp.limv).TOT.(['PP_p_', tmp.vstr])(loni,lati)=NaN;
                        data_all_NPP_rec_corr.(tmp.limv).TOT.(['sig_PP_', tmp.vstr])(loni,lati)=NaN;                        

                    end
                end
            end

            
% % % % %             %% correlation between NPP <-> NPP_rec in same dataset map --------------------------------------
% % % % %                 fig_cfg.c_lim = [-1 1];
% % % % %                 [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
% % % % %             
% % % % %                 tmp.X=grid.tlong([end, 1:end],:);
% % % % %                 tmp.Y=grid.tlat([end, 1:end],:);
% % % % % 
% % % % %                 %% data set
% % % % %                 tmp.C=data_NPP_rec_corr.(tmp.limv).TOT.(tmp.vstr) .* grid.mask_ocn;
% % % % % 
% % % % %                 tmp.C=tmp.C([end, 1:end],:); 
% % % % %             
% % % % %                 fig_cfg.fig_name='data_NPP_rec_corr_explained';
% % % % %                 fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
% % % % %                     'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
% % % % %                 %% map setting
% % % % %                 ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
% % % % %                     'fontname','freeserif'); 
% % % % %             
% % % % %                 axis off; 
% % % % %                 hold on;
% % % % %                 setm(ax_m,'origin',[0,180],'MapLatLimit',fig_cfg.y_lim);
% % % % %                 set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
% % % % %                 text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
% % % % %                 'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
% % % % %                 'fontsize',14,'fontname','freeserif','interpreter','none')
% % % % %             
% % % % %                 %% caxis & colorbar
% % % % %                 caxis(ax_m, fig_cfg.c_lim); 
% % % % %                 colormap(fig_cfg.c_map);
% % % % %                 cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
% % % % %                 set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
% % % % %                 title(cb,'R','fontsize',12);
% % % % %             
% % % % %                 %% draw on ax_m
% % % % %                 h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
% % % % %                 shading flat;
% % % % %                 geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
% % % % %             
% % % % %             
% % % % %                 %% frame and label setting
% % % % %                 setm(ax_m,'frame','on','FLineWidth',1);
% % % % %             
% % % % %                 label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
% % % % %                 label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
% % % % %                 mlabel; plabel;
% % % % %                 label_y=plabel; label_x=mlabel;
% % % % %                 for lxi=1:length(label_x)
% % % % %                     tmp.tmppos=label_x(lxi,1).Position;
% % % % %                     tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55;
% % % % %                     label_x(lxi,1).Position=tmp.tmppos;
% % % % %                     label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
% % % % %                 end
% % % % %                 for lyi=1:length(label_y)
% % % % %                     label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
% % % % %                 end
% % % % %             
% % % % %                 %% save
% % % % %                 dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, 'corr_explained', filesep, tmp.limv, filesep, 'TOT'];
% % % % %                 if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
% % % % %                 cfg.figname=[dirs.figdir, filesep, 'corr_exp_map_', tmp.limv, '_', 'TOT', '_', tmp.vstr, '.tif'];
% % % % %                 print(fig_h, cfg.figname, '-dpng');
% % % % %                 disp(cfg.figname)
% % % % %                 RemoveWhiteSpace([], 'file', cfg.figname);
% % % % %                 close all;


% % % % %         %% Potential Predictability of NPP, but by NPP_rec map --------------------------------------
% % % % %                 fig_cfg.c_lim = [-1 1];
% % % % %                 [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
% % % % %             
% % % % %                 tmp.X=grid.tlong([end, 1:end],:);
% % % % %                 tmp.Y=grid.tlat([end, 1:end],:);
% % % % % 
% % % % %                 %% data set
% % % % %                 tmp.C=data_NPP_rec_corr.(tmp.limv).TOT.(['PP_', tmp.vstr]) .* grid.mask_ocn;
% % % % % %                 switch tmp.limv
% % % % % %                     case {'Fe', 'N', 'SiO3', 'P'}
% % % % % %                         if strcmp(tmp.vstr, ['hcst_int_l',tmp.lyear_str])==1
% % % % % %                             tmp.mask=data_min_freq.([tmp.limv, '_lim_phyto']).TOT.(['model_l',tmp.lyear_str]).*grid.mask_ocn;
% % % % % %                         else
% % % % % %                             tmp.mask=data_min_freq.([tmp.limv, '_lim_phyto']).TOT.(tmp.vstr).*grid.mask_ocn;
% % % % % %                         end
% % % % % %                         tmp.mask(tmp.mask<0.8)=NaN;
% % % % % %                         tmp.C=tmp.C.*tmp.mask;
% % % % % %                 end
% % % % %                 tmp.C=tmp.C([end, 1:end],:); 
% % % % %             
% % % % %                 fig_cfg.fig_name='data_NPP_rec_PP';
% % % % %                 fig_h = figure('name',fig_cfg.fig_name,'PaperUnits','inches', ...
% % % % %                     'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
% % % % %                 %% map setting
% % % % %                 ax_m = axesm('MapProjection',fig_cfg.map_proj,'grid','on','fontsize',14, ...
% % % % %                     'fontname','freeserif'); 
% % % % %             
% % % % %                 axis off; 
% % % % %                 hold on;
% % % % %                 setm(ax_m,'origin',[0,180],'MapLatLimit',fig_cfg.y_lim);
% % % % %                 set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
% % % % %                 text(ax_m,fig_cfg.title_pos(1),fig_cfg.title_pos(2),fig_cfg.fig_name, ...
% % % % %                 'units','normalized', 'horizontalalignment','center', 'verticalalignment','middle', ...
% % % % %                 'fontsize',14,'fontname','freeserif','interpreter','none')
% % % % %             
% % % % %                 %% caxis & colorbar
% % % % %                 caxis(ax_m, fig_cfg.c_lim); 
% % % % %                 colormap(fig_cfg.c_map);
% % % % %                 cb = colorbar(ax_m,'units','inches','position',fig_cfg.cb_size);
% % % % %                 set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
% % % % %                 title(cb,'R','fontsize',12);
% % % % %             
% % % % %                 %% draw on ax_m
% % % % %                 h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
% % % % %                 shading flat;
% % % % %                 geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
% % % % %             
% % % % %             
% % % % %                 %% frame and label setting
% % % % %                 setm(ax_m,'frame','on','FLineWidth',1);
% % % % %             
% % % % %                 label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
% % % % %                 label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
% % % % %                 mlabel; plabel;
% % % % %                 label_y=plabel; label_x=mlabel;
% % % % %                 for lxi=1:length(label_x)
% % % % %                     tmp.tmppos=label_x(lxi,1).Position;
% % % % %                     tmp.tmppos(2)=-fig_cfg.ax_size(4)+1.55;
% % % % %                     label_x(lxi,1).Position=tmp.tmppos;
% % % % %                     label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
% % % % %                 end
% % % % %                 for lyi=1:length(label_y)
% % % % %                     label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
% % % % %                 end
% % % % %             
% % % % %                 %% save
% % % % %                 dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, 'corr_PP', filesep, tmp.limv, filesep, 'TOT'];
% % % % %                 if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
% % % % %                 cfg.figname=[dirs.figdir, filesep, 'corr_PP_map_', tmp.limv, '_', 'TOT', '_', tmp.vstr, '.tif'];
% % % % %                 print(fig_h, cfg.figname, '-dpng');
% % % % %                 disp(cfg.figname)
% % % % %                 RemoveWhiteSpace([], 'file', cfg.figname);
% % % % %                 close all;

        %% trend ratio of NPP, but by NPP_rec map --------------------------------------
                fig_cfg.c_lim = [-1 1];
                fig_cfg.c_lim = [-1 1];
                fig_cfg.name_rgn = 'Glob';
                fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
                fig_cfg.x_lim = [-180 180];
                fig_cfg.y_lim = [-80 89];
                fig_cfg.fig_size = [0,0,6,3.5];
                fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
                fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
                fig_cfg.title_pos = [0.5,0.93];
                fig_cfg.p_lim =0.1;
            
                tmp.X=grid.tlong([end, 1:end],:);
                tmp.Y=grid.tlat([end, 1:end],:);

                [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
            
                tmp.X=grid.tlong([end, 1:end],:);
                tmp.Y=grid.tlat([end, 1:end],:);

                %% data set
                tmp.C=data_all_NPP_rec.(tmp.limv).TOT.ratio.(tmp.vstr) .* grid.mask_ocn;

                maxabs=max(abs(tmp.C(:)-1));

%                 switch tmp.limv
%                     case {'Fe', 'N', 'SiO3', 'P'}
%                         if strcmp(tmp.vstr, ['hcst_int_l',tmp.lyear_str])==1
%                             tmp.mask=data_min_freq.([tmp.limv, '_lim_phyto']).TOT.(['model_l',tmp.lyear_str]).*grid.mask_ocn;
%                         else
%                             tmp.mask=data_min_freq.([tmp.limv, '_lim_phyto']).TOT.(tmp.vstr).*grid.mask_ocn;
%                         end
%                         tmp.mask(tmp.mask<0.8)=NaN;
%                         tmp.C=tmp.C.*tmp.mask;
%                 end
                tmp.C=tmp.C([end, 1:end],:); 
            
                fig_cfg.fig_name='NPP_tr_ratio';
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
                caxis(ax_m, [1-maxabs 1+maxabs]); 
                caxis(ax_m, [-2 2]); 

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
                end
            
                %% save
                dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, 'tr_ratio', filesep, tmp.limv, filesep, 'TOT'];
                if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
                cfg.figname=[dirs.figdir, filesep, 'tr_ratio_map_', tmp.limv, '_', 'TOT', '_', tmp.vstr, '.tif'];
                print(fig_h, cfg.figname, '-dpng');
                disp(cfg.figname)
                RemoveWhiteSpace([], 'file', cfg.figname);
                close all;


        %% trend of NPP-rec map --------------------------------------
                fig_cfg.c_lim = [-1 1];
                fig_cfg.c_lim = [-1 1];
                fig_cfg.name_rgn = 'Glob';
                fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
                fig_cfg.x_lim = [-180 180];
                fig_cfg.y_lim = [-80 89];
                fig_cfg.fig_size = [0,0,6,3.5];
                fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
                fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
                fig_cfg.title_pos = [0.5,0.93];
                fig_cfg.p_lim =0.1;
            
                tmp.X=grid.tlong([end, 1:end],:);
                tmp.Y=grid.tlat([end, 1:end],:);

                [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
            
                tmp.X=grid.tlong([end, 1:end],:);
                tmp.Y=grid.tlat([end, 1:end],:);

                %% data set
                tmp.C=data_all_NPP_rec.(tmp.limv).TOT.trend.(tmp.vstr) .* grid.mask_ocn;
                
                tmp.D1=data_all_NPP_rec.(tmp.limv).TOT.trend.(tmp.vstr) .* grid.mask_ocn;
                tmp.D2=data_all_NPP_phyto.TOT.trend.('photoC_TOT_zint_100m_assm') .* grid.mask_ocn;
                tmp.D=[tmp.D1 tmp.D2];
                maxabs=prctile(abs(tmp.D(:)), 90);


%                 switch tmp.limv
%                     case {'Fe', 'N', 'SiO3', 'P'}
%                         if strcmp(tmp.vstr, ['hcst_int_l',tmp.lyear_str])==1
%                             tmp.mask=data_min_freq.([tmp.limv, '_lim_phyto']).TOT.(['model_l',tmp.lyear_str]).*grid.mask_ocn;
%                         else
%                             tmp.mask=data_min_freq.([tmp.limv, '_lim_phyto']).TOT.(tmp.vstr).*grid.mask_ocn;
%                         end
%                         tmp.mask(tmp.mask<0.8)=NaN;
%                         tmp.C=tmp.C.*tmp.mask;
%                 end
                tmp.C=tmp.C([end, 1:end],:); 
            
                fig_cfg.fig_name='NPP_rec_trend';
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
                if maxabs~=0
                    caxis(ax_m, [-maxabs maxabs]); 
                end

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
                end
            
                %% save
                dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, 'tr_ratio', filesep, tmp.limv, filesep, 'TOT'];
                if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
                cfg.figname=[dirs.figdir, filesep, 'tr_npp-rec_map_', tmp.limv, '_', 'TOT', '_', tmp.vstr, '.tif'];
                print(fig_h, cfg.figname, '-dpng');
                disp(cfg.figname)
                RemoveWhiteSpace([], 'file', cfg.figname);
                close all;

        %% trend of NPP map --------------------------------------
                fig_cfg.c_lim = [-1 1];
                fig_cfg.c_lim = [-1 1];
                fig_cfg.name_rgn = 'Glob';
                fig_cfg.map_proj = 'eqdcylin';  % robinson, eqdcylin
                fig_cfg.x_lim = [-180 180];
                fig_cfg.y_lim = [-80 89];
                fig_cfg.fig_size = [0,0,6,3.5];
                fig_cfg.ax_size = [0.3,0.7,5.4,2.7];
                fig_cfg.cb_size = [5.15,0.8,0.15,2.3];
                fig_cfg.title_pos = [0.5,0.93];
                fig_cfg.p_lim =0.1;
            
                tmp.X=grid.tlong([end, 1:end],:);
                tmp.Y=grid.tlat([end, 1:end],:);

                [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
            
                tmp.X=grid.tlong([end, 1:end],:);
                tmp.Y=grid.tlat([end, 1:end],:);

                %% data set
                tmp.C=data_all_NPP_phyto.TOT.trend.('photoC_TOT_zint_100m_assm') .* grid.mask_ocn;

                tmp.D1=data_all_NPP_rec.(tmp.limv).TOT.trend.(tmp.vstr) .* grid.mask_ocn;
                tmp.D2=data_all_NPP_phyto.TOT.trend.('photoC_TOT_zint_100m_assm') .* grid.mask_ocn;
                tmp.D=[tmp.D1 tmp.D2];
                maxabs=prctile(abs(tmp.D(:)), 90);

%                 switch tmp.limv
%                     case {'Fe', 'N', 'SiO3', 'P'}
%                         if strcmp(tmp.vstr, ['hcst_int_l',tmp.lyear_str])==1
%                             tmp.mask=data_min_freq.([tmp.limv, '_lim_phyto']).TOT.(['model_l',tmp.lyear_str]).*grid.mask_ocn;
%                         else
%                             tmp.mask=data_min_freq.([tmp.limv, '_lim_phyto']).TOT.(tmp.vstr).*grid.mask_ocn;
%                         end
%                         tmp.mask(tmp.mask<0.8)=NaN;
%                         tmp.C=tmp.C.*tmp.mask;
%                 end
                tmp.C=tmp.C([end, 1:end],:); 
            
                fig_cfg.fig_name='NPP_trend';
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
                caxis(ax_m, [-maxabs maxabs]); 

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
                end
            
                %% save
                dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, 'tr_ratio', filesep, tmp.limv, filesep, 'TOT'];
                if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
                cfg.figname=[dirs.figdir, filesep, 'tr_npp_map_', tmp.limv, '_', 'TOT', '_', tmp.vstr, '.tif'];
                print(fig_h, cfg.figname, '-dpng');
                disp(cfg.figname)
                RemoveWhiteSpace([], 'file', cfg.figname);
                close all;

        end
    end
    toc;

    matname2=[tmp.rootpath, '/kimyy/Model/CESM2/ESP/HCST_EXP/mat/ocn/photoC_TOT_zint_100m/','all', '_NPP_l', tmp.lyear_str, '.mat'];
    save(matname2, 'data_all_NPP_rec', 'data_all_NPP_phyto', 'grid')
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


% save('/Users/kimyy/Model/CESM2/ESP/HCST_EXP/mat/ocn/photoC_TOT_zint_100m/tmp.mat')

% % % pcolor(data_mask_nut_dominant.ind.diat.model_l00'.*grid.mask_ocn'); shading flat; colorbar;
% % % figure; pcolor(data_mask_nut_dominant.val.diat.model_l00'.*grid.mask_ocn'); shading flat; colorbar;
% % % 
% % % pcolor(data_mask_nut_dominant.ind.diat.assm'.*grid.mask_ocn'); shading flat; colorbar;
% % % pcolor(data_mask_nut_dominant.ind.sp.assm'.*grid.mask_ocn'); shading flat; colorbar;
% % % pcolor(data_mask_nut_dominant.ind.diaz.assm'.*grid.mask_ocn'); shading flat; colorbar;


% data_NPP_diat_rec_assm = data_light.diat_light_lim_Cweight_avg_100m_assm .* ...
%                     lim_set .* ...
%                     1.7.^((data_TEMP.TEMP_assm-30.0)./10.0) .* ...
%                     mu_ref .* ...
%                     data_diatC.diatC_assm ...
%                     ./ 86400.0 .* 10000.0;


% pcolor(data_NPP_rec_corr.ALL.diat.assm'); shading flat; colorbar;
% pcolor(data_NPP_rec_corr.ALL.diat.lens2_l04'); shading flat; colorbar;
% pcolor(data_NPP_rec_corr.ALL.diat.model_l00'); shading flat; colorbar;
% 
% pcolor(data_NPP_rec_corr.L.diat.model_l00'); shading flat; colorbar;
% pcolor(data_NPP_rec_corr.T.diat.model_l00'); shading flat; colorbar;
% pcolor(data_NPP_rec_corr.V.diat.model_l00'); shading flat; colorbar;
% pcolor(data_NPP_rec_corr.Fe.diat.model_l00'); shading flat; colorbar;
% pcolor(data_NPP_rec_corr.N.diat.model_l00'); shading flat; colorbar;
% pcolor(data_NPP_rec_corr.P.diat.model_l00'); shading flat; colorbar;
% pcolor(data_NPP_rec_corr.SiO3.diat.model_l00'); shading flat; colorbar;
% 
% pcolor(data_NPP_rec_corr.L.sp.model_l00'); shading flat; colorbar;

% % % pcolor(data_NPP_rec_corr.L.diat.PP_model_l04'); shading flat; colorbar; colormap(fig_cfg.c_map); caxis([-1 1]);
% % % [fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_10', tmp.dropboxpath);
% % % 
% % % pcolor(data_NPP_rec_corr.T.diat.PP_model_l04'); shading flat; colorbar; colormap(fig_cfg.c_map); caxis([-1 1]);
% % % pcolor(data_NPP_rec_corr.V.diat.PP_model_l04'); shading flat; colorbar; colormap(fig_cfg.c_map); caxis([-1 1]);
% % % pcolor(data_NPP_rec_corr.ALL.diat.PP_model_l04'); shading flat; colorbar; colormap(fig_cfg.c_map); caxis([-1 1]);
% % % pcolor(data_NPP_rec_corr.P.diat.PP_model_l04'); shading flat; colorbar; colormap(fig_cfg.c_map); caxis([-1 1]);
% % % pcolor(data_NPP_rec_corr.SiO3.diat.PP_model_l04'); shading flat; colorbar; colormap(fig_cfg.c_map); caxis([-1 1]);
% % % pcolor(data_NPP_rec_corr.N.diat.PP_model_l04'); shading flat; colorbar; colormap(fig_cfg.c_map); caxis([-1 1]);
% % % pcolor(data_NPP_rec_corr.Fe.diat.PP_model_l04'); shading flat; colorbar; colormap(fig_cfg.c_map); caxis([-1 1]);
% % % 
% % % 
% % % pcolor(data_NPP_rec_corr.T.diat.PP_hcst_int_l04'.*grid.mask_ocn'); shading flat; colorbar; colormap(fig_cfg.c_map); caxis([-1 1]);
% % % pcolor(data_NPP_rec_corr.V.diat.PP_hcst_int_l04'.*grid.mask_ocn'); shading flat; colorbar; colormap(fig_cfg.c_map); caxis([-1 1]);
% % % pcolor(data_NPP_rec_corr.L.diat.PP_hcst_int_l04'.*grid.mask_ocn'); shading flat; colorbar; colormap(fig_cfg.c_map); caxis([-1 1]);
% % % pcolor(data_NPP_rec_corr.ALL.diat.PP_hcst_int_l04'.*grid.mask_ocn'); shading flat; colorbar; colormap(fig_cfg.c_map); caxis([-1 1]);


% plot(squeeze(data_NPP_phyto.diat.photoC_diat_zint_100m_model_l04(170,200,:)))
% hold on
% plot(squeeze(data_NPP_rec.ALL.diat.model_l04(170,200,:)))
% plot(squeeze(data_NPP_rec.T.diat.model_l04(170,200,:)))
% hold off
% plot(squeeze(data_NPP_rec.T.diat.lens2_l04(170,200,:)))

% 1. contribution of each component (light, N, P); L_driven NPP <-> NPP ...
% 2. PP of each component driven NPP (light, N, P); L_driven NPP(hcst) <-> NPP(assm) ...
% 3. PP of each nutrient driven NPP (light, N, P); Fe_driven NPP(hcst) <-> NPP(assm) ...

% nutrient contribution





%% subsampled test
% % % [grid.id_w, grid.id_e, grid.id_s, grid.id_n] = Func_0012_findind_Y(0.5, [180, 20], ...
% % %             grid.tlong, ...
% % %             grid.tlat, 'CESM2'); % find valid lon, lat index near station
% % % 
% % % tmp.NPP_ref=data_NPP.photoC_diat_zint_100m_assm(grid.id_w, grid.id_s, :);
% % % % mmol/m^3 cm/s
% % % tmp.NPP_recons=data_light.diat_light_lim_Cweight_avg_100m_assm(grid.id_w, grid.id_s,:).* ...
% % %     data_NO3.diat_N_lim_Cweight_avg_100m_assm(grid.id_w, grid.id_s,:) .* ...
% % %     1.7.^((data_TEMP.TEMP_assm(grid.id_w, grid.id_s,:)-30.0)./10.0) .* ...
% % %     mu_ref .* ...
% % %     data_diatC.diatC_assm(grid.id_w, grid.id_s,:) ...
% % %     ./86400.*10000; % /day -> /sec, 100m integration (cm -> 10000)
% % % % 	mmol/m^3
% % % %  mu_ref : /day
% % % 
% % % tmp.NPP_recons_N = mean(data_light.diat_light_lim_Cweight_avg_100m_assm(grid.id_w, grid.id_s,:)).* ...
% % %     data_NO3.diat_N_lim_Cweight_avg_100m_assm(grid.id_w, grid.id_s,:) .* ...
% % %     mean(1.7.^((data_TEMP.TEMP_assm(grid.id_w, grid.id_s,:)-30.0)./10.0)) .* ...
% % %     mu_ref .* ...
% % %     mean(data_diatC.diatC_assm(grid.id_w, grid.id_s,:)) ...
% % %     ./86400.*10000; % /day -> /sec, 100m integration (cm -> 10000)
% % % 
% % % tmp.NPP_recons_T = mean(data_light.diat_light_lim_Cweight_avg_100m_assm(grid.id_w, grid.id_s,:)).* ...
% % %     mean(data_NO3.diat_N_lim_Cweight_avg_100m_assm(grid.id_w, grid.id_s,:)) .* ...
% % %     1.7.^((data_TEMP.TEMP_assm(grid.id_w, grid.id_s,:)-30.0)./10.0) .* ...
% % %     mu_ref .* ...
% % %     mean(data_diatC.diatC_assm(grid.id_w, grid.id_s,:)) ...
% % %     ./86400.*10000; % /day -> /sec, 100m integration (cm -> 10000)
% % % 
% % % tmp.NPP_recons_L = data_light.diat_light_lim_Cweight_avg_100m_assm(grid.id_w, grid.id_s,:).* ...
% % %     mean(data_NO3.diat_N_lim_Cweight_avg_100m_assm(grid.id_w, grid.id_s,:)) .* ...
% % %     mean(1.7.^((data_TEMP.TEMP_assm(grid.id_w, grid.id_s,:)-30.0)./10.0)) .* ...
% % %     mu_ref .* ...
% % %     mean(data_diatC.diatC_assm(grid.id_w, grid.id_s,:)) ...
% % %     ./86400.*10000; % /day -> /sec, 100m integration (cm -> 10000)
% % % 
% % % tmp.NPP_recons_C = mean(data_light.diat_light_lim_Cweight_avg_100m_assm(grid.id_w, grid.id_s,:)).* ...
% % %     mean(data_NO3.diat_N_lim_Cweight_avg_100m_assm(grid.id_w, grid.id_s,:)) .* ...
% % %     mean(1.7.^((data_TEMP.TEMP_assm(grid.id_w, grid.id_s,:)-30.0)./10.0)) .* ...
% % %     mu_ref .* ...
% % %     data_diatC.diatC_assm(grid.id_w, grid.id_s,:) ...
% % %     ./86400.*10000; % /day -> /sec, 100m integration (cm -> 10000)
% % % 
% % % plot(squeeze(tmp.NPP_ref));
% % % hold on
% % % plot(squeeze(tmp.NPP_recons));
% % % plot(squeeze(tmp.NPP_recons_N));
% % % plot(squeeze(tmp.NPP_recons_T));
% % % 
% % % hold off
% % % legend
% % % std(tmp.NPP_recons)/std(tmp.NPP_ref)
% % % corrcoef(tmp.NPP_recons, tmp.NPP_ref)
% % % corrcoef(tmp.NPP_recons_N, tmp.NPP_ref)
% % % corrcoef(tmp.NPP_recons_T, tmp.NPP_ref)
% % % corrcoef(tmp.NPP_recons_L, tmp.NPP_ref)
% % % corrcoef(tmp.NPP_recons_C, tmp.NPP_ref)




% % %% ts check plots
% % xpoint=10;
% % ypoint=120;
% % plot(cfg.iyears,squeeze(data_NPP.photoC_diat_zint_100m_assm(xpoint,ypoint,:)));
% % hold on
% % plot(cfg.iyears,squeeze(data_NPP_diat_rec(xpoint,ypoint,:)));
% % legend({'NPP-diat', 'NPP-diat-rec'})
% % title(['xind=',num2str(xpoint),',','yind=',num2str(ypoint),',', ...
% %     'lon=',num2str(grid.tlong(xpoint,ypoint)), ',','lat=',num2str(grid.tlat(xpoint,ypoint))])
% % hold off
% % tmpcorr=corrcoef(data_NPP.photoC_diat_zint_100m_assm(xpoint,ypoint,:),data_NPP_diat_rec(xpoint,ypoint,:) );
% % text( ...
% %     min(cfg.iyears), max([squeeze(data_NPP.photoC_diat_zint_100m_assm(xpoint,ypoint,:)); ...
% %     squeeze(data_NPP_diat_rec(xpoint,ypoint,:))]), ...
% %     ['corr=',num2str(round(tmpcorr(1,2),2))], 'fontsize', 20);
% % set(gca,'fontsize', 20)


% %% contribution of nutrient limitation
% mask_dominant=NaN(size(grid.tlong,1), size(grid.tlat,2), 4);
% for loni=1:size(grid.tlong,1)
%     for lati=1:size(grid.tlat,2)
%         data_Fe.min_freq(loni,lati)=length(find(lim_set_ind(loni,lati,:)==1))/length(cfg.iyears);
%         data_NO3.min_freq(loni,lati)=length(find(lim_set_ind(loni,lati,:)==2))/length(cfg.iyears);
%         data_SiO3.min_freq(loni,lati)=length(find(lim_set_ind(loni,lati,:)==3))/length(cfg.iyears);
%         data_PO4.min_freq(loni,lati)=length(find(lim_set_ind(loni,lati,:)==4))/length(cfg.iyears);
%         tmp.min_freq = [data_Fe.min_freq(loni,lati); data_NO3.min_freq(loni,lati); ...
%             data_SiO3.min_freq(loni,lati); data_PO4.min_freq(loni,lati)];
%         [max_minfreq,maxind_minfreq]=max(tmp.min_freq);
%         mask_dominant(loni,lati,maxind_minfreq)=max_minfreq;
%     end
% end

% % pcolor(mask_dominant(:,:,1)'.*grid.mask_ocn'); shading flat; colorbar; colormap(jet(10)); set(gca,'fontsize', 20);
% % pcolor(mask_dominant(:,:,2)'.*grid.mask_ocn'); shading flat; colorbar; colormap(jet(10)); set(gca,'fontsize', 20);
% % pcolor(mask_dominant(:,:,4)'.*grid.mask_ocn'); shading flat; colorbar; colormap(jet(10)); set(gca,'fontsize', 20);
% % pcolor(mask_dominant(:,:,3)'.*grid.mask_ocn'); shading flat; colorbar; colormap(jet(10)); set(gca,'fontsize', 20);



% pcolor(data_Fe.min_freq'.*grid.mask_ocn'); shading flat; colorbar; colormap(jet(10)); set(gca,'fontsize', 20);
% pcolor(data_NO3.min_freq'.*grid.mask_ocn'); shading flat; colorbar; colormap(jet(10)); set(gca,'fontsize', 20);
% pcolor(data_PO4.min_freq'.*grid.mask_ocn'); shading flat; colorbar; colormap(jet(10)); set(gca,'fontsize', 20);
% pcolor(data_SiO3.min_freq'.*grid.mask_ocn'); shading flat; colorbar; colormap(jet(10)); set(gca,'fontsize', 20);
% pcolor((data_Fe.min_freq+data_NO3.min_freq+data_PO4.min_freq+data_SiO3.min_freq)'.*grid.mask_ocn'); shading flat; colorbar; colormap(jet(10)); set(gca,'fontsize', 20);

% % % %% reconstruction of NPP
% % % for loni=1:size(grid.tlong,1)
% % %     for lati=1:size(grid.tlat,2)
% % %         data_NPP_diat_rec = data_light.diat_light_lim_Cweight_avg_100m_assm(grid.id_w, grid.id_s,:).* ...
% % %         data_NO3.diat_N_lim_Cweight_avg_100m_assm(grid.id_w, grid.id_s,:) .* ...
% % %         mean(1.7.^((data_TEMP.TEMP_assm(grid.id_w, grid.id_s,:)-30.0)./10.0)) .* ...
% % %         mu_ref .* ...
% % %         mean(data_diatC.diatC_assm(grid.id_w, grid.id_s,:)) ...
% % %         ./86400.*10000; % /day -> /sec, 100m integration (cm -> 10000)
% % %     end
% % % end






% % % %% Fe
% % % for loni=1:size(grid.tlong,1)
% % %     for lati=1:size(grid.tlat,2)
% % %         if (data_Fe.(['diat_Fe_lim_Cweight_avg_100m_corr_assm_int_l', tmp.lyear_str])(loni,lati)>0 && ...
% % %                 data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_l', tmp.lyear_str])(loni,lati) >0 && ...
% % %                 data_Fe.(['diat_Fe_lim_Cweight_avg_100m_corr_assm_int_p_l', tmp.lyear_str])(loni,lati)<0.05 && ...
% % %                 data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_p_l', tmp.lyear_str])(loni,lati) <0.05)
% % %             
% % %             [tmp.corr, tmp.corr_p]= corrcoef(data_NPP.(['photoC_TOT_zint_100m_model_l',tmp.lyear_str])(loni,lati,:), ...
% % %                 data_Fe.(['diat_Fe_lim_Cweight_avg_100m_model_l', tmp.lyear_str])(loni,lati,:));
% % %             nut_corrset(loni,lati,1) = tmp.corr(1,2);
% % % %          [tmp.corr, tmp.corr_p]=corrcoef(tmp.data(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));
% % %         else
% % %             nut_corrset(loni,lati,1)=0;
% % %         end
% % %     end
% % % end
% % % 
% % % %% N
% % % for loni=1:size(grid.tlong,1)
% % %     for lati=1:size(grid.tlat,2)
% % %         if (data_NO3.(['diat_N_lim_Cweight_avg_100m_corr_assm_int_l', tmp.lyear_str])(loni,lati)>0 && ...
% % %                 data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_l', tmp.lyear_str])(loni,lati) >0 && ...
% % %                 data_NO3.(['diat_N_lim_Cweight_avg_100m_corr_assm_int_p_l', tmp.lyear_str])(loni,lati)<0.05 && ...
% % %                 data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_p_l', tmp.lyear_str])(loni,lati) <0.05)
% % %             
% % %             [tmp.corr, tmp.corr_p]= corrcoef(data_NPP.(['photoC_TOT_zint_100m_model_l',tmp.lyear_str])(loni,lati,:), ...
% % %                 data_NO3.(['diat_N_lim_Cweight_avg_100m_model_l', tmp.lyear_str])(loni,lati,:));
% % %             nut_corrset(loni,lati,2) = tmp.corr(1,2);
% % % %          [tmp.corr, tmp.corr_p]=corrcoef(tmp.data(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));
% % %         else
% % %             nut_corrset(loni,lati,2)=0;
% % %         end
% % %     end
% % % end
% % % 
% % % 
% % % %% SiO3
% % % for loni=1:size(grid.tlong,1)
% % %     for lati=1:size(grid.tlat,2)
% % %         if (data_SiO3.(['diat_SiO3_lim_Cweight_avg_100m_corr_assm_int_l', tmp.lyear_str])(loni,lati)>0 && ...
% % %                 data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_l', tmp.lyear_str])(loni,lati) >0 && ...
% % %                 data_SiO3.(['diat_SiO3_lim_Cweight_avg_100m_corr_assm_int_p_l', tmp.lyear_str])(loni,lati)<0.05 && ...
% % %                 data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_p_l', tmp.lyear_str])(loni,lati) <0.05)
% % %             
% % %             [tmp.corr, tmp.corr_p]= corrcoef(data_NPP.(['photoC_TOT_zint_100m_model_l',tmp.lyear_str])(loni,lati,:), ...
% % %                 data_SiO3.(['diat_SiO3_lim_Cweight_avg_100m_model_l', tmp.lyear_str])(loni,lati,:));
% % %             nut_corrset(loni,lati,3) = tmp.corr(1,2);
% % % %          [tmp.corr, tmp.corr_p]=corrcoef(tmp.data(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));
% % %         else
% % %             nut_corrset(loni,lati,3)=0;
% % %         end
% % %     end
% % % end
% % % 
% % %             
% % % %% P
% % % for loni=1:size(grid.tlong,1)
% % %     for lati=1:size(grid.tlat,2)
% % %         if (data_PO4.(['diat_P_lim_Cweight_avg_100m_corr_assm_int_l', tmp.lyear_str])(loni,lati)>0 && ...
% % %                 data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_l', tmp.lyear_str])(loni,lati) >0 && ...
% % %                 data_PO4.(['diat_P_lim_Cweight_avg_100m_corr_assm_int_p_l', tmp.lyear_str])(loni,lati)<0.05 && ...
% % %                 data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_p_l', tmp.lyear_str])(loni,lati) <0.05)
% % %             
% % %            [tmp.corr, tmp.corr_p]= corrcoef(data_NPP.(['photoC_TOT_zint_100m_model_l',tmp.lyear_str])(loni,lati,:), ...
% % %                 data_PO4.(['diat_P_lim_Cweight_avg_100m_model_l', tmp.lyear_str])(loni,lati,:));
% % %             nut_corrset(loni,lati,4) = tmp.corr(1,2);
% % % %          [tmp.corr, tmp.corr_p]=corrcoef(tmp.data(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));
% % %         else
% % %             nut_corrset(loni,lati,4)=0;
% % %         end
% % %     end
% % % end
% % % 
% % % %% Light
% % % for loni=1:size(grid.tlong,1)
% % %     for lati=1:size(grid.tlat,2)
% % %         if (data_light.(['diat_light_lim_Cweight_avg_100m_corr_assm_int_l', tmp.lyear_str])(loni,lati)>0 && ...
% % %                 data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_l', tmp.lyear_str])(loni,lati) >0 && ...
% % %                 data_light.(['diat_light_lim_Cweight_avg_100m_corr_assm_int_p_l', tmp.lyear_str])(loni,lati)<0.05 && ...
% % %                 data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_p_l', tmp.lyear_str])(loni,lati) <0.05)
% % %             
% % %            [tmp.corr, tmp.corr_p]= corrcoef(data_NPP.(['photoC_TOT_zint_100m_model_l',tmp.lyear_str])(loni,lati,:), ...
% % %                 data_light.(['diat_light_lim_Cweight_avg_100m_model_l', tmp.lyear_str])(loni,lati,:));
% % %             nut_corrset(loni,lati,5) = tmp.corr(1,2);
% % % %          [tmp.corr, tmp.corr_p]=corrcoef(tmp.data(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));
% % %         else
% % %             nut_corrset(loni,lati,5)=0;
% % %         end
% % %     end
% % % end
% % % 
% % % 
% % % %% TEMP
% % % for loni=1:size(grid.tlong,1)
% % %     for lati=1:size(grid.tlat,2)
% % %         if (data_TEMP.(['TEMP_corr_assm_int_l', tmp.lyear_str])(loni,lati)>0 && ...
% % %                 data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_l', tmp.lyear_str])(loni,lati) >0 && ...
% % %                 data_TEMP.(['TEMP_corr_assm_int_p_l', tmp.lyear_str])(loni,lati)<0.05 && ...
% % %                 data_NPP.(['photoC_TOT_zint_100m_corr_assm_int_p_l', tmp.lyear_str])(loni,lati) <0.05)
% % %             
% % %            [tmp.corr, tmp.corr_p]= corrcoef(data_NPP.(['photoC_TOT_zint_100m_model_l',tmp.lyear_str])(loni,lati,:), ...
% % %                 data_TEMP.(['TEMP_model_l', tmp.lyear_str])(loni,lati,:));
% % %             nut_corrset(loni,lati,6) = tmp.corr(1,2);
% % % %          [tmp.corr, tmp.corr_p]=corrcoef(tmp.data(isfinite(tmp.data_obs)), tmp.data_obs(isfinite(tmp.data_obs)));
% % %         else
% % %             nut_corrset(loni,lati,6)=0;
% % %         end
% % %     end
% % % end