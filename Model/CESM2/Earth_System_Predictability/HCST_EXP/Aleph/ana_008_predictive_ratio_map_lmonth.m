% %  Created 28-Nov-2023 by Yong-Yub Kim
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

% cfg.vars = {'photoC_TOT_zint_100m'};
cfg.vars = {'TS', 'SST', 'PSL', 'PRECT'};
% cfg.vars = {'photoC_TOT_zint_100m'};
cfg.lmonths=[12 13 24 25 36 37 48 49 60];
% cfg.lmonths=60;
% tmp.dimids= [1, 2, 4];
cfg.vlayer=1; % surf, vertical slice 

% cfg.vlayer=1:10; % 10layer. don't put more than 15

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
    % cfg.obs_iyears=f_obs_iyears(cfg.var);
    
    % dirs.hcstroot=['/mnt/lustre/proj/earth.system.predictability/HCST_EXP/archive/', cfg.comp, '/', cfg.var];
    % dirs.obsroot=['/Volumes/kyy_raid/kimyy/Observation/', cfg.obs_name, '/monthly_reg_', cfg.obs_fname_module(2:4)];
    % dirs.figroot=['/Volumes/kyy_raid/kimyy/Figure/CESM2/ESP/HCST_EXP/archive/', cfg.comp,'/', cfg.var];
    % dirs.lens2root=['/mnt/lustre/proj/earth.system.predictability/LENS2/archive_analysis/', cfg.comp, '/', cfg.var];
    % dirs.assmroot=['/mnt/lustre/proj/earth.system.predictability/ASSM_EXP/archive_analysis/', cfg.comp, '/', cfg.var];
    
    dirs.hcstroot=['/mnt/lustre/proj/earth.system.predictability/HCST_EXP/archive_transfer/', cfg.comp, '/', cfg.var];
    dirs.obsroot=['/mnt/lustre/proj/kimyy/Observation/', cfg.obs_name, '/monthly_reg_', cfg.obs_fname_module(2:4)];
    dirs.matroot=['/mnt/lustre/proj/kimyy/Model/CESM2/ESP/HCST_EXP/mat/', cfg.comp,'/', cfg.var];
    dirs.lens2root=['/mnt/lustre/proj/earth.system.predictability/LENS2/archive_transfer/', cfg.comp, '/', cfg.var, '/smbb'];
    dirs.assmroot=['/mnt/lustre/proj/earth.system.predictability/ASSM_EXP/archive_transfer/', cfg.comp, '/', cfg.var];
    
    
    
    cfg.iyears=cfg.obs_iyears;
    cfg.gnm='f09_g17';
%     cfg.proj_year=1;
    cfg.proj_year=5;

    cfg.max_iy=max(cfg.iyears);
    cfg.min_iy=min(cfg.iyears);
    
    cfg.len_t_m = length(cfg.lmonths);
    cfg.len_t_y = length(cfg.iyears);
    % cfg.len_t = cfg.len_t_y *cfg.len_t_m;
    
    % %% grid set(mask from model)
    % tmp.obsname=cfg.obsnames{1};
    % iyear=min(cfg.iyears);
    % cfg.casename_m=[cfg.gridname, '.hcst.', tmp.obsname, '-', cfg.assm_factor, 'p', cfg.ens_member];
    % cfg.casename=[cfg.casename_m, '_i', num2str(iyear)];
    % dirs.datadir= [dirs.hcstroot, filesep, cfg.casename_m, filesep, 'GMSV'];
    
    % f09_g17.hcst.en4.2_ba-10p1_i2021.pop.h.once.nc
    
    % [tmp.error_status, tmp.value]=system(['ls ', dirs.datadir, '/*once*']);  % b.e21.BHISTsmbb.f09_g17.assm.oras4_ba-10p1.pop.h.once.nc
    tmp.gridname = [dirs.hcstroot, tmp.fs, '../grid.nc'];
    tmp.maskname = '/mnt/lustre/proj/kimyy/Model/CESM2/ESP/grids/ocn/RECCAP2_region_masks_all_v20210412_POP2_grid.nc';
    
    % grid.region_mask=ncread(tmp.gridname, 'REGION_MASK'); 
    % grid.ocean_mask=NaN(size(grid.region_mask));
    % grid.ocean_mask(grid.region_mask>0)=1;
    % grid.tarea = ncread(tmp.gridname, 'TAREA');
    
    switch cfg.comp
        case {'ocn', 'ice'}
            grid.tlong=ncread(tmp.gridname, 'TLONG');
            grid.tlat=ncread(tmp.gridname, 'TLAT');
            grid.mask_ocn=ncread(tmp.maskname, 'open_ocean');
            grid.mask_ocn(grid.mask_ocn<-10e10)=NaN;
            grid.mask_ocn=grid.mask_ocn./grid.mask_ocn;
            grid.z_t=ncread(tmp.gridname, 'z_t')./100; % meter
            grid.dz=ncread(tmp.gridname, 'dz')./100; % meter
            grid.dz_res=reshape(grid.dz(cfg.vlayer), [1 1 length(grid.dz(cfg.vlayer))]);
        case {'atm', 'lnd'}
            grid.lon=ncread(tmp.gridname, 'lon');
            grid.lat=ncread(tmp.gridname, 'lat');
            [grid.tlat grid.tlong]=meshgrid(grid.lat, grid.lon);
    end
    
    grid.nlon=size(grid.tlong,1);
    grid.nlat=size(grid.tlat,2);
    % grid.ntime=cfg.proj_year.*12;
    
    % % model filename example
    % /mnt/lustre/proj/kimyy/Model/CESM2/ESP/HCST_EXP/archive/atm/SST/ens_all/ens_all_i2020
    % SST_f09_g17.hcst.ens_all_i2020.cam.h0.2024-12.nc
        
    tmp.varname=cfg.var;

    clear tmp.ydata tmp.ydata_lens2 tmp.ydata_obs tmp.ydata_assm

    tmp.cases=dir([dirs.hcstroot, '/f09*']);
    tmp.len_case=length(tmp.cases);

    tmp.cases_assm=dir([dirs.assmroot, '/f09*']);
    tmp.len_case_assm=length(tmp.cases_assm);

    tmp.cases_lens2=dir([dirs.lens2root, '/LE2*']);
    tmp.len_case_lens2=length(tmp.cases_lens2);
  %% variables initialization (initial only)
    data.([tmp.varname, '_model_ini'])=NaN(grid.nlon, grid.nlat, cfg.len_t_y,tmp.len_case);
    data.([tmp.varname, '_lens2_ini'])=NaN(grid.nlon, grid.nlat, cfg.len_t_y,tmp.len_case_lens2);
    data.([tmp.varname, '_assm_ini'])=NaN(grid.nlon, grid.nlat, cfg.len_t_y,tmp.len_case_lens2);

    for lmind=1:length(cfg.lmonths)
        lmon=cfg.lmonths(lmind);

        tmp.lmon_str=num2str(lmon, '%02i');        

        disp('HCST')

 %% variables initialization (tgt)
        data.([tmp.varname, '_model_tgt_', tmp.lmon_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y, tmp.len_case);
        data.([tmp.varname, '_model_inc_', tmp.lmon_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y, tmp.len_case);
        data2.([tmp.varname, '_model_spr_', tmp.lmon_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
        data2.([tmp.varname, '_model_inc_mean_', tmp.lmon_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
        data2.([tmp.varname, '_model_ratio_', tmp.lmon_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y); 

        for casei=1:tmp.len_case
            tmp.expn=tmp.cases(casei).name;
            disp([tmp.expn, ' lmon: ', num2str(lmon)]);
            for iyear=min(cfg.iyears):max(cfg.iyears)
                iy_ind=iyear-min(cfg.iyears)+1;
                iyear_str=num2str(iyear);
                if mod(lmon,12)==0
                    fyear=iyear+(lmon/12)-1;
                else
                    fyear=iyear+floor(lmon/12);
                end
                fyear_str=num2str(fyear);
                fmon=lmon-(fyear-iyear)*12;
                fmon_str=num2str(fmon,'%02i');
                tmp.fname_ini=[dirs.hcstroot, filesep, tmp.expn,  filesep, tmp.expn,'_i',iyear_str, ...
                    filesep, tmp.varname,'_', tmp.expn,'_i',iyear_str,cfg.obs_fname_module, ...
                    iyear_str, '-', '01', '.nc'];
                tmp.fname_tgt=[dirs.hcstroot, filesep, tmp.expn,  filesep, tmp.expn,'_i',iyear_str, ...
                    filesep, tmp.varname,'_', tmp.expn,'_i',iyear_str,cfg.obs_fname_module, ...
                    fyear_str, '-', fmon_str, '.nc'];
                
                data.([tmp.varname, '_model_ini_', tmp.lmon_str])(:,:,iy_ind,casei)= ...
                    ncread(tmp.fname_ini, tmp.varname);
                data.([tmp.varname, '_model_tgt_', tmp.lmon_str])(:,:,iy_ind,casei)= ...
                    ncread(tmp.fname_tgt, tmp.varname);
%                 ncinfo(tmp.fname)
                
            end
        end
        
                
        data.([tmp.varname, '_model_inc_', tmp.lmon_str])= ...
            data.([tmp.varname, '_model_tgt_', tmp.lmon_str]) - ...
            data.([tmp.varname, '_model_ini_', tmp.lmon_str]);
                data2.([tmp.varname, '_model_inc_mean_', tmp.lmon_str]) = ...
                    mean(data.([tmp.varname, '_model_inc_', tmp.lmon_str]), 4);
        data2.([tmp.varname, '_model_spr_', tmp.lmon_str])= ...
            std(data.([tmp.varname, '_model_tgt_', tmp.lmon_str]),1,4);
        data2.([tmp.varname, '_model_ratio_', tmp.lmon_str]) = ...
            abs(mean(data.([tmp.varname, '_model_inc_', tmp.lmon_str]), 4) ...
            ./ data2.([tmp.varname, '_model_spr_', tmp.lmon_str]));
        data2.([tmp.varname, '_model_ratio_mean_', tmp.lmon_str]) = ...
            mean(data2.([tmp.varname, '_model_ratio_', tmp.lmon_str]), 3);
        


%         tmp.data = squeeze(mean(data.([tmp.varname, '_model_inc_', tmp.lmon_str]), 4));
%         tmp.data = data.([tmp.varname, '_model_spr_', tmp.lmon_str]);
% 
%         pcolor(tmp.data(:,:,10)'); 
% 
%         pcolor(data.([tmp.varname, '_model_ratio_mean_', tmp.lmon_str])'); 
%         shading flat; colorbar; caxis([1 10]);
     
%         fprintf('%02d_%s_%s  ',lyear, ',', cfg.casename_m,'_', tmp.varname); lap_time = tic;
    

%% LENS2
disp('LENS2')

         %% variables initialization (tgt)
        data.([tmp.varname, '_lens2_tgt_', tmp.lmon_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y, tmp.len_case_lens2);
        data.([tmp.varname, '_lens2_inc_', tmp.lmon_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y, tmp.len_case_lens2);
        data2.([tmp.varname, '_lens2_spr_', tmp.lmon_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
        data2.([tmp.varname, '_lens2_inc_mean_', tmp.lmon_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
        data2.([tmp.varname, '_lens2_ratio_', tmp.lmon_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y); 

        for casei=1:tmp.len_case_lens2
            tmp.expn=tmp.cases_lens2(casei).name;
            disp([tmp.expn, ' lmon: ', num2str(lmon)]);
            for iyear=min(cfg.iyears):max(cfg.iyears)
                iy_ind=iyear-min(cfg.iyears)+1;
                iyear_str=num2str(iyear);
                if mod(lmon,12)==0
                    fyear=iyear+(lmon/12)-1;
                else
                    fyear=iyear+floor(lmon/12);
                end
                fyear_str=num2str(fyear);
                fmon=lmon-(fyear-iyear)*12;
                fmon_str=num2str(fmon,'%02i');
                tmp.fname_ini=[dirs.lens2root, filesep, tmp.expn,  filesep,  ...
                     tmp.varname,'_', tmp.expn,cfg.obs_fname_module, ...
                    iyear_str, '-', '01', '.nc'];
                tmp.fname_tgt=[dirs.lens2root, filesep, tmp.expn,  filesep,  ...
                     tmp.varname,'_', tmp.expn,cfg.obs_fname_module, ...
                    fyear_str, '-', fmon_str, '.nc'];
                
                data.([tmp.varname, '_lens2_ini_', tmp.lmon_str])(:,:,iy_ind,casei)= ...
                    ncread(tmp.fname_ini, tmp.varname);
                data.([tmp.varname, '_lens2_tgt_', tmp.lmon_str])(:,:,iy_ind,casei)= ...
                    ncread(tmp.fname_tgt, tmp.varname);
%                 ncinfo(tmp.fname)
                
            end
        end


                
        data.([tmp.varname, '_lens2_inc_', tmp.lmon_str])= ...
            data.([tmp.varname, '_lens2_tgt_', tmp.lmon_str]) - ...
            data.([tmp.varname, '_lens2_ini_', tmp.lmon_str]);
                data2.([tmp.varname, '_lens2_inc_mean_', tmp.lmon_str]) = ...
                    mean(data.([tmp.varname, '_lens2_inc_', tmp.lmon_str]), 4);
        data2.([tmp.varname, '_lens2_spr_', tmp.lmon_str])= ...
            std(data.([tmp.varname, '_lens2_tgt_', tmp.lmon_str]),1,4);
        data2.([tmp.varname, '_lens2_ratio_', tmp.lmon_str]) = ...
            abs(mean(data.([tmp.varname, '_lens2_inc_', tmp.lmon_str]), 4)) ...
            ./ data2.([tmp.varname, '_lens2_spr_', tmp.lmon_str]);
        data2.([tmp.varname, '_lens2_ratio_mean_', tmp.lmon_str]) = ...
            mean(data2.([tmp.varname, '_lens2_ratio_', tmp.lmon_str]), 3);
      
        

%% ASSM
disp('ASSM')
         %% variables initialization (tgt)
        data.([tmp.varname, '_assm_tgt_', tmp.lmon_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y, tmp.len_case_assm);
        data.([tmp.varname, '_assm_inc_', tmp.lmon_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y, tmp.len_case_assm);
        data2.([tmp.varname, '_assm_spr_', tmp.lmon_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
        data2.([tmp.varname, '_assm_inc_mean_', tmp.lmon_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y);
        data2.([tmp.varname, '_assm_ratio_', tmp.lmon_str])=NaN(grid.nlon, grid.nlat, cfg.len_t_y); 

        for casei=1:tmp.len_case_assm
            tmp.expn=tmp.cases_assm(casei).name;
            disp([tmp.expn, ' lmon: ', num2str(lmon)]);
            for iyear=min(cfg.iyears):max(cfg.iyears)
                iy_ind=iyear-min(cfg.iyears)+1;
                iyear_str=num2str(iyear);
                if mod(lmon,12)==0
                    fyear=iyear+(lmon/12)-1;
                else
                    fyear=iyear+floor(lmon/12);
                end
                fyear_str=num2str(fyear);
                fmon=lmon-(fyear-iyear)*12;
                fmon_str=num2str(fmon,'%02i');
                tmp.fname_ini=[dirs.assmroot, filesep, tmp.expn,  filesep,  ...
                     tmp.varname,'_', tmp.expn,cfg.obs_fname_module, ...
                    iyear_str, '-', '01', '.nc'];
                tmp.fname_tgt=[dirs.assmroot, filesep, tmp.expn,  filesep,  ...
                     tmp.varname,'_', tmp.expn,cfg.obs_fname_module, ...
                    fyear_str, '-', fmon_str, '.nc'];
                
                data.([tmp.varname, '_assm_ini_', tmp.lmon_str])(:,:,iy_ind,casei)= ...
                    ncread(tmp.fname_ini, tmp.varname);
                try
                    data.([tmp.varname, '_assm_tgt_', tmp.lmon_str])(:,:,iy_ind,casei)= ...
                    ncread(tmp.fname_tgt, tmp.varname);
                catch
                    data.([tmp.varname, '_assm_tgt_', tmp.lmon_str])(:,:,iy_ind,casei) = NaN;
                end
%                 ncinfo(tmp.fname)
                
            end
        end

        data.([tmp.varname, '_assm_inc_', tmp.lmon_str])= ...
            data.([tmp.varname, '_assm_tgt_', tmp.lmon_str]) - ...
            data.([tmp.varname, '_assm_ini_', tmp.lmon_str]);
        data2.([tmp.varname, '_assm_inc_mean_', tmp.lmon_str]) = ...
                mean(data.([tmp.varname, '_assm_inc_', tmp.lmon_str]), 4);
        data2.([tmp.varname, '_assm_spr_', tmp.lmon_str])= ...
            std(data.([tmp.varname, '_assm_tgt_', tmp.lmon_str]),1,4);
        data2.([tmp.varname, '_assm_ratio_', tmp.lmon_str]) = ...
            abs(mean(data.([tmp.varname, '_assm_inc_', tmp.lmon_str]), 4)) ...
            ./ data2.([tmp.varname, '_assm_spr_', tmp.lmon_str]);
        data2.([tmp.varname, '_assm_ratio_mean_', tmp.lmon_str]) = ...
            mean(data2.([tmp.varname, '_assm_ratio_', tmp.lmon_str]), 3);

        
        
        disp('abc')
        
        if ~exist(dirs.matroot,'dir'), mkdir(dirs.matroot); end
        fig_cfg.mat_name=[dirs.matroot, filesep, 'hcst_pred_ratio_assm_', tmp.varname, ...
            '_v', num2str(cfg.vlayer_1st, '%02i'), '_v', num2str(max(cfg.vlayer), '%02i'), ...
            '_l', tmp.lmon_str, 'm.mat'];
        save(fig_cfg.mat_name, 'data2')
        clear data2
%         fprintf('%7.1f sec\n', toc(lap_time) );

    end

    


% data.([tmp.varname, '_corr_obs', '_l3_4'])= ( data.([tmp.varname, '_corr_obs', '_l03']) + ...
%     data.([tmp.varname, '_corr_obs', '_l04']) ) / 2;
% data.([tmp.varname, '_corr_obs', '_l5_9'])= ( data.([tmp.varname, '_corr_obs', '_l05']) + ...
%     data.([tmp.varname, '_corr_obs', '_l06']) + ..._
%     data.([tmp.varname, '_corr_obs', '_l07']) + ...
%     data.([tmp.varname, '_corr_obs', '_l08']) + ...
%     data.([tmp.varname, '_corr_obs', '_l09']) ) / 5;
% 
% data.([tmp.varname, '_corr_assm', '_l3_4'])= ( data.([tmp.varname, '_corr_assm', '_l03']) + ...
%     data.([tmp.varname, '_corr_assm', '_l04']) ) / 2;
% data.([tmp.varname, '_corr_assm', '_l5_9'])= ( data.([tmp.varname, '_corr_assm', '_l05']) + ...
%     data.([tmp.varname, '_corr_assm', '_l06']) + ...
%     data.([tmp.varname, '_corr_assm', '_l07']) + ...
%     data.([tmp.varname, '_corr_assm', '_l08']) + ...
%     data.([tmp.varname, '_corr_assm', '_l09']) ) / 5;


    
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
