% %  Created 05-Oct-2022 by Yong-Yub Kim

clc; clear all; close all;

%% set path
tmp.dropboxpath = '/Volumes/kyy_raid/kimyy/Dropbox';
tmp.fs=filesep;
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);

%% model configuration
dirs.root='/Volumes/kyy_raid/earth.system.predictability/ASSM_EXP';
dirs.yoshi_root='/proj/yaoshi/DATA/CESM2_ODA';
dirs.archive=[dirs.root, filesep, 'archive'];
dirs.saveroot='/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/ASSM_EXP';

config.years=1960:2021;
config.months=1:12;
config.scenname='HIST';
config.gridname='f09_g17';
% config.obsnames={'oras4', 'projdv7.3', 'en4.2'};
config.obsnames={'projdv7.3', 'en4.2'};

% config.obsnames={'oras4'};

% config.ensnames={'ba-10p1'};

% config.ensnames={'ba-10p1', 'ba-10p2', 'ba-10p3', 'ba-10p4', 'ba-10p5'};

config.ensnames={'ba-10p1', 'ba-10p2', 'ba-10p3', 'ba-10p4', 'ba-10p5', ...
    'ba-20p1', 'ba-20p2', 'ba-20p3', 'ba-20p4', 'ba-20p5'};

config.component='ocn';
config.varnames={'SALT', 'TEMP', 'DIC', 'PO4', 'ALK', 'NO3', 'SiO3', 'PD'};
% config.varnames={'SALT', 'TEMP', 'DIC', 'PO4', 'ALK', 'NO3', 'SiO3', 'PD', 'NO2'};

% config.varnames={'PO4','NO3', 'SiO3', 'PD'};

config.len_t_y = length(config.years);
config.len_t_m = length(config.months);
config.len_obs= length(config.obsnames);
config.len_ens= length(config.ensnames);



% system(['ls ', dirs.datadir])
% b.e21.BHISTsmbb.f09_g17.assm.oras4_ba-10p1.pop.h.once.nc

% RCM_info.name={'test2127', 'test2128', 'test2129', 'test2130', 'test2131', 'ens2201'};

%% observation configuration
config_obs.staname = 'ALOHA';
config_obs.sta_lon = 360.0-158.0;
config_obs.sta_lat = 22.0+45.0/60.0;

dirs.obsroot = '/Volumes/kyy_raid/kimyy/Observation/HOT/ALOHA/Bottle';
dirs.obssavedir = [dirs.obsroot, filesep, 'mat'];
mkdir(dirs.obssavedir)
config_obs.filename = [dirs.obsroot, filesep, 'ALOHA_0_200m_depth.nc'];
config_obs.ctrl_var={'crn','stn', 'cast', 'ros', 'mdate', 'mtime', 'press'};

ncinfo(config_obs.filename)


for varind=1:length(config.varnames)
    tmp.varname=config.varnames{varind};

    tmp.varname_obs = f_varname_obs(tmp.varname);
    config.savename_obs=[dirs.obssavedir, filesep, 'HOT_', config_obs.staname, '_', ...
        tmp.varname_obs, '_', num2str(min(config.years)), '_', num2str(max(config.years)), '.mat'];
    config.savename_obs_cfg=[dirs.obssavedir,  filesep, 'cfg_HOT_', config_obs.staname, '_', ...
        tmp.varname_obs, '_', num2str(min(config.years)), '_', num2str(max(config.years)), '.mat'];
 
    load(config.savename_obs_cfg, 'OBS');


    config.savename=[dirs.saveroot, filesep, config_obs.staname, filesep, 'CESM2_', config_obs.staname, '_', ...
        tmp.varname, '_', num2str(min(config.years)), '_', num2str(max(config.years)), '.mat'];
    config.savename_cfg=[dirs.saveroot, filesep, config_obs.staname, filesep, 'cfg_CESM2_', config_obs.staname, '_', ...
        tmp.varname, '_', num2str(min(config.years)), '_', num2str(max(config.years)), '.mat'];

    load(config.savename_cfg, 'CESM2_data', 'CESM2_grid');

%     tmp.data=reshape(permute(CESM2_data.([tmp.varname,'_vm_200m']),[1 2 4 3]),[3 10 744]);
%     projd_ens10m=squeeze(mean(tmp.data(2,1:5,:),2));
%     projd_ens20m=squeeze(mean(tmp.data(2,6:10,:),2));
%     en4_ens10m=squeeze(mean(tmp.data(3,1:5,:),2));
%     en4_ens20m=squeeze(mean(tmp.data(3,6:10,:),2));
%     obs_data_m=OBS.data_m(:);
% 
%     plot(projd_ens10m, '-r^');
%     hold on
%     plot(projd_ens20m, '-r*');
%     hold off
   
    tmp.data=reshape(permute(CESM2_data.([tmp.varname,'_vm_200m']),[1 2 4 3]),[3 10 744]);
    switch tmp.varname_obs
        case 'llp' % meq/m^3
            tmp.data=tmp.data.*1000;
        case 'lln' %mmol/m^3
            tmp.data=tmp.data.*1000;
        case 'sigma' %g/cm^3
            tmp.data=tmp.data.*1000-1000;

    end 
    
    projd_ensm=squeeze(mean(tmp.data(2,:,:),2));
    
    en4_ensm=squeeze(mean(tmp.data(3,:,:),2));
    
    obs_data_m=OBS.data_m';
    obs_data_m=obs_data_m(:);

    [projd_ensm_m, projd_ensm_std, projd_ensm_ano, projd_ensm_ano_norm] = get_ano_norm(projd_ensm);
    [en4_ensm_m, en4_ensm_std, en4_ensm_ano, en4_ensm_ano_norm] = get_ano_norm(en4_ensm);
    [obs_data_m_m, obs_data_m_std, obs_data_m_ano, obs_data_m_ano_norm] = get_ano_norm(obs_data_m);

    time_m=min(config.years):1/12:max(config.years)+11/12;
    plot(time_m,projd_ensm_ano_norm, '-r', 'linewidth', 3);
    hold on
    plot(time_m,en4_ensm_ano_norm, '-g', 'linewidth', 3);
    plot(time_m,obs_data_m_ano_norm, 'k*');

    %%transparent each ensemble plot
    val_transparent = 0.5;
    rgb_tr_projd=get_transparent_rgb([1,0,0], val_transparent);
    rgb_tr_en4=get_transparent_rgb([0,1,0], val_transparent);
    for ensi=1:size(tmp.data,2)/2
        [projd_ens_m, projd_ens_std, projd_ens_ano, projd_ens_ano_norm] = get_ano_norm(squeeze(tmp.data(2,ensi,:)));
        plot(time_m, projd_ens_ano_norm, 'color', rgb_tr_projd);
        [en4_ens_m, en4_ens_std, en4_ens_ano, en4_ens_ano_norm] = get_ano_norm(squeeze(tmp.data(3,ensi,:)));
        plot(time_m, en4_ens_ano_norm, 'color', rgb_tr_en4);
    end
    val_transparent = 0.2;
    rgb_tr_projd=get_transparent_rgb([1,0,0], val_transparent);
    rgb_tr_en4=get_transparent_rgb([0,1,0], val_transparent);
    for ensi=size(tmp.data,2)/2+1:size(tmp.data,2)
        [projd_ens_m, projd_ens_std, projd_ens_ano, projd_ens_ano_norm] = get_ano_norm(squeeze(tmp.data(2,ensi,:)));
        plot(time_m, projd_ens_ano_norm, 'color', rgb_tr_projd);
        [en4_ens_m, en4_ens_std, en4_ens_ano, en4_ens_ano_norm] = get_ano_norm(squeeze(tmp.data(3,ensi,:)));
        plot(time_m, en4_ens_ano_norm, 'color', rgb_tr_en4);
    end
    plot(time_m,projd_ensm_ano_norm, '-r', 'linewidth', 3);
    plot(time_m,en4_ensm_ano_norm, '-g', 'linewidth', 3);
    ind_valid=find(isfinite(obs_data_m_ano_norm));
    tmp.coef=corrcoef(en4_ensm_ano_norm(ind_valid), obs_data_m_ano_norm(ind_valid));
    coef_en4=tmp.coef(1,2);
    tmp.coef=corrcoef(projd_ensm_ano_norm(ind_valid), obs_data_m_ano_norm(ind_valid));
    coef_projd=tmp.coef(1,2);

    xlabel('Year'); ylabel([tmp.varname, '(', tmp.varname_obs, ')']); 
    legend(['JMA7.3, R=', num2str(round(coef_projd,2))],...
        ['EN4.2, R=', num2str(round(coef_en4,2))],...
        'Obs(ALOHA)', 'Location', 'NorthWest');
    set(gca, 'fontsize', 20)
    grid minor
    hold off
    
    dirs.figdir='/Volumes/kyy_raid/kimyy/Figure/ESP/ASSM_EXP/ALOHA';
    system(['mkdir -p ', dirs.figdir])
    figname=[dirs.figdir, filesep, tmp.varname, '_ano_norm_', tmp.varname_obs, '_', ...
        num2str(min(config.years)), '_', num2str(max(config.years)), '.tif'];
    set(gcf, 'PaperPosition', [0, 0, 8, 4]);

   saveas(gcf,figname,'tif');

%% raw plot
    close all;
    plot(time_m,projd_ensm, '-r', 'linewidth', 3);
    hold on
    plot(time_m,en4_ensm, '-g', 'linewidth', 3);
    plot(time_m,obs_data_m, 'k*');

    %%transparent each ensemble plot
    val_transparent = 0.5;
    rgb_tr_en4=get_transparent_rgb([0,1,0], val_transparent);
    rgb_tr_projd=get_transparent_rgb([1,0,0], val_transparent);

    for ensi=1:size(tmp.data,2)/2
        [projd_ens_m, projd_ens_std, projd_ens_ano, projd_ens_ano_norm] = get_ano_norm(squeeze(tmp.data(2,ensi,:)));
        plot(time_m, squeeze(tmp.data(2,ensi,:)), 'color', rgb_tr_projd);
        [en4_ens_m, en4_ens_std, en4_ens_ano, en4_ens_ano_norm] = get_ano_norm(squeeze(tmp.data(3,ensi,:)));
        plot(time_m, squeeze(tmp.data(3,ensi,:)), 'color', rgb_tr_en4);
    end
    val_transparent = 0.2;
    rgb_tr_projd=get_transparent_rgb([1,0,0], val_transparent);
    rgb_tr_en4=get_transparent_rgb([0,1,0], val_transparent);
    for ensi=size(tmp.data,2)/2+1:size(tmp.data,2)
        [projd_ens_m, projd_ens_std, projd_ens_ano, projd_ens_ano_norm] = get_ano_norm(squeeze(tmp.data(2,ensi,:)));
        plot(time_m, squeeze(tmp.data(2,ensi,:)), 'color', rgb_tr_projd);
        [en4_ens_m, en4_ens_std, en4_ens_ano, en4_ens_ano_norm] = get_ano_norm(squeeze(tmp.data(3,ensi,:)));
        plot(time_m, squeeze(tmp.data(3,ensi,:)), 'color', rgb_tr_en4);
    end

    plot(time_m,projd_ensm, '-r', 'linewidth', 3);
    plot(time_m,en4_ensm, '-g', 'linewidth', 3);
    ind_valid=find(isfinite(obs_data_m_ano_norm));
    tmp.coef=corrcoef(en4_ensm_ano_norm(ind_valid), obs_data_m_ano_norm(ind_valid));
    coef_en4=tmp.coef(1,2);
    tmp.coef=corrcoef(projd_ensm_ano_norm(ind_valid), obs_data_m_ano_norm(ind_valid));
    coef_projd=tmp.coef(1,2);

    xlabel('Year'); ylabel([tmp.varname, '(', tmp.varname_obs, ')']); 
    legend(['JMA7.3, R=', num2str(round(coef_projd,2))],...
        ['EN4.2, R=', num2str(round(coef_en4,2))],...
        'Obs(ALOHA)', 'Location', 'NorthWest');
    set(gca, 'fontsize', 20)
    grid minor
    hold off

    figname=[dirs.figdir, filesep, tmp.varname, '_raw_', tmp.varname_obs, '_', ...
        num2str(min(config.years)), '_', num2str(max(config.years)), '.tif'];
    set(gcf, 'PaperPosition', [0, 0, 8, 4]);

   saveas(gcf,figname,'tif');



end


function var_obs = f_varname_obs(var_model)
%bsal; bottle salinity
%theta; potential temperature
%dic; dissolved inorganic carbon, DIC
%phos; phosphate
%alk; alkaliniy
%nit; Nitrate + Nitrite
%sil; Silicate
%sigma; potential density
    switch var_model
        case 'SALT'
            var_obs='bsal';  % PSS-78       
        case 'TEMP'
            var_obs='theta'; % ITS-90
        case 'DIC' % mmol/m^3
            var_obs='dic'; %umol/kg
        case 'PO4' %mmol/m^3ㅇㅅ
            var_obs='phos';  % standard method with low accuracy, umol/kg
%             var_obs='llp';  %high precision, nmol/kg
        case 'ALK' % meq/m^3
            var_obs='alk'; % 'ueq/kg'
        case 'NO3' %mmol/m^3
%             var_obs='nit'; %(NO2 + NO3 in ALOHA data) umol/kg
            var_obs='lln'; %(NO2 + NO3 in ALOHA data) nmol/kg
        case 'SiO3' %mmol/m^3
            var_obs='sil';  %umol/kg
        case 'PD' %g/cm^3
            var_obs='sigma';  %kg/m3
        case 'NO2'
            var_obs='no2'; %nmol/kg
    end
end


function [data_m, data_std, data_ano, data_ano_norm] = get_ano_norm(data)
    data_m=mean(data,'omitnan');
    data_std=std(data,'omitnan');
    data_ano=data-data_m;
    data_ano_norm=data_ano/data_std;
end

function rgb_transparent = get_transparent_rgb(rgb, val_transparent)
    %[1,0,0]; % red
    %[1, 165/255, 0]; orange
    %[0,1,0]; % lime green
    %[0,0,1]; % blue
    rgb_hsv = rgb2hsv(rgb);
    rgb_hsv(:,2) =  val_transparent;
    rgb_transparent = hsv2rgb(rgb_hsv);
end


% Lomb-Scargle transform?
% 
% tmptmp(1:12)=0;
% for i=1:744
%     if isfinite(obs_data_m(i))
%         tmptmp(mod(i-1,12)+1)=tmptmp(mod(i-1,12)+1)+1;
%     end
% end