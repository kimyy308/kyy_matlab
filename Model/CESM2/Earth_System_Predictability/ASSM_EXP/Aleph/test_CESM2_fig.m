% %  Created 13-Sep-2022 by Yong-Yub Kim
clc; clear all; close all;

%% set path
tmp.dropboxpath = '/mnt/lustre/proj/kimyy/Dropbox';
tmp.fs=filesep;
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);

%% model configuration
dirs.root='/mnt/lustre/proj/earth.system.predictability/ASSM_EXP';
dirs.yoshi_root='/proj/yaoshi/DATA/CESM2_ODA';
dirs.archive=[dirs.root, filesep, 'archive'];
dirs.saveroot='/mnt/lustre/proj/kimyy/Model/CESM2/ESP/ASSM_EXP';

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
config.varnames={'SALT', 'TEMP', 'DIC', 'PO4', 'ALK', 'NO3', 'SiO3', 'PD', 'NO2'};
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

dirs.obsroot = '/mnt/lustre/proj/kimyy/Observation/HOT/ALOHA/Bottle';
dirs.obssavedir = [dirs.obsroot, filesep, 'mat'];
mkdir(dirs.obssavedir)
config_obs.filename = [dirs.obsroot, filesep, 'ALOHA_0_200m_depth.nc'];
config_obs.ctrl_var={'crn','stn', 'cast', 'ros', 'mdate', 'mtime', 'press'};

ncinfo(config_obs.filename)

%bsal; bottle salinity
%theta; potential temperature
%dic; dissolved inorganic carbon, DIC
%phos; phosphate
%alk; alkaliniy
%nit; Nitrate + Nitrite
%sil; Silicate
%sigma; potential density
for varind=1:length(config.varnames)
    tmp.varname=config.varnames{varind};

    switch tmp.varname
        case 'SALT'
            tmp.varname_obs='bsal';  % PSS-78       
        case 'TEMP'
            tmp.varname_obs='theta'; % ITS-90
        case 'DIC' % mmol/m^3
            tmp.varname_obs='dic'; %umol/kg
        case 'PO4' %mmol/m^3
            tmp.varname_obs='phos';  % standard method with low accuracy, umol/kg
%             tmp.varname_obs='llp';  %high precision, nmol/kg
        case 'ALK' % meq/m^3
            tmp.varname_obs='alk'; % 'ueq/kg'
        case 'NO3' %mmol/m^3
%             tmp.varname_obs='nit'; %(NO2 + NO3 in ALOHA data) umol/kg
            tmp.varname_obs='lln'; %(NO2 + NO3 in ALOHA data) nmol/kg
        case 'SiO3' %mmol/m^3
            tmp.varname_obs='sil';  %umol/kg
        case 'PD' %g/cm^3
            tmp.varname_obs='sigma';  %kg/m3
        case 'NO2'
            tmp.varname_obs='no2'; %nmol/kg
    end
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
    projd_ensm=squeeze(mean(tmp.data(2,:,:),2));
    en4_ensm=squeeze(mean(tmp.data(3,:,:),2));
    projd_ensm_m=mean(projd_ensm,'omitnan');
    projd_ensm_std=std(projd_ensm,'omitnan');
    projd_ensm_ano=projd_ensm-projd_ensm_m;
    projd_ensm_ano_norm=projd_ensm_ano/projd_ensm_std;
    en4_ensm_m=mean(en4_ensm,'omitnan');
    en4_ensm_std=std(en4_ensm,'omitnan');
    en4_ensm_ano=en4_ensm-en4_ensm_m;
    en4_ensm_ano_norm=en4_ensm_ano/en4_ensm_std;
    obs_data_m=OBS.data_m';
    obs_data_m=obs_data_m(:);
    obs_data_m_m=mean(obs_data_m, 'omitnan');
    obs_data_m_std=std(obs_data_m, 'omitnan');
    obs_data_m_ano=obs_data_m-obs_data_m_m;
    obs_data_m_ano_norm=obs_data_m_ano/obs_data_m_std;
    time_m=min(config.years):1/12:max(config.years)+11/12;
    plot(time_m,projd_ensm_ano_norm, '-r', 'linewidth', 3);
    hold on
    plot(time_m,en4_ensm_ano_norm, '-g', 'linewidth', 3);
    plot(time_m,obs_data_m_ano_norm, 'k*');
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
    
    dirs.figdir='/mnt/lustre/proj/kimyy/Figure/ESMP/ASSM_EXP/ALOHA';
    figname=[dirs.figdir, filesep, tmp.varname, '_', tmp.varname_obs, '_', ...
        num2str(min(config.years)), '_', num2str(max(config.years)), '.tif'];
    set(gcf, 'PaperPosition', [0, 0, 8, 4]);

   saveas(gcf,figname,'tif');

end


% tmp.dd=squeeze(mean(CESM2_data.DIC_vm_200m,4));
% plot(tmp.dd)

%     



% 
% % Salinity
% tmp.data=squeeze(mean(SALT_vm_200m(1,1,:,:),4));
% plot(config.years,tmp.data-mean(tmp.data,'omitnan'), 'color', [0.8 0.8 0.8]); 
% config.savename_obs=[dirs.obssavedir, filesep, 'HOT_', config_obs.staname, '_', ...
%         'bsal', '_', num2str(min(config.years)), '_', num2str(max(config.years)), '.mat'];
% load(config.savename_obs);
% config.savename_obs_cfg=[dirs.obssavedir, filesep, 'cfg_HOT_', config_obs.staname, '_', ...
%         'bsal', '_', num2str(min(config.years)), '_', num2str(max(config.years)), '.mat'];
% load(config.savename_obs_cfg);
% 
% hold on
% tmp.data2=squeeze(mean(SALT_vm_200m(2,1,:,:),4));
% plot(config.years,tmp.data2-mean(tmp.data2,'omitnan'),'y');
% tmp.data3=squeeze(mean(SALT_vm_200m(3,1,:,:),4));
% plot(config.years,tmp.data3-mean(tmp.data3,'omitnan'), 'g');
% tmp.data_mean=(tmp.data+tmp.data2+tmp.data3)/3.0;
% plot(config.years,tmp.data_mean-mean(tmp.data_mean,'omitnan'), 'k');
% plot(config.years,OBS.data_y(:)-mean(OBS.data_y,'omitnan'), 'b'); 
% xlabel('year'); ylabel('SALT(bsal)'); legend('oras4', 'projd', 'en4', 'ens-mean', 'obs');
% set(gca, 'fontsize', 20)
% grid minor
% hold off
% 
% 
% 
% % Nitrate
% tmp.data=squeeze(mean(NO3_vm_200m(1,1,:,:),4));
% plot(config.years,tmp.data-mean(tmp.data,'omitnan'), 'color', [0.8 0.8 0.8]); 
% config.savename_obs=[dirs.obssavedir, filesep, 'HOT_', config_obs.staname, '_', ...
%         'lln', '_', num2str(min(config.years)), '_', num2str(max(config.years)), '.mat'];
% load(config.savename_obs);
% % OBS_nit=OBS;
% % config.savename_obs=[dirs.obssavedir, filesep, 'HOT_', config_obs.staname, '_', ...
% %         'no2', '_', num2str(min(config.years)), '_', num2str(max(config.years)), '.mat'];
% % load(config.savename_obs);
% % OBS_no2=OBS;
% config.savename_obs_cfg=[dirs.obssavedir, filesep, 'cfg_HOT_', config_obs.staname, '_', ...
%         'lln', '_', num2str(min(config.years)), '_', num2str(max(config.years)), '.mat'];
% load(config.savename_obs_cfg);
% % OBS.data_y=OBS_nit.data_y-OBS_no2.data_y;
% 
% hold on
% tmp.data2=squeeze(mean(NO3_vm_200m(2,1,:,:),4));
% plot(config.years,tmp.data2-mean(tmp.data2,'omitnan'),'y');
% tmp.data3=squeeze(mean(NO3_vm_200m(3,1,:,:),4));
% plot(config.years,tmp.data3-mean(tmp.data3,'omitnan'), 'g');
% tmp.data_mean=(tmp.data+tmp.data2+tmp.data3)/3.0;
% plot(config.years,tmp.data_mean-mean(tmp.data_mean,'omitnan'), 'k');
% plot(config.years,OBS.data_y(:)-mean(OBS.data_y,'omitnan'), 'b'); 
% xlabel('year'); ylabel('NO3(lln)'); legend('oras4', 'projd', 'en4', 'ens-mean', 'obs');
% set(gca, 'fontsize', 20)
% grid minor
% hold off
% 
% % Phosphate
% tmp.data=squeeze(mean(PO4_vm_200m(1,1,:,:),4));
% plot(config.years,tmp.data-mean(tmp.data,'omitnan'), 'color', [0.8 0.8 0.8]); 
% config.savename_obs=[dirs.obssavedir, filesep, 'HOT_', config_obs.staname, '_', ...
%         'llp', '_', num2str(min(config.years)), '_', num2str(max(config.years)), '.mat'];
% load(config.savename_obs);
% 
% config.savename_obs_cfg=[dirs.obssavedir, filesep, 'cfg_HOT_', config_obs.staname, '_', ...
%         'llp', '_', num2str(min(config.years)), '_', num2str(max(config.years)), '.mat'];
% load(config.savename_obs_cfg);
% 
% hold on
% tmp.data2=squeeze(mean(PO4_vm_200m(2,1,:,:),4));
% plot(config.years,tmp.data2-mean(tmp.data2,'omitnan'),'y');
% tmp.data3=squeeze(mean(PO4_vm_200m(3,1,:,:),4));
% plot(config.years,tmp.data3-mean(tmp.data3,'omitnan'), 'g');
% tmp.data_mean=(tmp.data+tmp.data2+tmp.data3)/3.0;
% plot(config.years,tmp.data_mean-mean(tmp.data_mean,'omitnan'), 'k');
% plot(config.years,(OBS.data_y(:)-mean(OBS.data_y,'omitnan'))/1000, 'b'); 
% xlabel('year'); ylabel('PO4(llp)'); legend('oras4', 'projd', 'en4', 'ens-mean', 'obs');
% set(gca, 'fontsize', 20)
% grid minors
% hold off