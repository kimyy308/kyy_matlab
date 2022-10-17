% %  Created 13-Sep-2022 by Yong-Yub Kim
clc; clear all; close all;

%% set path
tmp.dropboxpath = '/mnt/lustre/proj/kimyy/Dropbox';
tmp.fs=filesep;
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);

%% model configuration
dirs.root='/mnt/lustre/proj/earth.system.predictability/ASSM_EXP';
dirs.yoshi_root='/proj/yoshi/DATA/CESM2_ODA';
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
% config.varnames={'SALT', 'TEMP', 'DIC', 'PO4', 'ALK', 'NO3', 'SiO3', 'PD', 'NO2'};
config.varnames={'PO4','NO3', 'SiO3', 'PD'};

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
    tmp.varname_obs = f_varname_obs(tmp.varname);
    
    config.savename_obs=[dirs.obssavedir, filesep, 'HOT_', config_obs.staname, '_', ...
        tmp.varname_obs, '_', num2str(min(config.years)), '_', num2str(max(config.years)), '.mat'];
    config.savename_obs_cfg=[dirs.obssavedir,  filesep, 'cfg_HOT_', config_obs.staname, '_', ...
        tmp.varname_obs, '_', num2str(min(config.years)), '_', num2str(max(config.years)), '.mat'];


    OBS.data=ncread(config_obs.filename, tmp.varname_obs);
    OBS.data(OBS.data==-9)=NaN;
    for ctrlvarind=1:length(config_obs.ctrl_var)
        tmp.ctrlvarname=config_obs.ctrl_var{ctrlvarind};
        OBS.ctrl.(tmp.ctrlvarname)=ncread(config_obs.filename, tmp.ctrlvarname);
        tmp.nanind=find(OBS.ctrl.(tmp.ctrlvarname)==-9);
        OBS.ctrl.(tmp.ctrlvarname)(tmp.nanind)=NaN;
    end

%     OBS.ctrl.crn=ncread(config_obs.filename, 'crn');
%     OBS.ctrl.stn=ncread(config_obs.filename, 'stn');
%     OBS.ctrl.cast=ncread(config_obs.filename, 'cast');
%     OBS.ctrl.ros=ncread(config_obs.filename, 'ros');
%     OBS.ctrl.mdate=ncread(config_obs.filename, 'mdate');
%     OBS.ctrl.mtime=ncread(config_obs.filename, 'mtime');
%     OBS.ctrl.press=ncread(config_obs.filename, 'press');
    tmp.obsmask=ones(size(OBS.data));
    tmp.obsmask=tmp.obsmask.* OBS.ctrl.mdate .* OBS.ctrl.mtime .* OBS.data;
    tmp.obsmask(isfinite(tmp.obsmask))=1;
    OBS.ctrl.mask=tmp.obsmask;
    for ctrlvarind=1:length(config_obs.ctrl_var)
        tmp.ctrlvarname=config_obs.ctrl_var{ctrlvarind};
        OBS.ctrl.(tmp.ctrlvarname)=OBS.ctrl.(tmp.ctrlvarname) .* OBS.ctrl.mask;
        
    end
    OBS.data=OBS.data.*OBS.ctrl.mask;

    for yind=1:length(config.years)

        tmp.year=config.years(yind); tmp.yearstr=num2str(tmp.year, '%04i');
        disp(['OBS, ', tmp.varname, ', ', tmp.yearstr]);
    
        for mind=1:length(config.months)
            tmp.month=config.months(mind); tmp.monthstr=num2str(tmp.month, '%02i');
            tmp.mdate_min_str=[tmp.monthstr, '00', tmp.yearstr(end-1:end)];
            tmp.mdate_max_str=[tmp.monthstr, '99', tmp.yearstr(end-1:end)];
            tmp.mdate_min=str2num(tmp.mdate_min_str);
            tmp.mdate_max=str2num(tmp.mdate_max_str);

            tmp.valid_ind=find(mod(OBS.ctrl.mdate,100)==mod(tmp.year,100) ...
                & floor(OBS.ctrl.mdate/1e4)==tmp.month); %%find index for  valid year, month
            tmp.data_obs=OBS.data(tmp.valid_ind);
            if isempty(tmp.valid_ind)~=1
                OBS.data_m(yind,mind)=mean(tmp.data_obs);                
            else
                OBS.data_m(yind,mind)=NaN;
            end
        end
        OBS.data_y(yind)=mean(OBS.data_m(yind,:), 'omitnan');
    end
    save(config.savename_obs, '-struct', 'OBS');
    save(config.savename_obs_cfg, 'OBS', 'config', 'config_obs', 'dirs');
%     plot(OBS.data_m(:)); ylabel(tmp.varname_obs);
end



for varind=1:length(config.varnames)
    tmp.varname=config.varnames{varind};
  
    for obsind=1:length(config.obsnames)
        tmp.obsname=config.obsnames{obsind};
        config.savename=[dirs.saveroot, filesep, config_obs.staname, filesep, 'CESM2_', config_obs.staname, '_', ...
            tmp.obsname, '_', tmp.varname, '_', num2str(min(config.years)), '_', num2str(max(config.years)), '.mat'];
        config.savename_cfg=[dirs.saveroot, filesep, config_obs.staname, filesep, 'cfg_CESM2_', config_obs.staname, '_', ...
            tmp.obsname, '_', tmp.varname, '_', num2str(min(config.years)), '_', num2str(max(config.years)), '.mat'];
        if (exist(config.savename) ==0)
            for ensind=1:length(config.ensnames)
                tmp.ensname=config.ensnames{ensind};
                dirs.datadir=[dirs.archive, filesep, 'b.e21.B', config.scenname, 'smbb.', ...
                    config.gridname, '.assm.', tmp.obsname, '_', tmp.ensname, filesep, ...
                    config.component, filesep, 'hist'];
                dirs.yoshi_datadir=[dirs.yoshi_root, filesep, ...
                    'assm.', tmp.obsname, '_', tmp.ensname];
    
                %% read grid
                if ensind==1
                    [tmp.error_status, tmp.value]=system(['ls ', dirs.datadir, '| grep once']);
                    tmp.gridname = [dirs.datadir, filesep, tmp.value(1:end-1)];
                    CESM2_grid.lon_t=ncread(tmp.gridname, 'TLONG'); 
                    CESM2_grid.lat_t=ncread(tmp.gridname, 'TLAT'); 
                    CESM2_grid.region_mask=ncread(tmp.gridname, 'REGION_MASK'); 
                    CESM2_grid.ocean_mask=NaN(size(CESM2_grid.region_mask));
                    CESM2_grid.ocean_mask(CESM2_grid.region_mask>0)=1;
                    CESM2_grid.z_t=ncread(tmp.gridname, 'z_t');  % (cm)
                    CESM2_grid.dz=ncread(tmp.gridname, 'dz'); 
                    CESM2_grid.len_z_t=length(CESM2_grid.z_t);
                    
                    [CESM2_grid.obs_lon_ind, CESM2_grid.obs_lon_ind, CESM2_grid.obs_lat_ind, CESM2_grid.obs_lat_ind] = ...
                        Func_0012_findind_Y(1, [config_obs.sta_lon, config_obs.sta_lat], ...
                        CESM2_grid.lon_t.*CESM2_grid.ocean_mask, ...
                        CESM2_grid.lat_t.*CESM2_grid.ocean_mask, 'CESM2'); % find valid lon, lat index near station
                %% initialize
%                     CESM2_data.(tmp.varname)(1:config.len_obs, 1:config.len_ens, 1:config.len_t_y,1:config.len_t_m,1:CESM2_grid.len_z_t)=NaN;
%                     CESM2_data.([tmp.varname, '_vm'])(1:config.len_obs, 1:config.len_ens, 1:config.len_t_y,1:config.len_t_m)=NaN;
%                     CESM2_data.([tmp.varname, '_vm_200m'])(1:config.len_obs, 1:config.len_ens, 1:config.len_t_y,1:config.len_t_m)=NaN;
                    CESM2_data.(tmp.varname)(1:config.len_ens, 1:config.len_t_y,1:config.len_t_m,1:CESM2_grid.len_z_t)=NaN;
                    CESM2_data.([tmp.varname, '_vm'])(1:config.len_ens, 1:config.len_t_y,1:config.len_t_m)=NaN;
                    CESM2_data.([tmp.varname, '_vm_200m'])(1:config.len_ens, 1:config.len_t_y,1:config.len_t_m)=NaN;

                end
                %% read data
                for yind=1:length(config.years)
    
                    tmp.year=config.years(yind); tmp.yearstr=num2str(tmp.year, '%04i');
                    disp([tmp.varname, ', ', tmp.obsname, ', ', tmp.ensname, ', ', tmp.yearstr]);
    
                    for mind=1:length(config.months)
                        tmp.month=config.months(mind); tmp.monthstr=num2str(tmp.month, '%02i');
%                         tmp.filename=[dirs.datadir, filesep, 'b.e21.B', config.scenname, 'smbb.', ...
%                             config.gridname, '.assm.', config.obsname, '_', config.ensname, ...
%                             '.pop.h.', tmp.yearstr, '-', tmp.monthstr, '.nc'];
                        [tmp.error_status, tmp.value]=system(['ls ', dirs.yoshi_datadir, ...
                            '/*pop*', tmp.yearstr, '-', tmp.monthstr, '.* | grep ', 'h.', tmp.yearstr, '-', tmp.monthstr]);
%                         if ~isempty(tmp.value) && strcmp('NO2', tmp.varname)~=1
                        if strcmp(tmp.value(end-3:end-1), '.nc')==1 && strcmp('NO2', tmp.varname)~=1

%                             tmp.filename = [dirs.datadir, filesep, tmp.value(1:end-1)];
                            tmp.filename = tmp.value(1:end-1);
                    %         ncinfo([tmp.filename]);
    
                            tmp.data= ...
                                squeeze(ncread(tmp.filename, tmp.varname, ...
                                [CESM2_grid.obs_lon_ind, CESM2_grid.obs_lat_ind, 1, 1], ...
                                [1, 1, inf, inf]));
                            CESM2_data.(tmp.varname)(ensind,yind,mind,:)= tmp.data;
                            tmp.vmask=NaN(CESM2_grid.len_z_t,1);
                            tmp.vmask(isfinite(squeeze(tmp.data)))=1;
                            tmp.data_vm=sum(tmp.data.*CESM2_grid.dz.*tmp.vmask, 'omitnan') / sum(CESM2_grid.dz.*tmp.vmask, 'omitnan'); % vertical weighted mean
                            tmp.zind_200m=20;
                            tmp.data_vm_200m=sum(tmp.data(1:tmp.zind_200m).*CESM2_grid.dz(1:tmp.zind_200m).*tmp.vmask(1:tmp.zind_200m), 'omitnan') ...
                                / sum(CESM2_grid.dz(1:tmp.zind_200m).*tmp.vmask(1:tmp.zind_200m), 'omitnan'); % vertical weighted mean (~200m)
                            CESM2_data.([tmp.varname, '_vm'])(ensind,yind,mind) = tmp.data_vm;
                            CESM2_data.([tmp.varname, '_vm_200m'])(ensind,yind,mind) = tmp.data_vm_200m;
    
                        else
                            tmp.filename = NaN;
                            CESM2_data.(tmp.varname)(ensind,yind,mind,:)= NaN;
                            CESM2_data.([tmp.varname, '_vm'])(ensind,yind,mind) = NaN;
                            CESM2_data.([tmp.varname, '_vm_200m'])(ensind,yind,mind) = NaN;
                        end
                    end
                end
            end
            save(config.savename, '-struct', 'CESM2_data');
            save(config.savename_cfg, 'CESM2_data', 'CESM2_grid', 'config', 'config_obs', 'dirs');
            clear CESM2_data
        end
        
    end
    load(config.savename)

end





% tmp.dd=squeeze(mean(CESM2_data.DIC_vm_200m,4));
% plot(tmp.dd)

%     


% Salinity
tmp.data=squeeze(mean(SALT_vm_200m(1,1,:,:),4));
plot(config.years,tmp.data-mean(tmp.data,'omitnan'), 'color', [0.8 0.8 0.8]); 
config.savename_obs=[dirs.obssavedir, filesep, 'HOT_', config_obs.staname, '_', ...
        'bsal', '_', num2str(min(config.years)), '_', num2str(max(config.years)), '.mat'];
load(config.savename_obs);
config.savename_obs_cfg=[dirs.obssavedir, filesep, 'cfg_HOT_', config_obs.staname, '_', ...
        'bsal', '_', num2str(min(config.years)), '_', num2str(max(config.years)), '.mat'];
load(config.savename_obs_cfg);

hold on
tmp.data2=squeeze(mean(SALT_vm_200m(2,1,:,:),4));
plot(config.years,tmp.data2-mean(tmp.data2,'omitnan'),'y');
tmp.data3=squeeze(mean(SALT_vm_200m(3,1,:,:),4));
plot(config.years,tmp.data3-mean(tmp.data3,'omitnan'), 'g');
tmp.data_mean=(tmp.data+tmp.data2+tmp.data3)/3.0;
plot(config.years,tmp.data_mean-mean(tmp.data_mean,'omitnan'), 'k');
plot(config.years,OBS.data_y(:)-mean(OBS.data_y,'omitnan'), 'b'); 
xlabel('year'); ylabel('SALT(bsal)'); legend('oras4', 'projd', 'en4', 'ens-mean', 'obs');
set(gca, 'fontsize', 20)
grid minor
hold off



% Nitrate
tmp.data=squeeze(mean(NO3_vm_200m(1,1,:,:),4));
plot(config.years,tmp.data-mean(tmp.data,'omitnan'), 'color', [0.8 0.8 0.8]); 
config.savename_obs=[dirs.obssavedir, filesep, 'HOT_', config_obs.staname, '_', ...
        'lln', '_', num2str(min(config.years)), '_', num2str(max(config.years)), '.mat'];
load(config.savename_obs);
% OBS_nit=OBS;
% config.savename_obs=[dirs.obssavedir, filesep, 'HOT_', config_obs.staname, '_', ...
%         'no2', '_', num2str(min(config.years)), '_', num2str(max(config.years)), '.mat'];
% load(config.savename_obs);
% OBS_no2=OBS;
config.savename_obs_cfg=[dirs.obssavedir, filesep, 'cfg_HOT_', config_obs.staname, '_', ...
        'lln', '_', num2str(min(config.years)), '_', num2str(max(config.years)), '.mat'];
load(config.savename_obs_cfg);
% OBS.data_y=OBS_nit.data_y-OBS_no2.data_y;

hold on
tmp.data2=squeeze(mean(NO3_vm_200m(2,1,:,:),4));
plot(config.years,tmp.data2-mean(tmp.data2,'omitnan'),'y');
tmp.data3=squeeze(mean(NO3_vm_200m(3,1,:,:),4));
plot(config.years,tmp.data3-mean(tmp.data3,'omitnan'), 'g');
tmp.data_mean=(tmp.data+tmp.data2+tmp.data3)/3.0;
plot(config.years,tmp.data_mean-mean(tmp.data_mean,'omitnan'), 'k');
plot(config.years,OBS.data_y(:)-mean(OBS.data_y,'omitnan'), 'b'); 
xlabel('year'); ylabel('NO3(lln)'); legend('oras4', 'projd', 'en4', 'ens-mean', 'obs');
set(gca, 'fontsize', 20)
grid minor
hold off

% Phosphate
tmp.data=squeeze(mean(PO4_vm_200m(1,1,:,:),4));
plot(config.years,tmp.data-mean(tmp.data,'omitnan'), 'color', [0.8 0.8 0.8]); 
config.savename_obs=[dirs.obssavedir, filesep, 'HOT_', config_obs.staname, '_', ...
        'llp', '_', num2str(min(config.years)), '_', num2str(max(config.years)), '.mat'];
load(config.savename_obs);

config.savename_obs_cfg=[dirs.obssavedir, filesep, 'cfg_HOT_', config_obs.staname, '_', ...
        'llp', '_', num2str(min(config.years)), '_', num2str(max(config.years)), '.mat'];
load(config.savename_obs_cfg);

hold on
tmp.data2=squeeze(mean(PO4_vm_200m(2,1,:,:),4));
plot(config.years,tmp.data2-mean(tmp.data2,'omitnan'),'y');
tmp.data3=squeeze(mean(PO4_vm_200m(3,1,:,:),4));
plot(config.years,tmp.data3-mean(tmp.data3,'omitnan'), 'g');
tmp.data_mean=(tmp.data+tmp.data2+tmp.data3)/3.0;
plot(config.years,tmp.data_mean-mean(tmp.data_mean,'omitnan'), 'k');
plot(config.years,(OBS.data_y(:)-mean(OBS.data_y,'omitnan'))/1000, 'b'); 
xlabel('year'); ylabel('PO4(llp)'); legend('oras4', 'projd', 'en4', 'ens-mean', 'obs');
set(gca, 'fontsize', 20)
grid minors
hold off




function var_obs = f_varname_obs(var_model)
    switch var_model
        case 'SALT'
            var_obs='bsal';  % PSS-78       
        case 'TEMP'
            var_obs='theta'; % ITS-90
        case 'DIC' % mmol/m^3
            var_obs='dic'; %umol/kg
        case 'PO4' %mmol/m^3
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



