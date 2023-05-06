% %  Created 02-Feb-2023 by Yong-Yub Kim lead month
% %  Created 02-Feb-2023 by Yong-Yub Kim lead month (calculate climatology using evey initialized set)

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
dirs.hcstroot=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/HCST_EXP/archive/atm/', cfg.var];
dirs.obsroot=['/Volumes/kyy_raid/kimyy/Reanalysis/ERA5/', cfg.var, '/monthly_reg_cam'];
dirs.figroot=['/Volumes/kyy_raid/kimyy/Figure/CESM2/ESP/HCST_EXP/archive/atm/', cfg.var];
cfg.iyears=1970  :2020;
cfg.months=1:12;
% cfg.scenname='HIST';
cfg.gnm='f09_g17';
% cfg.assm_factor='10';
% cfg.ens_member='1';
cfg.proj_year=5;
% cfg.obsnames={'en4.2_ba'};
% cfg.ensnames={'ba-10p1'};
cfg.component='atm';
cfg.varnames={'PSL'};
cfg.len_t_y = length(cfg.iyears);
cfg.len_t_m = length(cfg.months);
cfg.len_t= cfg.len_t_y * cfg.len_t_m;
% cfg.len_obs= length(cfg.obsnames);
% cfg.len_ens= length(cfg.ensnames);
cfg.region1 = 'TBVAI';
cfg.region2 = 'TBVCP';
cfg.region = 'TBV';

cfg.max_lead_month = 36;

%% grid set(mask from model)
% tmp.gridname = [dirs.hcstroot, tmp.fs, 'grid.nc'];
% 
% grid.lon=ncread(tmp.gridname, 'lon');
% grid.lat=ncread(tmp.gridname, 'lat');
% [grid.tlat grid.tlong]=meshgrid(grid.lat, grid.lon);
% 
% grid.nlon=size(grid.tlong,1);
% grid.nlat=size(grid.tlat,2);
%% read & plot data
varn.varname=cfg.var;
for iyear=min(cfg.iyears):max(cfg.iyears)
    tmp.iyear_str=num2str(iyear, '%04i');
    %% set casename
    cfg.casename_m=['ens_all'];
    cfg.casename=[cfg.casename_m, '_i', tmp.iyear_str];

    %% read data
    for ind_r = 1:2
        tmp.ir= num2str(ind_r);
        dirs.(['datadir',tmp.ir]) = [dirs.hcstroot, filesep, ...
            cfg.casename_m, filesep, cfg.casename, tmp.fs, cfg.(['region',tmp.ir])];
        varn.(['modelvar',tmp.ir]) = [cfg.(['region',tmp.ir]), ...
            'm_', varn.varname, '_model_',  'i', tmp.iyear_str];
        varn.(['obsvar',tmp.ir]) = [cfg.(['region',tmp.ir]), ...
            'm_', varn.varname, '_obs_',  'i', tmp.iyear_str];

        for fy = iyear:iyear+cfg.proj_year-1
            tmp.fy=fy;
            tmp.fy_str=num2str(tmp.fy, '%04i');
            for mon=1:12
                tmp.mon_str=num2str(mon, '%02i');
                cfg.(['datafilename',tmp.ir])=[dirs.(['datadir',tmp.ir]), filesep, ...
                    'M_', cfg.(['region',tmp.ir]), '_', varn.varname, '_', cfg.gnm, '.hcst.', cfg.casename, '.cam.h0.', tmp.fy_str, '-', tmp.mon_str, '.nc'];

                data.(varn.(['modelvar',tmp.ir]))((fy-iyear)*12+mon) = ...
                    ncread(cfg.(['datafilename',tmp.ir]), varn.varname);
                data.time_tmp((fy-iyear)*12+mon)=ncread(cfg.datafilename1, 'time');
                cfg.(['obsfnm',tmp.ir]) = [dirs.obsroot, tmp.fs, cfg.(['region',tmp.ir]), tmp.fs, ...
                    'M_', cfg.(['region',tmp.ir]), '_', 'ERA5_', cfg.obs_var, '_reg_cesm2.',tmp.fy_str, tmp.mon_str, '.nc'];

                if  exist(cfg.(['obsfnm',tmp.ir]))==0
                    data.(varn.(['obsvar',tmp.ir]))((fy-iyear)*12+mon)=NaN;
                else
                    data.(varn.(['obsvar',tmp.ir]))((fy-iyear)*12+mon)=ncread(cfg.(['obsfnm',tmp.ir]), 'msl');
                end
            end
        end
    end


    %% get climatology
    varn.obsvar1_clim=[cfg.region1, 'm_', varn.varname, '_obs_clim'];
    varn.obsvar2_clim=[cfg.region2, 'm_', varn.varname, '_obs_clim'];
    varn.modelvar1_clim=[cfg.region1, 'm_', varn.varname, '_model_clim'];
    varn.modelvar2_clim=[cfg.region2, 'm_', varn.varname, '_model_clim'];
    for mon=1:12
        tmp.tind_clim=(0:12:48) + mon;
        data.(varn.modelvar1_clim)(mon,iyear-min(cfg.iyears)+1)=mean(data.(varn.modelvar1)(tmp.tind_clim), 'omitnan');
        data.(varn.modelvar2_clim)(mon,iyear-min(cfg.iyears)+1)=mean(data.(varn.modelvar2)(tmp.tind_clim), 'omitnan');
        data.(varn.obsvar1_clim)(mon,iyear-min(cfg.iyears)+1)=mean(data.(varn.obsvar1)(tmp.tind_clim), 'omitnan');
        data.(varn.obsvar2_clim)(mon,iyear-min(cfg.iyears)+1)=mean(data.(varn.obsvar2)(tmp.tind_clim), 'omitnan');  
    end
    
    %% read variables    
%     data.time_tmp=ncread(cfg.datafilename, 'time');
    data.time_tmp_leap=daynoleap2datenum(data.time_tmp, 0)-1; % reflection of leap year and -1 for month correction
    data.time_vec_tmp=datevec(data.time_tmp_leap);
    data.time_vec_tmp(:,1)=data.time_vec_tmp(:,1)+iyear; % y correction (0~4 + iyear)
    data.time_tmp_leap=datenum(data.time_vec_tmp); % recal time var
    %% assign variables according to lead year (model)
    for lmonth=0:cfg.max_lead_month
        tmp.lmonth_str=num2str(lmonth, '%03i');
        varn.modelvar1_lm=[cfg.region1, 'm_', varn.varname, '_model',  '_l', tmp.lmonth_str];
        varn.modelvar2_lm=[cfg.region2, 'm_', varn.varname, '_model',  '_l', tmp.lmonth_str];
        tmp.mdata1= data.(varn.modelvar1)(lmonth+1);
        tmp.mdata2= data.(varn.modelvar2)(lmonth+1);
        tmp.mind=(iyear-min(cfg.iyears))+1;
        data.(varn.modelvar1_lm)(tmp.mind)=tmp.mdata1;
        data.(varn.modelvar1_lm)(data.(varn.modelvar1_lm)==0)=NaN;
        data.(varn.modelvar2_lm)(tmp.mind)=tmp.mdata2;
        data.(varn.modelvar2_lm)(data.(varn.modelvar2_lm)==0)=NaN;
    end
    %% get merged time series for obs, model
    tmp.yind=(iyear-min(cfg.iyears))*12+1:(iyear-min(cfg.iyears))*12+12;
    tmp.ydata1= data.(varn.obsvar1)(1:12);
    tmp.ydata2= data.(varn.obsvar2)(1:12);
    tmp.ydata1_model= data.(varn.modelvar1)(1:12);
    tmp.ydata2_model= data.(varn.modelvar2)(1:12);
    data.([cfg.region1, 'm_', varn.varname, '_obs'])(tmp.yind)= tmp.ydata1;
    data.([cfg.region2, 'm_', varn.varname, '_obs'])(tmp.yind)= tmp.ydata2;
    if iyear~=max(cfg.iyears)
        data.([cfg.region1, 'm_', varn.varname, '_model'])(tmp.yind)= tmp.ydata1_model;
        data.([cfg.region2, 'm_', varn.varname, '_model'])(tmp.yind)= tmp.ydata2_model;
    else
        tmp.yind2=(iyear-min(cfg.iyears))*12+1:(iyear-min(cfg.iyears))*12+12*cfg.proj_year;
        data.([cfg.region1, 'm_', varn.varname, '_model'])(tmp.yind2) = data.(varn.modelvar1);
        data.([cfg.region2, 'm_', varn.varname, '_model'])(tmp.yind2) = data.(varn.modelvar2);
    end
    
%     tmp.ydata= ncread(cfg.datafilename, ['nino34m_assm_', tmp.varname], [1], 12);
%     data.([cfg.region, 'm_', tmp.varname, '_assm_', tmp.obsname_simple])(tmp.yind)= tmp.ydata;
    data.time_leap(tmp.yind)=data.time_tmp_leap(1:12);
    data.time_vec(tmp.yind,:)= data.time_vec_tmp(1:12,:);
%     [f_exp, gof_exp] = fit(trendtime_yearly',squeeze(cmems_sla_yearly(i,j,:)),'poly1');
%                             cmems_sla_yearly_poly1_fit(i,j,1:length(trendtime_yearly))=f_exp(trendtime_yearly);
end


disp('abc')
%% %% assign variables according to lead year (obs)
for iyear=min(cfg.iyears):max(cfg.iyears)
    tmp.iyear_str=num2str(iyear, '%04i');
    varn.varname=cfg.var;
    for lmonth=0:cfg.max_lead_month
        tmp.lmonth_str=num2str(lmonth, '%03i');
        tmp.yind=(iyear-min(cfg.iyears))+1;
        varn.obsvar1_lm=[cfg.region1, 'm_', varn.varname, '_obs_',  '_l', tmp.lmonth_str];
        varn.obsvar1=[cfg.region1, 'm_', varn.varname, '_obs'];
        varn.obsvar2_lm=[cfg.region2, 'm_', varn.varname, '_obs_',  '_l', tmp.lmonth_str];
        varn.obsvar2=[cfg.region2, 'm_', varn.varname, '_obs'];
%             tmp.assmvar_lm=[cfg.region, 'm_', tmp.varname, '_assm_',  '_l', tmp.lmonth_str];
%             tmp.assmvar=[cfg.region, 'm_', tmp.varname, '_assm_', tmp.obsname_simple];
        tmp.mind=(iyear-min(cfg.iyears))*12+lmonth+1;
        if tmp.mind<=length(data.(varn.obsvar1))
            tmp.ydata_obs1=data.(varn.obsvar1)(tmp.mind);
            tmp.ydata_obs2=data.(varn.obsvar2)(tmp.mind);
%             tmp.ydata_assm=data.(tmp.assmvar)(tmp.mind);
            data.(varn.obsvar1_lm)(tmp.yind)=tmp.ydata_obs1;
            data.(varn.obsvar2_lm)(tmp.yind)=tmp.ydata_obs2;
%             data.(tmp.assmvar_lm)(tmp.yind)=tmp.ydata_assm;
        else
            data.(varn.obsvar1_lm)(tmp.yind)=NaN;
            data.(varn.obsvar2_lm)(tmp.yind)=NaN;
%             data.(tmp.assmvar_lm)(tmp.yind)=NaN;
        end
    end
end

%% get fitted data and climatological data first to get anomaly
%% data_fit = data-data_det
% varn.obsvar1_all=[cfg.region1, 'm_', varn.varname, '_obs'];
% varn.obsvar2_all=[cfg.region2, 'm_', varn.varname, '_obs'];
% varn.modelvar1_all=[cfg.region1, 'm_', varn.varname, '_model'];
% varn.modelvar2_all=[cfg.region2, 'm_', varn.varname, '_model'];
% 
% varn.obsvar1_fit=[cfg.region1, 'm_', varn.varname, '_obs_fit'];
% varn.obsvar2_fit=[cfg.region2, 'm_', varn.varname, '_obs_fit'];
% varn.modelvar1_fit=[cfg.region1, 'm_', varn.varname, '_model_fit'];
% varn.modelvar2_fit=[cfg.region2, 'm_', varn.varname, '_model_fit'];
% 
% data.(varn.obsvar1_fit) = data.(varn.obsvar1_all) - Func_0028_detrend_linear_1d(data.(varn.obsvar1_all)');
% data.(varn.obsvar2_fit) = data.(varn.obsvar2_all) - Func_0028_detrend_linear_1d(data.(varn.obsvar2_all)');
% data.(varn.modelvar1_fit) = data.(varn.modelvar1_all) - Func_0028_detrend_linear_1d(data.(varn.modelvar1_all)');
% data.(varn.modelvar2_fit) = data.(varn.modelvar2_all) - Func_0028_detrend_linear_1d(data.(varn.modelvar2_all)');
% 
% tmp.tind1=length(data.(varn.obsvar1_all))+1;
% tmp.tind2=length(data.(varn.modelvar1_all));
% 
% data.(varn.obsvar1_fit)(tmp.tind1:tmp.tind2)=NaN;
% data.(varn.obsvar2_fit)(tmp.tind1:tmp.tind2)=NaN;


%% get climatology
% varn.obsvar1_clim=[cfg.region1, 'm_', varn.varname, '_obs_clim'];
% varn.obsvar2_clim=[cfg.region2, 'm_', varn.varname, '_obs_clim'];
% varn.modelvar1_clim=[cfg.region1, 'm_', varn.varname, '_model_clim'];
% varn.modelvar2_clim=[cfg.region2, 'm_', varn.varname, '_model_clim'];
% 
% tmp.ld_obs=length(data.(varn.obsvar1_all));
% tmp.ld=length(data.(varn.modelvar1_all));
% 
% data.(varn.obsvar1_clim)=reshape(data.(varn.obsvar1_all), [12, tmp.ld_obs/12]);
% data.(varn.obsvar1_clim)=mean(data.(varn.obsvar1_clim),2,'omitnan');
% 
% data.(varn.obsvar1_clim)=reshape(data.(varn.obsvar1_all)-data.(varn.obsvar1_fit)(tmp.ld_obs), [12, tmp.ld_obs/12]);
% data.(varn.obsvar2_clim)=reshape(data.(varn.obsvar2_all)-data.(varn.obsvar2_fit)(tmp.ld_obs), [12, tmp.ld_obs/12]);
% data.(varn.modelvar1_clim)=reshape(data.(varn.modelvar1_all)-data.(varn.modelvar1_fit), [12, tmp.ld/12]);
% data.(varn.modelvar2_clim)=reshape(data.(varn.modelvar2_all)-data.(varn.modelvar2_fit), [12, tmp.ld/12]);
% 
% data.(varn.obsvar1_clim)=mean(data.(varn.obsvar1_clim),2,'omitnan');
% data.(varn.obsvar2_clim)=mean(data.(varn.obsvar2_clim),2,'omitnan');
% data.(varn.modelvar1_clim)=mean(data.(varn.modelvar1_clim),2,'omitnan');
% data.(varn.modelvar2_clim)=mean(data.(varn.modelvar2_clim),2,'omitnan');





%% make anomaly for nino index
for lmonth=0:cfg.max_lead_month
    tmp.lmonth_str=num2str(lmonth, '%03i');
    tmp.avgwin_min=find(cfg.iyears==1996)-floor(lmonth/12);
    tmp.avgwin_max=find(cfg.iyears==2020)-floor(lmonth/12);
    
    tmp.tind_fit=(cfg.iyears-min(cfg.iyears))*12+1 + lmonth;
    tmp.tind_clim=mod(lmonth+1,12);
    if (tmp.tind_clim==0) tmp.tind_clim=12; end
    %% var define
    varn.modelvar1_lm=[cfg.region1, 'm_', varn.varname, '_model',  '_l', tmp.lmonth_str];
    varn.obsvar1_lm=[cfg.region1, 'm_', varn.varname, '_obs_',  '_l', tmp.lmonth_str];
    varn.modelvar2_lm=[cfg.region2, 'm_', varn.varname, '_model',  '_l', tmp.lmonth_str];
    varn.obsvar2_lm=[cfg.region2, 'm_', varn.varname, '_obs_',  '_l', tmp.lmonth_str];
%     tmp.assmvar_lm=[cfg.region, 'm_', tmp.varname, '_assm_',  '_l', tmp.lmonth_str];
    varn.modelvar_ano_lm=[cfg.region, 'm_', varn.varname, '_model_ano_',  '_l', tmp.lmonth_str];
    varn.obsvar_ano_lm=[cfg.region, 'm_', varn.varname, '_obs_ano_',  '_l', tmp.lmonth_str];
    varn.modelvar1_ano_lm=[cfg.region1, 'm_', varn.varname, '_model_ano_',  '_l', tmp.lmonth_str];
    varn.obsvar1_ano_lm=[cfg.region1, 'm_', varn.varname, '_obs_ano_',  '_l', tmp.lmonth_str];
    varn.modelvar2_ano_lm=[cfg.region2, 'm_', varn.varname, '_model_ano_',  '_l', tmp.lmonth_str];
    varn.obsvar2_ano_lm=[cfg.region2, 'm_', varn.varname, '_obs_ano_',  '_l', tmp.lmonth_str];
%     tmp.assmvar_ano_lm=[cfg.region, 'm_', tmp.varname, '_assm_ano_',  '_l', tmp.lmonth_str];

    %% remove climate signal
   
    data.(varn.modelvar1_ano_lm)=data.(varn.modelvar1_lm) ...
        - data.(varn.modelvar1_clim)(tmp.tind_clim,:);
    data.(varn.modelvar2_ano_lm)=data.(varn.modelvar2_lm) ...
        - data.(varn.modelvar2_clim)(tmp.tind_clim,:);
    data.(varn.obsvar1_ano_lm)=data.(varn.obsvar1_lm) ...
        - data.(varn.obsvar1_clim)(tmp.tind_clim,:);
    data.(varn.obsvar2_ano_lm)=data.(varn.obsvar2_lm) ...
        - data.(varn.obsvar2_clim)(tmp.tind_clim,:);

    %% detrending
    data.(varn.modelvar1_ano_lm)= Func_0028_detrend_linear_1d(data.(varn.modelvar1_ano_lm)','omitnan');
    data.(varn.modelvar2_ano_lm)= Func_0028_detrend_linear_1d(data.(varn.modelvar2_ano_lm)','omitnan');
    data.(varn.obsvar1_ano_lm)= Func_0028_detrend_linear_1d(data.(varn.obsvar1_ano_lm)','omitnan');
    data.(varn.obsvar2_ano_lm)= Func_0028_detrend_linear_1d(data.(varn.obsvar2_ano_lm)','omitnan');
    


%     data.(tmp.assmvar_ano_lm)=data.(tmp.assmvar_lm)-tmp.assm_ltm;
    data.(varn.modelvar_ano_lm)=data.(varn.modelvar1_ano_lm) - data.(varn.modelvar2_ano_lm);
    data.(varn.obsvar_ano_lm)=data.(varn.obsvar1_ano_lm) - data.(varn.obsvar2_ano_lm);
    varn.modelvar_ano_lm_rm=[cfg.region, 'm_', varn.varname, '_model_ano_',  '_l', tmp.lmonth_str, '_rm'];
    varn.obsvar_ano_lm_rm=[cfg.region, 'm_', varn.varname, '_obs_ano_',  '_l', tmp.lmonth_str, '_rm'];
    varn.modelvar1_ano_lm_rm=[cfg.region1, 'm_', varn.varname, '_model_ano_',  '_l', tmp.lmonth_str, '_rm'];
    varn.obsvar1_ano_lm_rm=[cfg.region1, 'm_', varn.varname, '_obs_ano_',  '_l', tmp.lmonth_str, '_rm'];
    varn.modelvar2_ano_lm_rm=[cfg.region2, 'm_', varn.varname, '_model_ano_',  '_l', tmp.lmonth_str, '_rm'];
    varn.obsvar2_ano_lm_rm=[cfg.region2, 'm_', varn.varname, '_obs_ano_',  '_l', tmp.lmonth_str, '_rm'];
%     tmp.assmvar_ano_lm_rm=[cfg.region, 'm_', tmp.varname, '_assm_ano_',  '_l', tmp.lmonth_str, '_rm'];
    
    tmp.rm_window1=1;
    tmp.rm_window=[0 tmp.rm_window1-1]; % forward
    data.(varn.modelvar_ano_lm_rm)=movmean(data.(varn.modelvar_ano_lm), tmp.rm_window, 'Endpoints', 'fill');
    data.(varn.obsvar_ano_lm_rm)=movmean(data.(varn.obsvar_ano_lm), tmp.rm_window, 'Endpoints', 'fill');
    data.(varn.modelvar1_ano_lm_rm)=movmean(data.(varn.modelvar1_ano_lm), tmp.rm_window, 'Endpoints', 'fill');
    data.(varn.obsvar1_ano_lm_rm)=movmean(data.(varn.obsvar1_ano_lm), tmp.rm_window, 'Endpoints', 'fill');
    data.(varn.modelvar2_ano_lm_rm)=movmean(data.(varn.modelvar2_ano_lm), tmp.rm_window, 'Endpoints', 'fill');
    data.(varn.obsvar2_ano_lm_rm)=movmean(data.(varn.obsvar2_ano_lm), tmp.rm_window, 'Endpoints', 'fill');
%     data.(tmp.assmvar_ano_lm_rm)=movmean(data.(tmp.assmvar_ano_lm), tmp.rm_window, 'Endpoints', 'fill');
end

% plot(data.('IODm_SST_model_ano__l000_rm'))
% hold on
% plot(data.('IODm_SST_model_ano__l000'))
% hold off
% 

% plot(data.('IODem_SST_obs__l008'))
% 
% plot(data.('IODm_SST_obs_ano__l000_rm'))
% hold on
% plot(data.('IODem_SST_obs_ano__l000'))
% hold off



%% set time value for figure
data.time_vec_extended=data.time_vec;
for lmonth=1:cfg.proj_year-1
    tmp.endind=size(data.time_vec_extended,1);
    data.time_vec_extended(tmp.endind+1:tmp.endind+12,1)=max(cfg.iyears)+lmonth;
    data.time_vec_extended(tmp.endind+1:tmp.endind+12,2)=data.time_vec_extended(1:12,2);
    data.time_vec_extended(tmp.endind+1:tmp.endind+12,3)=data.time_vec_extended(1:12,3);
end
data.time_leap_extended=datenum(data.time_vec_extended);



%% correlation coefficient (hcst&obs), region ------------------
% for varind=1:length(cfg.varnames)
    varn.varname=cfg.var;
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, cfg.region, filesep, 'Timeseries'];
    system(['mkdir -p ', dirs.figdir]);
%         for lmonth=0:cfg.max_lead_month
    for lmonth=0:cfg.max_lead_month
        tmp.lmonth_str=num2str(lmonth, '%03i');
        varn.modelvar_ano_lm_rm=[cfg.region, 'm_', varn.varname, '_model_ano_',  '_l', tmp.lmonth_str, '_rm'];
        varn.obsvar_ano_lm_rm=[cfg.region, 'm_', varn.varname, '_obs_ano_',  '_l', tmp.lmonth_str, '_rm'];
%         tmp.assmvar_ano_lm_rm=[cfg.region, 'm_', tmp.varname, '_assm_ano_',  '_l', tmp.lmonth_str, '_rm'];
        
        tmp.md=data.(varn.modelvar_ano_lm_rm);
        tmp.obs=data.(varn.obsvar_ano_lm_rm);
%         tmp.assm=data.(tmp.assmvar_ano_lm_rm);
        tmp.tind= (cfg.iyears-min(cfg.iyears))*12+1+lmonth;
        tmp.corrcoef=corrcoef(tmp.md(isfinite(tmp.obs)), tmp.obs(isfinite(tmp.obs)));
        data.([cfg.region, 'm_', varn.varname, '_corr'])(lmonth+1)=tmp.corrcoef(1,2);
    end
    %% correlation coefficient of NINO3.4 ind, function of lead year (3-mon running mean)
        cfg.figname=[dirs.figdir, filesep, cfg.region, 'm_hcst_obs_70-20m', '_corr_leadm_3rm_v2', '.tif'];
%             bar(0:cfg.max_lead_month,data.([cfg.region, 'm_', varn.varname, '_corr']),  'linewidth', 2)
        bar(1:cfg.max_lead_month,data.([cfg.region, 'm_', varn.varname, '_corr'])(1:cfg.max_lead_month),  'linewidth', 2)
        xlabel('Lead month'); ylabel(['corr. coef.,', cfg.region , ', obs']);
        set(gca, 'fontsize', 20)
        grid minor
        ylim([-0.4 1])
        saveas(gcf,cfg.figname,'tif');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;
    

%% correlation coefficient (hcst&obs), region1 ------------------
% for varind=1:length(cfg.varnames)
    varn.varname=cfg.var;
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, cfg.region, filesep, 'Timeseries'];
    system(['mkdir -p ', dirs.figdir]);
%         for lmonth=0:cfg.max_lead_month
    for lmonth=0:cfg.max_lead_month
        tmp.lmonth_str=num2str(lmonth, '%03i');
        varn.modelvar1_ano_lm_rm=[cfg.region1, 'm_', varn.varname, '_model_ano_',  '_l', tmp.lmonth_str, '_rm'];
        varn.obsvar1_ano_lm_rm=[cfg.region1, 'm_', varn.varname, '_obs_ano_',  '_l', tmp.lmonth_str, '_rm'];
%         tmp.assmvar_ano_lm_rm=[cfg.region, 'm_', tmp.varname, '_assm_ano_',  '_l', tmp.lmonth_str, '_rm'];
        
        tmp.md=data.(varn.modelvar1_ano_lm_rm);
        tmp.obs=data.(varn.obsvar1_ano_lm_rm);
%         tmp.assm=data.(tmp.assmvar_ano_lm_rm);
        tmp.tind= (cfg.iyears-min(cfg.iyears))*12+1+lmonth;
        tmp.corrcoef=corrcoef(tmp.md(isfinite(tmp.obs)), tmp.obs(isfinite(tmp.obs)));
        data.([cfg.region1, 'm_', varn.varname, '_corr'])(lmonth+1)=tmp.corrcoef(1,2);
    end
    %% correlation coefficient of NINO3.4 ind, function of lead year (3-mon running mean)
        cfg.figname=[dirs.figdir, filesep, cfg.region1, 'm_hcst_obs_70-20m', '_corr_leadm_3rm_v2', '.tif'];
%             bar(0:cfg.max_lead_month,data.([cfg.region, 'm_', varn.varname, '_corr']),  'linewidth', 2)
        bar(1:cfg.max_lead_month,data.([cfg.region1, 'm_', varn.varname, '_corr'])(1:cfg.max_lead_month),  'linewidth', 2)
        xlabel('Lead month'); ylabel(['corr. coef.,', cfg.region1 , ', obs']);
        set(gca, 'fontsize', 20)
        grid minor
        ylim([-0.4 1])
        saveas(gcf,cfg.figname,'tif');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;


%% correlation coefficient (hcst&obs), region2 ------------------
% for varind=1:length(cfg.varnames)
    varn.varname=cfg.var;
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, cfg.region, filesep, 'Timeseries'];
    system(['mkdir -p ', dirs.figdir]);
%         for lmonth=0:cfg.max_lead_month
    for lmonth=0:cfg.max_lead_month
        tmp.lmonth_str=num2str(lmonth, '%03i');
        varn.modelvar2_ano_lm_rm=[cfg.region2, 'm_', varn.varname, '_model_ano_',  '_l', tmp.lmonth_str, '_rm'];
        varn.obsvar2_ano_lm_rm=[cfg.region2, 'm_', varn.varname, '_obs_ano_',  '_l', tmp.lmonth_str, '_rm'];
%         tmp.assmvar_ano_lm_rm=[cfg.region, 'm_', tmp.varname, '_assm_ano_',  '_l', tmp.lmonth_str, '_rm'];
        
        tmp.md=data.(varn.modelvar2_ano_lm_rm);
        tmp.obs=data.(varn.obsvar2_ano_lm_rm);
%         tmp.assm=data.(tmp.assmvar_ano_lm_rm);
        tmp.tind= (cfg.iyears-min(cfg.iyears))*12+1+lmonth;
        tmp.corrcoef=corrcoef(tmp.md(isfinite(tmp.obs)), tmp.obs(isfinite(tmp.obs)));
        data.([cfg.region2, 'm_', varn.varname, '_corr'])(lmonth+1)=tmp.corrcoef(1,2);
    end
    %% correlation coefficient of NINO3.4 ind, function of lead year (3-mon running mean)
        cfg.figname=[dirs.figdir, filesep, cfg.region2, 'm_hcst_obs_70-20m', '_corr_leadm_3rm_v2', '.tif'];
%             bar(0:cfg.max_lead_month,data.([cfg.region, 'm_', varn.varname, '_corr']),  'linewidth', 2)
        bar(1:cfg.max_lead_month,data.([cfg.region2, 'm_', varn.varname, '_corr'])(1:cfg.max_lead_month),  'linewidth', 2)
        xlabel('Lead month'); ylabel(['corr. coef.,', cfg.region2 , ', obs']);
        set(gca, 'fontsize', 20)
        grid minor
        ylim([-0.4 1])
        saveas(gcf,cfg.figname,'tif');
        RemoveWhiteSpace([], 'file', cfg.figname);
        close all;



% end
%% spaghetti plot (region) by lead month -----------
system(['mkdir -p ', dirs.figdir]);
%         for lmonth=0:cfg.max_lead_month
for lmonth=0:cfg.max_lead_month
    tmp.lmonth_str=num2str(lmonth, '%03i');
    varn.modelvar_ano_lm_rm=[cfg.region, 'm_', varn.varname, '_model_ano_',  '_l', tmp.lmonth_str, '_rm'];
    varn.obsvar_ano_lm_rm=[cfg.region, 'm_', varn.varname, '_obs_ano_',  '_l', tmp.lmonth_str, '_rm'];
    
    tmp.md=data.(varn.modelvar_ano_lm_rm);
    tmp.obs=data.(varn.obsvar_ano_lm_rm);
%         tmp.tind= (cfg.iyears-min(cfg.iyears))*12+1+lmonth;
    tmp.plottime=cfg.iyears+lmonth*1/12;
    plot(tmp.plottime, tmp.md, 'b');
    hold on
    plot(tmp.plottime, tmp.obs, 'k', 'linewidth',2)
end
hold off
axis tight
    
cfg.figname=[dirs.figdir, filesep, 'spaghetti_', cfg.region, 'm_hcst_obs_70-20m', '_corr_leadm_3rm_v2', '.tif'];
xlabel('Year'); ylabel('()');
set(gca, 'fontsize', 20)
grid minor
set(gcf, 'PaperPosition', [0, 0, 8, 4]);
saveas(gcf,cfg.figname,'tif');
RemoveWhiteSpace([], 'file', cfg.figname);
close all;




%% spaghetti plot (region1) by lead month -----------
system(['mkdir -p ', dirs.figdir]);
%         for lmonth=0:cfg.max_lead_month
for lmonth=0:cfg.max_lead_month
    tmp.lmonth_str=num2str(lmonth, '%03i');
    varn.modelvar1_ano_lm_rm=[cfg.region1, 'm_', varn.varname, '_model_ano_',  '_l', tmp.lmonth_str, '_rm'];
    varn.obsvar1_ano_lm_rm=[cfg.region1, 'm_', varn.varname, '_obs_ano_',  '_l', tmp.lmonth_str, '_rm'];
    
    tmp.md=data.(varn.modelvar1_ano_lm_rm);
    tmp.obs=data.(varn.obsvar1_ano_lm_rm);
%         tmp.tind= (cfg.iyears-min(cfg.iyears))*12+1+lmonth;
    tmp.plottime=cfg.iyears+lmonth*1/12;
    plot(tmp.plottime, tmp.md, 'b');
    hold on
    plot(tmp.plottime, tmp.obs, 'k', 'linewidth',2)
end
hold off
axis tight
    
cfg.figname=[dirs.figdir, filesep, 'spaghetti_', cfg.region1, 'm_hcst_obs_70-20m', '_corr_leadm_3rm_v2', '.tif'];
xlabel('Year'); ylabel('()');
set(gca, 'fontsize', 20)
grid minor
set(gcf, 'PaperPosition', [0, 0, 8, 4]);
saveas(gcf,cfg.figname,'tif');
RemoveWhiteSpace([], 'file', cfg.figname);
close all;



%% spaghetti plot (region2) by lead month -----------
system(['mkdir -p ', dirs.figdir]);
%         for lmonth=0:cfg.max_lead_month
for lmonth=0:cfg.max_lead_month
    tmp.lmonth_str=num2str(lmonth, '%03i');
    varn.modelvar2_ano_lm_rm=[cfg.region2, 'm_', varn.varname, '_model_ano_',  '_l', tmp.lmonth_str, '_rm'];
    varn.obsvar2_ano_lm_rm=[cfg.region2, 'm_', varn.varname, '_obs_ano_',  '_l', tmp.lmonth_str, '_rm'];
    
    tmp.md=data.(varn.modelvar2_ano_lm_rm);
    tmp.obs=data.(varn.obsvar2_ano_lm_rm);
%         tmp.tind= (cfg.iyears-min(cfg.iyears))*12+1+lmonth;
    tmp.plottime=cfg.iyears+lmonth*1/12;
    plot(tmp.plottime, tmp.md, 'b');
    hold on
    plot(tmp.plottime, tmp.obs, 'k', 'linewidth',2)
end
hold off
axis tight
    
cfg.figname=[dirs.figdir, filesep, 'spaghetti_', cfg.region2, 'm_hcst_obs_70-20m', '_corr_leadm_3rm_v2', '.tif'];
xlabel('Year'); ylabel('()');
set(gca, 'fontsize', 20)
grid minor
set(gcf, 'PaperPosition', [0, 0, 8, 4]);
saveas(gcf,cfg.figname,'tif');
RemoveWhiteSpace([], 'file', cfg.figname);
close all;




%% spaghetti plot (SST raw, region1) by initial year
for iyear=min(cfg.iyears):max(cfg.iyears)
    tmp.iyear_str=num2str(iyear);
    tmp.plottime=iyear+1/24 : 1/12 : iyear+5 - 1/24;
    varn.modelvar=[cfg.region1, 'm_', varn.varname, '_model_',  'i', tmp.iyear_str];
    varn.obsvar=[cfg.region1, 'm_', varn.varname, '_obs_',  'i', tmp.iyear_str];
    
    plot(tmp.plottime,data.(varn.modelvar)-273.15,'b');
    hold on
    plot(tmp.plottime,data.(varn.obsvar),'k', 'linewidth',3);
end
axis tight
hold off
xlabel('Year'); ylabel('^oC');
set(gca, 'fontsize', 20)
grid minor
set(gcf, 'PaperPosition', [0, 0, 8, 4]);
cfg.figname=[dirs.figdir, filesep, 'spaghetti_', cfg.region1, '_SST_hcst_obs_70-20m_v2', '.tif'];
saveas(gcf,cfg.figname,'tif');
RemoveWhiteSpace([], 'file', cfg.figname);
close all;

%% spaghetti plot (SST raw, region2) by initial year
for iyear=min(cfg.iyears):max(cfg.iyears)
    tmp.iyear_str=num2str(iyear);
    tmp.plottime=iyear+1/24 : 1/12 : iyear+5 - 1/24;
    varn.modelvar=[cfg.region2, 'm_', varn.varname, '_model_',  'i', tmp.iyear_str];
    varn.obsvar=[cfg.region2, 'm_', varn.varname, '_obs_',  'i', tmp.iyear_str];
    
    plot(tmp.plottime,data.(varn.modelvar)-273.15,'b');
    hold on
    plot(tmp.plottime,data.(varn.obsvar),'k', 'linewidth',3);
end
axis tight
hold off
xlabel('Year'); ylabel('^oC');
set(gca, 'fontsize', 20)
grid minor
set(gcf, 'PaperPosition', [0, 0, 8, 4]);
cfg.figname=[dirs.figdir, filesep, 'spaghetti_', cfg.region2, '_SST_hcst_obs_70-20m_v2', '.tif'];
saveas(gcf,cfg.figname,'tif');
RemoveWhiteSpace([], 'file', cfg.figname);
close all;







function obsname_simple = f_obs_simple(obsname)
    switch obsname
        case 'en4.2_ba'
            obsname_simple='en4';
        case 'projdv7.3'
            obsname_simple='projd';
    end
end
function gmsst = f_nino34m_var(var_2d, area)
    mask_var=NaN(size(var_2d));
    mask_var(isfinite(var_2d))=1;
    area=area.*mask_var;
    var_2d=squeeze(var_2d);
    var_2d_sum=var_2d.*area;
    var_2d_sum=sum(var_2d_sum(:), 'omitnan');
    gmsst=var_2d_sum ./ sum(area(:), 'omitnan');
end
function dn = daynoleap2datenum(day, pivotyr)
%DAYNOLEAP2DATENUM Convert from days since, no leap to serial date number
%
% dn = daynoleap2datenum(day, pivotyr)
%
% A lot of model output is saved with a time scale of "days since
% YYYY-01-01 00:00:00", where every year has 365 day.  This function
% converts those dates into a serial date number.
%
% Input variables:
%
%   day:        number of days since pivot year
%
%   pivotyr:    pivot year, i.e. year that day count begins (on Jan 1)
% Determine which years in range are leap years
nyr = max(day./365);
yrs = pivotyr:(pivotyr + nyr);
isleap = @(x) (mod(x,4)==0 & mod(x,100)~=0) | mod(x,400) == 0;
leapyrs = yrs(isleap(yrs));
% Calculate date numbers
dayofyear = rem(day, 365);
yr = floor(day/365) + pivotyr;
dn = datenum(pivotyr, 1, 1) + day;
for ileap = 1:length(leapyrs)
  needsbump = yr > leapyrs(ileap) | (yr == leapyrs(ileap) & dayofyear > 59);
  dn(needsbump) = dn(needsbump) + 1;
end
end