% %  Created 28-Oct-2022 by Yong-Yub Kim
% %  Created 29-Oct-2022 by Yong-Yub Kim
% %  Created 29-Oct-2022 by Yong-Yub Kim addition of MSSS
% %  Created 29-Oct-2022 by Yong-Yub Kim lead month

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
% dirs.root='/mnt/lustre/proj/earth.system.predictability/HCST_EXP';
% dirs.yoshi_root='/proj/yoshi/DATA/CESM2_ODA';
% dirs.archive=[dirs.root, filesep, 'archive'];
% dirs.saveroot='/mnt/lustre/proj/kimyy/Model/CESM2/ESP/HCST_EXP';
dirs.hcstroot='/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/HCST_EXP';
dirs.figroot='/Volumes/kyy_raid/kimyy/Figure/CESM2/ESP/HCST_EXP';

config.iyears=1960:2021;
config.months=1:12;
config.scenname='HIST';
config.gridname='f09_g17';
config.assm_factor='10';
config.ens_member='1';
config.proj_year=10;

config.obsnames={'en4.2_ba'};
config.ensnames={'ba-10p1'};

config.component='ocn';
config.varnames={'temp'};
config.len_t_y = length(config.iyears);
config.len_t_m = length(config.months);
config.len_t= config.len_t_y * config.len_t_m;
config.len_obs= length(config.obsnames);
config.len_ens= length(config.ensnames);


%% grid set(mask from model)
% tmp.obsname=config.obsnames{1};
% iyear=min(config.iyears);
% config.casename_m=[config.gridname, '.hcst.', tmp.obsname, '-', config.assm_factor, 'p', config.ens_member];
% config.casename=[config.casename_m, '_i', num2str(iyear)];
% dirs.datadir= [dirs.archive, filesep, config.casename_m, filesep, config.casename, ...
%     filesep, 'ocn/hist'];
% 
% [tmp.error_status, tmp.value]=system(['ls ', dirs.datadir, '/*once*']);  % b.e21.BHISTsmbb.f09_g17.assm.oras4_ba-10p1.pop.h.once.nc
% tmp.gridname = [tmp.value(1:end-1)];
% grid.region_mask=ncread(tmp.gridname, 'REGION_MASK'); 
% grid.ocean_mask=NaN(size(grid.region_mask));
% grid.ocean_mask(grid.region_mask>0)=1;
% grid.tarea = ncread(tmp.gridname, 'TAREA');
% 
% grid.tlong=ncread(config.obs_en4_filename, 'TLONG');
% grid.tlat=ncread(config.obs_en4_filename, 'TLAT');
% grid.nlon=size(grid.tlong,1);
% grid.nlat=size(grid.tlong,2);
% grid.ntime=config.proj_year.*12;


% % model filename example
% /mnt/lustre/proj/earth.system.predictability/HCST_EXP/archive/
% f09_g17.hcst.en4.2_ba-10p1/f09_g17.hcst.en4.2_ba-10p1_i1993/ocn/hist/
% f09_g17.hcst.en4.2_ba-10p1_i1993.pop.h.1993-01.nc

%% read & plot data
for obsind=1:length(config.obsnames)
    tmp.obsname=config.obsnames{obsind};
    tmp.obsname_simple= f_obs_simple(tmp.obsname);
    for iyear=min(config.iyears):max(config.iyears)
        tmp.iyear_str=num2str(iyear, '%04i');
        config.casename_m=[config.gridname, '.hcst.', tmp.obsname, '-', config.assm_factor, 'p', config.ens_member];
        config.casename=[config.casename_m, '_i', num2str(iyear)];
        dirs.datadir= [dirs.hcstroot, filesep, config.casename_m, filesep, 'NINO34'];
        config.datafilename=[dirs.datadir, filesep, 'NINO34_', config.casename, '.nc'];

        for varind=1:length(config.varnames)
            tmp.varname=config.varnames{varind};
            %% read variables
            tmp.modelvar = ['nino34m_', tmp.varname, '_model_', tmp.obsname_simple, '_i', tmp.iyear_str];
            tmp.assmvar = ['nino34m_', tmp.varname, '_assm_', tmp.obsname_simple, '_i', tmp.iyear_str];
            tmp.obsvar = ['nino34m_', tmp.varname, '_obs_', tmp.obsname_simple, '_i', tmp.iyear_str];

            data.(tmp.modelvar)=ncread(config.datafilename, ['nino34m_', tmp.varname]);
            data.(tmp.assmvar)=ncread(config.datafilename, ['nino34m_assm_', tmp.varname]);
            data.(tmp.obsvar)=ncread(config.datafilename, ['nino34m_obs_', tmp.varname]);
            data.time_tmp=ncread(config.datafilename, 'time');
            data.time_tmp_leap=daynoleap2datenum(data.time_tmp, 0)-1; % reflection of leap year and -1 for month correction
            data.time_vec_tmp=datevec(data.time_tmp_leap);
            %% assign variables according to lead year
            for lmonth=0:config.proj_year*12-1
                tmp.lmonth_str=num2str(lmonth, '%03i');
                tmp.modelvar_lm=['nino34m_', tmp.varname, '_model_', tmp.obsname_simple, '_l', tmp.lmonth_str];
                tmp.ydata= ncread(config.datafilename, ['nino34m_', tmp.varname], lmonth+1, 1);
                tmp.yind=(iyear-min(config.iyears))+1;
                data.(tmp.modelvar_lm)(tmp.yind)=tmp.ydata;
                data.(tmp.modelvar_lm)(data.(tmp.modelvar_lm)==0)=NaN;
            end
            tmp.yind=(iyear-min(config.iyears))*12+1:(iyear-min(config.iyears))*12+12;
            tmp.ydata= ncread(config.datafilename, ['nino34m_obs_', tmp.varname], [1], 12);
            data.(['nino34m_', tmp.varname, '_obs_', tmp.obsname_simple])(tmp.yind)= tmp.ydata;
            tmp.ydata= ncread(config.datafilename, ['nino34m_assm_', tmp.varname], [1], 12);
            data.(['nino34m_', tmp.varname, '_assm_', tmp.obsname_simple])(tmp.yind)= tmp.ydata;
            data.time_leap(tmp.yind)=data.time_tmp_leap(1:12);
            data.time_vec(tmp.yind,:)= data.time_vec_tmp(1:12,:);
        end
%% save ncfile

%     [f_exp, gof_exp] = fit(trendtime_yearly',squeeze(cmems_sla_yearly(i,j,:)),'poly1');
%                             cmems_sla_yearly_poly1_fit(i,j,1:length(trendtime_yearly))=f_exp(trendtime_yearly);

    end
    disp('abc')

    for iyear=min(config.iyears):max(config.iyears)
        tmp.iyear_str=num2str(iyear, '%04i');
        for varind=1:length(config.varnames)
            tmp.varname=config.varnames{varind};
            for lmonth=0:config.proj_year*12-1
                tmp.lmonth_str=num2str(lmonth, '%03i');
                tmp.yind=(iyear-min(config.iyears))+1;
                tmp.obsvar_lm=['nino34m_', tmp.varname, '_obs_', tmp.obsname_simple, '_l', tmp.lmonth_str];
                tmp.obsvar=['nino34m_', tmp.varname, '_obs_', tmp.obsname_simple];
                tmp.assmvar_lm=['nino34m_', tmp.varname, '_assm_', tmp.obsname_simple, '_l', tmp.lmonth_str];
                tmp.assmvar=['nino34m_', tmp.varname, '_assm_', tmp.obsname_simple];
                tmp.mind=(iyear-min(config.iyears))*12+lmonth+1;
                if tmp.mind<=length(data.(tmp.obsvar))
                    tmp.ydata_obs=data.(tmp.obsvar)(tmp.mind);
                    tmp.ydata_assm=data.(tmp.assmvar)(tmp.mind);
                    data.(tmp.obsvar_lm)(tmp.yind)=tmp.ydata_obs;
                    data.(tmp.assmvar_lm)(tmp.yind)=tmp.ydata_assm;
                else
                    data.(tmp.obsvar_lm)(tmp.yind)=NaN;
                    data.(tmp.assmvar_lm)(tmp.yind)=NaN;
                end

            end
        end
    end
    
    %% make anomaly for nino index
    for lmonth=0:config.proj_year*12-1
        tmp.lmonth_str=num2str(lmonth, '%03i');
        tmp.avgwin_min=find(config.iyears==1991)-floor(lmonth/12);
        tmp.avgwin_max=find(config.iyears==2020)-floor(lmonth/12);
        tmp.modelvar_lm=['nino34m_', tmp.varname, '_model_', tmp.obsname_simple, '_l', tmp.lmonth_str];
        tmp.obsvar_lm=['nino34m_', tmp.varname, '_obs_', tmp.obsname_simple, '_l', tmp.lmonth_str];
        tmp.assmvar_lm=['nino34m_', tmp.varname, '_assm_', tmp.obsname_simple, '_l', tmp.lmonth_str];

        tmp.modelvar_ano_lm=['nino34m_', tmp.varname, '_model_ano_', tmp.obsname_simple, '_l', tmp.lmonth_str];
        tmp.obsvar_ano_lm=['nino34m_', tmp.varname, '_obs_ano_', tmp.obsname_simple, '_l', tmp.lmonth_str];
        tmp.assmvar_ano_lm=['nino34m_', tmp.varname, '_assm_ano_', tmp.obsname_simple, '_l', tmp.lmonth_str];

        tmp.model_ltm=mean(data.(tmp.modelvar_lm)(tmp.avgwin_min:tmp.avgwin_max));
        tmp.obs_ltm=mean(data.(tmp.obsvar_lm)(tmp.avgwin_min:tmp.avgwin_max));
        tmp.assm_ltm=mean(data.(tmp.assmvar_lm)(tmp.avgwin_min:tmp.avgwin_max));

        data.(tmp.modelvar_ano_lm)=data.(tmp.modelvar_lm)-tmp.model_ltm;
        data.(tmp.obsvar_ano_lm)=data.(tmp.obsvar_lm)-tmp.obs_ltm;
        data.(tmp.assmvar_ano_lm)=data.(tmp.assmvar_lm)-tmp.assm_ltm;

        tmp.modelvar_ano_lm_rm=['nino34m_', tmp.varname, '_model_ano_', tmp.obsname_simple, '_l', tmp.lmonth_str, '_rm'];
        tmp.obsvar_ano_lm_rm=['nino34m_', tmp.varname, '_obs_ano_', tmp.obsname_simple, '_l', tmp.lmonth_str, '_rm'];
        tmp.assmvar_ano_lm_rm=['nino34m_', tmp.varname, '_assm_ano_', tmp.obsname_simple, '_l', tmp.lmonth_str, '_rm'];
        
        tmp.rm_window=3;
        data.(tmp.modelvar_ano_lm_rm)=movmean(data.(tmp.modelvar_ano_lm), tmp.rm_window, 'Endpoints', 'fill');
        data.(tmp.obsvar_ano_lm_rm)=movmean(data.(tmp.obsvar_ano_lm), tmp.rm_window, 'Endpoints', 'fill');
        data.(tmp.assmvar_ano_lm_rm)=movmean(data.(tmp.assmvar_ano_lm), tmp.rm_window, 'Endpoints', 'fill');

    end

    %% set time value for figure
    data.time_vec_extended=data.time_vec;
    for lmonth=1:config.proj_year-1
        tmp.endind=size(data.time_vec_extended,1);
        data.time_vec_extended(tmp.endind+1:tmp.endind+12,1)=max(config.iyears)+lmonth;
        data.time_vec_extended(tmp.endind+1:tmp.endind+12,2)=data.time_vec_extended(1:12,2);
        data.time_vec_extended(tmp.endind+1:tmp.endind+12,3)=data.time_vec_extended(1:12,3);
    end
    data.time_leap_extended=datenum(data.time_vec_extended);

 %% normal time series and correlation coefficient (combined plot, hcst&assm)
    for varind=1:length(config.varnames)
        tmp.varname=config.varnames{varind};
        dirs.figdir= [dirs.figroot, filesep, config.casename_m, filesep, 'NINO34', filesep, 'Timeseries'];
        system(['mkdir -p ', dirs.figdir]);
%         for lmonth=0:config.proj_year*12-1
        for lmonth=0:35
            tmp.lmonth_str=num2str(lmonth, '%03i');
            tmp.modelvar_ano_lm_rm=['nino34m_', tmp.varname, '_model_ano_', tmp.obsname_simple, '_l', tmp.lmonth_str, '_rm'];
            tmp.obsvar_ano_lm_rm=['nino34m_', tmp.varname, '_obs_ano_', tmp.obsname_simple, '_l', tmp.lmonth_str, '_rm'];
            tmp.assmvar_ano_lm_rm=['nino34m_', tmp.varname, '_assm_ano_', tmp.obsname_simple, '_l', tmp.lmonth_str, '_rm'];

            tmp.md=data.(tmp.modelvar_ano_lm_rm);
            tmp.obs=data.(tmp.obsvar_ano_lm_rm);
            tmp.assm=data.(tmp.assmvar_ano_lm_rm);

            tmp.corrcoef=corrcoef(tmp.md(isfinite(tmp.assm)), tmp.assm(isfinite(tmp.assm)));
            data.(['nino34m_', tmp.varname, '_corr'])(lmonth+1)=tmp.corrcoef(1,2);

        end

        %% correlation coefficient of NINO3.4 ind, function of lead year (3-mon running mean)
            config.figname=[dirs.figdir, filesep, 'nino34m_hcst_assm_91-20m', '_corr_leadm_3rm', '.tif'];
%             bar(0:config.proj_year*12-1,data.(['nino34m_', tmp.varname, '_corr']),  'linewidth', 2)
            bar(1:36,data.(['nino34m_', tmp.varname, '_corr'])(1:36),  'linewidth', 2)

            xlabel('Lead month'); ylabel(['corr. coef.,', 'Nino3.4' , ', assm-ba']);
            set(gca, 'fontsize', 20)
            grid minor
            ylim([-0.4 1])
            saveas(gcf,config.figname,'tif');
            RemoveWhiteSpace([], 'file', config.figname);
            close all;
        
    end


end




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