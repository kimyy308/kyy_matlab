% %  Created 30-Jan-2023 by Yong-Yub Kim lead month

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

cfg.varnames={'RAIN'};
cfg.var='RAIN';

%% model configuration

dirs.hcstroot=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/HCST_EXP/archive/lnd/', cfg.var];
dirs.assmroot=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/ASSM_EXP/archive_transfer/lnd/', cfg.var];
dirs.lens2root=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/LENS2/archive_analysis/lnd/', cfg.var];
dirs.obsroot=['/Volumes/kyy_raid/kimyy/Observation/GPCC/monthly_reg_cam'];
dirs.figroot=['/Volumes/kyy_raid/kimyy/Figure/CESM2/ESP/HCST_EXP/archive/lnd/', cfg.var];

cfg.iyears=1960:2020;
cfg.months=1:12;
% cfg.scenname='HIST';
cfg.gnm='f09_g17';
% cfg.assm_factor='10';
% cfg.ens_member='1';
cfg.proj_year=5;

% cfg.obsnames={'en4.2_ba'};
% cfg.ensnames={'ba-10p1'};

cfg.component='lnd';

cfg.len_t_y = length(cfg.iyears);
cfg.len_t_m = length(cfg.months);
cfg.len_t= cfg.len_t_y * cfg.len_t_m;
% cfg.len_obs= length(cfg.obsnames);
% cfg.len_ens= length(cfg.ensnames);
cfg.region = 'GLO';

%% grid set(mask from model)
% tmp.gridname = [dirs.hcstroot, tmp.fs, 'grid.nc'];
% 
% grid.lon=ncread(tmp.gridname, 'lon');
% grid.lat=ncread(tmp.gridname, 'lat');
% [grid.tlat grid.tlong]=meshgrid(grid.lat, grid.lon);
% 
% grid.nlon=size(grid.tlong,1);
% grid.nlat=size(grid.tlat,2);


str.iy_min=num2str(min(cfg.iyears));
str.iy_max=num2str(max(cfg.iyears));

%% read & plot data
tmp.varname=cfg.var;
for iyear=min(cfg.iyears):max(cfg.iyears)
    tmp.iyear_str=num2str(iyear, '%04i');
    cfg.casename_m=['ens_all'];
    cfg.casename=[cfg.casename_m, '_i', tmp.iyear_str];
    dirs.datadir= [dirs.hcstroot, filesep, cfg.casename_m, filesep, cfg.casename, tmp.fs, cfg.region];
    dirs.assmdir= [dirs.assmroot, filesep, cfg.casename_m, filesep, cfg.region];
    dirs.lens2dir= [dirs.lens2root, filesep, 'ens_smbb', filesep, cfg.region];

    tmp.modelvar = [cfg.region, 'm_', tmp.varname, '_model_',  'i', tmp.iyear_str];
    tmp.modelvar_y = [cfg.region, 'y_', tmp.varname, '_model_',  'i', tmp.iyear_str];
    tmp.modelvar_inity = [cfg.region, 'y_', tmp.varname, '_model'];

    tmp.assmvar = [cfg.region, 'm_', tmp.varname, '_assm_',  'i', tmp.iyear_str];
    tmp.assmvar_y = [cfg.region, 'y_', tmp.varname, '_assm_',  'i', tmp.iyear_str];
    tmp.assmvar_y_raw = [cfg.region, 'y_', tmp.varname, '_assm'];

    tmp.lens2var = [cfg.region, 'm_', tmp.varname, '_lens2_',  'i', tmp.iyear_str];
    tmp.lens2var_y = [cfg.region, 'y_', tmp.varname, '_lens2_',  'i', tmp.iyear_str];
    tmp.lens2var_y_raw = [cfg.region, 'y_', tmp.varname, '_lens2'];

%     tmp.assmvar = [cfg.region, 'm_', tmp.varname, '_assm_',  '_i', tmp.iyear_str];
    tmp.obsvar = [cfg.region, 'm_', tmp.varname, '_obs_',  'i', tmp.iyear_str];
    tmp.obsvar_y = [cfg.region, 'y_', tmp.varname, '_obs_',  'i', tmp.iyear_str];
    tmp.obsvar_y_raw = [cfg.region, 'y_', tmp.varname, '_obs'];

    for fy = iyear:iyear+cfg.proj_year-1
        tmp.fy=fy;
        tmp.fy_str=num2str(tmp.fy, '%04i');
        for mon=1:12
            tmp.mon_str=num2str(mon, '%02i');
%             cfg.datafilename=[dirs.datadir, filesep, ...
%                 'M_', cfg.region, '_', tmp.varname, '_', cfg.gnm, '.hcst.', cfg.casename, '.cam.h0.', tmp.fy_str, '-', tmp.mon_str, '.nc'];
%             data.(tmp.modelvar)((fy-iyear)*12+mon)=ncread(cfg.datafilename, tmp.varname);
            cfg.datafilename=[dirs.datadir, filesep, ...
                'M_', cfg.region, '_', tmp.varname, '_', cfg.gnm, '.hcst.', 'ensmean_all_', 'i', tmp.iyear_str, '.clm2.h0.', tmp.fy_str, '-', tmp.mon_str, '.nc'];
            data.(tmp.modelvar)((fy-iyear)*12+mon)=ncread(cfg.datafilename, tmp.varname);
            data.time_tmp((fy-iyear)*12+mon)=ncread(cfg.datafilename, 'time');
            
            %% read obs
            cfg.obsfnm = [dirs.obsroot, tmp.fs, cfg.region, tmp.fs, 'M_', cfg.region, '_', 'GPCC_reg_cesm2.v5.',tmp.fy_str,tmp.mon_str, '.nc'];
            if exist(cfg.obsfnm)==0
                data.(tmp.obsvar)((fy-iyear)*12+mon)=NaN;
            else
                tmp.dd=ncread(cfg.obsfnm, 'precip');
                tmp.dd=tmp.dd./86400/eomday(tmp.fy,mon);
                data.(tmp.obsvar)((fy-iyear)*12+mon)=tmp.dd;
            end
            
            %% read assm
            cfg.assmfilename=[dirs.assmdir, filesep, ...
                'M_', cfg.region, '_', tmp.varname, '_ensmean_', tmp.fy_str, '-', tmp.mon_str, '.nc'];
            if exist(cfg.assmfilename)==0
                data.(tmp.assmvar)((fy-iyear)*12+mon)=NaN;
            else
                data.(tmp.assmvar)((fy-iyear)*12+mon)=ncread(cfg.assmfilename, tmp.varname);
            end

            %% read lens2
            cfg.lens2filename=[dirs.lens2dir, filesep, ...
                'M_', cfg.region, '_', tmp.varname, '_ensmean_', tmp.fy_str, '-', tmp.mon_str, '.nc'];
            if exist(cfg.lens2filename)==0
                data.(tmp.lens2var)((fy-iyear)*12+mon)=NaN;
            else
                data.(tmp.lens2var)((fy-iyear)*12+mon)=ncread(cfg.lens2filename, tmp.varname);
            end
            

        end
        tmp.ymean=mean(data.(tmp.modelvar)((fy-iyear)*12+1:(fy-iyear)*12+12));
        data.(tmp.modelvar_y)((fy-iyear+1))=tmp.ymean;
        tmp.ymean_obs=mean(data.(tmp.obsvar)((fy-iyear)*12+1:(fy-iyear)*12+12));
        data.(tmp.obsvar_y)((fy-iyear+1))=tmp.ymean_obs;
        tmp.ymean_assm=mean(data.(tmp.assmvar)((fy-iyear)*12+1:(fy-iyear)*12+12));
        data.(tmp.assmvar_y)((fy-iyear+1))=tmp.ymean_assm;
        tmp.ymean_lens2=mean(data.(tmp.lens2var)((fy-iyear)*12+1:(fy-iyear)*12+12));
        data.(tmp.lens2var_y)((fy-iyear+1))=tmp.ymean_lens2;
    end
    data.(tmp.obsvar_y_raw)(iyear-min(cfg.iyears)+1)=data.(tmp.obsvar_y)(1);
    data.(tmp.modelvar_inity)(iyear-min(cfg.iyears)+1)=data.(tmp.modelvar_y)(1);
    data.(tmp.assmvar_y_raw)(iyear-min(cfg.iyears)+1)=data.(tmp.assmvar_y)(1);
    data.(tmp.lens2var_y_raw)(iyear-min(cfg.iyears)+1)=data.(tmp.lens2var_y)(1);

end
% disp('abc')


dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, cfg.region, filesep, 'Timeseries'];


%% raw value spaghetti figure
plot(cfg.iyears,data.(tmp.obsvar_y_raw))
hold on
for iyear=min(cfg.iyears):max(cfg.iyears)
    tmp.iyear_str=num2str(iyear, '%04i');
    tmp.modelvar_y = [cfg.region, 'y_', tmp.varname, '_model_',  'i', tmp.iyear_str];
    plot(iyear:iyear+cfg.proj_year-1, data.(tmp.modelvar_y)-273.15)
end
hold off


%% anomaly figure
mval.model=mean(data.(tmp.modelvar_inity), 'omitnan');
mval.obs=mean(data.(tmp.obsvar_y_raw), 'omitnan');
mval.assm=mean(data.(tmp.assmvar_y_raw), 'omitnan');
mval.lens2=mean(data.(tmp.lens2var_y_raw), 'omitnan');

plot(cfg.iyears,data.(tmp.obsvar_y_raw)-mval.obs, 'k', 'linewidth', 3)
hold on
plot(cfg.iyears,data.(tmp.assmvar_y_raw)-mval.assm, 'r', 'linewidth', 3)
plot(cfg.iyears,data.(tmp.lens2var_y_raw)-mval.lens2, 'g', 'linewidth', 3)

for iyear=min(cfg.iyears):max(cfg.iyears)
    tmp.iyear_str=num2str(iyear, '%04i');
    tmp.modelvar_y = [cfg.region, 'y_', tmp.varname, '_model_',  'i', tmp.iyear_str];
    plot(iyear:iyear+cfg.proj_year-1, data.(tmp.modelvar_y)-mval.model, 'b')
end
hold off
axis tight
hold off
xlabel('Year'); ylabel('mm/s');
set(gca, 'fontsize', 15)
grid minor
axis tight
ylim([-2 2].*10^-6)

set(gcf, 'PaperPosition', [0, 0, 6, 4]);
% legend ('OBS', 'ASSM', 'LENS2', 'HCST', 'Location', 'Northwest')
legend ('OBS', 'ASSM', 'LENS2', 'HCST', 'Location', 'Southoutside', 'Orientation', 'Horizontal')

cfg.figname=[dirs.figdir, filesep, 'spaghetti_anom_', tmp.varname, '_hcst_all_', str.iy_min, '-', str.iy_max,'_proj_',num2str(cfg.proj_year),  'y.tif'];
mkdir(dirs.figdir)
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