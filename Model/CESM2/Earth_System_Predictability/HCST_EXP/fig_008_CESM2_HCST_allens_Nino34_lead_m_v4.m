% %  Created 27-Mar-2023 by Yong-Yub Kim lead month

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
dirs.obsroot=['/Volumes/kyy_raid/kimyy/Observation/ERSST/monthly_reg_cam'];
dirs.figroot=['/Volumes/kyy_raid/kimyy/Figure/CESM2/ESP/HCST_EXP/archive/atm/', cfg.var];

cfg.iyears=1970:2019;
% cfg.iyears=1982:2016;

cfg.months=1:12;
% cfg.scenname='HIST';
cfg.gnm='f09_g17';
% cfg.assm_factor='10';
% cfg.ens_member='1';
cfg.proj_year=3;

% cfg.obsnames={'en4.2_ba'};
% cfg.ensnames={'ba-10p1'};

cfg.component='atm';
cfg.varnames={'SST'};
cfg.len_t_y = length(cfg.iyears);
cfg.len_t_m = length(cfg.months);
cfg.len_t= cfg.len_t_y * cfg.len_t_m;
% cfg.len_obs= length(cfg.obsnames);
% cfg.len_ens= length(cfg.ensnames);
cfg.region = 'NINO34';

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

    tmp.modelvar = ['nino34m_', tmp.varname, '_model_',  'i', tmp.iyear_str];
%     tmp.assmvar = ['nino34m_', tmp.varname, '_assm_',  '_i', tmp.iyear_str];
    tmp.obsvar = ['nino34m_', tmp.varname, '_obs_',  'i', tmp.iyear_str];

    for fy = iyear:iyear+cfg.proj_year-1
        tmp.fy=fy;
        tmp.fy_str=num2str(tmp.fy, '%04i');
        for mon=1:12
            tmp.mon_str=num2str(mon, '%02i');
            cfg.datafilename=[dirs.datadir, filesep, ...
                'M_', cfg.region, '_', tmp.varname, '_', cfg.gnm, '.hcst.', cfg.casename, '.cam.h0.', tmp.fy_str, '-', tmp.mon_str, '.nc'];
            data.(tmp.modelvar)((fy-iyear)*12+mon)=ncread(cfg.datafilename, tmp.varname);
            data.time_tmp((fy-iyear)*12+mon)=ncread(cfg.datafilename, 'time');
            cfg.obsfnm = [dirs.obsroot, tmp.fs, cfg.region, tmp.fs, 'M_', cfg.region, '_', 'ersst_reg_cesm2.v5.',tmp.fy_str,tmp.mon_str, '.nc'];
            if exist(cfg.obsfnm)==0
                data.(tmp.obsvar)((fy-iyear)*12+mon)=NaN;
            else
                data.(tmp.obsvar)((fy-iyear)*12+mon)=ncread(cfg.obsfnm, 'sst');
            end
        end
    end

    %% read variables (assign variable to new one by lead month)
%     data.time_tmp=ncread(cfg.datafilename, 'time');
    data.time_tmp_leap=daynoleap2datenum(data.time_tmp, 0)-1; % reflection of leap year and -1 for month correction
    data.time_vec_tmp=datevec(data.time_tmp_leap);
    data.time_vec_tmp(:,1)=data.time_vec_tmp(:,1)+iyear; % y correction (0~4 + iyear)
    data.time_tmp_leap=datenum(data.time_vec_tmp); % recal time var
    %% assign variables according to lead year
    for lmonth=0:cfg.proj_year*12-1
        tmp.lmonth_str=num2str(lmonth, '%03i');
        tmp.modelvar_lm=['nino34m_', tmp.varname, '_model',  '_l', tmp.lmonth_str];
        tmp.mdata= data.(tmp.modelvar)(lmonth+1);
        tmp.mind=(iyear-min(cfg.iyears))+1;
        data.(tmp.modelvar_lm)(tmp.mind)=tmp.mdata;
        data.(tmp.modelvar_lm)(data.(tmp.modelvar_lm)==0)=NaN;
    end
    tmp.yind=(iyear-min(cfg.iyears))*12+1:(iyear-min(cfg.iyears))*12+12;
    tmp.ydata= data.(tmp.obsvar)(1:12);
    data.(['nino34m_', tmp.varname, '_obs'])(tmp.yind)= tmp.ydata;
%     tmp.ydata= ncread(cfg.datafilename, ['nino34m_assm_', tmp.varname], [1], 12);
%     data.(['nino34m_', tmp.varname, '_assm_', tmp.obsname_simple])(tmp.yind)= tmp.ydata;
    data.time_leap(tmp.yind)=data.time_tmp_leap(1:12);
    data.time_vec(tmp.yind,:)= data.time_vec_tmp(1:12,:);

%     [f_exp, gof_exp] = fit(trendtime_yearly',squeeze(cmems_sla_yearly(i,j,:)),'poly1');
%                             cmems_sla_yearly_poly1_fit(i,j,1:length(trendtime_yearly))=f_exp(trendtime_yearly);

end
% disp('abc')

for iyear=min(cfg.iyears):max(cfg.iyears)
    tmp.iyear_str=num2str(iyear, '%04i');
    tmp.varname=cfg.var;
    for lmonth=0:cfg.proj_year*12-1
        tmp.lmonth_str=num2str(lmonth, '%03i');
        tmp.yind=(iyear-min(cfg.iyears))+1;
        tmp.obsvar_lm=['nino34m_', tmp.varname, '_obs_',  '_l', tmp.lmonth_str];
        tmp.obsvar=['nino34m_', tmp.varname, '_obs'];
%             tmp.assmvar_lm=['nino34m_', tmp.varname, '_assm_',  '_l', tmp.lmonth_str];
%             tmp.assmvar=['nino34m_', tmp.varname, '_assm_', tmp.obsname_simple];
        tmp.mind=(iyear-min(cfg.iyears))*12+lmonth+1;
        if tmp.mind<=length(data.(tmp.obsvar))
            tmp.ydata_obs=data.(tmp.obsvar)(tmp.mind);
%             tmp.ydata_assm=data.(tmp.assmvar)(tmp.mind);
            data.(tmp.obsvar_lm)(tmp.yind)=tmp.ydata_obs;
%             data.(tmp.assmvar_lm)(tmp.yind)=tmp.ydata_assm;
        else
            data.(tmp.obsvar_lm)(tmp.yind)=NaN;
%             data.(tmp.assmvar_lm)(tmp.yind)=NaN;
        end
    end
end

%% make anomaly for nino index
for lmonth=0:cfg.proj_year*12-1
    tmp.lmonth_str=num2str(lmonth, '%03i');
% %     %% recent nino mean definition
% %     tmp.avgwin_min=find(cfg.iyears==1996)-floor(lmonth/12);
% %     tmp.avgwin_max=find(cfg.iyears==2020)-floor(lmonth/12);
    %% whole period nino mean definition
    %% ex) if iy= 1970 ~ 2019, 3y hindcast: 1972~2019 mean
    tmp.avgwin_min=find(cfg.iyears==min(cfg.iyears))+cfg.proj_year-1-floor(lmonth/12);
    tmp.avgwin_max=find(cfg.iyears==max(cfg.iyears))-floor(lmonth/12);
    tmp.avgwin=tmp.avgwin_min:tmp.avgwin_max;
    tmp.Nwin=length(tmp.avgwin);

    tmp.modelvar_lm=['nino34m_', tmp.varname, '_model',  '_l', tmp.lmonth_str];
    tmp.obsvar_lm=['nino34m_', tmp.varname, '_obs_',  '_l', tmp.lmonth_str];
%     tmp.assmvar_lm=['nino34m_', tmp.varname, '_assm_',  '_l', tmp.lmonth_str];

    tmp.modelvar_ano_lm=['nino34m_', tmp.varname, '_model_ano_',  '_l', tmp.lmonth_str];
    tmp.obsvar_ano_lm=['nino34m_', tmp.varname, '_obs_ano_',  '_l', tmp.lmonth_str];
%     tmp.assmvar_ano_lm=['nino34m_', tmp.varname, '_assm_ano_',  '_l', tmp.lmonth_str];

    tmp.model_ltm=mean(data.(['nino34m_', tmp.varname, '_model',  '_l', '000'])(tmp.avgwin_min:tmp.avgwin_max));
    tmp.obs_ltm=mean(data.(['nino34m_', tmp.varname, '_obs_',  '_l', '000'])(tmp.avgwin_min:tmp.avgwin_max));
    %% v4
%     tmp.assm_ltm=mean(data.(tmp.assmvar_lm)(tmp.avgwin_min:tmp.avgwin_max));

    data.(tmp.modelvar_ano_lm)=data.(tmp.modelvar_lm)-tmp.model_ltm;
    data.(tmp.obsvar_ano_lm)=data.(tmp.obsvar_lm)-tmp.obs_ltm;
%     data.(tmp.assmvar_ano_lm)=data.(tmp.assmvar_lm)-tmp.assm_ltm;

    tmp.modelvar_ano_lm_rm=['nino34m_', tmp.varname, '_model_ano_',  '_l', tmp.lmonth_str, '_rm'];
    tmp.obsvar_ano_lm_rm=['nino34m_', tmp.varname, '_obs_ano_',  '_l', tmp.lmonth_str, '_rm'];
%     tmp.assmvar_ano_lm_rm=['nino34m_', tmp.varname, '_assm_ano_',  '_l', tmp.lmonth_str, '_rm'];
    
%% movmean for averaged rm_window years -> wrong, so rm_window should set to 1
    tmp.rm_window=1;
    data.(tmp.modelvar_ano_lm_rm)=movmean(data.(tmp.modelvar_ano_lm), tmp.rm_window, 'Endpoints', 'fill');
    data.(tmp.obsvar_ano_lm_rm)=movmean(data.(tmp.obsvar_ano_lm), tmp.rm_window, 'Endpoints', 'fill');
%     data.(tmp.assmvar_ano_lm_rm)=movmean(data.(tmp.assmvar_ano_lm), tmp.rm_window, 'Endpoints', 'fill');
    
    %% separate El Nino(>=0.5), La Nina year (<=-0.5)
    if lmonth==0
        data.Y_Elnino=find(data.(tmp.obsvar_ano_lm_rm)>=0.5);
        data.Y_Elnino(data.Y_Elnino>max(tmp.avgwin)-floor((cfg.proj_year*12-1)/12))=[];
%         data.Y_Elnino(data.Y_Elnino<min(tmp.avgwin)+cfg.proj_year-1)=[];
        data.Y_Lanina=find(data.(tmp.obsvar_ano_lm_rm)<=-0.5);
        data.Y_Lanina(data.Y_Lanina>max(tmp.avgwin)-floor((cfg.proj_year*12-1)/12))=[];
%         data.Y_Lanina(data.Y_Lanina<min(tmp.avgwin)+cfg.proj_year-1)=[];
        data.Y_Neutral=find(data.(tmp.obsvar_ano_lm_rm)<0.5 & data.(tmp.obsvar_ano_lm_rm)>-0.5);
        data.Y_Neutral(data.Y_Neutral>max(tmp.avgwin)-floor((cfg.proj_year*12-1)/12))=[];
%         data.Y_Neutral(data.Y_Neutral<min(tmp.avgwin)+cfg.proj_year-1)=[];

        %% get the number of years
        tmp.NoY=min([length(data.Y_Elnino), length(data.Y_Lanina), length(data.Y_Neutral)]);
        [tmp.El_val, tmp.El_ind]=sort(data.(tmp.obsvar_ano_lm_rm)(data.Y_Elnino), 'descend'); % high -> low
        [tmp.Ln_val, tmp.Ln_ind]=sort(data.(tmp.obsvar_ano_lm_rm)(data.Y_Lanina), 'ascend'); % low -> high
        [tmp.N_val, tmp.N_ind]=sort(abs(data.(tmp.obsvar_ano_lm_rm)(data.Y_Neutral)), 'ascend'); % low -> high (abs)
        
        %% choose same number of years
        data.Y_Elnino = data.Y_Elnino(tmp.El_ind(1:tmp.NoY));
        data.Y_Lanina = data.Y_Lanina(tmp.Ln_ind(1:tmp.NoY));
        data.Y_Neutral = data.Y_Neutral(tmp.N_ind(1:tmp.NoY));
    end

    %% nRMSE get
    tmp.sig_0=sqrt( sum( data.(tmp.obsvar_ano_lm_rm)(tmp.avgwin).^2 ) ./ tmp.Nwin ); %denominator (sigma_0)
    tmp.RMSE=sqrt( sum( ( data.(tmp.modelvar_ano_lm_rm)(tmp.avgwin)-data.(tmp.obsvar_ano_lm_rm)(tmp.avgwin) ).^2 ) ./ tmp.Nwin );
    data.(['nino34m_', tmp.varname, '_RMSE'])(lmonth+1)=tmp.RMSE ./ tmp.sig_0;

    %% nRMSE get (Elnino start)
    tmp.Nwin_Elnino=length(data.Y_Elnino);
    tmp.sig_0=sqrt( sum( data.(tmp.obsvar_ano_lm_rm)(data.Y_Elnino).^2 ) ./ tmp.Nwin_Elnino ); %denominator (sigma_0)
    tmp.RMSE=sqrt( sum( ( data.(tmp.modelvar_ano_lm_rm)(data.Y_Elnino)-data.(tmp.obsvar_ano_lm_rm)(data.Y_Elnino) ).^2 ) ./ tmp.Nwin_Elnino);
    data.(['nino34m_', tmp.varname, '_RMSE_Elnino'])(lmonth+1)=tmp.RMSE ./ tmp.sig_0;

    %% nRMSE get (Lanina start)
    tmp.Nwin_Lanina=length(data.Y_Lanina);
    tmp.sig_0=sqrt( sum( data.(tmp.obsvar_ano_lm_rm)(data.Y_Lanina).^2 ) ./ tmp.Nwin_Lanina ); %denominator (sigma_0)
    tmp.RMSE=sqrt( sum( ( data.(tmp.modelvar_ano_lm_rm)(data.Y_Lanina)-data.(tmp.obsvar_ano_lm_rm)(data.Y_Lanina) ).^2 ) ./ tmp.Nwin_Lanina);
    data.(['nino34m_', tmp.varname, '_RMSE_Lanina'])(lmonth+1)=tmp.RMSE ./ tmp.sig_0;

    %% nRMSE get (Neutral start)
    tmp.Nwin_Neutral=length(data.Y_Neutral);
    tmp.sig_0=sqrt( sum( data.(tmp.obsvar_ano_lm_rm)(data.Y_Neutral).^2 ) ./ tmp.Nwin_Neutral ); %denominator (sigma_0)
    tmp.RMSE=sqrt( sum( ( data.(tmp.modelvar_ano_lm_rm)(data.Y_Neutral)-data.(tmp.obsvar_ano_lm_rm)(data.Y_Neutral) ).^2 ) ./ tmp.Nwin_Neutral);
    data.(['nino34m_', tmp.varname, '_RMSE_Neutral'])(lmonth+1)=tmp.RMSE ./ tmp.sig_0;


end

%% set time value for figure
data.time_vec_extended=data.time_vec;
for lmonth=1:cfg.proj_year-1
    tmp.endind=size(data.time_vec_extended,1);
    data.time_vec_extended(tmp.endind+1:tmp.endind+12,1)=max(cfg.iyears)+lmonth;
    data.time_vec_extended(tmp.endind+1:tmp.endind+12,2)=data.time_vec_extended(1:12,2);
    data.time_vec_extended(tmp.endind+1:tmp.endind+12,3)=data.time_vec_extended(1:12,3);
end
data.time_leap_extended=datenum(data.time_vec_extended);

%% normal time series and correlation coefficient (combined plot, hcst&obs)
% for varind=1:length(cfg.varnames)
    tmp.varname=cfg.var;
    dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, cfg.region, filesep, 'Timeseries'];
    system(['mkdir -p ', dirs.figdir]);
    for lmonth=0:cfg.proj_year*12-1
%     for lmonth=0:35

        tmp.lmonth_str=num2str(lmonth, '%03i');
        tmp.modelvar_ano_lm_rm=['nino34m_', tmp.varname, '_model_ano_',  '_l', tmp.lmonth_str, '_rm'];
        tmp.obsvar_ano_lm_rm=['nino34m_', tmp.varname, '_obs_ano_',  '_l', tmp.lmonth_str, '_rm'];
%         tmp.assmvar_ano_lm_rm=['nino34m_', tmp.varname, '_assm_ano_',  '_l', tmp.lmonth_str, '_rm'];

        tmp.md=data.(tmp.modelvar_ano_lm_rm);
        tmp.obs=data.(tmp.obsvar_ano_lm_rm);
%         tmp.assm=data.(tmp.assmvar_ano_lm_rm);
%         tmp.tind= (cfg.iyears-min(cfg.iyears))*12+1+lmonth;
        tmp.corrcoef=corrcoef(tmp.md(isfinite(tmp.obs)), tmp.obs(isfinite(tmp.obs)));
        data.(['nino34m_', tmp.varname, '_corr'])(lmonth+1)=tmp.corrcoef(1,2);

        %% Elnino start case
        tmp.md=data.(tmp.modelvar_ano_lm_rm)(data.Y_Elnino);
        tmp.obs=data.(tmp.obsvar_ano_lm_rm)(data.Y_Elnino);
        tmp.corrcoef=corrcoef(tmp.md(isfinite(tmp.obs)), tmp.obs(isfinite(tmp.obs)));
        data.(['nino34m_', tmp.varname, '_corr_Elnino'])(lmonth+1)=tmp.corrcoef(1,2);

        %% Lanina start case
        tmp.md=data.(tmp.modelvar_ano_lm_rm)(data.Y_Lanina);
        tmp.obs=data.(tmp.obsvar_ano_lm_rm)(data.Y_Lanina);
        tmp.corrcoef=corrcoef(tmp.md(isfinite(tmp.obs)), tmp.obs(isfinite(tmp.obs)));
        data.(['nino34m_', tmp.varname, '_corr_Lanina'])(lmonth+1)=tmp.corrcoef(1,2);

        %% Neutral start case
        tmp.md=data.(tmp.modelvar_ano_lm_rm)(data.Y_Neutral);
        tmp.obs=data.(tmp.obsvar_ano_lm_rm)(data.Y_Neutral);
        tmp.corrcoef=corrcoef(tmp.md(isfinite(tmp.obs)), tmp.obs(isfinite(tmp.obs)));
        data.(['nino34m_', tmp.varname, '_corr_Neutral'])(lmonth+1)=tmp.corrcoef(1,2);
    end
    
    %% lmonth=0 is not computed. computation from lmonth=1.

    %% correlation coefficient of NINO3.4 ind, function of lead year
    cfg.figname=[dirs.figdir, filesep, 'nino34m_v4_hcst_obs_',str.iy_min, '-', str.iy_max, 'm', '_corr_leadm', '.tif'];
%             bar(0:cfg.proj_year*12-1,data.(['nino34m_', tmp.varname, '_corr']),  'linewidth', 2)
    bar(1:cfg.proj_year*12, data.(['nino34m_', tmp.varname, '_corr'])(1:cfg.proj_year*12),  'linewidth', 2);

    xlabel('Lead month'); ylabel(['corr. coef.,', 'Nino3.4' , ', obs']);
    set(gca, 'fontsize', 20);
    grid minor
    ylim([-0.4 1]);
    saveas(gcf,cfg.figname,'tif');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;

    %% correlation coefficient of NINO3.4 ind, function of lead year (Elnino start only)
    cfg.figname=[dirs.figdir, filesep, 'nino34m_v4_hcst_obs_Elnino_',str.iy_min, '-', str.iy_max, 'm', '_corr_leadm', '.tif'];
    bar(1:cfg.proj_year*12, data.(['nino34m_', tmp.varname, '_corr_Elnino'])(1:cfg.proj_year*12),  'linewidth', 2);

    xlabel('Lead month'); ylabel(['corr. coef.,', 'Nino3.4 (E-init)']);
    set(gca, 'fontsize', 20);
    grid minor
    ylim([-0.4 1]);
    saveas(gcf,cfg.figname,'tif');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;

    %% correlation coefficient of NINO3.4 ind, function of lead year (Lanina start only)
    cfg.figname=[dirs.figdir, filesep, 'nino34m_v4_hcst_obs_Lanina_',str.iy_min, '-', str.iy_max, 'm', '_corr_leadm', '.tif'];
    bar(1:cfg.proj_year*12, data.(['nino34m_', tmp.varname, '_corr_Lanina'])(1:cfg.proj_year*12),  'linewidth', 2);

    xlabel('Lead month'); ylabel(['corr. coef.,', 'Nino3.4 (L-init)']);
    set(gca, 'fontsize', 20);
    grid minor
    ylim([-0.4 1]);
    saveas(gcf,cfg.figname,'tif');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;

    %% correlation coefficient of NINO3.4 ind, function of lead year (Neutral start only)
    cfg.figname=[dirs.figdir, filesep, 'nino34m_v4_hcst_obs_Neutral_',str.iy_min, '-', str.iy_max, 'm', '_corr_leadm', '.tif'];
    bar(1:cfg.proj_year*12, data.(['nino34m_', tmp.varname, '_corr_Neutral'])(1:cfg.proj_year*12),  'linewidth', 2);

    xlabel('Lead month'); ylabel(['corr. coef.,', 'Nino3.4 (L-init)']);
    set(gca, 'fontsize', 20);
    grid minor
    ylim([-0.4 1]);
    saveas(gcf,cfg.figname,'tif');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;

    %% correlation coefficient of NINO3.4 ind, function of lead year (All, Elnino, Lanina, Neutral start)
    cfg.figname=[dirs.figdir, filesep, 'nino34m_v4_hcst_obs_All_',str.iy_min, '-', str.iy_max, 'm', '_corr_leadm', '.tif'];
    plot(1:cfg.proj_year*12,data.(['nino34m_', tmp.varname, '_corr'])(1:cfg.proj_year*12),...
        '-o', 'color', 'k', 'linewidth', 2);
    hold on
    plot(1:cfg.proj_year*12,data.(['nino34m_', tmp.varname, '_corr_Elnino'])(1:cfg.proj_year*12),...
        '-o', 'color', 'r', 'linewidth', 2);
    plot(1:cfg.proj_year*12,data.(['nino34m_', tmp.varname, '_corr_Lanina'])(1:cfg.proj_year*12),...
        '-o', 'color', 'b', 'linewidth', 2);
    plot(1:cfg.proj_year*12,data.(['nino34m_', tmp.varname, '_corr_Neutral'])(1:cfg.proj_year*12),...
        '-o', 'color', 'g', 'linewidth', 2);
    xlabel('Lead month'); ylabel(['corr. coef.,', 'Nino3.4']);
    set(gca, 'fontsize', 20);
    grid minor
    ylim([-0.4 1]);
    legend({'All', 'E-start', 'L-start', 'N-start'})
    saveas(gcf,cfg.figname,'tif');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;
    

    %% nRMSE of NINO3.4 ind, function of lead year
    cfg.figname=[dirs.figdir, filesep, 'nino34m_v4_hcst_obs_',str.iy_min, '-', str.iy_max, 'm', '_nRMSE_leadm', '.tif'];
    bar(1:cfg.proj_year*12, data.(['nino34m_', tmp.varname, '_RMSE'])(1:cfg.proj_year*12),  'linewidth', 2);

    xlabel('Lead month'); ylabel(['nRMSE.,', 'Nino3.4']);
    set(gca, 'fontsize', 20);
    grid minor
    ylim([0.2 1.8]);
    saveas(gcf,cfg.figname,'tif');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;
 
    %% nRMSE of NINO3.4 ind, function of lead year (Elnino start only)
    cfg.figname=[dirs.figdir, filesep, 'nino34m_v4_hcst_obs_',str.iy_min, '-', str.iy_max, 'm', '_nRMSE_Elnino_leadm', '.tif'];
    bar(1:cfg.proj_year*12, data.(['nino34m_', tmp.varname, '_RMSE_Elnino'])(1:cfg.proj_year*12),  'linewidth', 2);

    xlabel('Lead month'); ylabel(['nRMSE.,', 'Nino3.4 (E-init)']);
    set(gca, 'fontsize', 20);
    grid minor
    ylim([0.2 1.8]);
    saveas(gcf,cfg.figname,'tif');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;

    %% nRMSE of NINO3.4 ind, function of lead year (Lanina start only)
    cfg.figname=[dirs.figdir, filesep, 'nino34m_v4_hcst_obs_',str.iy_min, '-', str.iy_max, 'm', '_nRMSE_Lanina_leadm', '.tif'];
    bar(1:cfg.proj_year*12, data.(['nino34m_', tmp.varname, '_RMSE_Lanina'])(1:cfg.proj_year*12),  'linewidth', 2);

    xlabel('Lead month'); ylabel(['nRMSE.,', 'Nino3.4 (L-init)']);
    set(gca, 'fontsize', 20);
    grid minor
    ylim([0.2 1.8]);
    saveas(gcf,cfg.figname,'tif');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;

    %% nRMSE of NINO3.4 ind, function of lead year (All, Elnino, Lanina, Neutral start)
    cfg.figname=[dirs.figdir, filesep, 'nino34m_v4_hcst_obs_',str.iy_min, '-', str.iy_max, 'm', '_nRMSE_All_leadm', '.tif'];
    plot(1:cfg.proj_year*12,data.(['nino34m_', tmp.varname, '_RMSE'])(1:cfg.proj_year*12),  ...
        '-o', 'color', 'k', 'linewidth', 2);
    hold on
    plot(1:cfg.proj_year*12,data.(['nino34m_', tmp.varname, '_RMSE_Elnino'])(1:cfg.proj_year*12),  ...
        '-o', 'color', 'r', 'linewidth', 2);
    plot(1:cfg.proj_year*12,data.(['nino34m_', tmp.varname, '_RMSE_Lanina'])(1:cfg.proj_year*12),  ...
        '-o', 'color', 'b', 'linewidth', 2);
    plot(1:cfg.proj_year*12,data.(['nino34m_', tmp.varname, '_RMSE_Neutral'])(1:cfg.proj_year*12),  ...
        '-o', 'color', 'g', 'linewidth', 2);
    xlabel('Lead month'); ylabel(['nRMSE.,', 'Nino3.4']);
    set(gca, 'fontsize', 20);
    grid minor
    ylim([0.2 1.8]);
    legend({'All', 'E-start', 'L-start', 'N-start'})
    saveas(gcf,cfg.figname,'tif');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;
% end

     %% correlation coefficient of NINO3.4 ind, function of lead year (to compare with SMYLE)
    cfg.figname=[dirs.figdir, filesep, 'SMYLE_nino34m_v4_hcst_obs_',str.iy_min, '-', str.iy_max, 'm', '_corr_leadm', '.tif'];
%             bar(0:cfg.proj_year*12-1,data.(['nino34m_', tmp.varname, '_corr']),  'linewidth', 2)
    plot(1:19,data.(['nino34m_', tmp.varname, '_corr'])(1:19),  ...
        '-o', 'color', 'b', 'linewidth', 2);
    hold on
    plot(1:3:19,data.(['nino34m_', tmp.varname, '_corr'])(1:3:19),  ...
        '-o', 'color', 'b', 'linewidth', 2, 'MarkerFaceColor', 'b');
    line([1 19], [0.5 0.5], 'color', 'k', 'linewidth', 2)
    xlabel('Lead month'); ylabel(['corr. coef.,', 'Nino3.4']);
    set(gca, 'fontsize', 20);
    grid minor
    ylim([0.1 1]);
    set(gcf, 'PaperPosition', [0, 0, 10, 4]);
    saveas(gcf,cfg.figname,'tif');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;
    

    %% nRMSE of NINO3.4 ind, function of lead year (to compare with SMYLE)
    cfg.figname=[dirs.figdir, filesep, 'SMYLE_nino34m_v4_hcst_obs_',str.iy_min, '-', str.iy_max, 'm', '_nRMSE_leadm', '.tif'];
%             bar(0:cfg.proj_year*12-1,data.(['nino34m_', tmp.varname, '_corr']),  'linewidth', 2)
    plot(1:19, data.(['nino34m_', tmp.varname, '_RMSE'])(1:19),  ...
        '-o', 'color', 'b','linewidth', 2);
    hold on
    plot(1:3:19, data.(['nino34m_', tmp.varname, '_RMSE'])(1:3:19),  ...
        '-o', 'color', 'b','linewidth', 2, 'MarkerFaceColor', 'b');
    line([1 19], [1.0 1.0], 'color', 'k', 'linewidth', 2)
    xlabel('Lead month'); ylabel(['nRMSE.,', 'Nino3.4']);
    set(gca, 'fontsize', 20);
    grid minor
    ylim([0.2 1.8]);
    set(gcf, 'PaperPosition', [0, 0, 10, 4]);
    saveas(gcf,cfg.figname,'tif');
    RemoveWhiteSpace([], 'file', cfg.figname);
    close all;


% % % % %% spaghetti plot (NINO34)
% % % % system(['mkdir -p ', dirs.figdir]);
% % % % %         for lmonth=0:cfg.proj_year*12-1
% % % %     for lmonth=0:35
% % % %         tmp.lmonth_str=num2str(lmonth, '%03i');
% % % %         tmp.modelvar_ano_lm_rm=['nino34m_', tmp.varname, '_model_ano_',  '_l', tmp.lmonth_str, '_rm'];
% % % %         tmp.obsvar_ano_lm_rm=['nino34m_', tmp.varname, '_obs_ano_',  '_l', tmp.lmonth_str, '_rm'];
% % % %         
% % % %         tmp.md=data.(tmp.modelvar_ano_lm_rm);
% % % %         tmp.obs=data.(tmp.obsvar_ano_lm_rm);
% % % % %         tmp.tind= (cfg.iyears-min(cfg.iyears))*12+1+lmonth;
% % % %         tmp.plottime=cfg.iyears+lmonth*1/12;
% % % %         plot(tmp.plottime, tmp.md, 'b');
% % % %         hold on
% % % %         plot(tmp.plottime, tmp.obs, 'k', 'linewidth',2)
% % % %     end
% % % %     hold off
% % % %     axis tight
% % % %     
% % % % %% correlation coefficient of NINO3.4 ind, function of lead year (3-mon running mean)
% % % % cfg.figname=[dirs.figdir, filesep, 'spaghetti_nino34m_hcst_obs_', str.iy_min, '-', str.iy_max, 'm', '_corr_leadm', '.tif'];
% % % % %             bar(0:cfg.proj_year*12-1,data.(['nino34m_', tmp.varname, '_corr']),  'linewidth', 2)
% % % % 
% % % % xlabel('Year'); ylabel('()');
% % % % set(gca, 'fontsize', 20)
% % % % grid minor
% % % % set(gcf, 'PaperPosition', [0, 0, 8, 4]);
% % % % saveas(gcf,cfg.figname,'tif');
% % % % RemoveWhiteSpace([], 'file', cfg.figname);
% % % % close all;


%% spaghetti plot (NINO34) (one observation line)
tmp.jet=jet(cfg.proj_year*12);
for lmonth=0:cfg.proj_year*12-1
    tmp.lmonth_str=num2str(lmonth, '%03i');
%     tmp.modelvar_ano_lm_rm=['nino34m_', tmp.varname, '_model_ano_',  '_l', tmp.lmonth_str, '_rm'];
    tmp.modelvar_ano_lm_rm=['nino34m_', tmp.varname, '_model_ano_',  '_l', tmp.lmonth_str];
    tmp.plottime=min(cfg.iyears)+1/24+lmonth*1/12 : 1 : max(cfg.iyears)-1/24+lmonth*1/12;
    tmp.md=data.(tmp.modelvar_ano_lm_rm);
    tmp.plottime=cfg.iyears+1/24+lmonth*(1/12);
%     plot(tmp.plottime, tmp.md, 'b');
    plot(tmp.plottime, tmp.md, 'color', tmp.jet(lmonth+1,:));
    hold on
end
cax=colorbar;
colormap(jet);
caxis([0 cfg.proj_year*12-1])
title(cax,'LM');

tmp=rmfield(tmp, 'obs');

for lmonth=0:cfg.proj_year*12-1
%     tmp.obsvar_ano_lm_rm=['nino34m_', tmp.varname, '_obs_ano_',  '_l', tmp.lmonth_str, '_rm'];
    tmp.lmonth_str=num2str(lmonth, '%03i');
    tmp.obsvar_ano_lm_rm=['nino34m_', tmp.varname, '_obs_ano_',  '_l', tmp.lmonth_str];
    tmp.obs(1+lmonth:12:(length(cfg.iyears)-1)*12+1+lmonth)=data.(tmp.obsvar_ano_lm_rm);
end

tmp.plottime=min(cfg.iyears)+1/24:1/12:max(cfg.iyears)-1/24+36*1/12;
plot(tmp.plottime, tmp.obs, 'k', 'linewidth',2)
hold off
axis tight
cfg.figname=[dirs.figdir, filesep, 'spaghetti_nino34m_v4_hcst_obs_oneline', '_', str.iy_min, '-', str.iy_max, '_leadm', '.tif'];
xlabel('Year'); ylabel('Nino 3.4 Index');
set(gca, 'fontsize', 20)
grid minor
set(gcf, 'PaperPosition', [0, 0, 8, 4]);
saveas(gcf,cfg.figname,'tif');
RemoveWhiteSpace([], 'file', cfg.figname);
close all;





%% spaghetti plot (SST raw)
for iyear=min(cfg.iyears):max(cfg.iyears)
    tmp.iyear_str=num2str(iyear);
    tmp.plottime=iyear+1/24 : 1/12 : iyear+cfg.proj_year - 1/24;
    tmp.modelvar=['nino34m_', tmp.varname, '_model_',  'i', tmp.iyear_str];
    tmp.obsvar=['nino34m_', tmp.varname, '_obs_',  'i', tmp.iyear_str];
    
    plot(tmp.plottime,data.(tmp.modelvar)-273.15,'b');
    hold on
    plot(tmp.plottime,data.(tmp.obsvar),'k', 'linewidth',3);
end
axis tight
hold off
xlabel('Year'); ylabel('^oC');
set(gca, 'fontsize', 20)
grid minor
set(gcf, 'PaperPosition', [0, 0, 8, 4]);
cfg.figname=[dirs.figdir, filesep, 'spaghetti_v4_', 'SST_hcst_obs_', str.iy_min, '-', str.iy_max,'_proj_',cfg.proj_year,  'y.tif'];
saveas(gcf,cfg.figname,'tif');
RemoveWhiteSpace([], 'file', cfg.figname);
close all;






% %% normal time series and correlation coefficient (combined plot, hcst&assm)
% for varind=1:length(cfg.varnames)
%     tmp.varname=cfg.var;
%     dirs.figdir= [dirs.figroot, filesep, cfg.casename_m, filesep, cfg.region, filesep, 'Timeseries'];
%     system(['mkdir -p ', dirs.figdir]);
% %         for lmonth=0:cfg.proj_year*12-1
%     for lmonth=0:35
%         tmp.lmonth_str=num2str(lmonth, '%03i');
%         tmp.modelvar_ano_lm_rm=['nino34m_', tmp.varname, '_model_ano_',  '_l', tmp.lmonth_str, '_rm'];
%         tmp.obsvar_ano_lm_rm=['nino34m_', tmp.varname, '_obs_ano_',  '_l', tmp.lmonth_str, '_rm'];
%         tmp.assmvar_ano_lm_rm=['nino34m_', tmp.varname, '_assm_ano_',  '_l', tmp.lmonth_str, '_rm'];
% 
%         tmp.md=data.(tmp.modelvar_ano_lm_rm);
%         tmp.obs=data.(tmp.obsvar_ano_lm_rm);
% %         tmp.assm=data.(tmp.assmvar_ano_lm_rm);
% 
%         tmp.corrcoef=corrcoef(tmp.md(isfinite(tmp.assm)), tmp.assm(isfinite(tmp.assm)));
%         data.(['nino34m_', tmp.varname, '_corr'])(lmonth+1)=tmp.corrcoef(1,2);
%     end
% 
%     %% correlation coefficient of NINO3.4 ind, function of lead year (3-mon running mean)
%         cfg.figname=[dirs.figdir, filesep, 'nino34m_hcst_assm_91-20m', '_corr_leadm', '.tif'];
% %             bar(0:cfg.proj_year*12-1,data.(['nino34m_', tmp.varname, '_corr']),  'linewidth', 2)
%         bar(1:36,data.(['nino34m_', tmp.varname, '_corr'])(1:36),  'linewidth', 2)
% 
%         xlabel('Lead month'); ylabel(['corr. coef.,', 'Nino3.4' , ', assm-ba']);
%         set(gca, 'fontsize', 20)
%         grid minor
%         ylim([-0.4 1])
%         saveas(gcf,cfg.figname,'tif');
%         RemoveWhiteSpace([], 'file', cfg.figname);
%         close all;
%     
% end




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