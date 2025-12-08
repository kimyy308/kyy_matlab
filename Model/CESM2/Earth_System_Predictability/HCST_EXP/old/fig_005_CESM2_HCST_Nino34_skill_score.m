% %  Created 28-Oct-2022 by Yong-Yub Kim
% %  Created 29-Oct-2022 by Yong-Yub Kim
% %  Created 29-Oct-2022 by Yong-Yub Kim addition of MSSS

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
            for lyear=0:config.proj_year-1
                tmp.lyear_str=num2str(lyear, '%02i');
                tmp.modelvar_ly=['nino34m_', tmp.varname, '_model_', tmp.obsname_simple, '_l', tmp.lyear_str];
                tmp.ydata= ncread(config.datafilename, ['nino34m_', tmp.varname], [(lyear)*12+1], 12);
                tmp.yind=(iyear-min(config.iyears))*12+lyear*12+1:(iyear-min(config.iyears))*12+lyear*12+12;
                data.(tmp.modelvar_ly)(tmp.yind)=tmp.ydata;
                data.(tmp.modelvar_ly)(data.(tmp.modelvar_ly)==0)=NaN;
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

    %% set time value for figure
    data.time_vec_extended=data.time_vec;
    for lyear=1:config.proj_year-1
        tmp.endind=size(data.time_vec_extended,1);
        data.time_vec_extended(tmp.endind+1:tmp.endind+12,1)=max(config.iyears)+lyear;
        data.time_vec_extended(tmp.endind+1:tmp.endind+12,2)=data.time_vec_extended(1:12,2);
        data.time_vec_extended(tmp.endind+1:tmp.endind+12,3)=data.time_vec_extended(1:12,3);
    end
    data.time_leap_extended=datenum(data.time_vec_extended);

    %% get tind of mean period to get anomaly
    tmp.tind_min=find(data.time_vec_extended(:,1)==1969, 1, 'first');
    tmp.tind_max=find(data.time_vec_extended(:,1)==2021, 1, 'last');
    for obsind=1:length(config.obsnames)
        tmp.obsname=config.obsnames{obsind};
        tmp.obsname_simple= f_obs_simple(tmp.obsname);
        
        tmp.obsvar=['nino34m_', tmp.varname, '_obs_', tmp.obsname_simple];
        tmp.obsvar_ltm_ly=['nino34m_', tmp.varname, '_obs_ltm_', tmp.obsname_simple]; %long-term mean name
        data.(tmp.obsvar_ltm_ly)=mean(data.(tmp.obsvar)(tmp.tind_min:tmp.tind_max));
        tmp.obsvaranom=['nino34m_', tmp.varname, '_obs_anom_', tmp.obsname_simple];
        data.(tmp.obsvaranom)=data.(tmp.obsvar)(tmp.tind_min:tmp.tind_max)-data.(tmp.obsvar_ltm_ly);
        tmp.mse_o_obs=['nino34m_', tmp.varname, '_mse_o_obs_', tmp.obsname_simple];
        data.(tmp.mse_o_obs)=mean(data.(tmp.obsvaranom).^2);

        tmp.assmvar=['nino34m_', tmp.varname, '_assm_', tmp.obsname_simple];
        tmp.assmvar_ltm_ly=['nino34m_', tmp.varname, '_assm_ltm_', tmp.obsname_simple]; %long-term mean name
        data.(tmp.assmvar_ltm_ly)=mean(data.(tmp.assmvar)(tmp.tind_min:tmp.tind_max));
        tmp.assmvaranom=['nino34m_', tmp.varname, '_assm_anom_', tmp.obsname_simple];
        data.(tmp.assmvaranom)=data.(tmp.assmvar)(tmp.tind_min:tmp.tind_max)-data.(tmp.assmvar_ltm_ly);
        tmp.mse_o_assm=['nino34m_', tmp.varname, '_mse_o_assm_', tmp.obsname_simple];
        data.(tmp.mse_o_assm)=mean(data.(tmp.assmvaranom).^2);

        for lyear=0:config.proj_year-1
            tmp.lyear_str=num2str(lyear, '%02i');
            tmp.modelvar_ly=['nino34m_', tmp.varname, '_model_', tmp.obsname_simple, '_l', tmp.lyear_str]; % var == temp fixed
            tmp.modelvar_ltm_ly=['nino34m_', tmp.varname, '_model_ltm_', tmp.obsname_simple, '_l', tmp.lyear_str]; %long-term mean name
            data.(tmp.modelvar_ltm_ly)=mean(data.(tmp.modelvar_ly)(tmp.tind_min:tmp.tind_max));
            tmp.modelvar_anom_ly=['nino34m_', tmp.varname, '_model_anom_', tmp.obsname_simple, '_l', tmp.lyear_str];
            data.(tmp.modelvar_anom_ly)=data.(tmp.modelvar_ly)(tmp.tind_min:tmp.tind_max) - data.(tmp.modelvar_ltm_ly);
            
            tmp.mse_h_obs_ly=['nino34m_', tmp.varname, '_mse_h_obs_', tmp.obsname_simple, '_l', tmp.lyear_str];
            data.(tmp.mse_h_obs_ly)=mean((data.(tmp.modelvar_anom_ly) - data.(tmp.obsvaranom)).^2);
            tmp.msss_obs_ly=['nino34m_', tmp.varname, '_msss_obs_', tmp.obsname_simple, '_l', tmp.lyear_str];
            data.(tmp.msss_obs_ly)=1- data.(tmp.mse_h_obs_ly)/data.(tmp.mse_o_obs);
            
            tmp.mse_h_assm_ly=['nino34m_', tmp.varname, '_mse_h_assm_', tmp.obsname_simple, '_l', tmp.lyear_str];
            data.(tmp.mse_h_assm_ly)=mean((data.(tmp.modelvar_anom_ly) - data.(tmp.assmvaranom)).^2);
            tmp.msss_assm_ly=['nino34m_', tmp.varname, '_msss_assm_', tmp.obsname_simple, '_l', tmp.lyear_str];
            data.(tmp.msss_assm_ly)=1- data.(tmp.mse_h_assm_ly)/data.(tmp.mse_o_obs);
        end


    end

%     data.time_y=config.iyears;
%     data.time_y_extended=[config.iyears, max(config.iyears)+1:max(config.iyears)+config.proj_year-1];

 %% normal time series and correlation coefficient (combined plot, hcst&obs)
    for varind=1:length(config.varnames)
        tmp.varname=config.varnames{varind};
        dirs.figdir= [dirs.figroot, filesep, config.casename_m, filesep, 'NINO34', filesep, 'Timeseries'];
        system(['mkdir -p ', dirs.figdir]);
        for lyear=0:config.proj_year-1
            tmp.lyear_str=num2str(lyear, '%02i');
            tmp.md=data.(['nino34m_', tmp.varname, '_model_', tmp.obsname_simple, '_l', tmp.lyear_str]) ...
                - data.(['nino34m_', tmp.varname, '_model_ltm_', tmp.obsname_simple, '_l', tmp.lyear_str]);
            tmp.obs=data.(['nino34m_', tmp.varname, '_obs_', tmp.obsname_simple]) ...
                - data.(['nino34m_', tmp.varname, '_obs_ltm_', tmp.obsname_simple]);
            
            subplot(5,2, lyear+1)
            plot(data.time_leap_extended(1:length(tmp.md)),tmp.md,  '-o', 'color', 'g', 'linewidth', 2)
            xlim([min(data.time_leap_extended), max(data.time_leap_extended)])
            datetick('x', 'yymmm', 'keeplimits')
%             switch tmp.varname
%                 case 'temp'
%                     ylim([18.2 19.8])
%                 case 'salt'
%                     ylim([34.3 34.4])
%             end
            hold on;
            plot(data.time_leap, tmp.obs,  '-^', 'color', 'k', 'linewidth', 2)

            tmp.corrcoef=corrcoef(tmp.md(isfinite(tmp.md(1:length(data.time_leap)))), tmp.obs(isfinite(tmp.md(1:length(data.time_leap)))));
            data.(['nino34m_', tmp.varname, '_corr'])(lyear+1)=tmp.corrcoef(1,2);

%             xlabel('Year'); 
            ylabel([tmp.varname, '-l', tmp.lyear_str]); 
            legend(['HCST-l', tmp.lyear_str, ', R=', num2str(round(tmp.corrcoef(1,2),2))],...
                [tmp.obsname_simple, '(obs-ba)'],...
                'Location', 'NorthWest');
%             set(gca, 'fontsize', 20)
            grid minor
            hold off
                    
            set(gcf, 'PaperPosition', [0, 0, 10, 15]);
        
%             close all;
        end
        config.figname=[dirs.figdir, filesep, 'nino34m_hcst_obs_', tmp.varname, '.tif'];
        saveas(gcf,config.figname,'tif');
        RemoveWhiteSpace([], 'file', config.figname);
        close all;

        %% correlation coefficient of NINO3.4 ind, function of lead year
            config.figname=[dirs.figdir, filesep, 'nino34m_hcst_obs_', tmp.varname, '_corr', '.tif'];
%             plot(0:config.proj_year-1,data.(['nino34m_', tmp.varname, '_corr']),  '-o', 'linewidth', 2)
            bar(0:config.proj_year-1,data.(['nino34m_', tmp.varname, '_corr']),  'linewidth', 2)

            xlabel('Lead year'); ylabel(['corr. coef.,', tmp.varname, ', obs-ba']);
            set(gca, 'fontsize', 20)
            grid minor
            ylim([-1 1])
            saveas(gcf,config.figname,'tif');
            RemoveWhiteSpace([], 'file', config.figname);
            close all;

        %% correlation coefficient of SST, function of lead year (1, 2, 3-4, 5-9)
            config.figname=[dirs.figdir, filesep, 'nino34m_hcst_obs_', tmp.varname, '_corr_2', '.tif'];
            tmp.timelabel={'0','1','2', '3~4', '5~9'};
            tmp.data(1)=data.(['nino34m_', tmp.varname, '_corr'])(1);
            tmp.data(2)=data.(['nino34m_', tmp.varname, '_corr'])(2);
            tmp.data(3)=data.(['nino34m_', tmp.varname, '_corr'])(3);
            tmp.data(4)=mean(data.(['nino34m_', tmp.varname, '_corr'])(4:5));
            tmp.data(5)=mean(data.(['nino34m_', tmp.varname, '_corr'])(6:10));
            bar(1:5,tmp.data,  'linewidth', 2)
            xticklabels(tmp.timelabel);
            xlabel('Lead year'); 
%             ylabel(['ACC,', tmp.varname, ', obs-ba']);
            ylabel(['ACC']);
            set(gca, 'fontsize', 20)
            grid minor
            ylim([0 1])
            saveas(gcf,config.figname,'tif');
            RemoveWhiteSpace([], 'file', config.figname);
            close all;
        
    end


    %% normal time series and correlation coefficient (combined plot, hcst&assm)
    for varind=1:length(config.varnames)
        tmp.varname=config.varnames{varind};
        dirs.figdir= [dirs.figroot, filesep, config.casename_m, filesep, 'NINO34', filesep, 'Timeseries'];
        system(['mkdir -p ', dirs.figdir]);
        for lyear=0:config.proj_year-1
            tmp.lyear_str=num2str(lyear, '%02i');
            tmp.md=data.(['nino34m_', tmp.varname, '_model_', tmp.obsname_simple, '_l', tmp.lyear_str]) ...
                - data.(['nino34m_', tmp.varname, '_model_ltm_', tmp.obsname_simple, '_l', tmp.lyear_str]);
            tmp.assm=data.(['nino34m_', tmp.varname, '_assm_', tmp.obsname_simple]) ...
                - data.(['nino34m_', tmp.varname, '_assm_ltm_', tmp.obsname_simple]);
            
            subplot(5,2, lyear+1)
            plot(data.time_leap_extended(1:length(tmp.md)),tmp.md,  '-o', 'color', 'g', 'linewidth', 2)
            xlim([min(data.time_leap_extended), max(data.time_leap_extended)])
            datetick('x', 'yymmm', 'keeplimits')
%             switch tmp.varname
%                 case 'temp'
%                     ylim([18.2 19.8])
%                 case 'salt'
%                     ylim([34.3 34.4])
%             end
            hold on;
            plot(data.time_leap, tmp.assm,  '-^', 'color', 'k', 'linewidth', 2)

            tmp.corrcoef=corrcoef(tmp.md(isfinite(tmp.md(1:length(data.time_leap)))), tmp.assm(isfinite(tmp.md(1:length(data.time_leap)))));
            data.(['nino34m_', tmp.varname, '_corr'])(lyear+1)=tmp.corrcoef(1,2);

%             xlabel('Year'); 
            ylabel([tmp.varname, '-l', tmp.lyear_str]); 
            legend(['HCST-l', tmp.lyear_str, ', R=', num2str(round(tmp.corrcoef(1,2),2))],...
                [tmp.obsname_simple, '(assm-ba)'],...
                'Location', 'NorthWest');
%             set(gca, 'fontsize', 20)
            grid minor
            hold off
                    
            set(gcf, 'PaperPosition', [0, 0, 10, 15]);
        
%             close all;
        end
        config.figname=[dirs.figdir, filesep, 'nino34m_hcst_assm_', tmp.varname, '.tif'];
        saveas(gcf,config.figname,'tif');
        RemoveWhiteSpace([], 'file', config.figname);
        close all;

        %% correlation coefficient of NINO3.4 ind, function of lead year
            config.figname=[dirs.figdir, filesep, 'nino34m_hcst_assm_', tmp.varname, '_corr', '.tif'];
%             plot(0:config.proj_year-1,data.(['nino34m_', tmp.varname, '_corr']),  '-o', 'linewidth', 2)
            bar(0:config.proj_year-1,data.(['nino34m_', tmp.varname, '_corr']),  'linewidth', 2)

            xlabel('Lead year'); ylabel(['corr. coef.,', tmp.varname, ', assm-ba']);
            set(gca, 'fontsize', 20)
            grid minor
            ylim([-1 1])
            saveas(gcf,config.figname,'tif');
            RemoveWhiteSpace([], 'file', config.figname);
            close all;

        %% correlation coefficient of SST, function of lead year (0, 1, 2, 3-4, 5-9)
            config.figname=[dirs.figdir, filesep, 'nino34m_hcst_assm_', tmp.varname, '_corr_2', '.tif'];
            tmp.timelabel={'0','1','2', '3~4', '5~9'};
            tmp.data(1)=data.(['nino34m_', tmp.varname, '_corr'])(1);
            tmp.data(2)=data.(['nino34m_', tmp.varname, '_corr'])(2);
            tmp.data(3)=data.(['nino34m_', tmp.varname, '_corr'])(3);
            tmp.data(4)=mean(data.(['nino34m_', tmp.varname, '_corr'])(4:5));
            tmp.data(5)=mean(data.(['nino34m_', tmp.varname, '_corr'])(6:10));
            bar(1:5,tmp.data,  'linewidth', 2)
            xticklabels(tmp.timelabel);
            xlabel('Lead year'); 
%             ylabel(['ACC,', tmp.varname, ', assm-ba']);
            ylabel(['ACC']);
            set(gca, 'fontsize', 20)
            grid minor
            ylim([0 1])
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





% %                 % % %         cmems poly1 fitting yearly
% %                 tic;
% %                 for i=1:cmems_lonsize_cut
% %                     for j=1:cmems_latsize_cut
% %                         if isfinite(cmems_sla_yearly(i,j,1))
% %                             [f_exp, gof_exp] = fit(trendtime_yearly',squeeze(cmems_sla_yearly(i,j,:)),'poly1');
% %                             cmems_sla_yearly_poly1_fit(i,j,1:length(trendtime_yearly))=f_exp(trendtime_yearly);
% %                             cmems_sla_yearly_poly1_fit_rsquare(i,j)=gof_exp.rsquare;
% %                             cmems_sla_yearly_poly1_fit_rmse(i,j)=gof_exp.rmse;
% %                         else
% %                             cmems_sla_yearly_poly1_fit(i,j,1:length(trendtime_yearly))=NaN;
% %                             cmems_sla_yearly_poly1_fit_rsquare(i,j)=NaN;
% %                             cmems_sla_yearly_poly1_fit_rmse(i,j)=NaN;
% %                         end
% % %                         p=polyfit(trendtime_yearly,squeeze(cmems_sla_yearly(i,j,:))',1);
% % %                         cmems_trend_yearly(i,j)=p(1) * 1000.0 ;
% %                     end
% %                 end
% %                 disp('cmems poly1 fitting yearly complete') 
% %                 toc;
% %                 
% %                 for i=1:cmems_lonsize_cut
% %                     for j=1:cmems_latsize_cut
% %                         for k=1:size(cmems_sla_yearly,3)
% %                             cmems_sla_yearly_poly1_fit_detrended(i,j,k)=cmems_sla_yearly(i,j,k)-cmems_sla_yearly_poly1_fit(i,j,k);
% %                         end
% %                     end
% %                 end



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