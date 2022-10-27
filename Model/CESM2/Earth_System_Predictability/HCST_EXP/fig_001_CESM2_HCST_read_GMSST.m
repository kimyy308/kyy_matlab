% %  Created 05-Oct-2022 by Yong-Yub Kim
% %  Created 17-Oct-2022 by Yong-Yub Kim

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
% config.varnames={'TEMP','SALT'};
config.varnames={'temp', 'salt'};
config.len_t_y = length(config.iyears);
config.len_t_m = length(config.months);
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
        dirs.datadir= [dirs.hcstroot, filesep, config.casename_m, filesep, 'GMSV'];
        config.datafilename=[dirs.datadir, filesep, 'GMSV_', config.casename, '.nc'];

        for varind=1:length(config.varnames)
            tmp.varname=config.varnames{varind};
%             grid.(['time_i', tmp.iyear_str])=NaN(1,grid.ntime);
            %% variables initialization
%             data.([tmp.varname, '_model_', tmp.obsname_simple, '_i', tmp.iyear_str])=NaN(grid.nlon, grid.nlat, grid.ntime);
%             data.([tmp.varname, '_bias_', tmp.obsname_simple, '_i', tmp.iyear_str])=NaN(grid.nlon, grid.nlat, grid.ntime);
%             data.([tmp.varname, '_obs_', tmp.obsname_simple])=NaN(grid.nlon, grid.nlat, grid.ntime);
%             data.([tmp.varname, '_sq_err_', tmp.obsname_simple, '_i', tmp.iyear_str])=NaN(grid.nlon, grid.nlat, grid.ntime);
%             data.([tmp.varname, '_rmse_', tmp.obsname_simple, '_i', tmp.iyear_str])=NaN(grid.nlon, grid.nlat);

%             data.(['gm_', tmp.varname, '_model_', tmp.obsname_simple, '_i', tmp.iyear_str])=ncread(config.datafilename, ['gm_', tmp.varname]);
%             data.(['gm_', tmp.varname, '_bias_', tmp.obsname_simple, '_i', tmp.iyear_str])=ncread(config.datafilename, ['gm_bias_', tmp.varname]);
%             data.(['gm_', tmp.varname, '_obs_', tmp.obsname_simple, '_i', tmp.iyear_str])=ncread(config.datafilename, ['gm_obs_', tmp.varname]);
%             data.(['gm_', tmp.varname, '_sq_err_', tmp.obsname_simple, '_i', tmp.iyear_str])=ncread(config.datafilename, ['gm_sq_err_', tmp.varname]);
            data.time=ncread(config.datafilename, 'time');

%             tmp.corrcoef= ...
%                 corrcoef(data.(['gm_', tmp.varname, '_model_', tmp.obsname_simple, '_i', tmp.iyear_str]), ...
%                 data.(['gm_', tmp.varname, '_obs_', tmp.obsname_simple, '_i', tmp.iyear_str]));
%             data.(['corr_gm_', tmp.varname])(iyear-min(config.iyears)+1)= tmp.corrcoef(1,2);
%             data.(['gm_', tmp.varname, '_rmse_', tmp.obsname_simple, '_i', tmp.iyear_str])= ncread(config.datafilename, ['gm_rmse_', tmp.varname]);

            for lyear=0:config.proj_year-1
                tmp.lyear_str=num2str(lyear, '%02i');
                tmp.ymean= mean(ncread(config.datafilename, ['gm_', tmp.varname], [(lyear)*12+1], 12));
                data.(['gm_', tmp.varname, '_model_', tmp.obsname_simple, '_l', tmp.lyear_str])(iyear-min(config.iyears)+1+lyear)= tmp.ymean;
                data.(['gm_', tmp.varname, '_model_', tmp.obsname_simple, '_l', tmp.lyear_str])(data.(['gm_', tmp.varname, '_model_', tmp.obsname_simple, '_l', tmp.lyear_str])==0)=NaN;
            end
            tmp.ymean= mean(ncread(config.datafilename, ['gm_obs_', tmp.varname], [1], 12));
            data.(['gm_', tmp.varname, '_obs_', tmp.obsname_simple])(iyear-min(config.iyears)+1)= tmp.ymean;
            tmp.ymean= mean(ncread(config.datafilename, ['gm_assm_', tmp.varname], [1], 12));
            data.(['gm_', tmp.varname, '_assm_', tmp.obsname_simple])(iyear-min(config.iyears)+1)= tmp.ymean;            
        end
%% save ncfile

%     [f_exp, gof_exp] = fit(trendtime_yearly',squeeze(cmems_sla_yearly(i,j,:)),'poly1');
%                             cmems_sla_yearly_poly1_fit(i,j,1:length(trendtime_yearly))=f_exp(trendtime_yearly);

    end
    disp('abc')
    
    data.time_y=config.iyears;
    data.time_y_extended=[config.iyears, max(config.iyears)+1:max(config.iyears)+config.proj_year-1];


%     %% normal time series and correlation coefficient (for each plot)
%     for varind=1:length(config.varnames)
%         tmp.varname=config.varnames{varind};
%         dirs.figdir= [dirs.figroot, filesep, config.casename_m, filesep, 'GMSV', filesep, 'Timeseries'];
%         system(['mkdir -p ', dirs.figdir]);
%         for lyear=0:config.proj_year-1
%             tmp.lyear_str=num2str(lyear, '%02i');
%             config.figname=[dirs.figdir, filesep, 'gm_', tmp.varname, '_l', tmp.lyear_str, '.tif'];
%             tmp.md=data.(['gm_', tmp.varname, '_model_', tmp.obsname_simple, '_l', tmp.lyear_str]);
%             tmp.obs=data.(['gm_', tmp.varname, '_obs_', tmp.obsname_simple]);
%             
%             plot(data.time_y_extended(1:length(tmp.md)),tmp.md,  '-o', 'color', 'g', 'linewidth', 2)
%             xlim([min(data.time_y_extended), max(data.time_y_extended)])
%             switch tmp.varname
%                 case 'temp'
%                     ylim([18.2 19.8])
%                 case 'salt'
%                     ylim([34.3 34.4])
%             end
%             hold on;
%             plot(data.time_y, tmp.obs,  '-^', 'color', 'k', 'linewidth', 2)
% 
%             tmp.corrcoef=corrcoef(tmp.md(isfinite(tmp.md(1:length(data.time_y)))), tmp.obs(isfinite(tmp.md(1:length(data.time_y)))));
%             data.(['gm_', tmp.varname, '_corr'])(lyear+1)=tmp.corrcoef(1,2);
% 
%             xlabel('Year'); ylabel([tmp.varname]); 
%             legend(['HCST-l', tmp.lyear_str, ', R=', num2str(round(tmp.corrcoef(1,2),2))],...
%                 tmp.obsname_simple,...
%                 'Location', 'NorthWest');
%             set(gca, 'fontsize', 20)
%             grid minor
%             hold off
%                     
%             set(gcf, 'PaperPosition', [0, 0, 8, 4]);
%         
%             saveas(gcf,config.figname,'tif');
%             close all;
%         end
% 
%         %% correlation coefficient of SST, function of lead year
%             config.figname=[dirs.figdir, filesep, 'gm_', tmp.varname, '_corr', '.tif'];
% %             plot(0:config.proj_year-1,data.(['gm_', tmp.varname, '_corr']),  '-o', 'linewidth', 2)
%             bar(0:config.proj_year-1,data.(['gm_', tmp.varname, '_corr']),  'linewidth', 2)
% 
%             xlabel('Lead year'); ylabel(['corr. coef.,', tmp.varname]);
%             set(gca, 'fontsize', 20)
%             grid minor
%             ylim([-1 1])
%             saveas(gcf,config.figname,'tif');
%             close all;
% 
%     end

 %% normal time series and correlation coefficient (combined plot, hcst&obs)
    for varind=1:length(config.varnames)
        tmp.varname=config.varnames{varind};
        dirs.figdir= [dirs.figroot, filesep, config.casename_m, filesep, 'GMSV', filesep, 'Timeseries'];
        system(['mkdir -p ', dirs.figdir]);
        for lyear=0:config.proj_year-1
            tmp.lyear_str=num2str(lyear, '%02i');
            tmp.md=data.(['gm_', tmp.varname, '_model_', tmp.obsname_simple, '_l', tmp.lyear_str]);
            tmp.obs=data.(['gm_', tmp.varname, '_obs_', tmp.obsname_simple]);
            
            subplot(5,2, lyear+1)
            plot(data.time_y_extended(1:length(tmp.md)),tmp.md,  '-o', 'color', 'g', 'linewidth', 2)
            xlim([min(data.time_y_extended), max(data.time_y_extended)])
            switch tmp.varname
                case 'temp'
                    ylim([18.2 19.8])
                case 'salt'
                    ylim([34.3 34.4])
            end
            hold on;
            plot(data.time_y, tmp.obs,  '-^', 'color', 'k', 'linewidth', 2)

            tmp.corrcoef=corrcoef(tmp.md(isfinite(tmp.md(1:length(data.time_y)))), tmp.obs(isfinite(tmp.md(1:length(data.time_y)))));
            data.(['gm_', tmp.varname, '_corr'])(lyear+1)=tmp.corrcoef(1,2);

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
        config.figname=[dirs.figdir, filesep, 'gm_hcst_obs_', tmp.varname, '.tif'];
        saveas(gcf,config.figname,'tif');
        RemoveWhiteSpace([], 'file', config.figname);
        close all;

        %% correlation coefficient of SST, function of lead year
            config.figname=[dirs.figdir, filesep, 'gm_hcst_obs_', tmp.varname, '_corr', '.tif'];
%             plot(0:config.proj_year-1,data.(['gm_', tmp.varname, '_corr']),  '-o', 'linewidth', 2)
            bar(0:config.proj_year-1,data.(['gm_', tmp.varname, '_corr']),  'linewidth', 2)

            xlabel('Lead year'); ylabel(['corr. coef.,', tmp.varname, ', obs-ba']);
            set(gca, 'fontsize', 20)
            grid minor
            ylim([-1 1])
            saveas(gcf,config.figname,'tif');
            RemoveWhiteSpace([], 'file', config.figname);
            close all;

        %% correlation coefficient of SST, function of lead year (1, 2, 3-4, 5-9)
            config.figname=[dirs.figdir, filesep, 'gm_hcst_obs_', tmp.varname, '_corr_2', '.tif'];
            tmp.timelabel={'1','2', '3~4', '5~9'};
            tmp.data(1)=data.(['gm_', tmp.varname, '_corr'])(2);
            tmp.data(2)=data.(['gm_', tmp.varname, '_corr'])(3);
            tmp.data(3)=mean(data.(['gm_', tmp.varname, '_corr'])(4:5));
            tmp.data(4)=mean(data.(['gm_', tmp.varname, '_corr'])(6:10));
            bar(1:4,tmp.data,  'linewidth', 2)
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
        dirs.figdir= [dirs.figroot, filesep, config.casename_m, filesep, 'GMSV', filesep, 'Timeseries'];
        system(['mkdir -p ', dirs.figdir]);
        for lyear=0:config.proj_year-1
            tmp.lyear_str=num2str(lyear, '%02i');
            tmp.md=data.(['gm_', tmp.varname, '_model_', tmp.obsname_simple, '_l', tmp.lyear_str]);
            tmp.obs=data.(['gm_', tmp.varname, '_assm_', tmp.obsname_simple]);
            
            subplot(5,2, lyear+1)
            plot(data.time_y_extended(1:length(tmp.md)),tmp.md,  '-o', 'color', 'g', 'linewidth', 2)
            xlim([min(data.time_y_extended), max(data.time_y_extended)])
            switch tmp.varname
                case 'temp'
                    ylim([18.2 19.8])
                case 'salt'
                    ylim([34.3 34.4])
            end
            hold on;
            plot(data.time_y, tmp.obs,  '-^', 'color', 'k', 'linewidth', 2)

            tmp.corrcoef=corrcoef(tmp.md(isfinite(tmp.md(1:length(data.time_y)))), tmp.obs(isfinite(tmp.md(1:length(data.time_y)))));
            data.(['gm_', tmp.varname, '_corr'])(lyear+1)=tmp.corrcoef(1,2);

%             xlabel('Year'); 
            ylabel([tmp.varname, '-l', tmp.lyear_str]); 
            legend(['HCST-l', tmp.lyear_str, ', R=', num2str(round(tmp.corrcoef(1,2),2))],...
                [tmp.obsname_simple, '(assm)'],...
                'Location', 'NorthWest');
%             set(gca, 'fontsize', 20)
            grid minor
            hold off
                    
            set(gcf, 'PaperPosition', [0, 0, 10, 15]);
        
%             close all;
        end
        config.figname=[dirs.figdir, filesep, 'gm_hcst_assm_', tmp.varname, '.tif'];
        saveas(gcf,config.figname,'tif');
        RemoveWhiteSpace([], 'file', config.figname);
        close all;

        %% correlation coefficient of SST, function of lead year
            config.figname=[dirs.figdir, filesep, 'gm_hcst_assm_', tmp.varname, '_corr', '.tif'];
%             plot(0:config.proj_year-1,data.(['gm_', tmp.varname, '_corr']),  '-o', 'linewidth', 2)
            bar(0:config.proj_year-1,data.(['gm_', tmp.varname, '_corr']),  'linewidth', 2)

            xlabel('Lead year'); ylabel(['corr. coef.,', tmp.varname, ', assm']);
            set(gca, 'fontsize', 20)
            grid minor
            ylim([-1 1])
            saveas(gcf,config.figname,'tif');
            RemoveWhiteSpace([], 'file', config.figname);
            close all;

        %% correlation coefficient of SST, function of lead year (1, 2, 3-4, 5-9)
            config.figname=[dirs.figdir, filesep, 'gm_hcst_assm_', tmp.varname, '_corr_2', '.tif'];
            tmp.timelabel={'1','2', '3~4', '5~9'};
            tmp.data(1)=data.(['gm_', tmp.varname, '_corr'])(2);
            tmp.data(2)=data.(['gm_', tmp.varname, '_corr'])(3);
            tmp.data(3)=mean(data.(['gm_', tmp.varname, '_corr'])(4:5));
            tmp.data(4)=mean(data.(['gm_', tmp.varname, '_corr'])(6:10));
            bar(1:4,tmp.data,  'linewidth', 2)
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


%     %% detrended time series and correlation coefficient
%     for varind=1:length(config.varnames)
%         tmp.varname=config.varnames{varind};
%         dirs.figdir= [dirs.figroot, filesep, config.casename_m, filesep, 'GMSV', filesep, 'Timeseries_det'];
%         system(['mkdir -p ', dirs.figdir]);
%         for lyear=0:config.proj_year-1
%             tmp.lyear_str=num2str(lyear, '%02i');
%             config.figname=[dirs.figdir, filesep, 'gm_', tmp.varname, '_det_l', tmp.lyear_str, '.tif'];
%             tmp.md=data.(['gm_', tmp.varname, '_model_', tmp.obsname_simple, '_l', tmp.lyear_str]);
%             tmp.obs=data.(['gm_', tmp.varname, '_obs_', tmp.obsname_simple]);
%             
%             tmp.md_finite=tmp.md(isfinite(tmp.md));
%             [f_exp, gof_exp] = fit(data.time_y_extended(isfinite(tmp.md))',squeeze(tmp.md_finite)','poly1');
%             tmp.md_det_finite=f_exp(data.time_y_extended(isfinite(tmp.md)));
%             tmp.md_det=tmp.md;
%             tmp.md_det(isfinite(tmp.md_det))=tmp.md_det(isfinite(tmp.md_det))' - tmp.md_det_finite;
%             
%             tmp.obs_finite=tmp.obs(isfinite(tmp.obs));
%             [f_exp, gof_exp] = fit(data.time_y_extended(isfinite(tmp.obs))',squeeze(tmp.obs_finite)','poly1');
%             tmp.obs_det_finite=f_exp(data.time_y_extended(isfinite(tmp.obs)));
%             tmp.obs_det=tmp.obs;
%             tmp.obs_det(isfinite(tmp.obs_det))=tmp.obs_det(isfinite(tmp.obs_det))' - tmp.obs_det_finite;
% 
% 
%             plot(data.time_y_extended(1:length(tmp.md_det)),tmp.md_det, '-o', 'color', 'g', 'linewidth', 2)
%             xlim([min(data.time_y_extended), max(data.time_y_extended)])
%             switch tmp.varname
%                 case 'temp'
%                     ylim([-0.4 0.4])
%                 case 'salt'
%                     ylim([-0.05 0.05])
%             end
%             hold on;
%             plot(data.time_y, tmp.obs_det,  '-^', 'color', 'k', 'linewidth', 2)
%             tmp.corrcoef=corrcoef(tmp.md_det(isfinite(tmp.md_det(1:length(data.time_y)))), tmp.obs_det(isfinite(tmp.md_det(1:length(data.time_y)))));
%             data.(['gm_', tmp.varname, '_det_corr'])(lyear+1)=tmp.corrcoef(1,2);
%             
%             xlabel('Year'); ylabel([tmp.varname, '-detrended']); 
%             legend(['HCST-det-l', tmp.lyear_str, ', R=', num2str(round(tmp.corrcoef(1,2),2))],...
%                 [tmp.obsname_simple, '-det'],...
%                 'Location', 'NorthWest');
%             set(gca, 'fontsize', 20)
%             grid minor
%             hold off
%                     
%             set(gcf, 'PaperPosition', [0, 0, 8, 4]);
%         
%             saveas(gcf,config.figname,'tif');
%             close all;
%         end
%         
%         %% correlation coefficient of SST, function of lead year
%         config.figname=[dirs.figdir, filesep, 'gm_', tmp.varname, '_det_corr', '.tif'];
% %         plot(0:config.proj_year-1,data.(['gm_', tmp.varname, '_det_corr']),  '-o', 'linewidth', 2)
%         bar(0:config.proj_year-1,data.(['gm_', tmp.varname, '_det_corr']),  'linewidth', 2)
%         
%         xlabel('Lead year'); ylabel(['corr. coef.,', tmp.varname]);
%         set(gca, 'fontsize', 20)
%         grid minor
%         ylim([-1 1])
%         saveas(gcf,config.figname,'tif');
%         close all;
% 
%     end


 %% detrended time series and correlation coefficient (combined plot, hcst&obs)
     for varind=1:length(config.varnames)
        tmp.varname=config.varnames{varind};
        dirs.figdir= [dirs.figroot, filesep, config.casename_m, filesep, 'GMSV', filesep, 'Timeseries_det'];
        system(['mkdir -p ', dirs.figdir]);
        for lyear=0:config.proj_year-1
            tmp.lyear_str=num2str(lyear, '%02i');
            tmp.md=data.(['gm_', tmp.varname, '_model_', tmp.obsname_simple, '_l', tmp.lyear_str]);
            tmp.obs=data.(['gm_', tmp.varname, '_obs_', tmp.obsname_simple]);
            
            tmp.md_finite=tmp.md(isfinite(tmp.md));
            [f_exp, gof_exp] = fit(data.time_y_extended(isfinite(tmp.md))',squeeze(tmp.md_finite)','poly1');
            tmp.md_det_finite=f_exp(data.time_y_extended(isfinite(tmp.md)));
            tmp.md_det=tmp.md;
            tmp.md_det(isfinite(tmp.md_det))=tmp.md_det(isfinite(tmp.md_det))' - tmp.md_det_finite;
            
            tmp.obs_finite=tmp.obs(isfinite(tmp.obs));
            [f_exp, gof_exp] = fit(data.time_y_extended(isfinite(tmp.obs))',squeeze(tmp.obs_finite)','poly1');
            tmp.obs_det_finite=f_exp(data.time_y_extended(isfinite(tmp.obs)));
            tmp.obs_det=tmp.obs;
            tmp.obs_det(isfinite(tmp.obs_det))=tmp.obs_det(isfinite(tmp.obs_det))' - tmp.obs_det_finite;

            subplot(5,2, lyear+1)
            plot(data.time_y_extended(1:length(tmp.md_det)),tmp.md_det, '-o', 'color', 'g', 'linewidth', 2)
            xlim([min(data.time_y_extended), max(data.time_y_extended)])
            switch tmp.varname
                case 'temp'
                    ylim([-0.4 0.4])
                case 'salt'
                    ylim([-0.05 0.05])
            end
            hold on;
            plot(data.time_y, tmp.obs_det,  '-^', 'color', 'k', 'linewidth', 2)
            tmp.corrcoef=corrcoef(tmp.md_det(isfinite(tmp.md_det(1:length(data.time_y)))), tmp.obs_det(isfinite(tmp.md_det(1:length(data.time_y)))));
            data.(['gm_', tmp.varname, '_det_corr'])(lyear+1)=tmp.corrcoef(1,2);
            
%             xlabel('Year'); 
%             ylabel([tmp.varname, '-detrended']); 
            ylabel([tmp.varname, '-l', tmp.lyear_str]); 

            legend(['HCST-det-l', tmp.lyear_str, ', R=', num2str(round(tmp.corrcoef(1,2),2))],...
                [tmp.obsname_simple, '-det(obs-ba)'],...
                'Location', 'NorthWest');
%             set(gca, 'fontsize', 20)
            grid minor
            hold off
                    
        end
        config.figname=[dirs.figdir, filesep, 'gm_hcst_obs_', tmp.varname, '_det.tif'];
        set(gcf, 'PaperPosition', [0, 0, 10, 15]);
        saveas(gcf,config.figname,'tif');
        RemoveWhiteSpace([], 'file', config.figname);
        close all;

        %% correlation coefficient of SST, function of lead year
        config.figname=[dirs.figdir, filesep, 'gm_hcst_obs_', tmp.varname, '_det_corr', '.tif'];
%         plot(0:config.proj_year-1,data.(['gm_', tmp.varname, '_det_corr']),  '-o', 'linewidth', 2)
        bar(0:config.proj_year-1,data.(['gm_', tmp.varname, '_det_corr']),  'linewidth', 2)
        
        xlabel('Lead year'); ylabel(['corr. coef.,', tmp.varname, '-det, obs-ba']);
        set(gca, 'fontsize', 20)
        grid minor
        ylim([-1 1])
        saveas(gcf,config.figname,'tif');
        RemoveWhiteSpace([], 'file', config.figname);
        close all;

    end
    
     %% detrended time series and correlation coefficient (combined plot, hcst&assm)
     for varind=1:length(config.varnames)
        tmp.varname=config.varnames{varind};
        dirs.figdir= [dirs.figroot, filesep, config.casename_m, filesep, 'GMSV', filesep, 'Timeseries_det'];
        system(['mkdir -p ', dirs.figdir]);
        for lyear=0:config.proj_year-1
            tmp.lyear_str=num2str(lyear, '%02i');
            tmp.md=data.(['gm_', tmp.varname, '_model_', tmp.obsname_simple, '_l', tmp.lyear_str]);
            tmp.obs=data.(['gm_', tmp.varname, '_assm_', tmp.obsname_simple]);
            
            tmp.md_finite=tmp.md(isfinite(tmp.md));
            [f_exp, gof_exp] = fit(data.time_y_extended(isfinite(tmp.md))',squeeze(tmp.md_finite)','poly1');
            tmp.md_det_finite=f_exp(data.time_y_extended(isfinite(tmp.md)));
            tmp.md_det=tmp.md;
            tmp.md_det(isfinite(tmp.md_det))=tmp.md_det(isfinite(tmp.md_det))' - tmp.md_det_finite;
            
            tmp.obs_finite=tmp.obs(isfinite(tmp.obs));
            [f_exp, gof_exp] = fit(data.time_y_extended(isfinite(tmp.obs))',squeeze(tmp.obs_finite)','poly1');
            tmp.obs_det_finite=f_exp(data.time_y_extended(isfinite(tmp.obs)));
            tmp.obs_det=tmp.obs;
            tmp.obs_det(isfinite(tmp.obs_det))=tmp.obs_det(isfinite(tmp.obs_det))' - tmp.obs_det_finite;

            subplot(5,2, lyear+1)
            plot(data.time_y_extended(1:length(tmp.md_det)),tmp.md_det, '-o', 'color', 'g', 'linewidth', 2)
            xlim([min(data.time_y_extended), max(data.time_y_extended)])
            switch tmp.varname
                case 'temp'
                    ylim([-0.4 0.4])
                case 'salt'
                    ylim([-0.05 0.05])
            end
            hold on;
            plot(data.time_y, tmp.obs_det,  '-^', 'color', 'k', 'linewidth', 2)
            tmp.corrcoef=corrcoef(tmp.md_det(isfinite(tmp.md_det(1:length(data.time_y)))), tmp.obs_det(isfinite(tmp.md_det(1:length(data.time_y)))));
            data.(['gm_', tmp.varname, '_det_corr'])(lyear+1)=tmp.corrcoef(1,2);
            
%             xlabel('Year'); 
%             ylabel([tmp.varname, '-detrended']); 
            ylabel([tmp.varname, '-l', tmp.lyear_str]); 

            legend(['HCST-det-l', tmp.lyear_str, ', R=', num2str(round(tmp.corrcoef(1,2),2))],...
                [tmp.obsname_simple, '-det(assm)'],...
                'Location', 'NorthWest');
%             set(gca, 'fontsize', 20)
            grid minor
            hold off
                    
        end
        config.figname=[dirs.figdir, filesep, 'gm_hcst_assm_', tmp.varname, '_det.tif'];
        set(gcf, 'PaperPosition', [0, 0, 10, 15]);
        saveas(gcf,config.figname,'tif');
        RemoveWhiteSpace([], 'file', config.figname);
        close all;

        %% correlation coefficient of SST, function of lead year
        config.figname=[dirs.figdir, filesep, 'gm_hcst_assm_', tmp.varname, '_det_corr', '.tif'];
%         plot(0:config.proj_year-1,data.(['gm_', tmp.varname, '_det_corr']),  '-o', 'linewidth', 2)
        bar(0:config.proj_year-1,data.(['gm_', tmp.varname, '_det_corr']),  'linewidth', 2)
        
        xlabel('Lead year'); ylabel(['corr. coef.,', tmp.varname, '-det, assm']);
        set(gca, 'fontsize', 20)
        grid minor
        ylim([-1 1])
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

function gmsst = f_gm_var(var_2d, area)
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