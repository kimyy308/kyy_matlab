% %  Created 25-Oct-2022 by Yong-Yub Kim
% %  Created 09-Nov-2022 by Yong-Yub Kim

clc; clear all; close all;

%% set path
tmp.dropboxpath = '/Volumes/kyy_raid/kimyy/Dropbox';
tmp.fs=filesep;
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);

%% model configuration
dirs.root='/Volumes/kyy_raid/earth.system.predictability/ASSM_EXP';
dirs.yoshi_root='/proj/yoshi/DATA/CESM2_ODA';
dirs.archive=[dirs.root, filesep, 'archive'];
dirs.saveroot='/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/ASSM_EXP';

config.years=1999:2009;
config.months=1:12;
config.scenname='HIST';
config.gridname='f09_g17';
% config.obsnames={'oras4', 'projdv7.3', 'en4.2'};
config.obsnames={'projdv7.3', 'en4.2'};

% config.obsnames={'en4.2'};

% config.ensnames={'ba-10p1'};

% config.ensnames={'ba-10p1', 'ba-10p2', 'ba-10p3', 'ba-10p4', 'ba-10p5'};

config.ensnames={'ba-10p1', 'ba-10p2', 'ba-10p3', 'ba-10p4', 'ba-10p5', ...
    'ba-20p1', 'ba-20p2', 'ba-20p3', 'ba-20p4', 'ba-20p5'};

config.component='ocn';
config.varnames={'TEMP', 'SALT', 'DIC', 'ALK', 'NO3', 'PO4', 'SiO3'};
% config.varnames={'SALT', 'TEMP', 'DIC', 'PO4', 'ALK', 'NO3', 'SiO3', 'PD', 'NO2'};

% config.varnames={'PO4','NO3', 'SiO3', 'PD'};

config.len_t_y = length(config.years);
config.len_t_m = length(config.months);
config.len_obs= length(config.obsnames);
config.len_ens= length(config.ensnames);



% system(['ls ', dirs.datadir])
% b.e21.BHISTsmbb.f09_g17.assm.oras4_ba-10p1.pop.h.once.nc


%% observation configuration (K2, Text; .csv)
config_obs.staname = 'K2';
config_obs.sta_lon = 160;
config_obs.sta_lat = 47;
config_obs.avgdepth=[100, 150, 200];

dirs.obsroot = '/Volumes/kyy_raid/kimyy/Observation';
dirs.obssavedir = [dirs.obsroot, filesep, 'mat'];
mkdir(dirs.obssavedir)
config_obs.filename = [dirs.obsroot, filesep, 'K2/www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0100115/K2_bottle_data.csv'];
%    filename: /Volumes/kyy_raid/kimyy/Observation/K2/www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0100115/K2_bottle_data.csv

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 31);

% Specify range and delimiter
opts.DataLines = [3, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["EXPOCODE", "SECT_ID", "DATE", "TIME", "LATITUDE", "LONGITUDE", "BOT_DEPTH", "CTDPRS", "CTDTMP", "CTDSAL", "CTDSAL_FLAG_W", "THETA", "SIGMA_THETA", "SALNTY", "SALNTY_FLAG_W", "OXYGEN", "OXYGEN_FLAG_W", "NITRAT", "NITRAT_FLAG_W", "NITRIT", "NITRIT_FLAG_W", "SILCAT", "SILCAT_FLAG_W", "PHSPHT", "PHSPHT_FLAG_W", "AMMONIUM", "AMMONIUM_FLAG_W", "ALKALI", "ALKALI_FLAG_W", "TCARBN", "TCARBN_FLAG_W"];
opts.VariableTypes = ["categorical", "double", "datetime", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["EXPOCODE", "TIME"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "DATE", "InputFormat", "MM/dd/yyyy");
opts = setvaropts(opts, "SECT_ID", "TrimNonNumeric", true);
opts = setvaropts(opts, "SECT_ID", "ThousandsSeparator", ",");

% Import the data
OBS_bdat = readtable(config_obs.filename, opts);

opts.DataLines = [2, 2];
opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];

OBS_unit = readmatrix(config_obs.filename, opts);

%% Clear temporary variables
clear opts

OBS_datenum=datenum(table2array(OBS_bdat(:,3)));
[OBS_year, OBS_month, OBS_day] = datevec(OBS_datenum);

% OBS_bdat.Properties.VariableNames  %% variable names

tmp.NaNval=-999;
tmp.depth_ind=find(strcmp(OBS_bdat.Properties.VariableNames, 'CTDPRS')==1);
OBS.depth=f_standard_depth_obs(1);
OBS.d_layer=zeros(1,length(OBS.depth));
OBS.d_layer(1:end-1)=diff(OBS.depth)/2;
OBS.d_layer(2:end)=OBS.d_layer(2:end)+diff(OBS.depth)/2;
tmp.time_ind=find(strcmp(OBS_bdat.Properties.VariableNames, 'TIME')==1);
tmp.time=table2array(OBS_bdat(:,tmp.time_ind));

%% get time (str -> double, only hour), datenum+time
for timei=1:length(tmp.time)
    tmp.time_char=char(tmp.time(timei));
    if strcmp(tmp.time_char,'-999')==1
        OBS_time(timei)=0;
    else
        tmp.time_str_split=split(tmp.time_char, ':');
        OBS_time(timei)=str2num(tmp.time_str_split{1}) / 24.0;
    end
end
OBS_datenum=OBS_datenum+OBS_time';

%% get observation data
for varind=1:length(config.varnames)
    tmp.varname=config.varnames{varind}; %varname
    [tmp.varname_obs, tmp.varname_flag] = f_varname_obs(tmp.varname); %variable name of observation
    tmp.var_ind=find(strcmp(OBS_bdat.Properties.VariableNames, tmp.varname_obs)==1);
    tmp.data_len=size(OBS_bdat,1);
    tmp.datenum_uniq=unique(OBS_datenum);
    tmp.data_len_uniq=length(tmp.datenum_uniq);
  
    for ti=1:tmp.data_len_uniq
        tmp.data_ind=find(OBS_datenum==tmp.datenum_uniq(ti));
        tmp.data=table2array(OBS_bdat(tmp.data_ind, tmp.var_ind));
        tmp.depth=table2array(OBS_bdat(tmp.data_ind, tmp.depth_ind));
        
        tmp.invalid_ind=find(tmp.data==tmp.NaNval);
        tmp.data(tmp.invalid_ind)=NaN;
        tmp.valid_ind=isfinite(tmp.data);
        if sum(isnan(tmp.data))==length(tmp.data) % all data is nan
            tmp.data_interped(1:length(OBS.depth)) = NaN;
        else
            tmp.data_interped=interp1(tmp.depth(tmp.valid_ind), tmp.data(tmp.valid_ind), OBS.depth);
        end
        

        % put interpolated data
        OBS.(tmp.varname_obs)(:,ti)=tmp.data_interped;
        OBS.datenum(ti)=tmp.datenum_uniq(ti);
        % get depth averaged data
        for avgi=1:length(config_obs.avgdepth)
            tmp.varname_obs_avg=[tmp.varname_obs, '_', num2str(config_obs.avgdepth(avgi), '%04i')];
            [tmp.val,tmp.avgdep_ind]=find(OBS.depth==config_obs.avgdepth(avgi));
            tmp.data=OBS.(tmp.varname_obs)(1:tmp.avgdep_ind,ti);
            tmp.d_layer=OBS.d_layer(1:tmp.avgdep_ind)';
            tmp.d_layer(isnan(tmp.data))=NaN;
            tmp.allsum_d_layer=sum(OBS.d_layer(1:tmp.avgdep_ind));
            tmp.sum_d_layer=sum(tmp.d_layer,'omitnan');
            tmp.data_cover_rate=tmp.sum_d_layer/tmp.allsum_d_layer;
            if tmp.data_cover_rate>0.8
                OBS.(tmp.varname_obs_avg)(ti)= sum(tmp.data.*tmp.d_layer, 'omitnan') ./ tmp.sum_d_layer;
            else
                OBS.(tmp.varname_obs_avg)(ti)=NaN;
            end
        end
    end
    
%     plot(OBS.datenum+0.1,OBS.THETA_0100)
%     datetick('x','yymmm')


%% get model data
    for obsi=1:config.len_obs
        tmp.obsname_assm=config.obsnames{obsi};
        config.savename=[dirs.saveroot, filesep, config_obs.staname, filesep, 'CESM2_', config_obs.staname, '_', ...
            tmp.obsname_assm, '_', tmp.varname, '_', num2str(min(config.years)), '_', num2str(max(config.years)), '.mat'];
        config.savename_cfg=[dirs.saveroot, filesep, config_obs.staname, filesep, 'cfg_CESM2_', config_obs.staname, '_', ...
            tmp.obsname_assm, '_', tmp.varname, '_', num2str(min(config.years)), '_', num2str(max(config.years)), '.mat'];
        load(config.savename_cfg, 'CESM2_data', 'CESM2_grid');
        
        tmp.obs_abb=tmp.obsname_assm(1:3);
        CESM2_alldata.(tmp.obs_abb)=CESM2_data;
        CESM2_allgrid.(tmp.obs_abb)=CESM2_grid;

        % get depth averaged data
        for avgi=1:length(config_obs.avgdepth)
            for ensi=1:config.len_ens
                for tyi=1:config.len_t_y
                    for tmi=1:config.len_t_m
                        tmp.varname_model_avg=[tmp.varname, '_', num2str(config_obs.avgdepth(avgi), '%04i')]; % varname_save
                        tmp.val=find(CESM2_allgrid.(tmp.obs_abb).z_t/100.0 > config_obs.avgdepth(avgi)); % get depths > avg depth ind
                        tmp.avgdep_ind=min(tmp.val); % closest depth to avg depth ind
                        tmp.data=squeeze(CESM2_alldata.(tmp.obs_abb).(tmp.varname)(ensi,tyi,tmi,1:tmp.avgdep_ind)); % get temporary data
                        tmp.d_layer=CESM2_allgrid.(tmp.obs_abb).dz(1:tmp.avgdep_ind); % get dz
                        tmp.d_layer(isnan(tmp.data))=NaN;  % remove invalid layer
                        tmp.allsum_d_layer=sum(CESM2_allgrid.(tmp.obs_abb).dz(1:tmp.avgdep_ind)); % total sum of dz
                        tmp.sum_d_layer=sum(tmp.d_layer,'omitnan'); % total sum of valid dz
                        tmp.data_cover_rate=tmp.sum_d_layer/tmp.allsum_d_layer; % data cover rate
                        if tmp.data_cover_rate>0.8
                            CESM2_alldata.(tmp.obs_abb).(tmp.varname_model_avg)(ensi,(tyi-1)*12+tmi)= ...
                                sum(tmp.data.*tmp.d_layer, 'omitnan') ./ tmp.sum_d_layer;
                        else
                            CESM2_alldata.(tmp.obs_abb).(tmp.varname_model_avg)(ensi,(tyi-1)*12+tmi)=NaN;
                        end
                    end
                end
            end
        end
    end
    
    %% get subsequent datenum of model
    for tyi=1:config.len_t_y
        for tmi=1:config.len_t_m
            CESM2_grid.datenum((tyi-1)*12+tmi)=datenum(config.years(tyi), config.months(tmi), 15);
        end
    end


%     plot(OBS.datenum+0.1,OBS.([tmp.varname_obs,'_0100']), 'k*')
%     datetick('x','yymmm')
%     hold on
%     for obsi=1:config.len_obs
%         for ensi=1:config.len_ens
%             tmp.varname_model_avg=[tmp.varname, '_', '0100']; % varname_save
%             tmp.data=CESM2_alldata.(tmp.obs_abb).(tmp.varname_model_avg)(ensi,:);
%             plot(CESM2_grid.datenum, tmp.data, 'b')
%         end
%     end
%     hold off
    
    disp(tmp.varname)

    %% get ensemble mean, anomaly, normalization, mean, ensemble spread
    for avgi=1:length(config_obs.avgdepth)
        tmp.varname_obs_avg=[tmp.varname_obs, '_', num2str(config_obs.avgdepth(avgi), '%04i')];
        tmp.varname_ano_obs_avg=[tmp.varname_obs, '_ano_', num2str(config_obs.avgdepth(avgi), '%04i')];
        tmp.varname_norm_obs_avg=[tmp.varname_obs, '_norm_', num2str(config_obs.avgdepth(avgi), '%04i')];
        tmp.varname_m_obs_avg=[tmp.varname_obs, '_m_', num2str(config_obs.avgdepth(avgi), '%04i')];
        tmp.varname_std_obs_avg=[tmp.varname_obs, '_std_', num2str(config_obs.avgdepth(avgi), '%04i')];

        [OBS.(tmp.varname_m_obs_avg), OBS.(tmp.varname_std_obs_avg), ...
            OBS.(tmp.varname_ano_obs_avg), OBS.(tmp.varname_norm_obs_avg)] = ...
        get_ano_norm(OBS.(tmp.varname_obs_avg));
        
        tmp.varname_model_avg=[tmp.varname, '_', num2str(config_obs.avgdepth(avgi), '%04i')];
        tmp.varname_ano_model_avg=[tmp.varname, '_ano_', num2str(config_obs.avgdepth(avgi), '%04i')];
        tmp.varname_norm_model_avg=[tmp.varname, '_norm_', num2str(config_obs.avgdepth(avgi), '%04i')];
        tmp.varname_m_model_avg=[tmp.varname, '_m_', num2str(config_obs.avgdepth(avgi), '%04i')];
        tmp.varname_std_model_avg=[tmp.varname, '_std_', num2str(config_obs.avgdepth(avgi), '%04i')];
        for obsi=1:config.len_obs
            tmp.obsname_assm=config.obsnames{obsi};
            tmp.obs_abb=tmp.obsname_assm(1:3);
            tmp.data=CESM2_alldata.(tmp.obs_abb).(tmp.varname_model_avg);
            tmp.data_ensmean=mean(tmp.data,1);
            CESM2_alldata.([tmp.obs_abb,'_ensm']).(tmp.varname_model_avg)=tmp.data_ensmean;
            [CESM2_alldata.([tmp.obs_abb,'_ensm']).(tmp.varname_m_model_avg), ...
                CESM2_alldata.([tmp.obs_abb,'_ensm']).(tmp.varname_std_model_avg), ...
                CESM2_alldata.([tmp.obs_abb,'_ensm']).(tmp.varname_ano_model_avg), ...
                CESM2_alldata.([tmp.obs_abb,'_ensm']).(tmp.varname_norm_model_avg)] = ...
                get_ano_norm(CESM2_alldata.([tmp.obs_abb,'_ensm']).(tmp.varname_model_avg));
        end
    end

    
    %% normalized anomaly plot
    
    time_m=min(config.years):1/12:max(config.years)+11/12;
    tmp.plot_color(1,1:3)=[1,0,0];  % red, projd
    tmp.plot_color(2,1:3)=[0,1,0];  % green, en4 
    for avgi=1:length(config_obs.avgdepth)
        tmp.varname_obs_avg=[tmp.varname_obs, '_', num2str(config_obs.avgdepth(avgi), '%04i')];
        tmp.varname_ano_obs_avg=[tmp.varname_obs, '_ano_', num2str(config_obs.avgdepth(avgi), '%04i')];
        tmp.varname_norm_obs_avg=[tmp.varname_obs, '_norm_', num2str(config_obs.avgdepth(avgi), '%04i')];
        tmp.varname_m_obs_avg=[tmp.varname_obs, '_m_', num2str(config_obs.avgdepth(avgi), '%04i')];
        tmp.varname_std_obs_avg=[tmp.varname_obs, '_std_', num2str(config_obs.avgdepth(avgi), '%04i')];

        tmp.varname_model_avg=[tmp.varname, '_', num2str(config_obs.avgdepth(avgi), '%04i')];
        tmp.varname_ano_model_avg=[tmp.varname, '_ano_', num2str(config_obs.avgdepth(avgi), '%04i')];
        tmp.varname_norm_model_avg=[tmp.varname, '_norm_', num2str(config_obs.avgdepth(avgi), '%04i')];
        tmp.varname_m_model_avg=[tmp.varname, '_m_', num2str(config_obs.avgdepth(avgi), '%04i')];
        tmp.varname_std_model_avg=[tmp.varname, '_std_', num2str(config_obs.avgdepth(avgi), '%04i')];

        hold on
        for obsi=1:config.len_obs
            %% ensemble mean & observation plot
            tmp.obsname_assm=config.obsnames{obsi};
            tmp.obs_abb=tmp.obsname_assm(1:3);
            tmp.data=CESM2_alldata.([tmp.obs_abb,'_ensm']).(tmp.varname_norm_model_avg);
            plot(CESM2_grid.datenum,tmp.data, 'color', tmp.plot_color(obsi,:), 'linewidth', 3);

            %% ensemble member plot (10p)
            tmp.val_transparent = 0.5;
            tmp.plot_color_tr_10p(obsi,:)=get_transparent_rgb(tmp.plot_color(obsi,:), tmp.val_transparent);
            tmp.size_ens=size(CESM2_alldata.(tmp.obs_abb).(tmp.varname_model_avg));
            for ensi=1:tmp.size_ens/2
                [CESM2_alldata.([tmp.obs_abb]).(tmp.varname_m_model_avg)(ensi), ...
                CESM2_alldata.([tmp.obs_abb]).(tmp.varname_std_model_avg)(ensi), ...
                CESM2_alldata.([tmp.obs_abb]).(tmp.varname_ano_model_avg)(ensi,:), ...
                CESM2_alldata.([tmp.obs_abb]).(tmp.varname_norm_model_avg)(ensi,:)] = ...
                    get_ano_norm(CESM2_alldata.([tmp.obs_abb]).(tmp.varname_model_avg)(ensi,:));
                tmp.data=CESM2_alldata.([tmp.obs_abb]).(tmp.varname_norm_model_avg)(ensi,:); %data
                tsplot{ensi}=plot(CESM2_grid.datenum, tmp.data, 'color', tmp.plot_color_tr_10p(obsi,:));
                set(get(get(tsplot{ensi},'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % turn off legend informatino
            end

            %% ensemble member plot (20p)
            tmp.val_transparent = 0.2;
            tmp.plot_color_tr_20p(obsi,:)=get_transparent_rgb(tmp.plot_color(obsi,:), tmp.val_transparent);
            tmp.size_ens=size(CESM2_alldata.(tmp.obs_abb).(tmp.varname_model_avg));
            for ensi=tmp.size_ens/2+1:tmp.size_ens
                [CESM2_alldata.([tmp.obs_abb]).(tmp.varname_m_model_avg)(ensi), ...
                CESM2_alldata.([tmp.obs_abb]).(tmp.varname_std_model_avg)(ensi), ...
                CESM2_alldata.([tmp.obs_abb]).(tmp.varname_ano_model_avg)(ensi,:), ...
                CESM2_alldata.([tmp.obs_abb]).(tmp.varname_norm_model_avg)(ensi,:)] = ...
                    get_ano_norm(CESM2_alldata.([tmp.obs_abb]).(tmp.varname_model_avg)(ensi,:));
                tmp.data=CESM2_alldata.([tmp.obs_abb]).(tmp.varname_norm_model_avg)(ensi,:); %data
                tsplot{ensi}=plot(CESM2_grid.datenum, tmp.data, 'color', tmp.plot_color_tr_20p(obsi,:));
                set(get(get(tsplot{ensi},'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % turn off legend informatino
            end

            datetick('x', 'yymmm')
        end
        plot(OBS.datenum,OBS.(tmp.varname_norm_obs_avg), 'k*');
        xlabel('Year'); ylabel([tmp.varname, '(', tmp.varname_obs, ')']);
%         lgd=legend([config.obsnames{1}], [config.obsnames{2}], ['OBS(', config_obs.staname,')'], 'Location', 'NorthWest');
        set(gca, 'fontsize', 20)
%         set(lgd, 'fontsize',5)
        grid minor
        hold off
        
        dirs.figdir=['/Volumes/kyy_raid/kimyy/Figure/CESM2/ESP/ASSM_EXP/', ...
            config_obs.staname, '/', num2str(config_obs.avgdepth(avgi), '%04i'), 'm'];
        system(['mkdir -p ', dirs.figdir])
        tmp.figname=[dirs.figdir, filesep, tmp.varname, '_ano_norm_', tmp.varname_obs, '_', ...
            num2str(config_obs.avgdepth(avgi), '%04i'), 'm_avg_', ...
            num2str(min(config.years)), '_', num2str(max(config.years)), '.tif'];
        set(gcf, 'PaperPosition', [0, 0, 8, 4]);
        saveas(gcf,tmp.figname,'tif');
        RemoveWhiteSpace([], 'file', tmp.figname);
        close all;
        %% corrcoef
%         ind_valid=find(isfinite(obs_data_m_ano_norm));
%         tmp.coef=corrcoef(en4_ensm_ano_norm(ind_valid), obs_data_m_ano_norm(ind_valid));
%         coef_en4=tmp.coef(1,2);
%         tmp.coef=corrcoef(projd_ensm_ano_norm(ind_valid), obs_data_m_ano_norm(ind_valid));
%         coef_projd=tmp.coef(1,2);

    end

    
%% raw plot
    time_m=min(config.years):1/12:max(config.years)+11/12;
    tmp.plot_color(1,1:3)=[1,0,0];  % red, projd
    tmp.plot_color(2,1:3)=[0,1,0];  % green, en4 
    for avgi=1:length(config_obs.avgdepth)
        tmp.varname_obs_avg=[tmp.varname_obs, '_', num2str(config_obs.avgdepth(avgi), '%04i')];
        tmp.varname_model_avg=[tmp.varname, '_', num2str(config_obs.avgdepth(avgi), '%04i')];
        
        hold on
        for obsi=1:config.len_obs
            %% ensemble mean & observation plot
            tmp.obsname_assm=config.obsnames{obsi};
            tmp.obs_abb=tmp.obsname_assm(1:3);
            tmp.data=CESM2_alldata.([tmp.obs_abb,'_ensm']).(tmp.varname_model_avg); %data
            plot(CESM2_grid.datenum,tmp.data, 'color', tmp.plot_color(obsi,:), 'linewidth', 3);

            %% ensemble member plot (10p)
            tmp.val_transparent = 0.5;
            tmp.plot_color_tr_10p(obsi,:)=get_transparent_rgb(tmp.plot_color(obsi,:), tmp.val_transparent);
            tmp.size_ens=size(CESM2_alldata.(tmp.obs_abb).(tmp.varname_model_avg));
            for ensi=1:tmp.size_ens/2
                tmp.data=CESM2_alldata.([tmp.obs_abb]).(tmp.varname_model_avg)(ensi,:); %data
                tsplot{ensi}=plot(CESM2_grid.datenum, tmp.data, 'color', tmp.plot_color_tr_10p(obsi,:));
                set(get(get(tsplot{ensi},'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % turn off legend informatino
            end

            %% ensemble member plot (20p)
            tmp.val_transparent = 0.2;
            tmp.plot_color_tr_20p(obsi,:)=get_transparent_rgb(tmp.plot_color(obsi,:), tmp.val_transparent);
            tmp.size_ens=size(CESM2_alldata.(tmp.obs_abb).(tmp.varname_model_avg));
            for ensi=tmp.size_ens/2+1:tmp.size_ens
                tmp.data=CESM2_alldata.([tmp.obs_abb]).(tmp.varname_model_avg)(ensi,:); %data
                tsplot{ensi}=plot(CESM2_grid.datenum, tmp.data, 'color', tmp.plot_color_tr_20p(obsi,:));
                set(get(get(tsplot{ensi},'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % turn off legend informatino
            end

            datetick('x', 'yymmm')
        end
        plot(OBS.datenum,OBS.(tmp.varname_obs_avg), 'k*');
        xlabel('Year'); ylabel([tmp.varname, '(', tmp.varname_obs, ')']);
%         lgd=legend([config.obsnames{1}], [config.obsnames{2}], ['OBS(', config_obs.staname,')'], 'Location', 'northoutside');
        set(gca, 'fontsize', 20)
%         set(lgd, 'fontsize',5)
        grid minor
        hold off
        
        dirs.figdir=['/Volumes/kyy_raid/kimyy/Figure/CESM2/ESP/ASSM_EXP/', ...
            config_obs.staname, '/', num2str(config_obs.avgdepth(avgi), '%04i'), 'm'];
        system(['mkdir -p ', dirs.figdir])
        tmp.figname=[dirs.figdir, filesep, tmp.varname, '_raw_', tmp.varname_obs, '_', ...
            num2str(config_obs.avgdepth(avgi), '%04i'), 'm_avg_', ...
            num2str(min(config.years)), '_', num2str(max(config.years)), '.tif'];
        set(gcf, 'PaperPosition', [0, 0, 8, 4]);
        saveas(gcf,tmp.figname,'tif');
        RemoveWhiteSpace([], 'file', tmp.figname);
        close all;
    end

end


function [var_obs, var_flag] = f_varname_obs(var_model)
% metadata
% https://www.ncei.noaa.gov/data/oceans/ncei/ocads/metadata/0100115.html
var_flag=NaN;
var_obs=NaN;
    switch var_model
        case 'lon'
            var_obs='LONGITUDE';
        case 'lat'
            var_obs='LATITUDE';
        case 'time'
            var_obs='TIME';
        case 'depth'
            var_obs='CTDPRS';
        case 'SALT'
            var_obs='SALNTY'; 
            var_flag='SALNTY_FLAG_W';
        case 'TEMP'
            var_obs='THETA'; % Potential temperature
%             var_obs='CTDTMP'; % Water temperature           
        case 'DIC' % mmol/m^3
            var_obs='TCARBN'; % DISSOLVED INORGANIC CARBON, umol/kg
            var_flag='TCARBN_FLAG_W';
        case 'PO4' %mmol/m^3
            var_obs='PHSPHT';  % phosphate, umol/kg
            var_flag='PHSPHT_FLAG_W';
        case 'ALK' % meq/m^3
            var_obs='ALKALI'; % total alkalinity, umol/kg
            var_flag='ALKALI_FLAG_W';
        case 'NO3' %mmol/m^3
            var_obs='NITRAT'; % NITRATE
            var_flag='NITRAT_FLAG_W';
        case 'SiO3' %mmol/m^3
            var_obs='SILCAT';  %silicate, umol/kg
            var_flag='SILCAT_FLAG_W';            
        case 'NO2'
            var_obs='NITRIT'; %NITRITE umol /kg
            var_flag='NITRIT_FLAG_W';            
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

function std_dep=f_standard_depth_obs(flag)
    switch flag
        case 1
            std_dep=[0:10:200, 225:25:300, 350:50:500, 600:100:1500, 1750:250:6000];
    end
end