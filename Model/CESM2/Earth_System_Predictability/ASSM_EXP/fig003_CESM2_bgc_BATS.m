% %  Created 25-Oct-2022 by Yong-Yub Kim
% %  Created 05-Nov-2022 by Yong-Yub Kim
% %  Created 16-Nov-2022 by Yong-Yub Kim

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

config.years=1988:2021;
config.months=1:12;
config.scenname='HIST';
config.gridname='f09_g17';
config.obsnames={'projdv7.3', 'en4.2'};
config.ensnames={'ba-10p1', 'ba-10p2', 'ba-10p3', 'ba-10p4', 'ba-10p5', ...
    'ba-20p1', 'ba-20p2', 'ba-20p3', 'ba-20p4', 'ba-20p5'};

config.component='ocn';
config.varnames={'TEMP', 'SALT', 'DIC', 'ALK', 'NO3', 'PO4', 'SiO3'};
% config.varnames={'SALT'};

config.len_t_y = length(config.years);
config.len_t_m = length(config.months);
config.len_obs= length(config.obsnames);
config.len_ens= length(config.ensnames);


%% observation configuration (BATS, Excel; .xlsx)

config_obs.staname = 'BATS';
% config_obs.sta_lon = [360-66.1690, 360-60.4470]; % degrees west to 0~360 degrees east
% config_obs.sta_lat = [24.7590, 35.6670]; % degrees north
config_obs.sta_lon = 360-64.1725; % degrees west to 0~360 degrees east
config_obs.sta_lat = 31.6840; % degrees north

config_obs.avgdepth=[100, 140, 200];

dirs.obsroot = '/Volumes/kyy_raid/kimyy/Observation';
dirs.obssavedir = [dirs.obsroot, filesep, 'mat'];
mkdir(dirs.obssavedir)
config_obs.filename = [dirs.obsroot, filesep, 'BATS/batsftp.bios.edu/BATS/bottle/bats_bottle.xlsx'];



%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 34);

% Specify sheet and range
opts.Sheet = "bats_bottle";
opts.DataRange = "A61:AH64841";

% Specify column names and types
opts.VariableNames = ["ID", "yyyymmdd", "decy", "time", "latN", "lonW", "Depth", "Temp", "CTD_S", "salt", "Sig_th", "O2", "OxFix", "Anom1", "CO2", "alk", "NO3NO2", "NO2", "PO4", "SiO2", "POC", "PON", "TOC", "TN", "Bact", "POP", "TDP", "SRP", "Bsi", "Lsi", "Pro", "Syn", "Piceu", "Naneu"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
OBS_bdat = readtable(config_obs.filename, opts, "UseExcel", false);

%% Clear temporary variables
clear opts

%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 1);

% Specify sheet and range
opts.Sheet = "bats_bottle";
opts.DataRange = "A23:A55";

% Specify column names and types
opts.VariableNames = "ID";
opts.VariableTypes = "double";

% Import the data
OBS_unit = readtable(config_obs.filename, opts, "UseExcel", false);

%% Clear temporary variables
clear opts

% OBS_config.lat_min=min(table2array(BATS_bdat(:,5)));
% OBS_config.lat_max=max(table2array(BATS_bdat(:,5)));
% OBS_config.lon_min=min(table2array(BATS_bdat(:,6)));
% OBS_config.lon_max=max(table2array(BATS_bdat(:,6)));


% K2_datenum=datenum(table2array(K2_bdat(:,3)));
% [K2_year, K2_month, K2_day] = datevec(K2_datenum);

tmp.NaNval=-999;
OBS.depth_ind=find(strcmp(OBS_bdat.Properties.VariableNames, 'Depth')==1);
OBS.depth_interp=f_standard_depth_obs(1);
OBS.d_layer=zeros(1,length(OBS.depth_interp));
OBS.d_layer(1:end-1)=diff(OBS.depth_interp)/2;
OBS.d_layer(2:end)=OBS.d_layer(2:end)+diff(OBS.depth_interp)/2;
tmp.yymmdd_ind=find(strcmp(OBS_bdat.Properties.VariableNames, 'yyyymmdd')==1);
tmp.yymmdd=table2array(OBS_bdat(:,tmp.yymmdd_ind));
tmp.hour_ind=find(strcmp(OBS_bdat.Properties.VariableNames, 'time')==1);
tmp.hour=table2array(OBS_bdat(:,tmp.hour_ind));
% tmp.time=tmp.yymmdd+tmp.hour;
OBS_datenum=datenum(tmp.yymmdd/10000.0, mod(tmp.yymmdd,10000)/100, mod(tmp.yymmdd, 100), ...
    floor(tmp.hour/100.0), mod(tmp.hour,100), 0);
tmp.decy_ind=find(strcmp(OBS_bdat.Properties.VariableNames, 'decy')==1);
tmp.decy=table2array(OBS_bdat(:,tmp.decy_ind));

OBS.lat_ind=find(strcmp(OBS_bdat.Properties.VariableNames, 'latN')==1);
OBS.lat=table2array(OBS_bdat(:,OBS.lat_ind));
OBS.latmean=mean(OBS.lat);
OBS.latstd=std(OBS.lat);
OBS.lon_ind=find(strcmp(OBS_bdat.Properties.VariableNames, 'lonW')==1);
OBS.lon=table2array(OBS_bdat(:,OBS.lon_ind));
OBS.lonmean=(mean(OBS.lon));
OBS.lonstd=(std(OBS.lon));

OBS_datenum=OBS_datenum+OBS.lon/1000.0;
OBS_datenum=OBS_datenum+OBS.lat/10000.0;




% find(tmp.lat==min(tmp.lat))
% histogram(tmp.lon)
abc=1;

%% get observation data
for varind=1:length(config.varnames)
    tmp.varname=config.varnames{varind}; %varname
    [tmp.varname_obs, tmp.varname_flag] = f_varname_obs(tmp.varname); %variable name of observation
    OBS.var_ind=find(strcmp(OBS_bdat.Properties.VariableNames, tmp.varname_obs)==1);
    OBS.data_len=size(OBS_bdat,1);
    OBS.datenum_uniq=unique(OBS_datenum);
    OBS.data_len_uniq=length(OBS.datenum_uniq);
    
    for ti=1:OBS.data_len_uniq
        OBS.data_ind=find(OBS_datenum==OBS.datenum_uniq(ti));
        tmp.data=table2array(OBS_bdat(OBS.data_ind, OBS.var_ind));
        tmp.depth=table2array(OBS_bdat(OBS.data_ind, OBS.depth_ind));
        OBS.lon=table2array(OBS_bdat(OBS.data_ind, OBS.lon_ind));
        OBS.lat=table2array(OBS_bdat(OBS.data_ind, OBS.lat_ind));

        
        tmp.invalid_ind=find(tmp.data==tmp.NaNval);
        tmp.data(tmp.invalid_ind)=NaN;
        tmp.valid_ind=isfinite(tmp.data);
        if (sum(isnan(tmp.data))==length(tmp.data) || ... % all data is nan
                unique(OBS.lon)<OBS.lonmean-3*OBS.lonstd || ... % Too far from center station
                unique(OBS.lon)>OBS.lonmean+3*OBS.lonstd || ...
                unique(OBS.lat)<OBS.latmean-3*OBS.latstd || ...
                unique(OBS.lat)>OBS.latmean+3*OBS.latstd)
            tmp.data_interped(1:length(OBS.depth_interp)) = NaN;
        else
            tmp.depth_valid=tmp.depth(tmp.valid_ind);
            tmp.data_valid=tmp.data(tmp.valid_ind);
            [tmp.uniqval, tmp.uniqind] = unique(tmp.depth_valid);
            if length(tmp.depth_valid(tmp.uniqind))==1
                tmp.data_interped(1:length(OBS.depth_interp)) = NaN;
            else
                tmp.data_interped=interp1(tmp.depth_valid(tmp.uniqind), tmp.data_valid(tmp.uniqind), OBS.depth_interp);
            end
        end
        

        % put interpolated data
        OBS.(tmp.varname_obs)(:,ti)=tmp.data_interped;
        OBS.datenum(ti)=OBS.datenum_uniq(ti);
        % get depth averaged data
        for avgi=1:length(config_obs.avgdepth)
            tmp.varname_obs_avg=[tmp.varname_obs, '_', num2str(config_obs.avgdepth(avgi), '%04i')];
            [tmp.val,tmp.avgdep_ind]=find(OBS.depth_interp==config_obs.avgdepth(avgi));
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
%     plot(OBS.datenum,OBS.Temp_0140, 'k-*')
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
        CESM2_grid.lon_cut=CESM2_grid.lon_t(CESM2_grid.obs_lon_ind_w:CESM2_grid.obs_lon_ind_e, CESM2_grid.obs_lat_ind_s:CESM2_grid.obs_lat_ind_n);
        CESM2_grid.lat_cut=CESM2_grid.lat_t(CESM2_grid.obs_lon_ind_w:CESM2_grid.obs_lon_ind_e, CESM2_grid.obs_lat_ind_s:CESM2_grid.obs_lat_ind_n);

        CESM2_alldata.(tmp.obs_abb)=CESM2_data;  % (ens, x, y, t_y, t_m, z)
        CESM2_allgrid.(tmp.obs_abb)=CESM2_grid;
       
        [tmp.lon_ind_w, tmp.lon_ind_e, tmp.lat_ind_s, tmp.lat_ind_n] = ...
                        Func_0012_findind_Y(0.5, [360-OBS.lonmean, OBS.latmean], ...
                        CESM2_grid.lon_cut, ...
                        CESM2_grid.lat_cut, 'CESM2'); % find valid lon, lat index near station
        tmp.lon_ind=tmp.lon_ind_w;
        tmp.lat_ind=tmp.lat_ind_s;
        % get depth averaged data
        for avgi=1:length(config_obs.avgdepth)
            for ensi=1:config.len_ens
                for tyi=1:config.len_t_y
                    for tmi=1:config.len_t_m
                        tmp.varname_model_avg=[tmp.varname, '_', num2str(config_obs.avgdepth(avgi), '%04i')]; % varname_save
                        tmp.val=find(CESM2_allgrid.(tmp.obs_abb).z_t/100.0 > config_obs.avgdepth(avgi)); % get depths > avg depth ind
                        tmp.avgdep_ind=min(tmp.val); % closest depth to avg depth ind
                        tmp.data=squeeze(CESM2_alldata.(tmp.obs_abb).(tmp.varname)(ensi,tmp.lon_ind,tmp.lat_ind,tyi,tmi,1:tmp.avgdep_ind)); % get temporary data
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

    plot(OBS.datenum+0.1,OBS.([tmp.varname_obs,'_0100']), 'k*')
    datetick('x','yymmm')
    hold on
    for obsi=1:config.len_obs
        for ensi=1:config.len_ens
            tmp.varname_model_avg=[tmp.varname, '_', '0100']; % varname_save
            tmp.data=CESM2_alldata.(tmp.obs_abb).(tmp.varname_model_avg)(ensi,:);
            plot(CESM2_grid.datenum, tmp.data, 'b')
        end
    end
    hold off
    

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



%loop

%CESM2 data get

%rearrange CESM2 data according to BATS observation

%~140m integration

% time series (NO2, PO4, SiO2)


%loop










function [var_obs, var_flag] = f_varname_obs(var_model)
% metadata
% http://batsftp.bios.edu/BATS/bottle/bats_bottle.txt

% yyyymmdd = Year Month Day   
% decy   = Decimal Year     
% time   = Time (hhmm)      
% latN   = Latitude (Deg N) 
% lonW   = Longitude (Deg W)
% Depth  = Depth (m)                  
% Temp   = Temperature ITS-90 (C)    
% CTD_S  = CTD Salinity (PSS-78)      
% Sal1   = Salinity-1 (PSS-78)        
% Sig-th = Sigma-Theta (kg/m^3)       
% O2(1)  = Oxygen-1 (umol/kg)          
% OxFixT = Oxygen Fix Temp (C)        
% Anom1  = Oxy Anomaly-1 (umol/kg)    
% CO2    = dissolved inorganic carbon (umol/kg)              
% Alk    = Alkalinity (uequiv)        
% NO3    = Nitrate+Nitrite-1 (umol/kg)
% NO2    = Nitrite-1 (umol/kg)        
% PO4    = Phosphate-1 (umol/kg)      
% Si     = Silicate-1 (umol/kg)       
% POC    = Particulate Organic Carbon; POC (ug/kg)                
% PON    = PON (ug/kg)                
% TOC    = Total Organic Carbon; TOC (umol/kg)                
% TN     = TN (umol/kg)  NOTE: Prior to BATS 121, DON is reported instead of TON
% Bact   = Bacteria enumeration (cells*10^8/kg)
% POP    = POP (umol/kg)
% TDP    = Total dissolved Phosphorus (nmol/kg)
% SRP    = Low-level phosphorus (nmol/kg)
% BSi    = Particulate biogenic silica (umol/kg)
% LSi    = Particulate lithogenic silica  (umol/kg)
% Pro    = Prochlorococcus (cells/ml)
% Syn    = Synechococcus (cells/ml)
% Piceu  = Picoeukaryotes (cells/ml)
% Naneu  = Nanoeukaryotes (cells/ml)   
% /Quality flags
% -999 = Bad or missing data
%  0 = Less than detection limit


var_flag=NaN;
var_obs=NaN;
    switch var_model
        case 'lon'
            var_obs='lonW';
        case 'lat'
            var_obs='latN';
        case 'time'
            var_obs='TIME';
        case 'depth'
            var_obs='Depth';
        case 'SALT'
            var_obs='salt'; 
%             var_flag='SALNTY_FLAG_W';
        case 'TEMP'
%             var_obs='Temp'; % Potential temperature
            var_obs='Temp'; % Water temperature           
        case 'DIC' % mmol/m^3
            var_obs='CO2'; % DISSOLVED INORGANIC CARBON, umol/kg
%             var_flag='TCARBN_FLAG_W';
        case 'PO4' %mmol/m^3
            var_obs='PO4';  % phosphate, umol/kg
%             var_flag='PHSPHT_FLAG_W';
        case 'ALK' % meq/m^3
            var_obs='alk'; % total alkalinity, umol/kg
%             var_flag='ALKALI_FLAG_W';
        case 'NO3' %mmol/m^3
            var_obs='NO3NO2'; % NITRATE
%             var_flag='NITRAT_FLAG_W';
        case 'SiO3' %mmol/m^3
            var_obs='SiO2';  %silicate, umol/kg
%             var_flag='SILCAT_FLAG_W';            
        case 'NO2'
            var_obs='NO2'; %NITRITE umol /kg
%             var_flag='NITRIT_FLAG_W';            
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