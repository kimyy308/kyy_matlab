% %  Created 30-Nov-2022 by Yong-Yub Kim  

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

% config.years=1967:2021;
config.years=1960:2021;

config.months=1:12;
config.scenname='HIST';
config.gridname='f09_g17';
config.obsnames={'projdv7.3', 'en4.2'};
config.ensnames={'ba-10p1', 'ba-10p2', 'ba-10p3', 'ba-10p4', 'ba-10p5', ...
    'ba-20p1', 'ba-20p2', 'ba-20p3', 'ba-20p4', 'ba-20p5'};

config.component='ocn';
config.varnames={'TEMP', 'SALT', 'DIC', 'ALK', 'NO3', 'PO4', 'SiO3'};
% config.varnames={'NO3'};

config.len_t_y = length(config.years);
config.len_t_m = length(config.months);
config.len_t = config.len_t_y * config.len_t_m;
config.len_obs= length(config.obsnames);
config.len_ens= length(config.ensnames);
% 
% 
% 
%% observation configuration (BATS, Excel; .xlsx)

config_obs.staname = 'CalCOFI';
% config_obs.sta_lon = [360-66.1690, 360-60.4470]; % degrees west to 0~360 degrees east
% config_obs.sta_lat = [24.7590, 35.6670]; % degrees north
% config_obs.sta_lon = 360-64.1725; % degrees west to 0~360 degrees east
% config_obs.sta_lat = 31.6840; % degrees north

config_obs.avgdepth=[10, 100, 150];

dirs.obsroot = '/Volumes/kyy_raid/kimyy/Observation/CalCOFI';
dirs.obssavedir = [dirs.obsroot, filesep, 'mat'];
mkdir(dirs.obssavedir)
OBS.tmp=0;


%% file read & save
% % % %% Import data from CalCOFI csv file
% % % %% Set up the Import Options and import the data

tmp.tmpsavename_1st=[dirs.obssavedir, '/obs.mat'];
if (exist(tmp.tmpsavename_1st) ==0)
    opts = delimitedTextImportOptions("NumVariables", 62);
    % Specify range and delimiter
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";
    % Specify column names and types
    opts.VariableNames = ["Cst_Cnt", "Btl_Cnt", "Sta_ID", "Depth_ID", "Depthm", ...
        "T_degC", "Salnty", "O2ml_L", "STheta", "O2Sat", "Oxy_molKg", "BtlNum", "RecInd", ...
        "T_prec", "T_qual", "S_prec", "S_qual", "P_qual", "O_qual", "SThtaq", "O2Satq", ...
        "ChlorA", "Chlqua", "Phaeop", "Phaqua", "PO4uM", "PO4q", "SiO3uM", "SiO3qu", "NO2uM", ...
        "NO2q", "NO3uM", "NO3q", "NH3uM", "NH3q", "C14As1", "C14A1p", "C14A1q", "C14As2", ...
        "C14A2p", "C14A2q", "DarkAs", "DarkAp", "DarkAq", "MeanAs", "MeanAp", "MeanAq", ...
        "IncTim", "LightP", "R_Depth", "R_TEMP", "R_Sal", "R_DYNHT", "R_Nuts", ...
        "R_Oxy_molKg", "DIC1", "DIC2", "TA1", "TA2", "pH1", "pH2", "DICQualityComment"];
    opts.VariableTypes = ["double", "double", "categorical", "string", "double", ...
        "double", "double", "double", "double", "double", "string", "string", "double", ...
        "double", "double", "double", "double", "double", "double", "double", "double", ...
        "string", "double", "string", "double", "double", "double", "double", "double", ...
        "double", "double", "double", "double", "double", "double", "string", "string", ...
        "double", "string", "string", "double", "string", "string", "double", "string", ...
        "string", "double", "string", "string", "double", "double", "double", "double", ...
        "string", "string", "double", "double", "double", "double", "double", "double", "string"];
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    % Specify variable properties
    opts = setvaropts(opts, ["Depth_ID", "Oxy_molKg", "BtlNum", "ChlorA", "Phaeop", "C14As1", "C14A1p", "C14As2", "C14A2p", "DarkAs", "DarkAp", "MeanAs", "MeanAp", "IncTim", "LightP", "R_Nuts", "R_Oxy_molKg", "DICQualityComment"], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["Sta_ID", "Depth_ID", "Oxy_molKg", "BtlNum", "ChlorA", "Phaeop", "C14As1", "C14A1p", "C14As2", "C14A2p", "DarkAs", "DarkAp", "MeanAs", "MeanAp", "IncTim", "LightP", "R_Nuts", "R_Oxy_molKg", "DICQualityComment"], "EmptyFieldRule", "auto");
    opts = setvaropts(opts, ["O2ml_L", "O2Sat", "PO4uM", "SiO3uM", "NO2uM", "NO3uM", "NH3uM", "DIC1", "DIC2", "TA1", "TA2", "pH1", "pH2"], "TrimNonNumeric", true);
    opts = setvaropts(opts, ["O2ml_L", "O2Sat", "PO4uM", "SiO3uM", "NO2uM", "NO3uM", "NH3uM", "DIC1", "DIC2", "TA1", "TA2", "pH1", "pH2"], "ThousandsSeparator", ",");
    % Import the data
    OBS_bdat = readtable([dirs.obsroot,'/CalCOFI_Database_194903-202001_csv_22Sep2021/194903-202001_Bottle.csv'], opts);
    %% Clear temporary variables
    clear opts
    
    tmp.NaNval=-999;
    %% interpolated depth configuration
    OBS.depth_ind=find(strcmp(OBS_bdat.Properties.VariableNames, 'Depthm')==1);
    OBS.depth_interp=f_standard_depth_obs(1);
    OBS.d_layer=zeros(1,length(OBS.depth_interp));
    OBS.d_layer(1:end-1)=diff(OBS.depth_interp)/2;
    OBS.d_layer(2:end)=OBS.d_layer(2:end)+diff(OBS.depth_interp)/2;
    
    OBS.Sta_ID_ind=find(strcmp(OBS_bdat.Properties.VariableNames, 'Sta_ID')==1);
    OBS.Sta_ID=table2array(OBS_bdat(:,OBS.Sta_ID_ind));
    OBS.Depth_ID_ind=find(strcmp(OBS_bdat.Properties.VariableNames, 'Depth_ID')==1);
    OBS.Depth_ID=table2array(OBS_bdat(:,OBS.Depth_ID_ind));
    tmp.data=char(OBS.Depth_ID);
    OBS.Depth_ID_char=string(tmp.data(:,1:21));
    OBS.len_data = length(OBS.Sta_ID);
    
    for di=1:OBS.len_data  % it takes 111 sec
    %     for di=1:1
        tmp.sta=split(char(OBS.Sta_ID(di)));
        OBS.line(di)=str2num(tmp.sta{1});
        OBS.sta(di)=str2num(tmp.sta{2});
        tmp.depid=split(OBS.Depth_ID(di), '-');
        tmp.yms=char(tmp.depid(2));
        OBS.year(di)=str2num([char(tmp.depid(1)), tmp.yms(1:2)]);
        OBS.month(di)=str2num(tmp.yms(3:4));
        [OBS.lat(di), OBS.lon(di)] = cc2lat(OBS.line(di),OBS.line(di));
    end
    save([dirs.obssavedir, '/obs.mat'], 'OBS', 'OBS_bdat')
else
    load([dirs.obssavedir, '/obs.mat'], 'OBS', 'OBS_bdat')
    %% load from saved file
    OBS.depth_ind=find(strcmp(OBS_bdat.Properties.VariableNames, 'Depthm')==1);
    tmp.data=char(OBS.Depth_ID);
    OBS.Depth_ID_char=string(tmp.data(:,1:21));
    tmp.NaNval=-999;
end

%% get observation data
for varind=1:length(config.varnames)
    tmp.varname=config.varnames{varind}; %varname
    [tmp.varname_obs, tmp.varname_flag] = f_varname_obs(tmp.varname); %variable name of observation
    OBS.var_ind=find(strcmp(OBS_bdat.Properties.VariableNames, tmp.varname_obs)==1);
    OBS.Depth_ID_uniq=unique(OBS.Depth_ID_char);
    OBS.len_Depth_ID_uniq=length(OBS.Depth_ID_uniq);
    
    %% remove >3std data for quality control
    tmp.alldata=table2array(OBS_bdat(:, OBS.var_ind));
%     tmp.alldata(tmp.alldata==-999)=NaN;
    tmp.alldat_med=median(tmp.alldata, 'omitnan');
    tmp.alldat_3std=3*std(tmp.alldata, 'omitnan');
    tmp.alldat_upper_lim = tmp.alldat_med + tmp.alldat_3std;
    tmp.alldat_lower_lim = tmp.alldat_med - tmp.alldat_3std;


    %%Initialization
    OBS.len_z=length(OBS.depth_interp);
    OBS.(tmp.varname_obs)=NaN(OBS.len_z, OBS.len_Depth_ID_uniq);
    OBS.year_ti=NaN(1,OBS.len_Depth_ID_uniq);
    OBS.month_ti=NaN(1,OBS.len_Depth_ID_uniq);


    
    tmp.tmpsavename_2nd=[dirs.obssavedir, '/obs_2', '_', tmp.varname, '.mat'];
    if (exist(tmp.tmpsavename_2nd) ==0)
        for ti=1:OBS.len_Depth_ID_uniq   %% 470 sec, 512?
            OBS.data_ind=find(strcmp(OBS.Depth_ID_char,OBS.Depth_ID_uniq(ti)));
            tmp.data=table2array(OBS_bdat(OBS.data_ind, OBS.var_ind));
            tmp.depth=table2array(OBS_bdat(OBS.data_ind, OBS.depth_ind));
            %% quality control
            tmp.data(tmp.data < tmp.alldat_lower_lim) = NaN;
            tmp.data(tmp.data > tmp.alldat_upper_lim) = NaN;
    
            tmp.invalid_ind=find(tmp.data==tmp.NaNval);
            tmp.data(tmp.invalid_ind)=NaN;
            tmp.valid_ind=isfinite(tmp.data);
            
            tmp.depth_valid=tmp.depth(tmp.valid_ind);
            tmp.data_valid=tmp.data(tmp.valid_ind);
            [tmp.uniqval, tmp.uniqind] = unique(tmp.depth_valid);
            if length(tmp.depth_valid(tmp.uniqind))<=1
                tmp.data_interped(1:length(OBS.depth_interp)) = NaN;
            else
                tmp.data_interped=interp1(tmp.depth_valid(tmp.uniqind), tmp.data_valid(tmp.uniqind), OBS.depth_interp);
            end
    
            OBS.(tmp.varname_obs)(:,ti)=tmp.data_interped;
            OBS.year_ti(ti)=unique(OBS.year(OBS.data_ind));
            OBS.month_ti(ti)=unique(OBS.month(OBS.data_ind));

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
        save(tmp.tmpsavename_2nd, 'OBS')
    else
        load(tmp.tmpsavename_2nd, 'OBS')
    end
%     OBS.lonmean=mean(OBS.lon);
%     OBS.latmean=mean(OBS.lat);
    
    OBS.lonmean=-117;
    OBS.latmean=32;
    
%     histogram(OBS.lat)
    disp('abc')
    tic;
    for yi=1:length(config.years)
        tmp.year=config.years(yi);
        for mi=1:length(config.months)
            tmp.month=config.months(mi);
            tmp.ind_tmp=find(OBS.year_ti==tmp.year & OBS.month_ti==tmp.month);
            tmp.data=OBS.(tmp.varname_obs)(:,tmp.ind_tmp);
            OBS.([tmp.varname_obs, '_ym'])(:,(yi-1)*12+mi)=mean(tmp.data,2,'omitnan');
            tmp.datenum=datenum(tmp.year, tmp.month, 28);
            OBS.datenum((yi-1)*12+mi)=tmp.datenum;
%             if (isnan(OBS.([tmp.varname_obs, '_ym'])(1,(yi-1)*12+mi)) && tmp.year==2018)
%                 disp('abc')
%             end
        end
    end
    toc;
    a=1;
%     plot(OBS.T_degC(1,:))
%     plot(OBS.datenum, OBS.T_degC_ym(1,:), '*'); datetick;
%     plot(OBS.datenum, OBS.T_degC_ym(1,:)); datetick;

    %% get model data
    for obsi=1:config.len_obs
        tmp.obsname_assm=config.obsnames{obsi};
        config.savename=[dirs.saveroot, filesep, config_obs.staname, filesep, 'CESM2_', config_obs.staname, '_', ...
            tmp.obsname_assm, '_', tmp.varname, '_', num2str(min(config.years)), '_', num2str(max(config.years)), '.mat'];
        config.savename_cfg=[dirs.saveroot, filesep, config_obs.staname, filesep, 'cfg_CESM2_', config_obs.staname, '_', ...
            tmp.obsname_assm, '_', tmp.varname, '_', num2str(min(config.years)), '_', num2str(max(config.years)), '.mat'];
        load(config.savename_cfg, 'CESM2_data', 'CESM2_grid');
        
        %% get subsequent datenum of model
        for tyi=1:config.len_t_y
            for tmi=1:config.len_t_m
                CESM2_grid.datenum((tyi-1)*12+tmi)=datenum(config.years(tyi), config.months(tmi), 15);
            end
        end

        tmp.obs_abb=tmp.obsname_assm(1:3);
        CESM2_grid.lon_cut=CESM2_grid.lon_t(CESM2_grid.obs_lon_ind_w:CESM2_grid.obs_lon_ind_e, CESM2_grid.obs_lat_ind_s:CESM2_grid.obs_lat_ind_n);
        CESM2_grid.lat_cut=CESM2_grid.lat_t(CESM2_grid.obs_lon_ind_w:CESM2_grid.obs_lon_ind_e, CESM2_grid.obs_lat_ind_s:CESM2_grid.obs_lat_ind_n);
        CESM2_grid.ocean_mask_cut=CESM2_grid.ocean_mask(CESM2_grid.obs_lon_ind_w:CESM2_grid.obs_lon_ind_e, CESM2_grid.obs_lat_ind_s:CESM2_grid.obs_lat_ind_n);

        CESM2_alldata.(tmp.obs_abb)=CESM2_data;  % (ens, x, y, t_y, t_m, z)
        CESM2_allgrid.(tmp.obs_abb)=CESM2_grid;
        

        tmp.lon_nan=CESM2_grid.lon_cut.*CESM2_grid.ocean_mask_cut;
        tmp.lat_nan=CESM2_grid.lat_cut.*CESM2_grid.ocean_mask_cut;

        [tmp.lon_ind_w, tmp.lon_ind_e, tmp.lat_ind_s, tmp.lat_ind_n] = ...
                        Func_0012_findind_Y(1, [360+min(OBS.lon), 360+max(OBS.lon), ...
                        median(OBS.lat)-std(OBS.lat), median(OBS.lat)+std(OBS.lat)], ...
                        CESM2_grid.lon_cut, ...
                        CESM2_grid.lat_cut, 'CESM2'); % find valid lon, lat index near station
%         [tmp.lon_ind_w, tmp.lon_ind_e, tmp.lat_ind_s, tmp.lat_ind_n] = ...
%                         Func_0012_findind_Y(2, [360+OBS.lonmean, OBS.latmean], ...
%                         tmp.lon_nan, ...
%                         tmp.lat_nan, 'CESM2'); % find valid lon, lat index near station

%         tmp.lon_ind=tmp.lon_ind_w;
%         tmp.lat_ind=tmp.lat_ind_s;
        tmp.lon_ind=tmp.lon_ind_w:tmp.lon_ind_e;
        tmp.lat_ind=tmp.lat_ind_s:tmp.lat_ind_n;
        %% get depth averaged data
        for avgi=1:length(config_obs.avgdepth)
            for ensi=1:config.len_ens
                for tyi=1:config.len_t_y
                    for tmi=1:config.len_t_m
                        tmp.varname_model_avg=[tmp.varname, '_', num2str(config_obs.avgdepth(avgi), '%04i')]; % varname_save
                        tmp.val=find(CESM2_allgrid.(tmp.obs_abb).z_t/100.0 > config_obs.avgdepth(avgi)); % get depths > avg depth ind
                        tmp.avgdep_ind=min(tmp.val); % closest depth to avg depth ind
                        tmp.data=squeeze(CESM2_alldata.(tmp.obs_abb).(tmp.varname)(ensi,tmp.lon_ind,tmp.lat_ind,tyi,tmi,1:tmp.avgdep_ind)); % get temporary data
                        tmp.data=squeeze(mean(mean(tmp.data,2,'omitnan'),1,'omitnan'));
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
                %% vertical average
                tmp.varname_m_model_avg=[tmp.varname, '_m_', num2str(config_obs.avgdepth(avgi), '%04i')]; % varname_save
                tmp.varname_std_model_avg=[tmp.varname, '_std_', num2str(config_obs.avgdepth(avgi), '%04i')]; % varname_save
                tmp.varname_ano_model_avg=[tmp.varname, '_ano_', num2str(config_obs.avgdepth(avgi), '%04i')]; % varname_save
                tmp.varname_norm_model_avg=[tmp.varname, '_norm_', num2str(config_obs.avgdepth(avgi), '%04i')]; % varname_save

                [CESM2_alldata.([tmp.obs_abb]).(tmp.varname_m_model_avg)(ensi), ...
                    CESM2_alldata.([tmp.obs_abb]).(tmp.varname_std_model_avg)(ensi), ...
                    CESM2_alldata.([tmp.obs_abb]).(tmp.varname_ano_model_avg)(ensi,1:config.len_t), ...
                    CESM2_alldata.([tmp.obs_abb]).(tmp.varname_norm_model_avg)(ensi,1:config.len_t)] = ...
                    get_ano_norm(CESM2_alldata.(tmp.obs_abb).(tmp.varname_model_avg)(ensi,:));
            end
            %% ensemble mean
            tmp.varname_ano_ensm_model_avg=[tmp.varname, '_ano_ensm_', num2str(config_obs.avgdepth(avgi), '%04i')]; % varname_save
            CESM2_alldata.([tmp.obs_abb]).(tmp.varname_ano_ensm_model_avg) = mean(CESM2_alldata.([tmp.obs_abb]).(tmp.varname_ano_model_avg),1);
             
            %% 12-month running mean filter - model
            tmp.varname_model_avg_12rm=[tmp.varname, '_', num2str(config_obs.avgdepth(avgi), '%04i'), '_12rm'];
            tmp.varname_model_avg_12rm_clim=[tmp.varname, '_', num2str(config_obs.avgdepth(avgi), '%04i'), '_12rm_clim'];
            tmp.varname_model_avg_12rm_sm=[tmp.varname, '_', num2str(config_obs.avgdepth(avgi), '%04i'), '_12rm_sm'];       
            for tyi=1:config.len_t_y
                tmp.year=config.years(tyi);
                for tmi=1:config.len_t_m
                    tmp.month=config.months(tmi);
                    tmp.datenum=datenum(tmp.year, tmp.month, 28);
                    if (tmp.datenum > min(CESM2_grid.datenum)+365.25/2 && tmp.datenum < max(CESM2_grid.datenum)-365.25/2)
                        tmp.ind_12runmean = ...
                            find(CESM2_grid.datenum>tmp.datenum-365.25/2 & CESM2_grid.datenum<tmp.datenum+365.25/2);
                        CESM2_alldata.(tmp.obs_abb).(tmp.varname_model_avg_12rm)(1:config.len_ens,(tyi-1)*12+tmi) = ...
                            mean(CESM2_alldata.(tmp.obs_abb).(tmp.varname_model_avg)(1:config.len_ens,tmp.ind_12runmean), 2, 'omitnan');
                    else
                        CESM2_alldata.(tmp.obs_abb).(tmp.varname_model_avg_12rm)(1:config.len_ens, (tyi-1)*12+tmi) = NaN;
                    end
                end
            end
            %% get get climatology from 12-month running mean filtered data
            tmp.data=reshape(CESM2_alldata.(tmp.obs_abb).(tmp.varname_model_avg_12rm), [config.len_ens, config.len_t_m, config.len_t_y]);
            CESM2_alldata.(tmp.obs_abb).(tmp.varname_model_avg_12rm_clim)=mean(tmp.data,3, 'omitnan');
            for tmi=1:12
                CESM2_alldata.(tmp.obs_abb).(tmp.varname_model_avg_12rm_sm)(1:config.len_ens,(0:config.len_t_y-1)*12+tmi) = ...
                    CESM2_alldata.(tmp.obs_abb).(tmp.varname_model_avg_12rm)(1:config.len_ens,(0:config.len_t_y-1)*12+tmi) ...
                    - CESM2_alldata.(tmp.obs_abb).(tmp.varname_model_avg_12rm_clim)(1:config.len_ens, tmi);
            end
            %% ensemble mean
            tmp.varname_ensm_model_avg=[tmp.varname, '_ensm_', num2str(config_obs.avgdepth(avgi), '%04i')]; % varname_save
            CESM2_alldata.([tmp.obs_abb]).(tmp.varname_ensm_model_avg) = mean(CESM2_alldata.([tmp.obs_abb]).(tmp.varname_model_avg),1);
            
            tmp.varname_ensm_model_avg_12rm_sm=[tmp.varname, '_ensm_', num2str(config_obs.avgdepth(avgi), '%04i'), '_12rm_sm']; % varname_save
            CESM2_alldata.([tmp.obs_abb]).(tmp.varname_ensm_model_avg_12rm_sm) = mean(CESM2_alldata.([tmp.obs_abb]).(tmp.varname_model_avg_12rm_sm),1);

        end

        %% get depth slice data
        for avgi=1:length(config_obs.avgdepth)
            for ensi=1:config.len_ens
                for tyi=1:config.len_t_y
                    for tmi=1:config.len_t_m
                        tmp.varname_model_avg=[tmp.varname, '_slice_', num2str(config_obs.avgdepth(avgi), '%04i')]; % varname_save
                        tmp.val=find(CESM2_allgrid.(tmp.obs_abb).z_t/100.0 > config_obs.avgdepth(avgi)); % slice depth
                        tmp.avgdep_ind=min(tmp.val); % closest depth to slice depth ind
                        tmp.data=squeeze(CESM2_alldata.(tmp.obs_abb).(tmp.varname)(ensi,tmp.lon_ind,tmp.lat_ind,tyi,tmi,tmp.avgdep_ind)); % get temporary data
                        tmp.data=squeeze(mean(mean(tmp.data,2,'omitnan'),1,'omitnan'));
                        
                        CESM2_alldata.(tmp.obs_abb).(tmp.varname_model_avg)(ensi,(tyi-1)*12+tmi)= tmp.data;
                    end
                end
            end
            %% ensemble mean
            tmp.varname_ensm_model_avg=[tmp.varname, '_slice_ensm_', num2str(config_obs.avgdepth(avgi), '%04i')]; % varname_save
            CESM2_alldata.([tmp.obs_abb]).(tmp.varname_ensm_model_avg) = mean(CESM2_alldata.([tmp.obs_abb]).(tmp.varname_model_avg),1);
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
    
% plot(OBS.datenum, CESM2_alldata.en4.([tmp.varname,'_0010'])(1,:)-mean(CESM2_alldata.en4.([tmp.varname,'_0010'])(1,:)))
% hold on
% plot(OBS.datenum, OBS.([tmp.varname_obs,'_ym'])(2,:)-mean(OBS.([tmp.varname_obs,'_ym'])(2,:), 'omitnan'), '*'); datetick;
% hold off

% plot(OBS.datenum, CESM2_alldata.en4.([tmp.varname,'_slice_0150'])(1,:)-mean(CESM2_alldata.en4.([tmp.varname,'_slice_0150'])(1,:)))
% hold on
% plot(OBS.datenum, OBS.([tmp.varname_obs,'_ym'])(16,:)-mean(OBS.([tmp.varname_obs,'_ym'])(16,:), 'omitnan'), '*'); datetick;
% hold off


% pcolor(CESM2_allgrid.en4.lon_cut, CESM2_allgrid.en4.lat_cut, squeeze(CESM2_data.SiO3_vm(1,:,:,62,1))); shading flat; colorbar;

    %% time series plot (anomaly) - slice
    
    time_m=min(config.years):1/12:max(config.years)+11/12;
    tmp.plot_color(1,1:3)=[1,0,0];  % red, projd
    tmp.plot_color(2,1:3)=[0,1,0];  % green, en4 
    for avgi=1:length(config_obs.avgdepth)
        tmp.depth=config_obs.avgdepth(avgi);

        tmp.varname_obs_avg=[tmp.varname_obs,'_ym'];
%         tmp.varname_obs_avg=[tmp.varname_obs, '_', num2str(config_obs.avgdepth(avgi), '%04i')];
%         tmp.varname_obs_avg_12rm_sm=[tmp.varname_obs, '_', num2str(config_obs.avgdepth(avgi), '%04i'), '_12rm_sm'];

%         tmp.varname_ensm_model_avg_12rm_sm=[tmp.varname, '_ensm_', num2str(config_obs.avgdepth(avgi), '%04i'), '_12rm_sm']; % varname_save
        tmp.varname_ensm_model_avg=[tmp.varname, '_slice_ensm_', num2str(config_obs.avgdepth(avgi), '%04i')]; % varname_save

        tmp.varname_model_avg=[tmp.varname, '_slice_', num2str(config_obs.avgdepth(avgi), '%04i')];
%         tmp.varname_model_avg_12rm_sm=[tmp.varname, '_', num2str(config_obs.avgdepth(avgi), '%04i'), '_12rm_sm'];       

        hold on
        for obsi=1:config.len_obs
            %% ensemble mean & observation plot
            tmp.obsname_assm=config.obsnames{obsi};
            tmp.obs_abb=tmp.obsname_assm(1:3);
            tmp.depthind=find(OBS.depth_interp==tmp.depth);
%             tmp.data=CESM2_alldata.([tmp.obs_abb]).(tmp.varname_ensm_model_avg_12rm_sm);
            tmp.data=OBS.([tmp.varname_obs,'_ym'])(tmp.depthind,:)-mean(OBS.([tmp.varname_obs,'_ym'])(tmp.depthind,:), 'omitnan');
%             plot(OBS.datenum,tmp.data, 'color', tmp.plot_color(obsi,:), 'linewidth', 3);
            plot(OBS.datenum,tmp.data, 'k*');
            
            tmp.llim=-std(tmp.data,'omitnan')*5;
            tmp.ulim=std(tmp.data,'omitnan')*5;

            %% ensemble member plot (10p)
            tmp.val_transparent = 0.7;
            tmp.plot_color_tr_10p(obsi,:)=get_transparent_rgb(tmp.plot_color(obsi,:), tmp.val_transparent);
            tmp.size_ens=size(CESM2_alldata.(tmp.obs_abb).(tmp.varname_model_avg));
            for ensi=1:tmp.size_ens/2
                tmp.data=CESM2_alldata.([tmp.obs_abb]).(tmp.varname_model_avg)(ensi,:); %data
                tmp.data=tmp.data-mean(tmp.data,'omitnan');
                tsplot{ensi}=plot(OBS.datenum, tmp.data, 'color', tmp.plot_color_tr_10p(obsi,:));
                set(get(get(tsplot{ensi},'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % turn off legend informatino
            end

            %% ensemble member plot (20p)
            tmp.val_transparent = 0.5;
            tmp.plot_color_tr_20p(obsi,:)=get_transparent_rgb(tmp.plot_color(obsi,:), tmp.val_transparent);
            tmp.size_ens=size(CESM2_alldata.(tmp.obs_abb).(tmp.varname_model_avg));
            for ensi=tmp.size_ens/2+1:tmp.size_ens
                tmp.data=CESM2_alldata.([tmp.obs_abb]).(tmp.varname_model_avg)(ensi,:); %data
                tmp.data=tmp.data-mean(tmp.data,'omitnan');
                tsplot{ensi}=plot(OBS.datenum, tmp.data, 'color', tmp.plot_color_tr_20p(obsi,:));
                set(get(get(tsplot{ensi},'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % turn off legend informatino
            end

            datetick('x', 'yyyy')
        end
%         plot(OBS.datenum_12rm,OBS.(tmp.varname_obs_avg_12rm_sm), 'k*');
%         plot(OBS.datenum,OBS.(tmp.varname_obs_avg), 'k', 'linewidth', 3);

        xlabel('Year'); ylabel([tmp.varname, '(', tmp.varname_obs, ')']);
%         lgd=legend([config.obsnames{1}], [config.obsnames{2}], ['OBS(', config_obs.staname,')'], 'Location', 'NorthWest');
        set(gca, 'fontsize', 20)
%         set(lgd, 'fontsize',5)
        grid minor
        hold off
        
        dirs.figdir=['/Volumes/kyy_raid/kimyy/Figure/CESM2/ESP/ASSM_EXP/', ...
            config_obs.staname, '/', num2str(config_obs.avgdepth(avgi), '%04i'), 'm_slice'];
        system(['mkdir -p ', dirs.figdir]);
        tmp.figname=[dirs.figdir, filesep, tmp.varname, '_', tmp.varname_obs, '_', ...
            num2str(config_obs.avgdepth(avgi), '%04i'), 'm_slice_', ...
            num2str(min(config.years)), '_', num2str(max(config.years)), '.tif'];
        ylim([tmp.llim, tmp.ulim]);
        set(gcf, 'PaperPosition', [0, 0, 8, 4]);
        saveas(gcf,tmp.figname,'tif');
        RemoveWhiteSpace([], 'file', tmp.figname);
        close all;
     end

end





% Depth_ID	n.a.	Uses the Cast_ID prefix ([Century]-[Year][Month][ShipCode]-[CastType][Julian Day]-[CastTime]-[Line][Sta]) 
% but adds three additional variables: [Depth][Bottle]-[Rec_Ind]
% Sta_ID	n.a.	Line and Station [Line] [Station]

% % Sta_Code (Station Codes - found in Cast table)
% % "ST" - Standard CalCOFI Station
% % "SCO" - SCCOOS nearshore/20m Station
% % "NRO" - Not Regularly Occupied Original CalCOFI Station
% % "OCO" - Occasionally CalCOFI Occupied
% % "IMX" - IMECOCAL Occupied
% % "NST" - Non-Standard Station
% % "MBR" - MBARI Occupied Station
% % Data_Type (Data Type - found in Cast table)
% % "PR" - Productivity Cast
% % "HY" - Hydrographic Cast
% % "10" - Ten-meter Cast
% % "CT" - Compressed CTD Cast (Low Resolution)
% % "MX" - Mixed CTD and Bottle Data
% % _qual, qua, or q (Quality Code - found in Bottle table; associated with discrete samples; examples: "O_qual", "Chlqua", "PO4q")
% %  Blank - Data OK
% % "4" - Value zeroed due to value below detection limit
% % "6" - Data taken from CTD sensor
% % "8" - Technician thinks value is suspect
% % "9" - Missing Data
% % RecInd (Record Indicator - found in Bottle table)
% % "3" - Observed Data
% % "4" - Educated office guess (ghost)
% % "5" - Data from STD or CTD sensor
% % "6" - Duplicate Depth
% % "7" - Interpolated to a standard depth
% 
% 

function [var_obs, var_flag] = f_varname_obs(var_model)
% metadata

%mmol = 10^-3 mol
%umol = 10^-6 mol
% 1L ~ 1kg
%1000L = 1m^3
% mmol/m^3 ~= umol/kg

var_flag=NaN;
var_obs=NaN;
    switch var_model
        case 'depth'
            var_obs='Depthm';
        case 'SALT'
            var_obs='Salnty'; 
%             var_flag='SALNTY_FLAG_W';
        case 'TEMP'
%             var_obs='Temp'; % Potential temperature
            var_obs='T_degC'; % Water temperature in degrees Celsius        
        case 'DIC' % mmol/m^3
            var_obs='DIC1'; % DISSOLVED INORGANIC CARBON, (umol/kg)
%             var_flag='TCARBN_FLAG_W';
        case 'PO4' %mmol/m^3
            var_obs='PO4uM';  % phosphate, (umol/liter)
%             var_flag='PHSPHT_FLAG_W';
        case 'ALK' % meq/m^3
            var_obs='TA1'; % total alkalinity, umol/kg
%             var_flag='ALKALI_FLAG_W';
        case 'NO3' %mmol/m^3
            var_obs='NO3uM'; % NITRATE (umol/liter)
%             var_flag='NITRAT_FLAG_W';
        case 'SiO3' %mmol/m^3
            var_obs='SiO3uM';  %silicate, (umol/liter)
%             var_flag='SILCAT_FLAG_W';            
        case 'NO2'
            var_obs='NO2uM'; %NITRITE umol /liter
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


function [lat, lon] = cc2lat(li,st)
% DESCRIPTION: Use this function to convert from calcofi grid coordinates 
% to latitude and longitude. This uses a slightly modified (because the 
% original is wrong) version of the CalCOFI gridding algorithm 
% (Eber and Hewitt 1979) (Also see errata 2006).  
% Basically a little trig and some reference point defining.
% 
% INPUT: The calcofi line and station values
% OUTPUT: This function outputs a degree decimal latitude and longitude
% 
% ASSUMPTIONS: All lat/long and station values are from the Northwestern 
%               hemisphere and within the middle latitudes.
% REFERENCE: Eber and Hewitt 1979, Conversion Algorithms for the CalCOFI
% Station Grid
% WRITTEN BY: Robert Thombley (2006), Scripps Institution of Oceanography,
% CalCOFI
% MODIFIED BY: Augusto Valencia (2014), Universidad de Baja California-UABC,
% based on Weber & Moore 2013.

	rlat = 34.15 - .2 * (li - 80)*cos(cRad(30));
	lat = rlat - (1/15)*(st - 60)*sin(cRad(30));
	l1 = (mctr(lat) - mctr(34.15))*tan(cRad(30));
	l2 = (mctr(rlat) - mctr(lat))/(cos(cRad(30))*sin(cRad(30)));
	lon = l1 + l2 + 121.15;
    if lon<=180
        lon = -1*lon;
    else
        lon = (-1*lon)+360;   %Obtain a positive number greater than 180�
    end
end

function rad = cRad(deg)
% DESCRIPTION: A simple helper function that converts degrees to radians
% INPUT: angle in degrees
% OUPUT: angle in radians
% WRITTEN BY: Robert Thombley SIO 2006
	rad = deg * pi/180;
end

function val = mctr(t1)
% DESCRIPTION: Mercator transform function.  
% INPUT: Latitude value in decimal degrees
% OUTPUT: value in mercator units
% WRITTEN BY: Robert Thombley SIO 2006
	val = 180/pi* (log(tan(cRad(45 + t1/2))) - 0.00676866 * sin(cRad(t1)));
end

% % %%read SSH & Tarea
% % SSH=ncread('b.e21.BSSP370cmip6.f09_g17.LE2-1301.009.pop.h.SSH.201501-202412.nc', 'SSH');
% % TAREA=ncread('b.e21.BSSP370cmip6.f09_g17.LE2-1301.009.pop.h.SSH.201501-202412.nc', 'TAREA');
% % pcolor(SSH(:,:,1)'); shading flat; colorbar;
% % pcolor(TAREA'); shading flat; colorbar;
% % 
% % tmpd=squeeze(SSH(:,:,1));
% % pcolor(tmpd'); shading flat; colorbar;
% % 
% % TAREA_valid=TAREA;
% % TAREA_valid(isnan(tmpd))=NaN;
% % pcolor(TAREA_valid'); shading flat; colorbar;
% % 
% % for ti=1:72
% %     tti(ti)=sum(SSH(:,:,ti).*TAREA_valid, 'omitnan')/sum(TAREA_valid, 'omitnan');
% % end
% % plot(tti)
% % 
% % SSH=ncread('b.e21.BSSP370cmip6.f09_g17.LE2-1301.009.pop.h.SSH.209501-210012.nc', 'SSH');
% % 
% % for ti=1:72
% %     tti(ti+72)=sum(SSH(:,:,ti).*TAREA_valid, 'omitnan')/sum(TAREA_valid, 'omitnan');
% % end
% % plot(tti)