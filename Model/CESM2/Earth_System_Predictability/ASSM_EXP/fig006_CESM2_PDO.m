% %  Created 30-Nov-2022 by Yong-Yub Kim  


%The SST anomalies are obtained by ...
%removing both the climatological annual cycle 
% and the global-mean SST anomaly from the data at each gridpoint.

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

% cfg.years=1967:2021;
cfg.years=1960:2020;

cfg.months=1:12;
cfg.scenname='HIST';
cfg.gridname='f09_g17';
cfg.obsnames={'en4.2_ba','projdv7.3_ba' };
cfg.ensnames={'10p1', '10p2', '10p3', '10p4', '10p5', ...
    '20p1', '20p2', '20p3', '20p4', '20p5'};

cfg.component='ocn';
cfg.varnames={'TEMP'};

cfg.len_t_y = length(cfg.years);
cfg.len_t_m = length(cfg.months);
cfg.len_t = cfg.len_t_y * cfg.len_t_m;
cfg.len_obs= length(cfg.obsnames);
cfg.len_ens= length(cfg.ensnames);
% 

for varind=1:length(cfg.varnames);     tmp.varname=cfg.varnames{varind}; %varname
    for obsi=1:cfg.len_obs;     tmp.obsname_assm=cfg.obsnames{obsi}; %obsname
        for ensi=1:cfg.len_ens;     tmp.ensname=cfg.ensnames{ensi}; %ensname
            cfg.PDO_ncname=[dirs.saveroot, filesep, 'PDO', filesep, 'CESM2_', 'PDO', ...
                '_', tmp.obsname_assm, '-', tmp.ensname, '_', tmp.varname, '_', ...
                '1960_2021', '.nc'];
            cfg.GLO_ncname=[dirs.saveroot, filesep, 'GLO', filesep, 'CESM2_', 'GLO', ...
                '_', tmp.obsname_assm, '-', tmp.ensname, '_', tmp.varname, '_', ...
                '1960_2021', '.nc'];
        
            CESM2_grid.tlong=ncread(cfg.PDO_ncname, 'TLONG');
            CESM2_grid.tlat=ncread(cfg.PDO_ncname, 'TLAT');
            %% get subsequent datenum of model
            for tyi=1:cfg.len_t_y
                for tmi=1:cfg.len_t_m
                    CESM2_grid.datenum((tyi-1)*12+tmi)=datenum(cfg.years(tyi), cfg.months(tmi), 15);
                end
            end
    
            %% get global mean SST, SST in PDO region
            CESM2_data.gm_temp=ncread(cfg.GLO_ncname, 'gm_temp', 1, cfg.len_t);
            CESM2_data.pdo_temp=ncread(cfg.PDO_ncname, 'TEMP', [1 1 1], [inf inf cfg.len_t]);
            CESM2_data.pdom_temp=ncread(cfg.PDO_ncname, 'gm_temp', [1], [cfg.len_t]);

            %% get matrix size
            [cfg.len_x, cfg.len_y, cfg.len_t]= size(CESM2_data.pdo_temp);
    
            %% get climatology (GMSST)
            tmp.data=reshape(CESM2_data.gm_temp, [cfg.len_t_m cfg.len_t_y]);
            CESM2_data.clim_gm_temp=squeeze(mean(tmp.data, 2));
            tmp.data=reshape(CESM2_data.pdom_temp, [cfg.len_t_m cfg.len_t_y]);
            CESM2_data.clim_pdom_temp=squeeze(mean(tmp.data, 2));


            %% get GMSST-A
            for tyi=1:cfg.len_t_y
                for tmi=1:cfg.len_t_m
                    tmp.ti=(tyi-1)*cfg.len_t_m+tmi;
                    CESM2_data.gm_temp_ano(tmp.ti)=CESM2_data.gm_temp(tmp.ti) - CESM2_data.clim_gm_temp(tmi);
                    CESM2_data.pdom_temp_ano(tmp.ti)=CESM2_data.pdom_temp(tmp.ti) - CESM2_data.clim_pdom_temp(tmi);
                end
            end
            plot(CESM2_data.gm_temp_ano);
            hold on;
            plot(CESM2_data.pdom_temp_ano);
            hold off;

            
            %% subtract GMSST-A from PDO SST
            for ti=1:cfg.len_t
                CESM2_data.pdo_temp_det(:,:,ti)=CESM2_data.pdo_temp(:,:,ti)-CESM2_data.gm_temp_ano(ti);
            end
%             CESM2_data.pdo_temp_det=CESM2_data.pdo_temp;
            
            %% get climatology(PDO SST)
            tmp.data=reshape(CESM2_data.pdo_temp_det, [cfg.len_x cfg.len_y cfg.len_t_m cfg.len_t_y]);
            CESM2_data.clim_pdo_temp=squeeze(mean(tmp.data,4));
    
            %% get PDO SST-A
            for tyi=1:cfg.len_t_y
                for tmi=1:cfg.len_t_m
                    tmp.ti=(tyi-1)*cfg.len_t_m+tmi;
                    CESM2_data.pdo_temp_ano(:,:,tmp.ti)=CESM2_data.pdo_temp_det(:,:,tmp.ti) - CESM2_data.clim_pdo_temp(:,:,tmi);
                end
            end
    
            [CESM2_data.lv(obsi,ensi,:,:,:), CESM2_data.pc(obsi,ensi,:,:), CESM2_data.var_exp(obsi,ensi,:)] = Func_0024_EOF_3d(CESM2_data.pdo_temp_ano,4);
            
            a=squeeze(CESM2_data.pc(obsi,ensi,:,:));
            b=squeeze(CESM2_data.lv(obsi,ensi,:,:,:));
            save('a.mat', 'a', 'b', 'CESM2_grid');
%             plot(CESM2_data.gm_temp)
%             plot(CESM2_data.pdom_temp)
%             hold on
%             plot(CESM2_grid.datenum, squeeze(CESM2_data.pc(obsi,ensi,:,1))/std(CESM2_data.pc(obsi,ensi,:,1))); datetick;
%             plot(CESM2_grid.datenum, -squeeze(CESM2_data.pc(obsi,ensi,:,2))/std(CESM2_data.pc(obsi,ensi,:,2))); datetick;

    %         pcolor(-lv(:,:,1)'); shading flat; colorbar; colormap jet
            load('/Volumes/kyy_raid/Data/Observation/ERSST/EOF.mat')
            for modei=1:2
            %% PC plot
                tmp.time=CESM2_grid.datenum;
                tmp.data=squeeze(CESM2_data.pc(obsi,ensi,:,modei))/std(CESM2_data.pc(obsi,ensi,:,modei));
                tmp.data_obs=-ERSST.pc(:,modei) / std(ERSST.pc(:,modei));
                tmp.corrcoef=corrcoef(tmp.data,tmp.data_obs);
                tmp.corrcoef_org=tmp.corrcoef;
                if tmp.corrcoef_org(1,2)<0
                    tmp.data=-tmp.data;
                    tmp.corrcoef=-tmp.corrcoef;
                end
                
                plot(tmp.time,tmp.data, 'linewidth', 2);
                hold on
                plot(tmp.time,tmp.data_obs, 'linewidth', 0.5, 'color', 'k')
                hold off
                xlabel('Year');
                ylabel(['(normalized PC ', num2str(modei),')']);
                datetick('x', 'yyyy');
                
                set(gca, 'fontsize', 20);
                grid minor
                dirs.figdir=['/Volumes/kyy_raid/kimyy/Figure/CESM2/ESP/ASSM_EXP/PDO_SST_EOF/', ...
                    'mode', num2str(modei)];
                system(['mkdir -p ', dirs.figdir]);
                tmp.figname=[dirs.figdir, filesep, ...
                    tmp.obsname_assm, '_', tmp.ensname, '_', 'PC_', num2str(modei), '_', ...
                    num2str(min(cfg.years)), '_', num2str(max(cfg.years)), '.tif'];
                set(gcf, 'PaperPosition', [0, 0, 8, 4]);
                text(tmp.time(5), min(tmp.data), num2str(round(tmp.corrcoef(1,2),2)))
%                 'units','normalized','horizontalalignment','center','verticalalignment','middle','fontsize',14,'fontname','freeserif','interpreter','none')


                saveas(gcf,tmp.figname,'tif');
                RemoveWhiteSpace([], 'file', tmp.figname);
                close all;

            %% LV pcolor
                fig_config.name_rgn = 'Glob';
%                 fig_config.map_proj = 'eqdcylin';  % robinson, eqdcylin
                fig_config.map_proj = 'robinson';  % robinson, eqdcylin

                fig_config.x_lim = [min(CESM2_grid.tlong(:)) max(CESM2_grid.tlong(:))];
                fig_config.y_lim = [min(CESM2_grid.tlat(:)) max(CESM2_grid.tlat(:))];
                fig_config.fig_size = [0,0,6,3.5];
                fig_config.ax_size = [0.3,0.7,5.4,2.7];
%                 fig_config.cb_size = [5.15,0.8,0.15,2.3];
                fig_config.cb_size = [6.85,1.7,0.15,2.3];
                fig_config.title_pos = [0.5,0.93];
                fig_config.p_lim =0.1;
                fig_config.c_lim = [-1 1];
%                 [fig_config.c_map, tmp.err_stat] = Func_0009_get_colormaps('jet', tmp.dropboxpath);
                [fig_config.c_map, tmp.err_stat] = Func_0009_get_colormaps('byr', tmp.dropboxpath);

                %% map setting
                ax_m = axesm('MapProjection',fig_config.map_proj,'grid','on','fontsize',14, ...
                    'fontname','freeserif'); 
                axis off; 
                hold on;
                setm(ax_m,'origin',[0,180],'MapLatLimit',fig_config.y_lim,'MapLonLimit',fig_config.x_lim);
                
                %% caxis & colorbar
                if tmp.corrcoef_org(1,2)<0
                    tmp.C=-squeeze(CESM2_data.lv(obsi,ensi,:,:,modei));                    
                else
                    tmp.C=squeeze(CESM2_data.lv(obsi,ensi,:,:,modei));
                end
                fig_config.c_lim = [-max(abs(tmp.C(:))) max(abs(tmp.C(:)))];
%                 fig_config.c_lim = [-0.08 0.08];

                caxis(ax_m, fig_config.c_lim); 
                colormap(fig_config.c_map);
                cb = colorbar(ax_m,'units','inches','position',fig_config.cb_size);
                set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');

                tmp.X=CESM2_grid.tlong;
                tmp.Y=CESM2_grid.tlat;
                if tmp.corrcoef_org(1,2)<0
                    tmp.C=-squeeze(CESM2_data.lv(obsi,ensi,:,:,modei));                    
                else
                    tmp.C=squeeze(CESM2_data.lv(obsi,ensi,:,:,modei));
                end
                S = shaperead('landareas.shp');
                h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
                shading flat;
                geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
                setm(ax_m,'frame','on','FLineWidth',1);
                label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
                label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
                mlabel; plabel;
                label_y=plabel; label_x=mlabel;

                tmp.figname=[dirs.figdir, filesep, ...
                    tmp.obsname_assm, '_', tmp.ensname, '_', 'LV_', num2str(modei), '_', ...
                    num2str(min(cfg.years)), '_', num2str(max(cfg.years)), '.tif'];
                saveas(gcf,tmp.figname,'tif');
                RemoveWhiteSpace([], 'file', tmp.figname);
                close all;

            %% ERSST LV pcolor

                tmp.figname=[dirs.figdir, filesep, ...
                    'ERSST', '_', 'LV_', num2str(modei), '_', ...
                    num2str(min(cfg.years)), '_', num2str(max(cfg.years)), '.tif'];
                if ~exist(tmp.figname)
                        
                    fig_config.name_rgn = 'Glob';
    %                 fig_config.map_proj = 'eqdcylin';  % robinson, eqdcylin
                    fig_config.map_proj = 'robinson';  % robinson, eqdcylin
    
                    fig_config.x_lim = [min(ERSST.lon2_pdo(:)) max(ERSST.lon2_pdo(:))];
                    fig_config.y_lim = [min(ERSST.lat2_pdo(:)) max(ERSST.lat2_pdo(:))];
                    fig_config.fig_size = [0,0,6,3.5];
                    fig_config.ax_size = [0.3,0.7,5.4,2.7];
    %                 fig_config.cb_size = [5.15,0.8,0.15,2.3];
                    fig_config.cb_size = [6.85,1.7,0.15,2.3];
                    fig_config.title_pos = [0.5,0.93];
                    fig_config.p_lim =0.1;
                    fig_config.c_lim = [-1 1];
                    [fig_config.c_map, tmp.err_stat] = Func_0009_get_colormaps('byr', tmp.dropboxpath);
    
                    %% map setting
                    ax_m = axesm('MapProjection',fig_config.map_proj,'grid','on','fontsize',14, ...
                        'fontname','freeserif'); 
                    axis off; 
                    hold on;
                    setm(ax_m,'origin',[0,180],'MapLatLimit',fig_config.y_lim,'MapLonLimit',fig_config.x_lim);
                    
                    %% caxis & colorbar
                    tmp.C=-squeeze(ERSST.lv(:,:,modei));                    
                    fig_config.c_lim = [-max(abs(tmp.C(:))) max(abs(tmp.C(:)))];
                    caxis(ax_m, fig_config.c_lim); 

                    colormap(fig_config.c_map);
                    cb = colorbar(ax_m,'units','inches','position',fig_config.cb_size);
                    set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
    
                    tmp.X=ERSST.lon2_pdo;
                    tmp.Y=ERSST.lat2_pdo;
                    tmp.C=-squeeze(ERSST.lv(:,:,modei));                    
    
                    S = shaperead('landareas.shp');
                    h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
                    shading flat;
                    geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);
                    setm(ax_m,'frame','on','FLineWidth',1);
                    label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
                    label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
                    mlabel; plabel;
                    label_y=plabel; label_x=mlabel;
    
                    
                    saveas(gcf,tmp.figname,'tif');
                    RemoveWhiteSpace([], 'file', tmp.figname);
                    close all;
                end
            end

        end

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