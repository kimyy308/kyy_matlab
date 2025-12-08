% %  Updated 21-Apr-2021 by Yong-Yub Kim, 

close all; clear all;  clc;
warning off;

OBS_info.dataroot = ['/Volumes/kyy_raid/Data/Observation/CMEMS/', filesep];

tmp.regionname = 'pollock_egg3'; % NWP, AKP4, ES_KHOA, YS, ...

OBS_info.years=1995:2014; %Jan-Feb mean
OBS_info.season='JF-';
tmp.testname='CMEMS';
tmp.fs=filesep;  

%%     set dropbox path
    addpath(genpath('/Volumes/kyy_raid/kimyy/Dropbox/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_20/run/MICT_pollack/2022_future_pollock/subroutine/'))
tmp.dropboxpath = '/Volumes/kyy_raid/kimyy/Dropbox';
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
[tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Model' ...
    tmp.fs, 'ROMS', tmp.fs, 'Analysis', tmp.fs, 'Figure', tmp.fs, 'nwp_1_20', tmp.fs ...
    'run', tmp.fs, 'SSH', tmp.fs, '2phase_2nd', tmp.fs, 'subroutine']));

%% get colormaps
[cmaps.byrmap3, tmp.error_status] = Func_0009_get_colormaps('byr3', tmp.dropboxpath);
[cmaps.byrmap, tmp.error_status] = Func_0009_get_colormaps('byr2', tmp.dropboxpath);        
[cmaps.yrmap, tmp.error_status] = Func_0009_get_colormaps('yr', tmp.dropboxpath);    

[OBS_grid.refpolygon, OBS_grid.domain, tmp.error_status] = Func_0007_get_polygon_data_from_regionname(tmp.regionname);

tmp.variable ='temp';

dirs.figrawdir =strcat('/Users/kimyy/Desktop/backup/Research/Ph_D_course/2022_pollock_future/figure',filesep, tmp.testname, filesep); % % where figure files will be saved            
tmp.param_script =['/Volumes/kyy_raid/kimyy/Dropbox/source/matlab/Model/ROMS/Analysis/Figure/', 'nwp_1_20', '/run/fig_param/fig_param2_kyy_', tmp.regionname, '.m'];
dirs.filedir = OBS_info.dataroot; % % where data files are                      

tmp.param_script =[tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Model', tmp.fs, ...
'ROMS', tmp.fs, 'Analysis', tmp.fs, 'Figure', tmp.fs, 'nwp_1_20', tmp.fs, 'run', tmp.fs, ...
'fig_param', tmp.fs, 'fig_param_kyy_cmems_', tmp.regionname, '.m'];




% start-------------------- earlier decadal SST, SSS plot
%% set variable name & figure directory
    tmp.variable='temp';
    dirs.figdir=[dirs.figrawdir, tmp.regionname, tmp.fs, ...
        num2str(min(OBS_info.years)), '_', num2str(max(OBS_info.years)), tmp.fs];
    if (exist(strcat(dirs.figdir) , 'dir') ~= 7)
        mkdir(strcat(dirs.figdir));
    end 

%% set figure file name
    tmp.tifname=strcat(dirs.figdir, tmp.testname, '_', ...
        num2str(min(OBS_info.years),'%04i'), '_',num2str(max(OBS_info.years),'%04i'), ...
        '_', OBS_info.season,'.tif'); %% ~_year_month.jpg

        OBS_info.filename=[OBS_info.dataroot,'AKP4cmems_gos_1995_2014_JanFeb.nc'];
        OBS_grid.lon=ncread(OBS_info.filename, 'lon');
        OBS_grid.lat=ncread(OBS_info.filename, 'lat');
        OBS_data.raw_u=ncread(OBS_info.filename, 'cmems_u');
        OBS_data.raw_v=ncread(OBS_info.filename, 'cmems_v');
        ncinfo(OBS_info.filename);
        
        [OBS_grid.lat_rho, OBS_grid.lon_rho] = meshgrid(OBS_grid.lat, OBS_grid.lon);
        [OBS_grid.lon_min, OBS_grid.lon_max, OBS_grid.lat_min, OBS_grid.lat_max] = ...
                            findind_Y(1/20, OBS_grid.domain(1:4), OBS_grid.lon_rho, OBS_grid.lat_rho);
        OBS_grid.cut_lon_rho = ...
            OBS_grid.lon_rho(OBS_grid.lon_min(1):OBS_grid.lon_max(1), OBS_grid.lat_min(1):OBS_grid.lat_max(1));
        OBS_grid.cut_lat_rho = ...
            OBS_grid.lat_rho(OBS_grid.lon_min(1):OBS_grid.lon_max(1), OBS_grid.lat_min(1):OBS_grid.lat_max(1));
        OBS_data.cut_raw_u=OBS_data.raw_u(OBS_grid.lon_min(1):OBS_grid.lon_max(1), OBS_grid.lat_min(1):OBS_grid.lat_max(1),:);
        OBS_data.cut_raw_v=OBS_data.raw_v(OBS_grid.lon_min(1):OBS_grid.lon_max(1), OBS_grid.lat_min(1):OBS_grid.lat_max(1),:);

        OBS_grid.mask_model = double(inpolygon(OBS_grid.cut_lon_rho,OBS_grid.cut_lat_rho,OBS_grid.refpolygon(:,1),OBS_grid.refpolygon(:,2)));
        OBS_grid.mask_model(OBS_grid.mask_model==0)=NaN;
        OBS_data.mean_u=mean(OBS_data.cut_raw_u, 3).*OBS_grid.mask_model;
        OBS_data.mean_v=mean(OBS_data.cut_raw_v, 3).*OBS_grid.mask_model;

        RCM_grid=OBS_grid;
        run(tmp.param_script);

        %% reference vector
        if (isfield(tmp, 'ref_vec_x_range') ~= 1)
            tmp.ref_vec_x_ind = find(abs(OBS_grid.cut_lon_rho(:,1)-param.m_quiver_ref_text_x_location) ...
                == min(abs(OBS_grid.cut_lon_rho(:,1)-param.m_quiver_ref_text_x_location)))+1;
            tmp.ref_vec_y_ind = find(abs(OBS_grid.cut_lat_rho(1,:)-param.m_quiver_ref_text_y_location) ...
                == min(abs(OBS_grid.cut_lat_rho(1,:)-param.m_quiver_ref_text_y_location)))+param.m_quiver_y_interval*2;
%             tmp.ref_vec_x_range = round(tmp.ref_vec_x_ind-(param.m_quiver_x_interval/2)) : ...
%                 round(tmp.ref_vec_x_ind-(param.m_quiver_x_interval/2))+param.m_quiver_x_interval*1;
%             tmp.ref_vec_y_range = round(tmp.ref_vec_y_ind-(param.m_quiver_y_interval/2)) : ...
%                 round(tmp.ref_vec_y_ind-(param.m_quiver_y_interval/2))+param.m_quiver_y_interval*1;
            tmp.ref_vec_x_range = round(tmp.ref_vec_x_ind-(param.m_quiver_x_interval/2)) : ...
                round(tmp.ref_vec_x_ind-(param.m_quiver_x_interval/2));
            tmp.ref_vec_y_range = round(tmp.ref_vec_y_ind-(param.m_quiver_y_interval/2)) : ...
                round(tmp.ref_vec_y_ind-(param.m_quiver_y_interval/2));

        end
        OBS_data.mean_u(tmp.ref_vec_x_range,tmp.ref_vec_y_range)=param.m_quiver_ref_u_value;
        OBS_data.mean_v(tmp.ref_vec_x_range,tmp.ref_vec_y_range)=param.m_quiver_ref_v_value;   



        m_proj(param.m_proj_name,'lon',[OBS_grid.domain(1) OBS_grid.domain(2)],'lat',[OBS_grid.domain(3) OBS_grid.domain(4)]);
            hold on;
       m_quiver(OBS_grid.cut_lon_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)', ...
                    OBS_grid.cut_lat_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)', ...
                    OBS_data.mean_u(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)' * param.m_quiver_vector_size, ...
                    OBS_data.mean_v(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)' * param.m_quiver_vector_size, ...
                    'color', param.m_quiver_vector_color, 'AutoScale','off','LineWidth', param.m_quiver_LineWidth);


        m_gshhs_i('color',param.m_gshhs_line_color)  
        m_gshhs_i('patch',param.m_gshhs_land_color);   % gray colored land
        
        m_text(param.m_quiver_ref_text_x_location, param.m_quiver_ref_text_y_location, param.m_quiver_ref_text, 'FontSize', param.m_quiver_ref_text_fontsize); 
    
        m_grid('fontsize', param.m_grid_fontsize, 'box', param.m_grid_box_type, 'tickdir', param.m_grid_tickdir_type);
%         if min(OBS_info.years) == max(OBS_info.years)
%             tmp.titlename = strcat('vec', ', ', OBS_info.season(1:3), ', ', tmp.abb,',(',num2str(max(OBS_info.years),'%04i'),') ');                        
%         else
%             tmp.titlename = strcat('vec', ', ', OBS_info.season(1:3), ', ', tmp.abb,',(',num2str(min(OBS_info.years),'%04i'),'-',num2str(max(OBS_info.years),'%04i'),') ');
%         end
%         title(tmp.titlename,'fontsize',param.m_pcolor_title_fontsize);  %%title 

    
        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [param.hor_paper_size_x, param.hor_paper_size_y]);
        set(gcf,'PaperPosition', [param.paper_position_hor param.paper_position_ver param.paper_position_width param.paper_position_height]) 
        saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
        close all;
        clear OBS_data.mean
        OBS_grid=rmfield(OBS_grid, 'lon_rho');
