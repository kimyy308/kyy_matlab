% %  Updated 21-Apr-2021 by Yong-Yub Kim, 

close all; clear all;  clc;
warning off;

OBS_info.dataroot = ['/Volumes/kyy_raid/Data/Observation/OISST/monthly_kimyy/', filesep];

tmp.regionname = 'pollock_egg3'; % NWP, AKP4, ES_KHOA, YS, ...

OBS_info.years=1995:2014; %Jan-Feb mean
OBS_info.season='JF-';
tmp.testname='OISST';
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
'fig_param', tmp.fs, 'fig_param2_kyy_', tmp.regionname, '.m'];




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

        OBS_info.filename=[OBS_info.dataroot,'avhrr_only_monthly_v2_1995-2014_JanFeb.nc'];
        OBS_grid.lon=ncread(OBS_info.filename, 'lon');
        OBS_grid.lat=ncread(OBS_info.filename, 'lat');
        OBS_data.raw=ncread(OBS_info.filename, 'temp');
        ncinfo(OBS_info.filename);
        
        [OBS_grid.lat_rho, OBS_grid.lon_rho] = meshgrid(OBS_grid.lat, OBS_grid.lon);
        [OBS_grid.lon_min, OBS_grid.lon_max, OBS_grid.lat_min, OBS_grid.lat_max] = ...
                            findind_Y(1/20, OBS_grid.domain(1:4), OBS_grid.lon_rho, OBS_grid.lat_rho);
        OBS_grid.cut_lon_rho = ...
            OBS_grid.lon_rho(OBS_grid.lon_min(1):OBS_grid.lon_max(1), OBS_grid.lat_min(1):OBS_grid.lat_max(1));
        OBS_grid.cut_lat_rho = ...
            OBS_grid.lat_rho(OBS_grid.lon_min(1):OBS_grid.lon_max(1), OBS_grid.lat_min(1):OBS_grid.lat_max(1));
        OBS_data.cut_raw=OBS_data.raw(OBS_grid.lon_min(1):OBS_grid.lon_max(1), OBS_grid.lat_min(1):OBS_grid.lat_max(1),:);

        OBS_grid.mask_model = double(inpolygon(OBS_grid.cut_lon_rho,OBS_grid.cut_lat_rho,OBS_grid.refpolygon(:,1),OBS_grid.refpolygon(:,2)));
        OBS_grid.mask_model(OBS_grid.mask_model==0)=NaN;
        OBS_data.mean=mean(OBS_data.cut_raw, 3).*OBS_grid.mask_model;

        RCM_grid=OBS_grid;
        run(tmp.param_script);


        m_proj(param.m_proj_name,'lon',[OBS_grid.domain(1) OBS_grid.domain(2)],'lat',[OBS_grid.domain(3) OBS_grid.domain(4)]);
        hold on;

        [tmp.m_value, tmp.error_status] = Func_0011_get_area_weighted_mean(OBS_data.mean, OBS_grid.cut_lon_rho, OBS_grid.cut_lat_rho);
        m_pcolor(OBS_grid.cut_lon_rho',OBS_grid.cut_lat_rho',OBS_data.mean');
        shading(gca,param.m_pcolor_shading_method);   
        [C,h2]=m_contour(OBS_grid.cut_lon_rho',OBS_grid.cut_lat_rho', OBS_data.mean', [2, 5, 10], 'color','k', ...
                'linewidth', 1.5, 'linestyle', '-');
            clabel(C,h2,'FontSize',13,'Color','k', ...
                'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
 
        m_gshhs_i('color',param.m_gshhs_line_color)  
        m_gshhs_i('patch',param.m_gshhs_land_color);   % gray colored land

        m_grid('fontsize', param.m_grid_fontsize, 'box', param.m_grid_box_type, 'tickdir', param.m_grid_tickdir_type);
        % set colorbar 
        h = colorbar;
%                     if (strcmp(tmp.variable,'SSH')==1)
%                         colormap(flip(cool));
%                     else
            colormap(jet);
%                     end

        set(h,'fontsize',param.colorbar_fontsize);
        title(h,'(^oC)','fontsize',param.colorbar_title_fontsize);
        caxis([-2,20])

        disp(['M = ', num2str(tmp.m_value)]);
%         m_text(param.m_pcolor_ref_text_x_location, param.m_pcolor_ref_text_y_location, ['M = ', num2str(tmp.m_value)], 'FontSize', param.m_quiver_ref_text_fontsize); 

        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [param.hor_paper_size_x, param.hor_paper_size_y]);
        set(gcf,'PaperPosition', [param.paper_position_hor param.paper_position_ver param.paper_position_width param.paper_position_height]) 
        saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
        close all;
        clear OBS_data.mean
        OBS_grid=rmfield(OBS_grid, 'lon_rho');
