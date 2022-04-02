% %  Updated 09-Oct-2021 by Yong-Yub Kim, structure


% start-------------------- earlier decadal current plot

dirs.figdir=[dirs.figrawdir,'surface', tmp.fs, tmp.regionname, tmp.fs, 'wind', tmp.fs, ...
    num2str(min(RCM_info.years)), '_', num2str(max(RCM_info.years)), tmp.fs];
if (exist(strcat(dirs.figdir) , 'dir') ~= 7)
    mkdir(strcat(dirs.figdir));
end 
%             outfile = strcat(figdir,regionname,);

tmp.tifname=strcat(dirs.figdir, tmp.testname, '_clim_wind_',num2str(min(RCM_info.years),'%04i'), ...
    '_',num2str(max(RCM_info.years),'%04i'), ...
                    '_', num2str(RCM_info.months(1),'%04i'), '-', num2str(RCM_info.months(end),'%04i'),'.tif'); %% ~_year_month.jpg
tmp.variable='wind';
if (exist(tmp.tifname , 'file') ~= 2 || flags.fig_switch(7)==2)      
    tmp.matname = [dirs.matdir, tmp.testname, '_', tmp.regionname, '_', tmp.variable, ...
        '_mean_', num2str(min(RCM_info.years),'%04i'), '-', num2str(max(RCM_info.years),'%04i'), ...
                    '_', num2str(RCM_info.months(1),'%04i'), '-', num2str(RCM_info.months(end),'%04i'),'.mat'];
    if (exist(tmp.matname , 'file') ~= 2 || flags.fig_switch(7)==2)
        for yearij=1:length(RCM_info.years)
            tmp.tempyear=RCM_info.years(yearij);
            tmp.yearstr=num2str(tmp.tempyear, '%04i');
            for monthij=1:length(RCM_info.months)
                tmp.tempmonth=RCM_info.months(monthij);
                tmp.monthstr=num2str(tmp.tempmonth, '%02i');
%                         D:\Data\Model\ROMS\nwp_1_20\backup_surf\test2102\run\v\1985\pck_test2102_v_monthly_1985_01.nc
                switch tmp.testname
                    case {'test2102', 'test2103', 'test2104', 'test2105', 'test2106'}
                        tmp.Uwindfilename=[dirs.filedir, 'Uwind', tmp.fs, tmp.yearstr, filesep,'pck_', tmp.testname, '_Uwind_monthly_', tmp.yearstr, '_', tmp.monthstr, '.nc'];
                        tmp.Vwindfilename=[dirs.filedir, 'Vwind', tmp.fs, tmp.yearstr, filesep,'pck_', tmp.testname, '_Vwind_monthly_', tmp.yearstr, '_', tmp.monthstr, '.nc'];
                    otherwise
                        tmp.Uwindfilename=[dirs.filedir, 'Uwind', tmp.fs, tmp.yearstr, filesep,'NWP_pck_', tmp.testname, '_Uwind_monthly_', tmp.yearstr, '_', tmp.monthstr, '.nc'];
                        tmp.Vwindfilename=[dirs.filedir, 'Vwind', tmp.fs, tmp.yearstr, filesep,'NWP_pck_', tmp.testname, '_Vwind_monthly_', tmp.yearstr, '_', tmp.monthstr, '.nc'];
                end
                
               
                if (isfield(RCM_grid, 'lon_rho') ~= 1)
                    for gridi=1:length(RCM_grid.gridname)
                        RCM_grid.(RCM_grid.gridname{gridi})=ncread(RCM_grid.(['filename_', RCM_grid.gridname{gridi}]), RCM_grid.gridname{gridi});
                    end
%                             lat_rho=ncread(filename, 'lat_rho');
%                             lon_u=ncread(filename, 'lon_u');
%                             lat_u=ncread(filename, 'lat_u');
%                             lon_v=ncread(filename, 'lon_v');
%                             lat_v=ncread(filename, 'lat_v');

                    [RCM_grid.lon_min, RCM_grid.lon_max, RCM_grid.lat_min, RCM_grid.lat_max] = ...
                        findind_Y(1/20, RCM_grid.domain(1:4), RCM_grid.lon_rho, RCM_grid.lat_rho);
%                             [lon_u_min, lon_u_max, lat_u_min, lat_u_max] = findind_Y(1/20, lonlat(1:4), lon_u, lat_u);
%                             [lon_v_min, lon_v_max, lat_v_min, lat_v_max] = findind_Y(1/20, lonlat(1:4), lon_v, lat_v);
                    cut_lon_rho = ...
                        RCM_grid.lon_rho(RCM_grid.lon_min(1):RCM_grid.lon_max(1), RCM_grid.lat_min(1):RCM_grid.lat_max(1));
                    cut_lat_rho = ...
                        RCM_grid.lat_rho(RCM_grid.lon_min(1):RCM_grid.lon_max(1), RCM_grid.lat_min(1):RCM_grid.lat_max(1));
                end
                tmp.data_info = ncinfo(tmp.Uwindfilename, 'Uwind'); 
                tmp.Uwind = ncread(tmp.Uwindfilename,'Uwind', ...
                    [RCM_grid.lon_min(1) RCM_grid.lat_min(1) 1], ...
                    [RCM_grid.lon_max(1)-RCM_grid.lon_min(1)+1 RCM_grid.lat_max(1)-RCM_grid.lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                tmp.Vwind = ncread(tmp.Vwindfilename,'Vwind', ...
                    [RCM_grid.lon_min(1) RCM_grid.lat_min(1) 1], ...
                    [RCM_grid.lon_max(1)-RCM_grid.lon_min(1)+1 RCM_grid.lat_max(1)-RCM_grid.lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                if (isfield(tmp, 'mean_Uwind') ~= 1)
                    tmp.mean_Uwind=zeros(size(tmp.Uwind));
                    tmp.mean_Vwind=zeros(size(tmp.Vwind));
                end
                tmp.mean_Uwind=tmp.mean_Uwind + (tmp.Uwind / (length(RCM_info.years) * length(RCM_info.months)));
                tmp.mean_Vwind=tmp.mean_Vwind + (tmp.Vwind / (length(RCM_info.years) * length(RCM_info.months)));
            end
        end

        Uwind= tmp.mean_Uwind;
        Vwind= tmp.mean_Vwind;
        if (isfield(tmp, 'ref_vec_x_range') ~= 1)
            tmp.ref_vec_x_ind = find(abs(cut_lon_rho(:,1)-param.m_quiver_ref_text_x_location) ...
                == min(abs(cut_lon_rho(:,1)-param.m_quiver_ref_text_x_location)));
            tmp.ref_vec_y_ind = find(abs(cut_lat_rho(1,:)-param.m_quiver_ref_text_y_location) ...
                == min(abs(cut_lat_rho(1,:)-param.m_quiver_ref_text_y_location)))+param.m_quiver_wind_y_interval*2;
            tmp.ref_vec_x_range = round(tmp.ref_vec_x_ind-(param.m_quiver_x_interval/2)) : ...
                round(tmp.ref_vec_x_ind-(param.m_quiver_x_interval/2))+param.m_quiver_wind_x_interval;
            tmp.ref_vec_y_range = round(tmp.ref_vec_y_ind-(param.m_quiver_y_interval/2)) : ...
                round(tmp.ref_vec_y_ind-(param.m_quiver_y_interval/2))+param.m_quiver_wind_y_interval;
        end
        Uwind(tmp.ref_vec_x_range,tmp.ref_vec_y_range)=param.m_quiver_ref_uwind_value;
        Vwind(tmp.ref_vec_x_range,tmp.ref_vec_y_range)=param.m_quiver_ref_vwind_value;     
        save(tmp.matname, 'Uwind','Vwind', 'cut_lon_rho', 'cut_lat_rho');
    else
        load(tmp.matname);
    end

    RCM_grid.mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,RCM_grid.refpolygon(:,1),RCM_grid.refpolygon(:,2)));
    RCM_grid.mask_model(RCM_grid.mask_model==0)=NaN;
    Uwind=Uwind(1:size(RCM_grid.mask_model,1),:).*RCM_grid.mask_model;
    Vwind=Vwind(:,1:size(RCM_grid.mask_model,2)).*RCM_grid.mask_model;  
    
    m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);
    hold on;

    m_gshhs_i('color',param.m_gshhs_line_color)  
    m_gshhs_i('patch',param.m_gshhs_land_color);   % gray colored land

    uvplot=m_quiver(cut_lon_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_wind_x_interval:end)', ...
                    cut_lat_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_wind_y_interval:end)', ...
                    Uwind(1:param.m_quiver_x_interval:end, 1:param.m_quiver_wind_x_interval:end)' * param.m_quiver_wind_vector_size, ...
                    Vwind(1:param.m_quiver_x_interval:end, 1:param.m_quiver_wind_y_interval:end)' * param.m_quiver_wind_vector_size, ...
                    'color', param.m_quiver_vector_color, 'AutoScale','off','LineWidth', param.m_quiver_LineWidth);

    m_text(param.m_quiver_ref_text_x_location, param.m_quiver_ref_text_y_location, param.m_quiver_wind_ref_text, 'FontSize', param.m_quiver_ref_text_fontsize); 
    m_grid('fontsize', param.m_grid_fontsize, 'box', param.m_grid_box_type, 'tickdir', param.m_grid_tickdir_type);
    if min(RCM_info.years) == max(RCM_info.years)
        param.titlename = strcat('Wind, ',  RCM_info.season(1:3), ', ', tmp.abb,',(',num2str(min(RCM_info.years),'%04i'),') ');  %% + glacier contribution
    else
        param.titlename = strcat('Wind, ',  RCM_info.season(1:3), ', ', tmp.abb,',(',num2str(min(RCM_info.years),'%04i'),'-',num2str(max(RCM_info.years),'%04i'),') ');  %% + glacier contribution
    end
    title(param.titlename,'fontsize',param.m_pcolor_title_fontsize);  %%title

    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperSize', [param.hor_paper_size_x, param.hor_paper_size_y]);
    set(gcf,'PaperPosition', [param.paper_position_hor param.paper_position_ver param.paper_position_width param.paper_position_height]) 
    saveas(gcf, tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
%                 system(['magick ', tifname, ' -trim ', tifname]);

    close all;
%             clear RCM_grid.lon_rho tmp.mean_Uwind tmp.ref_vec_x_range
    RCM_grid=rmfield(RCM_grid, 'lon_rho');
    tmp=rmfield(tmp, 'mean_Uwind');
    tmp=rmfield(tmp, 'ref_vec_x_range');
end
