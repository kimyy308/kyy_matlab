function status=vec_HYCOM_monthly_surface_UV(testname, outfile, filedir, lonlat, tempyear, inputmonth, shadlev, conlev, var, regionflag, param_script)
   run(param_script);
   for monthij = 1:length(inputmonth)
        tempmonth = inputmonth(monthij);
        % ex : rootdir\test37\data\2001\test37_monthly_2001_01.nc
        filename = strcat(filedir, num2str(tempyear,'%04i'), '\', ...
                'HYCOM_', num2str(tempyear,'%04i'), num2str(tempmonth,'%02i'), '.nc');

        % read data
        if (exist('lon')==0)
            lon = ncread(filename,'lon');
            lat = ncread(filename,'lat');

            lon_west = abs(lon - (lonlat(1)-1));
            min_lon_west=min(lon_west);
            lon_east = abs(lon - (lonlat(2)+1));
            min_lon_east=min(lon_east);
            lat_south = abs(lat - (lonlat(3)-1));
            min_lat_south=min(lat_south);
            lat_north = abs(lat - (lonlat(4)+1));
            min_lat_north=min(lat_north);

            lon_min = find(lon_west == min_lon_west);
            lon_max = find(lon_east == min_lon_east);
            lat_min = find(lat_south == min_lat_south);
            if (lonlat(4)>47)
                lat_max=length(lat);
            else
                lat_max = find(lat_north == min_lat_north);
            end
            
%             lon = lon(lon_min(1):lon_max(1));
%             lat = lat(lat_min(1):lat_max(1));
            lon = ncread(filename,'longitude', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
            lat = ncread(filename,'latitude', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
        end
        
        
        data_info_u = ncinfo(filename, 'u');
        data_info_v = ncinfo(filename, 'v');
        
        u_rho = ncread(filename,'u',[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);
        v_rho = ncread(filename,'v',[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);
        u_rho(find(abs(u_rho)>10000))=NaN;
        v_rho(find(abs(v_rho)>10000))=NaN;
        
        data_info = ncinfo(filename, varname);  %% [lon lat depth time] -> [1601 1201 33 1]
        if (strcmp(var,'SSH'))
            data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
%         elseif (strcmp(var,'YSBCW'))
%             data = ncread(filename,varname,[1 1 1 1], [data_info.Size(1) data_info.Size(2) 1 1]);
        else
            data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);
        end
        data(find(abs(data)>10000))=NaN;
        
        % Reference vector value setting
        
        u_rho(m_quiver_ref_vec_x_range,m_quiver_ref_vec_y_range) = m_quiver_ref_u_value;
        v_rho(m_quiver_ref_vec_x_range,m_quiver_ref_vec_y_range) = m_quiver_ref_v_value; 
        
        % plot
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        hold on;
        m_pcolor(lon,lat,data);
        shading(gca,m_pcolor_shading_method);
        uvplot=m_quiver(lon(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
                        lat(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
                        u_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
                        v_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
                        m_quiver_vector_size,'color',m_quiver_vector_color, 'LineWidth', m_quiver_LineWidth);
        m_gshhs_i('color',m_gshhs_line_color)  
        m_gshhs_i('patch',m_gshhs_land_color);
        titlename = strcat('Surface current', ' (',char(calendarname(tempmonth)), ', ', num2str(tempyear),')');
        title(titlename,'fontsize', m_quiver_title_fontsize);
        m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_text, 'FontSize', m_quiver_ref_text_fontsize); 
        
        [C,h2]=m_contour(lon,lat, data, conlev, m_contour_color, 'linewidth', m_contour_linewidth);
        
        % set colorbar 
        h = colorbar;
        colormap(colormap_style);
        set(h,'fontsize',colorbar_fontsize);
        title(h,colorbar_title,'fontsize',colorbar_title_fontsize);
        caxis(shadlev);
        
        % set grid
        m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
        
        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
        
        jpgname=strcat(outfile, '_', testname, '_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg
        saveas(gcf,jpgname,'jpg');

        disp(' ')
        disp([num2str(tempyear), '_', num2str(tempmonth), '_', 'UV', ' plot is created.'])
        disp(' ')
        disp([' File path is : ',jpgname])
        disp(' ')
        
        hold off
        close all;
    end
    status=1;
return