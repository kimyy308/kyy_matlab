function status=plot_ROMS_monthly_surface_UV_data(testname, outfile, filedir, lonlat, tempyear, inputmonth, shadlev, conlev, var, regionflag, param_script)
   run(param_script);
   for monthij = 1:length(inputmonth)
        tempmonth = inputmonth(monthij);
        % ex : rootdir\test37\data\2001\test37_monthly_2001_01.nc
        filename = strcat(filedir, num2str(tempyear,'%04i'), '\', ...
                testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');

        % read data
        if (exist('lon_u')==0)
            lon_u = ncread(filename,'lon_u');
            lat_u = ncread(filename,'lat_u');
            lon_v = ncread(filename,'lon_v');
            lat_v = ncread(filename,'lat_v');
            lon = ncread(filename,'lon_rho');
            lat = ncread(filename,'lat_rho');
        end
        data_info_u = ncinfo(filename, 'u');
        data_info_v = ncinfo(filename, 'v');
        u = ncread(filename,'u',[1 1 data_info_u.Size(3) 1], [data_info_u.Size(1) data_info_u.Size(2) 1 1]);
        v = ncread(filename,'v',[1 1 data_info_v.Size(3) 1], [data_info_v.Size(1) data_info_v.Size(2) 1 1]);
        u_rho=griddata(lon_u, lat_u, u, lon, lat);
        v_rho=griddata(lon_v, lat_v, v, lon, lat);
        
        varname = 'temp';
        data_info = ncinfo(filename, varname);
        data = ncread(filename,varname,[1 1 data_info.Size(3) 1], [data_info.Size(1) data_info.Size(2) 1 1]);

        
        
        % Reference vector value setting
        u_rho(m_quiver_ref_vec_x_range,m_quiver_ref_vec_y_range) = m_quiver_ref_u_value;
        v_rho(m_quiver_ref_vec_x_range,m_quiver_ref_vec_y_range) = m_quiver_ref_v_value; 
        
        % plot
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        hold on;
        m_pcolor(lon,lat,data);
        shading(gca,m_pcolor_shading_method);
        hold on;
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
%         fig=gcf;
%         print(jpgname,'-djpeg');
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