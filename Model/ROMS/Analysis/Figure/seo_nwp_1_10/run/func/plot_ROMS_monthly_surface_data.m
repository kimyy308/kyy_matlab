function status=plot_ROMS_monthly_surface_data(testname, outfile, filedir, lonlat, tempyear, inputmonth, shadlev, conlev, var, param_script)
   run(param_script);
   for monthij = 1:length(inputmonth)
        tempmonth = inputmonth(monthij);
        % ex : rootdir\test37\data\2001\test37_monthly_2001_01.nc
        filename = strcat(filedir, num2str(tempyear,'%04i'), '\', ...
                testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');

        % read data
        lon = ncread(filename,'lon_rho');
        lat = ncread(filename,'lat_rho');
        data_info = ncinfo(filename, varname);
        if (strcmp(var,'SSH'))
            data = ncread(filename,varname,[1 1 1], [data_info.Size(1) data_info.Size(2) 1]);
        elseif (strcmp(var,'YSBCW'))
            data = ncread(filename,varname,[1 1 1 1], [data_info.Size(1) data_info.Size(2) 1 1]);
        else
            data = ncread(filename,varname,[1 1 data_info.Size(3) 1], [data_info.Size(1) data_info.Size(2) 1 1]);
        end
        % plot
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        hold on;
        m_pcolor(lon,lat,data);
        shading(gca,m_pcolor_shading_method);
        m_gshhs_i('color',m_gshhs_line_color)  
        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        titlename = strcat(var, ' (',char(calendarname(tempmonth)), ', ', num2str(tempyear),')');
        title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

        % contour
        [C,h2]=m_contour(lon,lat, data, conlev, m_contour_color, 'linewidth', m_contour_linewidth);
        clabel(C,h2,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
            'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);
        
        % set colorbar 
        h = colorbar;
        colormap(colormap_style);
        set(h,'fontsize',colorbar_fontsize);
        title(h,colorbar_title,'fontsize',colorbar_title_fontsize);
        caxis(shadlev);
        
        % set grid
        m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
        
% % note
% % Starting in R2014b, you can use dot notation to query and set properties.
% % If you are using an earlier release, use the get and set functions instead.
%         ax = gca;
%         outerpos = ax.OuterPosition;
%         ti = ax.TightInset; 
%         left = outerpos(1) + ti(1);
%         bottom = outerpos(2) + ti(2);
%         ax_width = outerpos(3) - ti(1) - ti(3);
%         ax_height = outerpos(4) - ti(2) - ti(4);
%         ax.Position = [left bottom ax_width ax_height];
%         
%         fig = gcf;
%         fig.PaperPositionMode = 'auto'
%         fig_pos = fig.PaperPosition;
%         fig.PaperSize = [fig_pos(3) fig_pos(4)];
%         print(fig,jpgname,'-djpg');

        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
        
        jpgname=strcat(outfile, '_', testname, '_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg
        saveas(gcf,jpgname,'jpg');

        disp(' ')
        disp([num2str(tempyear), '_', num2str(tempmonth), '_', var, ' plot is created.'])
        disp(' ')
        disp([' File path is : ',jpgname])
        disp(' ')
        
        hold off
        close all;
    end
    status=1;
return