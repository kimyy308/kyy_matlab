function status=plot_OISST_monthly_hor_t(testname, outfile, filedir, section, tempyear,inputmonth, shadlev, conlev,  param_script)
% clear all;close all;
%==========================================================================
% %  Updated 01-Oct-2018 by Yong-Yub Kim 

lonlat = section(1:4);
var = 'temp';
varname = 'temp';
run(param_script);
%==========================================================================

for monthij=1:length(inputmonth)
    tempmonth = inputmonth(monthij);
    filename = strcat(filedir, testname, num2str(tempyear,'%04i'), '.nc');

    if (exist('lon' , 'var') ~= 1)
        lon  = ncread(filename, 'lon');
        lat  = ncread(filename, 'lat');
        Nlon=length(lon);
        Nlat=length(lat);
        dl=1/4;
        [lon_rho lat_rho]=meshgrid(lon, lat);
        [lon_min, lon_max, lat_min, lat_max] = findind_Y(dl, section(1:4), lon_rho', lat_rho');
    end
    
    data_info = ncinfo(filename, varname); 
    data = ncread(filename,varname,[lon_min(1) lat_min(1) tempmonth], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
    data=squeeze(data);

    cut_data = permute(data, [2 1]);  %% permute [x y z] -> [z y x]
    
    cut_lon_rho = lon_rho(lat_min(1):lat_max(1), lon_min(1):lon_max(1));
    cut_lat_rho = lat_rho(lat_min(1):lat_max(1), lon_min(1):lon_max(1));

    data=squeeze(cut_data);
    data(data<-100)=NaN;
    data(data>100)=NaN;
    
    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperSize', [500, vert_paper_size_y]);
    set(gcf,'PaperPosition', [paper_position_hor, paper_position_ver, 500, vert_paper_size_y]) 
    
% % %     plot
    m_proj(m_proj_name,'lon',[section(1) section(2)],'lat',[section(3) section(4)]);
    hold on;
    m_pcolor(cut_lon_rho,cut_lat_rho,data);
    shading(gca,m_pcolor_shading_method);
    m_gshhs_i('color',m_gshhs_line_color)  
    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
                            
    titlename = strcat('OISST, surf (',char(calendarname(tempmonth)), ', ', num2str(tempyear),')');

    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
    % contour
    [C,h2]=m_contour(cut_lon_rho,cut_lat_rho, data, conlev, m_contour_color, 'linewidth', m_contour_linewidth);
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
     
    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

    jpgname=strcat(outfile, '_', testname, '_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg
    grid off;
    drawnow;
    saveas(gcf,jpgname,'jpg');

    disp(' ')
    disp([num2str(tempyear), '_', num2str(tempmonth), '_', var, ' plot is created.'])
    disp(' ')
    disp([' File path is : ',jpgname])
    disp(' ')

    hold off
    close all;
end
status = 1;
end