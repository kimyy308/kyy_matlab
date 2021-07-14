function status=plot_seo_ROMS_monthly_hor_t_uv2(testname, outfile, filedir, section, tempyear, tempyear2,inputmonth, shadlev, conlev, var,  param_script)
% clear all;close all;
%==========================================================================
% % % vert_param;
% %  Updated 09-Apr-2018 by Yong-Yub Kim
% %  Updated 22-Aug-2018 by Yong-Yub Kim 

lonlat = section(1:4);
refdepth=(section(5)+section(6))/2.0   %% section(5) and section (6) should be equal, but prepare for wrong case

Vstretching=1;
Vtransform=1;
hc=5;
theta_s=5;
theta_b=0.4;
N=20;
grdname = strcat(filedir,'roms_grid_final.nc');
run(param_script);
%==========================================================================

for monthij=1:length(inputmonth)
    tempmonth = inputmonth(monthij);
    filename = strcat(filedir, testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempyear2,'%04i'),'_',num2str(tempmonth,'%02i'), '.nc');
    
    if (exist('lon_rho' , 'var') ~= 1)
        gd = read_grid(grdname,Vtransform,Vstretching,theta_s,theta_b,hc,N);
        lon_rho  = gd.lon_rho;
        lat_rho  = gd.lat_rho; 
        lon_u  = gd.lon_u;
        lat_u  = gd.lat_u; 
        lon_v  = gd.lon_v;
        lat_v  = gd.lat_v; 
        mask_rho = gd.mask_rho;
        h = gd.h;
        N = gd.N;
        depth=gd.z_r;
        [Nlat, Nlon, Nz]=size(lon_rho);
        dl=(max(max(lon_rho))-min(min(lon_rho)))/Nlon;
        depth_u=0.5*(depth(:,:,1:Nlon-1)+depth(:,:,2:Nlon));
        depth_v=0.5*(depth(:,1:Nlat-1,:)+depth(:,2:Nlat,:));
        
        [lon_min, lon_max, lat_min, lat_max] = findind_Y(dl, section(1:4), lon_rho, lat_rho);
        [lon_u_min, lon_u_max, lat_u_min, lat_u_max] = findind_Y(dl, section(1:4), lon_u, lat_u);
        [lon_v_min, lon_v_max, lat_v_min, lat_v_max] = findind_Y(dl, section(1:4), lon_v, lat_v);
    end
    

    data_info = ncinfo(filename, varname); 
    if (length(data_info.Dimensions)==4)
        data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 data_info.Size(3) 1]);  %% cut horizontal area [x,y,z] (wider than target area)
        u = ncread(filename,'u',[lon_u_min(1) lat_u_min(1) 1 1], [lon_u_max(1)-lon_u_min(1)+1 lat_u_max(1)-lat_u_min(1)+1 data_info.Size(3) 1]);  %% cut horizontal area [x,y,z] (wider than target area)
        v = ncread(filename,'v',[lon_v_min(1) lat_v_min(1) 1 1], [lon_v_max(1)-lon_v_min(1)+1 lat_v_max(1)-lat_v_min(1)+1 data_info.Size(3) 1]);  %% cut horizontal area [x,y,z] (wider than target area)
    else
        data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 data_info.Size(3)]);  %% cut horizontal area [x,y,z] (wider than target area)
        u = ncread(filename,'u',[lon_u_min(1) lat_u_min(1) 1], [lon_u_max(1)-lon_u_min(1)+1 lat_u_max(1)-lat_u_min(1)+1 data_info.Size(3)]);  %% cut horizontal area [x,y,z] (wider than target area)
        v = ncread(filename,'v',[lon_v_min(1) lat_v_min(1) 1], [lon_v_max(1)-lon_v_min(1)+1 lat_v_max(1)-lat_v_min(1)+1 data_info.Size(3)]);  %% cut horizontal area [x,y,z] (wider than target area)
    end
    data=squeeze(data);
    u=squeeze(u);
    v=squeeze(v);

    cut_data = permute(data, [3 2 1]);  %% permute [x y z] -> [z y x]
    cut_u = permute(u, [3 2 1]);  %% permute [x y z] -> [z y x]
    cut_v = permute(v, [3 2 1]);  %% permute [x y z] -> [z y x]
    cut_depth = depth(:,lat_min(1):lat_max(1), lon_min(1):lon_max(1));
    cut_depth_u = depth(:,lat_u_min(1):lat_u_max(1), lon_u_min(1):lon_u_max(1));
    cut_depth_v = depth(:,lat_v_min(1):lat_v_max(1), lon_v_min(1):lon_v_max(1));
    
    cut_lon_rho = lon_rho(lat_min(1):lat_max(1), lon_min(1):lon_max(1));
    cut_lat_rho = lat_rho(lat_min(1):lat_max(1), lon_min(1):lon_max(1));
    cut_lon_u = lon_u(lat_u_min(1):lat_u_max(1), lon_u_min(1):lon_u_max(1));
    cut_lat_u = lat_u(lat_u_min(1):lat_u_max(1), lon_u_min(1):lon_u_max(1));
    cut_lon_v = lon_v(lat_v_min(1):lat_v_max(1), lon_v_min(1):lon_v_max(1));
    cut_lat_v = lat_v(lat_v_min(1):lat_v_max(1), lon_v_min(1):lon_v_max(1));
    
    cut_mask_rho = mask_rho(lat_min(1):lat_max(1), lon_min(1):lon_max(1));
    
    cut_data2 = roms_vinterp(cut_data,cut_depth,refdepth);
    cut_u2 = roms_vinterp(cut_u,cut_depth_u,refdepth);
    cut_v2 = roms_vinterp(cut_v,cut_depth_v,refdepth);
    
    cut_u3=griddata(double(cut_lon_u), double(cut_lat_u), double(cut_u2(:,:)), double(cut_lon_rho), double(cut_lat_rho));
    cut_v3=griddata(double(cut_lon_v), double(cut_lat_v), double(cut_v2(:,:)), double(cut_lon_rho), double(cut_lat_rho));
    u_rho=cut_u3';
    v_rho=cut_v3';
    u_rho(u_rho<-100)=NaN;
    u_rho(u_rho>100)=NaN;
    v_rho(v_rho<-100)=NaN;
    u_rho(v_rho>100)=NaN;
    data=cut_data2;
    data(data<-100)=NaN;
    data(data>100)=NaN;
    if (exist('ref_vec_x_range' , 'var') ~= 1)
        ref_vec_x_ind = find(abs(cut_lon_rho(1,:)-m_quiver_ref_text_x_location) == min(abs(cut_lon_rho(1,:)-m_quiver_ref_text_x_location)));
        ref_vec_y_ind = find(abs(cut_lat_rho(:,1)-m_quiver_ref_text_y_location) == min(abs(cut_lat_rho(:,1)-m_quiver_ref_text_y_location)))+m_quiver_y_interval*2;
        ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval/2;
        ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2))+m_quiver_y_interval/2;
    end
    u_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_u_value;
    v_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_v_value; 
    
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
    
    uvplot=m_quiver(cut_lon_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
                        cut_lat_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
                        u_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
                        v_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
                        'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth);
                        
if refdepth==0
    titlename = strcat('T, surf (',char(calendarname(tempmonth)), ', ', num2str(tempyear),'-',num2str(tempyear2),')');
else  
    titlename = strcat('T, ', num2str(refdepth),'m (',char(calendarname(tempmonth)), ', ', num2str(tempyear),'-',num2str(tempyear2),')');
end
    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
    m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_text, 'FontSize', m_quiver_ref_text_fontsize); 

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