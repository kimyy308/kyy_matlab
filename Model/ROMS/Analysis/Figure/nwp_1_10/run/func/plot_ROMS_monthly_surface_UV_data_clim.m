function status=plot_ROMS_monthly_surface_UV_data(testname, outfile, outputdir, section, inputyear, inputmonth, shadlev, conlev, var, param_script, shadi, uvi)
   % %  Updated 08-Jun-2018 by Yong-Yub Kim

lonlat = section(1:4);   
run(param_script);
% tempmonth=inputmonth(1);
% filename = strcat(filedir, ...
%                 testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');

for monthij = 1:length(inputmonth)
    tempmonth = inputmonth(monthij);
    clear data u v
    for yearij=1:length(inputyear)
        tempyear = inputyear(yearij);
        filedir= [outputdir,'\', num2str(tempyear,'%04i'),'\'];

        % ex : rootdir\test37\data\2001\test37_monthly_2001_01.nc
        filename = strcat(filedir, ...
                testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');
        if (exist('lon_rho' , 'var') ~= 1)
            Vstretching = ncread(filename,'Vstretching')';
            Vtransform = ncread(filename,'Vtransform')';
            theta_s = ncread(filename,'theta_s')';
            theta_b = ncread(filename,'theta_b')';
            s_rho = ncread(filename,'s_rho')';
            N=length(s_rho);
            hc = ncread(filename,'hc')';    

            gd = read_grid(filename,Vtransform,Vstretching,theta_s,theta_b,hc,N);
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

            lon_west = abs(lon_rho - (section(1)-1));
            min_lon_west=min(lon_west(1,:));
            lon_east = abs(lon_rho - (section(2)+1));
            min_lon_east=min(lon_east(1,:));
            lat_south = abs(lat_rho - (section(3)-1));
            min_lat_south=min(lat_south(:,1));
            lat_north = abs(lat_rho - (section(4)+1));
            min_lat_north=min(lat_north(:,1));

            lon_u_west = abs(lon_u - (section(1)-1));
            min_lon_u_west=min(lon_u_west(1,:));
            lon_u_east = abs(lon_u - (section(2)+1));
            min_lon_u_east=min(lon_u_east(1,:));
            lat_u_south = abs(lat_u - (section(3)-1));
            min_lat_u_south=min(lat_u_south(:,1));
            lat_u_north = abs(lat_u - (section(4)+1));
            min_lat_u_north=min(lat_u_north(:,1));

            lon_v_west = abs(lon_v - (section(1)-1));
            min_lon_v_west=min(lon_v_west(1,:));
            lon_v_east = abs(lon_v - (section(2)+1));
            min_lon_v_east=min(lon_v_east(1,:));
            lat_v_south = abs(lat_v - (section(3)-1));
            min_lat_v_south=min(lat_v_south(:,1));
            lat_v_north = abs(lat_v - (section(4)+1));
            min_lat_v_north=min(lat_v_north(:,1));


            lon_min = find(lon_west(1,:) == min_lon_west);
            lon_max = find(lon_east(1,:) == min_lon_east);
            lat_min = find(lat_south(:,1) == min_lat_south);
            lat_max = find(lat_north(:,1) == min_lat_north);

            lon_u_min = find(lon_u_west(1,:) == min_lon_u_west);
            lon_u_max = find(lon_u_east(1,:) == min_lon_u_east);
            lat_u_min = find(lat_u_south(:,1) == min_lat_u_south);
            lat_u_max = find(lat_u_north(:,1) == min_lat_u_north);

            lon_v_min = find(lon_v_west(1,:) == min_lon_v_west);
            lon_v_max = find(lon_v_east(1,:) == min_lon_v_east);
            lat_v_min = find(lat_v_south(:,1) == min_lat_v_south);
            lat_v_max = find(lat_v_north(:,1) == min_lat_v_north);
        end
        data_info = ncinfo(filename, varname); 
        if shadi==1
            tempdata = ncread(filename,varname,[lon_min(1) lat_min(1) data_info.Size(3) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
            if (exist('data')==0)
%                 data=zeros([size(tempdata),12]);
                data=zeros(size(tempdata));
            end
            data=data + tempdata / (length(inputyear));
        end
        if uvi==1
            tempu = ncread(filename,'u',[lon_u_min(1) lat_u_min(1) data_info.Size(3) 1], [lon_u_max(1)-lon_u_min(1)+1 lat_u_max(1)-lat_u_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
            tempv = ncread(filename,'v',[lon_v_min(1) lat_v_min(1) data_info.Size(3) 1], [lon_v_max(1)-lon_v_min(1)+1 lat_v_max(1)-lat_v_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
            if (exist('u')==0)
                u=zeros(size(tempu));
            end
            if (exist('v')==0)
                v=zeros(size(tempv));
            end
            u = u + tempu / (length(inputyear));
            v = v + tempv / (length(inputyear));
        end
    end

    
%     if (length(data_info.Dimensions)==4)
%         data = ncread(filename,varname,[lon_min(1) lat_min(1) data_info.Size(3) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
%         u = ncread(filename,'u',[lon_u_min(1) lat_u_min(1) data_info.Size(3) 1], [lon_u_max(1)-lon_u_min(1)+1 lat_u_max(1)-lat_u_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
%         v = ncread(filename,'v',[lon_v_min(1) lat_v_min(1) data_info.Size(3) 1], [lon_v_max(1)-lon_v_min(1)+1 lat_v_max(1)-lat_v_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
%     else
%         data = ncread(filename,varname,[lon_min(1) lat_min(1) data_info.Size(3)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
%         u = ncread(filename,'u',[lon_u_min(1) lat_u_min(1) data_info.Size(3)], [lon_u_max(1)-lon_u_min(1)+1 lat_u_max(1)-lat_u_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
%         v = ncread(filename,'v',[lon_v_min(1) lat_v_min(1) data_info.Size(3)], [lon_v_max(1)-lon_v_min(1)+1 lat_v_max(1)-lat_v_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
%     end
%     data=squeeze(data);
%     u=squeeze(u);
%     v=squeeze(v);

    if shadi==1
        cut_data = permute(data, [3 2 1]);  %% permute [x y z] -> [z y x]
        cut_data = squeeze(cut_data);
    end
    cut_lon_rho = lon_rho(lat_min(1):lat_max(1), lon_min(1):lon_max(1));
    cut_lat_rho = lat_rho(lat_min(1):lat_max(1), lon_min(1):lon_max(1));
    cut_lon_u = lon_u(lat_u_min(1):lat_u_max(1), lon_u_min(1):lon_u_max(1));
    cut_lat_u = lat_u(lat_u_min(1):lat_u_max(1), lon_u_min(1):lon_u_max(1));
    cut_lon_v = lon_v(lat_v_min(1):lat_v_max(1), lon_v_min(1):lon_v_max(1));
    cut_lat_v = lat_v(lat_v_min(1):lat_v_max(1), lon_v_min(1):lon_v_max(1));
    
    cut_mask_rho = mask_rho(lat_min(1):lat_max(1), lon_min(1):lon_max(1));
    cut_depth = depth(:,lat_min(1):lat_max(1), lon_min(1):lon_max(1));
%     refdepth=(section(5)+section(6))/2.0   %% section(5) and section (6) should be equal, but prepare for wrong case
%     vert_dist=abs(cut_depth-refdepth);
    if uvi==1
        clear cut_u cut_v;
        for i=1:1
            cut_u(:,:,i)=griddata(double(cut_lon_u), double(cut_lat_u), double(u(:,:,i))', double(cut_lon_rho), double(cut_lat_rho));
            cut_v(:,:,i)=griddata(double(cut_lon_v), double(cut_lat_v), double(v(:,:,i))', double(cut_lon_rho), double(cut_lat_rho));
        end

        cut_u= permute(cut_u, [3 1 2]);  %% permute [y x z] -> [z y x]
        cut_v = permute(cut_v, [3 1 2]);  %% permute [y x z] -> [z y x]
        u_rho=squeeze(cut_u);
        v_rho=squeeze(cut_v);
        
        if (exist('ref_vec_x_range' , 'var') ~= 1)
            ref_vec_x_ind = find(abs(cut_lon_rho(1,:)-m_quiver_ref_text_x_location) == min(abs(cut_lon_rho(1,:)-m_quiver_ref_text_x_location)));
            ref_vec_y_ind = find(abs(cut_lat_rho(:,1)-m_quiver_ref_text_y_location) == min(abs(cut_lat_rho(:,1)-m_quiver_ref_text_y_location)))+m_quiver_y_interval;
            ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2));
            ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2));
        end
        u_rho(ref_vec_y_range,ref_vec_x_range)=m_quiver_ref_u_value;
        v_rho(ref_vec_y_range,ref_vec_x_range)=m_quiver_ref_v_value;    
    end
    

    
    
        
            
% %         % read data
% %         if (exist('lon_u')==0)
% %             lon_u = ncread(filename,'lon_u');
% %             lat_u = ncread(filename,'lat_u');
% %             lon_v = ncread(filename,'lon_v');
% %             lat_v = ncread(filename,'lat_v');
% %             lon = ncread(filename,'lon_rho');
% %             lat = ncread(filename,'lat_rho');
% %         end
% %         data_info_u = ncinfo(filename, 'u');
% %         data_info_v = ncinfo(filename, 'v');
% %         u = ncread(filename,'u',[1 1 data_info_u.Size(3) 1], [data_info_u.Size(1) data_info_u.Size(2) 1 1]);
% %         v = ncread(filename,'v',[1 1 data_info_v.Size(3) 1], [data_info_v.Size(1) data_info_v.Size(2) 1 1]);
% %         u_rho=griddata(lon_u, lat_u, u, lon, lat);
% %         v_rho=griddata(lon_v, lat_v, v, lon, lat);
% %         
% %         varname = 'temp';
% %         data_info = ncinfo(filename, varname);
% %         data = ncread(filename,varname,[1 1 data_info.Size(3) 1], [data_info.Size(1) data_info.Size(2) 1 1]);
% % 
% %         
% %         
% %         % Reference vector value setting
% %         u_rho(m_quiver_ref_vec_x_range,m_quiver_ref_vec_y_range) = m_quiver_ref_u_value;
% %         v_rho(m_quiver_ref_vec_x_range,m_quiver_ref_vec_y_range) = m_quiver_ref_v_value; 
        

    % % %     plot
    
    m_proj(m_proj_name,'lon',[section(1) section(2)],'lat',[section(3) section(4)]);
    hold on;
    if shadi==1
        m_pcolor(cut_lon_rho,cut_lat_rho,cut_data);
    end
    shading(gca,m_pcolor_shading_method);
    m_gshhs_i('color',m_gshhs_line_color)  
    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    
    if uvi==1
    uvplot=m_quiver(cut_lon_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
                        cut_lat_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
                        u_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end) * m_quiver_vector_size, ...
                        v_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end) * m_quiver_vector_size, ...
                        'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth);
    end           

%     titlename = strcat('temp- ', num2str(0),'m (',char(calendarname(tempmonth)), ', ', num2str(tempyear),')');
    titlename = strcat('surf (',char(calendarname(tempmonth)), ')');

    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
    if uvi==1
        m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_text, 'FontSize', m_quiver_ref_text_fontsize); 
    end
    % contour
    if shadi==1
        [C,h2]=m_contour(cut_lon_rho,cut_lat_rho, cut_data, conlev, m_contour_color, 'linewidth', m_contour_linewidth);
        clabel(C,h2,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
            'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);
    
        % set colorbar 
        h = colorbar;
        colormap(colormap_style);
        set(h,'fontsize',colorbar_fontsize);
        title(h,colorbar_title,'fontsize',colorbar_title_fontsize);
        caxis(shadlev);
    end
    % set grid
    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
     
    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

%     jpgname=strcat(outfile, '_', testname, '_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg
    jpgname=strcat(outfile, '_', testname, '_', 'CLIM', '_', num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg

%     pause(2); %% prevent too short time between previous command and save command
    drawnow; %% prevent too short time between previous command and save command
    saveas(gcf,jpgname,'jpg');

    disp(' ')
    disp([num2str(tempyear), '_', num2str(tempmonth), '_', var, ' plot is created.'])
    disp(' ')
    disp([' File path is : ',jpgname])
    disp(' ')

    hold off
    close all;



%         % plot
%         m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
%         hold on;
%         m_pcolor(lon,lat,data);
%         shading(gca,m_pcolor_shading_method);
%         hold on;
%         uvplot=m_quiver(lon(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
%                         lat(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
%                         u_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
%                         v_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
%                         m_quiver_vector_size,'color',m_quiver_vector_color, 'LineWidth', m_quiver_LineWidth);
%         m_gshhs_i('color',m_gshhs_line_color)  
%         m_gshhs_i('patch',m_gshhs_land_color);
%         titlename = strcat('Surface current', ' (',char(calendarname(tempmonth)), ', ', num2str(tempyear),')');
%         title(titlename,'fontsize', m_quiver_title_fontsize);
%         m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_text, 'FontSize', m_quiver_ref_text_fontsize); 
%         
%         % set colorbar 
%         h = colorbar;
%         colormap(colormap_style);
%         set(h,'fontsize',colorbar_fontsize);
%         title(h,colorbar_title,'fontsize',colorbar_title_fontsize);
%         caxis(shadlev);
%         
%         % set grid
%         m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
%         
%         
%         set(gcf, 'PaperUnits', 'points');
%         set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%         set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
%         
%         jpgname=strcat(outfile, '_', testname, '_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg
%         saveas(gcf,jpgname,'jpg');
% %         fig=gcf;
% %         print(jpgname,'-djpeg');
%         disp(' ')
%         disp([num2str(tempyear), '_', num2str(tempmonth), '_', 'UV', ' plot is created.'])
%         disp(' ')
%         disp([' File path is : ',jpgname])
%         disp(' ')
%         
%         hold off
%         close all;

    status=1;
end

return