function status=plot_MOM_monthly_horizontal_t_uv(testname, outfile, filedir, section, tempyear, inputmonth, shadlev, conlev, var,  param_script)
% clear all;close all;
%==========================================================================
% % % vert_param;
% % % Updated 06-Apr-2018 by Yong-Yub Kim.
% % % Updated 09-Apr-2018 by Yong-Yub Kim.

lonlat = section(1:4);
run(param_script);

%==========================================================================



% calendar{1} = ' Jan'; calendar{2} = ' Feb'; calendar{3} = ' Mar'; 
% calendar{4} = ' Apr'; calendar{5} = ' May'; calendar{6} = ' Jun';
% calendar{7} = ' Jul'; calendar{8} = ' Aug'; calendar{9} = ' Sep'; 
% calendar{10} = ' Oct'; calendar{11} = ' Nov'; calendar{12} = ' Dec';

for monthij=1:length(inputmonth)
    tempmonth = inputmonth(monthij);
    filename = strcat(filedir, 'ocean_snap_', num2str(tempyear,'%04i'), num2str(tempmonth,'%02i'), '.nc');    
    if (exist('lon' , 'var') ~= 1)
        lon  = ncread(filename,'xt_ocean');;
        lat  = ncread(filename,'yt_ocean'); 
        lon_u  = ncread(filename,'xu_ocean');
        lat_u  = ncread(filename,'yu_ocean'); 
        lon_v  = ncread(filename,'xu_ocean');
        lat_v  = ncread(filename,'yu_ocean'); 
%         h = gd.h;
        depth=ncread(filename, 'st_ocean');
        N = length(depth);
        
        
        lon_west = abs(lon - (section(1)-1));
        min_lon_west=min(lon_west(:));
        lon_east = abs(lon - (section(2)+1));
        min_lon_east=min(lon_east(:));
        lat_south = abs(lat - (section(3)-1));
        min_lat_south=min(lat_south(:));
        lat_north = abs(lat - (section(4)+1));
        min_lat_north=min(lat_north(:));

        lon_u_west = abs(lon_u - (section(1)-1));
        min_lon_u_west=min(lon_u_west(:));
        lon_u_east = abs(lon_u - (section(2)+1));
        min_lon_u_east=min(lon_u_east(:));
        lat_u_south = abs(lat_u - (section(3)-1));
        min_lat_u_south=min(lat_u_south(:));
        lat_u_north = abs(lat_u - (section(4)+1));
        min_lat_u_north=min(lat_u_north(:));

        lon_v_west = abs(lon_v - (section(1)-1));
        min_lon_v_west=min(lon_v_west(:));
        lon_v_east = abs(lon_v - (section(2)+1));
        min_lon_v_east=min(lon_v_east(:));
        lat_v_south = abs(lat_v - (section(3)-1));
        min_lat_v_south=min(lat_v_south(:));
        lat_v_north = abs(lat_v - (section(4)+1));
        min_lat_v_north=min(lat_v_north(:));


        lon_min = find(lon_west(:) == min_lon_west);
        lon_max = find(lon_east(:) == min_lon_east);
        lat_min = find(lat_south(:) == min_lat_south);
        lat_max = find(lat_north(:) == min_lat_north);

        lon_u_min = find(lon_u_west(:) == min_lon_u_west);
        lon_u_max = find(lon_u_east(:) == min_lon_u_east);
        lat_u_min = find(lat_u_south(:) == min_lat_u_south);
        lat_u_max = find(lat_u_north(:) == min_lat_u_north);

        lon_v_min = find(lon_v_west(:) == min_lon_v_west);
        lon_v_max = find(lon_v_east(:) == min_lon_v_east);
        lat_v_min = find(lat_v_south(:) == min_lat_v_south);
        lat_v_max = find(lat_v_north(:) == min_lat_v_north);
    
        [lon_rho lat_rho]=meshgrid(lon,lat);
        [lon_u lat_u]=meshgrid(lon_u,lat_u);
        [lon_v lat_v]=meshgrid(lon_v,lat_v);
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
    
    cut_lon_rho = lon_rho(lat_min(1):lat_max(1), lon_min(1):lon_max(1));
    cut_lat_rho = lat_rho(lat_min(1):lat_max(1), lon_min(1):lon_max(1));
    cut_lon_u = lon_u(lat_u_min(1):lat_u_max(1), lon_u_min(1):lon_u_max(1));
    cut_lat_u = lat_u(lat_u_min(1):lat_u_max(1), lon_u_min(1):lon_u_max(1));
    cut_lon_v = lon_v(lat_v_min(1):lat_v_max(1), lon_v_min(1):lon_v_max(1));
    cut_lat_v = lat_v(lat_v_min(1):lat_v_max(1), lon_v_min(1):lon_v_max(1));
    cut_depth = -repmat(depth,1,lat_max(1)-lat_min(1)+1,lon_max(1)-lon_min(1)+1);

    refdepth=(section(5)+section(6))/2.0   %% section(5) and section (6) should be equal, but prepare for wrong case
    vert_dist=abs(-depth-refdepth);
    for i=1:N
        cut_u(:,:,i)=griddata(double(cut_lon_u), double(cut_lat_u), double(u(:,:,i))', double(cut_lon_rho), double(cut_lat_rho));
        cut_v(:,:,i)=griddata(double(cut_lon_v), double(cut_lat_v), double(v(:,:,i))', double(cut_lon_rho), double(cut_lat_rho));
    end
    cut_u = permute(cut_u, [3 1 2]);  %% permute [y x z] -> [z y x]
    cut_v = permute(cut_v, [3 1 2]);  %% permute [y x z] -> [z y x]
    cut_data(find(cut_data<-100))=NaN;
    cut_data(find(cut_data>100))=NaN;
    cut_u(find(cut_u<-100))=NaN;
    cut_u(find(cut_u>100))=NaN;
    cut_v(find(cut_v<-100))=NaN;
    cut_v(find(cut_v>100))=NaN;

    if (exist('z_ind' , 'var') ~= 1)
        z_ind=find(squeeze(vert_dist(:))==min(squeeze(vert_dist(:))));  %% find closest depth from refdepth
    end
    
    for i=1:length(cut_data(1,:,1))
        for j=1:length(cut_data(1,1,:))
            if (z_ind==1)
%                 cut_depth2(1:3)=cut_depth(z_ind:z_ind+2,i,j);  %% if closest cell is first vertical cell, get 3 cells from first to third cell 
%                 cut_data2(1:3)=cut_data(z_ind:z_ind+2,i,j);
%                 cut_u2(1:3)=cut_u(z_ind:z_ind+2,i,j);
%                 cut_v2(1:3)=cut_v(z_ind:z_ind+2,i,j);
                cut_data3(i,j)=cut_data(z_ind,i,j); %% if closest cell is first vertical cell, get data from surface layer
                cut_u3(i,j)=cut_u(z_ind,i,j);
                cut_v3(i,j)=cut_v(z_ind,i,j);
            elseif (z_ind==N)
                cut_depth2(1:3)=cut_depth(z_ind-2:z_ind,i,j);  %% if closest cell is bottom vertical cell, get 3 cells from bottom-2 to bottom cell 
                cut_data2(1:3)=cut_data(z_ind-2:z_ind,i,j);
                cut_u2(1:3)=cut_u(z_ind-2:z_ind,i,j);
                cut_v2(1:3)=cut_v(z_ind-2:z_ind,i,j);
                cut_data3(i,j)=interpn(cut_depth2,cut_data2,refdepth,'linear');  %% vertical interpolation of all horizontal columns
                cut_u3(i,j)=interpn(cut_depth2,cut_u2,refdepth,'linear');  %% vertical interpolation of all horizontal columns
                cut_v3(i,j)=interpn(cut_depth2,cut_v2,refdepth,'linear');  %% vertical interpolation of all horizontal columns
            else
                cut_depth2(1:3)=cut_depth(z_ind-1:z_ind+1,i,j); %% get 3 cells from near reference depth
                cut_data2(1:3)=cut_data(z_ind-1:z_ind+1,i,j);
                cut_u2(1:3)=cut_u(z_ind-1:z_ind+1,i,j);
                cut_v2(1:3)=cut_v(z_ind-1:z_ind+1,i,j);
                cut_data3(i,j)=interpn(cut_depth2,cut_data2,refdepth,'linear');  %% vertical interpolation of all horizontal columns
                cut_u3(i,j)=interpn(cut_depth2,cut_u2,refdepth,'linear');  %% vertical interpolation of all horizontal columns
                cut_v3(i,j)=interpn(cut_depth2,cut_v2,refdepth,'linear');  %% vertical interpolation of all horizontal columns
            end
            
        end
    end
    clear data u v;
    data=cut_data3;
    u_rho=cut_u3;
    v_rho=cut_v3;
    u_rho(m_quiver_ref_vec_y_range,m_quiver_ref_vec_x_range)=m_quiver_ref_u_value;
    v_rho(m_quiver_ref_vec_y_range,m_quiver_ref_vec_x_range)=m_quiver_ref_v_value;
    
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
                        u_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end) * m_quiver_vector_size, ...
                        v_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end) * m_quiver_vector_size, ...
                        'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth);
    

    titlename = strcat('T ', num2str(refdepth),'m (',char(calendarname(tempmonth)), ', ', num2str(tempyear),')');
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