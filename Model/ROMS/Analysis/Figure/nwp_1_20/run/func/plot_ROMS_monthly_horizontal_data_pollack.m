function status=plot_seo_ROMS_monthly_horizontal_data_pollack(testname, outfile, filedir, section, tempyear, inputmonth, shadlev, conlev, var,  param_script)
% clear all;close all;
% %  Updated 08-Jun-2018 by Yong-Yub Kim


lonlat = section(1:4);   
run(param_script);
tempmonth=inputmonth(1);
filename = strcat(filedir, ...
                testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');
Vstretching = ncread(filename,'Vstretching')';
Vtransform = ncread(filename,'Vtransform')';
theta_s = ncread(filename,'theta_s')';
theta_b = ncread(filename,'theta_b')';
s_rho = ncread(filename,'s_rho')';
N=length(s_rho);
hc = ncread(filename,'hc')';


for monthij=1:length(inputmonth)
    tempmonth = inputmonth(monthij);
    filename = strcat(filedir, ...
                testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');
    if (exist('lon_rho' , 'var') ~= 1)
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
    
    lon_west = abs(lon_rho - (section(1)-1));
    min_lon_west=min(lon_west(1,:));
    lon_east = abs(lon_rho - (section(2)+1));
    min_lon_east=min(lon_east(1,:));
    lat_south = abs(lat_rho - (section(3)-1));
    min_lat_south=min(lat_south(:,1));
    lat_north = abs(lat_rho - (section(4)+1));
    min_lat_north=min(lat_north(:,1));
    
    lon_min = find(lon_west(1,:) == min_lon_west);
    lon_max = find(lon_east(1,:) == min_lon_east);
    lat_min = find(lat_south(:,1) == min_lat_south);
    lat_max = find(lat_north(:,1) == min_lat_north);
    
    data_info = ncinfo(filename, varname); 
    if (length(data_info.Dimensions)==4)
        data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 data_info.Size(3) 1]);  %% cut horizontal area [x,y,z] (wider than target area)
    else
        data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 data_info.Size(3)]);  %% cut horizontal area [x,y,z] (wider than target area)
    end
    data=squeeze(data);
    
    cut_data = permute(data, [3 2 1]);  %% permute [x y z] -> [z y x]
    
    cut_lon_rho = lon_rho(lat_min(1):lat_max(1), lon_min(1):lon_max(1));
    cut_lat_rho = lat_rho(lat_min(1):lat_max(1), lon_min(1):lon_max(1));
    
    cut_mask_rho = mask_rho(lat_min(1):lat_max(1), lon_min(1):lon_max(1));
    cut_depth = depth(:,lat_min(1):lat_max(1), lon_min(1):lon_max(1));
    refdepth=(section(5)+section(6))/2.0   %% section(5) and section (6) should be equal, but prepare for wrong case
    vert_dist=abs(cut_depth-refdepth);
    
    cut_h = h(lat_min(1):lat_max(1), lon_min(1):lon_max(1));
    cut_depth = depth(:,lat_min(1):lat_max(1), lon_min(1):lon_max(1));
    
    for i=1:length(cut_data(1,:,1))
        for j=1:length(cut_data(1,1,:))
            if (cut_mask_rho(i,j)==0)
                cut_data3(i,j)=NaN;  %% land -> NaN
            else
                if (cut_depth(1,i,j)>refdepth)
                    cut_data3(i,j)=NaN;  %% If deepest water depth is shallower than refdepth -> NaN
                else
                    z_ind=find(squeeze(vert_dist(:,i,j))==min(squeeze(vert_dist(:,i,j))));  %% find closest depth from refdepth
                    if (z_ind==N)
                        cut_data3(i,j)=cut_data(z_ind,i,j); %% if closest cell is first vertical cell, get data from surface layer
                    elseif (z_ind==1)
                        cut_depth2(1:3)=cut_depth(z_ind:z_ind+2,i,j);  %% if closest cell is first vertical cell, get 3 cells from first to third cell 
                        cut_data2(1:3)=cut_data(z_ind:z_ind+2,i,j);
                        cut_data3(i,j)=interpn(cut_depth2,cut_data2,refdepth,'linear');  %% vertical interpolation of all horizontal columns
                    else
                        cut_depth2(1:3)=cut_depth(z_ind-1:z_ind+1,i,j); %% get 3 cells from near reference depth
                        cut_data2(1:3)=cut_data(z_ind-1:z_ind+1,i,j);
                        cut_data3(i,j)=interpn(cut_depth2,cut_data2,refdepth,'linear');  %% vertical interpolation of all horizontal columns
                    end
                end
            end
        end
    end
    clear data;
    cut_data3(find(cut_h>500))=NaN; %%pollacks don't lay below the 500m depth
    cut_data3(find(cut_data3>5))=NaN; %%pollacks don't lay hotter than 5 degree
    cut_data3(find(cut_data3<2))=NaN; %%pollacks don't lay colder than 2 degree
    data=cut_data3;
    
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
    titlename = strcat(var, ' (',char(calendarname(tempmonth)), ', ', num2str(tempyear),')');
    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
   
%     % contour
%     [C,h2]=m_contour(cut_lon_rho,cut_lat_rho, data, conlev, m_contour_color, 'linewidth', m_contour_linewidth);
%     clabel(C,h2,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
%         'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);

    % set colorbar 
    hc = colorbar;
%     colormap(colormap_style);
    colormap([1,0,0]);   %% colormap == red
    set(hc,'fontsize',colorbar_fontsize);
    title(hc,colorbar_title,'fontsize',colorbar_title_fontsize);
    caxis(shadlev);

    % set grid
    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);
     
    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

    jpgname=strcat(outfile, '_spawn_', testname, '_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg
    saveas(gcf,jpgname,'jpg');

    disp(' ')
    disp([num2str(tempyear), '_', num2str(tempmonth), '_', var, ' plot is created.'])
    disp(' ')
    disp([' File path is : ',jpgname])
    disp(' ')

    hold off
    close all;
    status = 1;
end