function status=plot_MOM_monthly_horizontal_data(testname, outfile, filedir, section, tempyear, inputmonth, shadlev, conlev, var,  param_script)
% clear all;close all;
%==========================================================================
% % Updated 05-Apr-2018 by Yong-Yub Kim

% % % vert_param;

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
        lon  = ncread(filename,'xt_ocean');
        lat  = ncread(filename,'yt_ocean');
        depth = ncread(filename, 'st_ocean');
    
    
    lon_west = abs(lon - (section(1)-2));
    min_lon_west=min(lon_west);
    lon_east = abs(lon - (section(2)+2));
    min_lon_east=min(lon_east);
    lat_south = abs(lat - (section(3)-2));
    min_lat_south=min(lat_south);
    lat_north = abs(lat - (section(4)+2));
    min_lat_north=min(lat_north);

    lon_min = find(lon_west == min_lon_west);
    lon_max = find(lon_east == min_lon_east);
    lat_min = find(lat_south == min_lat_south);
    lat_max = find(lat_north == min_lat_north);
    
    data_info = ncinfo(filename, varname); 
    data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 data_info.Size(3) 1]);  %% cut horizontal area [x,y,z] (wider than target area)
    data=squeeze(data);
    cut_data = permute(data, [3 2 1]);  %% permute [x y z] -> [z y x]
    [lon_rho lat_rho]=meshgrid(lon,lat);
    cut_lon_rho = lon_rho(lat_min(1):lat_max(1), lon_min(1):lon_max(1));
    cut_lat_rho = lat_rho(lat_min(1):lat_max(1), lon_min(1):lon_max(1));
    cut_depth = -repmat(depth,1,lat_max(1)-lat_min(1)+1,lon_max(1)-lon_min(1)+1);
% %     cut_depth = depth(:,lat_min(1):lat_max(1), lon_min(1):lon_max(1));
    refdepth=(section(5)+section(6))/2.0   %% section(5) and section (6) should be equal, but prepare for wrong case
    vert_dist=abs(cut_depth-refdepth);  %% Vertical distance from reference depth
    end
%     dist=sqrt((cut_lon_rho-section(1)).^2+(cut_lat_rho-section(3)).^2); %% get distance from station 1
%     min_dist=min(min(dist));
%     dist2=sqrt((cut_lon_rho-section(2)).^2+(cut_lat_rho-section(4)).^2);  %% get distance from station 2
%     min_dist2=min(min(dist2));                
%     [x1,y1]=find(dist==min_dist);  %% find closest point from station 1. [lat lon]
%     [x2,y2]=find(dist2==min_dist2); %% find closest point from station 2. [lat lon]
    
    cut_data(find(cut_data<-100))=NaN;
    cut_data(find(cut_data>100))=NaN;
    if (exist('z_ind' , 'var') ~= 1)
        z_ind=find(squeeze(vert_dist(:))==min(squeeze(vert_dist(:))));  %% find closest depth from refdepth
    end
    for i=1:length(cut_data(1,:,1))
        for j=1:length(cut_data(1,1,:))
            if (z_ind==1)
                cut_depth2(1:3)=cut_depth(z_ind:z_ind+2,i,j);  %% if closest cell is first vertical cell, get 3 cells from first to third cell 
                cut_data2(1:3)=cut_data(z_ind:z_ind+2,i,j);
            else
                cut_depth2(1:3)=cut_depth(z_ind-1:z_ind+1,i,j); %% get 3 cells from near reference depth
                cut_data2(1:3)=cut_data(z_ind-1:z_ind+1,i,j);
            end
            cut_data3(i,j)=interpn(cut_depth2,cut_data2,refdepth,'linear');  %% vertical interpolation of all horizontal columns
        end
    end
    clear data;
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