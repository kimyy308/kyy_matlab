function status=plot_seo_monthly_vertical_data(testname, outfile, filedir, section, tempyear, inputmonth, shadlev, conlev, var,  param_script)
% clear all;close all;
%==========================================================================
% % vert_param;

run(param_script);
%==========================================================================



% calendar{1} = ' Jan'; calendar{2} = ' Feb'; calendar{3} = ' Mar'; 
% calendar{4} = ' Apr'; calendar{5} = ' May'; calendar{6} = ' Jun';
% calendar{7} = ' Jul'; calendar{8} = ' Aug'; calendar{9} = ' Sep'; 
% calendar{10} = ' Oct'; calendar{11} = ' Nov'; calendar{12} = ' Dec';

for monthij=1:length(inputmonth)
    tempmonth = inputmonth(monthij);
    filename = strcat(filedir, num2str(tempyear,'%04i'), '\', ...
                testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');
    
    if (exist('lon_rho' , 'var') ~= 1)
        gd = read_grid(filename,Vtransform,Vstretching,theta_s,theta_b,hc,N);
        lon_rho  = gd.lon_rho;
        lat_rho  = gd.lat_rho; 
        mask_rho = gd.mask_rho;
        h = gd.h;
        N = gd.N;
        depth=gd.z_r;
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
    data = ncread(filename,varname,[lon_min lat_min 1 1], [lon_max-lon_min+1 lat_max-lat_min+1 data_info.Size(3) 1]);
    data = permute(data, [3 2 1]);
    cut_lon_rho = lon_rho(lat_min:lat_max, lon_min:lon_max);
    cut_lat_rho = lat_rho(lat_min:lat_max, lon_min:lon_max);
    cut_depth = depth(:,lat_min:lat_max, lon_min:lon_max);
    
    dist=sqrt((cut_lon_rho-section(1)).^2+(cut_lat_rho-section(3)).^2); %% get distance from station 1
    min_dist=min(min(dist));
    dist2=sqrt((cut_lon_rho-section(2)).^2+(cut_lat_rho-section(4)).^2);  %% get distance from station 2
    min_dist2=min(min(dist2));                
    [x1,y1]=find(dist==min_dist);  %% find closest point from station 1. [lat lon]
    [x2,y2]=find(dist2==min_dist2); %% find closest point from station 2. [lat lon]

    lat1=cut_lat_rho(x1(1),y1(1));  lon1=cut_lon_rho(x1(1),y1(1));
    lat2=cut_lat_rho(x2(1),y2(1));  lon2=cut_lon_rho(x2(1),y2(1));
    if (lon2-lon1) >= (lat2-lat1)
        lon_line = lon1:mean(gradient(cut_lon_rho(1,:))):lon2;  %% for 1/20^o horizontal resolution
        lat_line = (lon_line-lon1)/((lon2-lon1)/(lat2-lat1))+lat1; %% weight for latitude index (lat1 : 0.05/((lon2-lon1) * (lat2-lat1) : lat2 
        x=repmat(lon_line,gd.N,1);  %% copy lon_line (gd.N times) to make matrix 
        x_label='Longitude(^oE)';
%                 domaxis=[domaxis(3) domaxis(4) domaxis(5) domaxis(6) domaxis(5) domaxis(6)];
        Temp=zeros(gd.N,length(lon_line)); %% initialize temp matrix that size is same with x
    else
        lat_line=[min(lat1,lat2):mean(gradient(cut_lat_rho(:,1))):max(lat1,lat2)];
        lon_line=(lat_line-lat1)*((lon2-lon1)/(lat2-lat1))+lon1;
        x=repmat(lat_line,gd.N,1);
        x_label='Latitude(^oN)';
        Temp=zeros(gd.N,length(lat_line)); %% initialize temp matrix that size is same with x
    end

    if (x1(1)==x2(1) || y1(1)==y2(1)) %% fixed lon or lat
        Temp(:,:) = squeeze(data(:,min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1))));
        Yi(:,:)= squeeze(cut_depth(:,min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1))));
    else
        for k=1:1:gd.N
            lon_range=cut_lon_rho(min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1)));
            lat_range=cut_lat_rho(min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1)));
            data_range=squeeze(data(k,min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1))));  %%data(zlevel, latmin : latmax, lonmin : lonmax)
            depth_range=squeeze(cut_depth(k,min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1)))); 
            Temp(k,:)=griddata(lon_range,lat_range,data_range,lon_line,lat_line); %%get interpolated data section
            Yi(k,:)=griddata(lon_range,lat_range,depth_range,lon_line,lat_line); %%get interpolated depth section
        end
    end
            
    data=Temp;

    max(section(2)-section(1), section(4)-section(3));

%     set(gcf, 'PaperUnits', 'points');
%     set(gcf, 'PaperSize', [vert_paper_size_y*max(section(2)-section(1),section(4)-section(3))/2.0, vert_paper_size_y]);
%     set(gcf,'PaperPosition', [paper_position_hor, paper_position_ver, vert_paper_size_y*max(section(2)-section(1),section(4)-section(3))/2.0, vert_paper_size_y]) 
    
    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperSize', [500, vert_paper_size_y]);
    set(gcf,'PaperPosition', [paper_position_hor, paper_position_ver, 500, vert_paper_size_y]) 
    
%         text_posi_x=(section(4)-section(3))/20+section(3);
%         text_posi_y1=(section(6)-section(5))/20+section(5);
%         text_posi_y2=2*(section(6)-section(5))/20+section(5);
%         text_posi_y3=3*(section(6)-section(5))/20+section(5);

    hold on
    pcolor(x,Yi,data);

%     data2 = imresize(data, 100);
%     x2 = imresize(x, 100);
%     Yi2 = imresize(Yi, 100);
%     pcolor(x2,Yi2,data2)
    if (lon2-lon1) >= (lat2-lat1)        
        axis([section(1) section(2) section(5) section(6)]);
    else
        axis([section(3) section(4) section(5) section(6)]);
    end
    shading(gca,m_pcolor_shading_method);
    set(gca,'box', vert_grid_box, 'linewidth',vert_grid_box_linewidth, 'tickdir',vert_grid_tickdir_type,'fontsize',vert_grid_fontsize)
    xlabel(x_label,'FontSize',vert_grid_fontsize,'fontweight','bold')
    ylabel('Depth(m)','FontSize',vert_grid_fontsize,'fontweight','bold')

    titlename = strcat(varname, ' (',char(calendarname(tempmonth)), ', ', num2str(tempyear),')');
    title(titlename,'fontsize',vert_pcolor_title_fontsize); % stlee
    hold on
    [C,h]=contour(x,Yi,data,conlev,vert_contour_color,'linewidth',vert_contour_linewidth);              
    if strcmp(param_script,'fig_param_kyy_schematic')
%         no contour label
    else
        clabel(C,h,'FontSize',vert_contour_label_fontsize,'Color',vert_contour_label_color,'labelspacing',vert_contour_labelspacing,'Rotation',vert_contour_rotation,'fontweight', vert_contour_fontweight);
    end
    h = colorbar;
    colormap(colormap_style);
    set(h,'fontsize',colorbar_fontsize);
    title(h,colorbar_title,'fontsize',colorbar_title_fontsize);
    caxis(shadlev);
            
    jpgname=strcat(outfile, '_', testname, '_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg
    saveas(gcf,jpgname,'jpg');
    close all
end
status = 1;
end