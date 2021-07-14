function status=plot_vertical_var(testname, outfile, filedir, section, tempyear, inputmonth, shadlev, conlev, var, param_script)
% clear all;close all;
%==========================================================================
% % This function needs 
% % parameter : vert_param.m (Vtransform, Vstretching, theta_s, theta_b, hc(Tcline), N(the number of vertical level);
% % function : 'read_grid.m(similar to grd.m)', 'stretching.m', 'zlevs.m'.
% % library : 'netcdf_old'

% % Updated 27-Apr-2018 Yong-Yub Kim
% %  Updated 09-May-2018 by Yong-Yub Kim
%==========================================================================

lonlat = section(1:4);
run(param_script);
%==========================================================================

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
    
    clear cut_data cut_mask_rho cut_depth
    cut_data = permute(data, [3 2 1]);  %% permute [x y z] -> [z y x]
    
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
clear cut_u cut_v
    for i=1:N
        cut_u(:,:,i)=griddata(double(cut_lon_u), double(cut_lat_u), double(u(:,:,i))', double(cut_lon_rho), double(cut_lat_rho));
        cut_v(:,:,i)=griddata(double(cut_lon_v), double(cut_lat_v), double(v(:,:,i))', double(cut_lon_rho), double(cut_lat_rho));
    end
    cut_u= permute(cut_u, [3 1 2]);  %% permute [y x z] -> [z y x]
    cut_v = permute(cut_v, [3 1 2]);  %% permute [y x z] -> [z y x]
    switch(var)
        case('vert_u')
        clear cut_data
        cut_data=cut_u;
        case('vert_v')
        clear cut_data
        cut_data=cut_v;
    end
    
    dist=sqrt((cut_lon_rho-section(1)).^2+(cut_lat_rho-section(3)).^2); %% get distance from station 1
    min_dist=min(min(dist));
    dist2=sqrt((cut_lon_rho-section(2)).^2+(cut_lat_rho-section(4)).^2);  %% get distance from station 2
    min_dist2=min(min(dist2));                
    [x1,y1]=find(dist==min_dist);  %% find closest point from station 1. [lat lon]
    [x2,y2]=find(dist2==min_dist2); %% find closest point from station 2. [lat lon]

    lat1=cut_lat_rho(x1(1),y1(1));  lon1=cut_lon_rho(x1(1),y1(1));
    lat2=cut_lat_rho(x2(1),y2(1));  lon2=cut_lon_rho(x2(1),y2(1));
% % lon1=cut_lon_rho(lat_min(1),lon_min(1));
% % lon2=cut_lon_rho(lat_max(1),lon_max(1));
% % lat1=cut_lat_rho(lat_min(1),lon_min(1));
% % lat2=cut_lat_rho(lat_max(1),lon_max(1));
% % x1=lon_min; x2=lon_max; y1=lat_min; y2=lat_max;

    if (lon2-lon1) >= (lat2-lat1)
        lon_line = cut_lon_rho(round((x1(1)+x2(1))/2),min(y1(1),y2(1)):max(y1(1),y2(1)));  %% for 1/20^o horizontal resolution
        lat_line = (lon_line-lon1)/((lon2-lon1)/(lat2-lat1))+lat1; %% weight for latitude index (lat1 : 0.05/((lon2-lon1) * (lat2-lat1) : lat2 
        x=repmat(lon_line,gd.N,1);  %% copy lon_line (gd.N times) to make matrix 
        x_label='Longitude(^oE)';
%                 domaxis=[domaxis(3) domaxis(4) domaxis(5) domaxis(6) domaxis(5) domaxis(6)];
%         Temp=zeros(gd.N,length(lon_line)); %% initialize temp matrix that size is same with x
        titlename = strcat(varname,', ',num2str(mean(section(3:4))),'N (',char(calendarname(tempmonth)), ', ', num2str(tempyear),')');
    else
        lat_line=cut_lat_rho(min(x1(1),x2(1)):max(x1(1),x2(1)),round((y1(1)+y2(1))/2));
        lon_line=(lat_line-lat1)*((lon2-lon1)/(lat2-lat1))+lon1;
        x=repmat(lat_line',gd.N,1);
        x_label='Latitude(^oN)';
%         Temp=zeros(gd.N,length(lat_line)); %% initialize temp matrix that size is same with x
        titlename = strcat(varname,', ',num2str(mean(section(1:2))),'E (',char(calendarname(tempmonth)), ', ', num2str(tempyear),')');
    end

    if (x1(1)==x2(1) || y1(1)==y2(1)) %% fixed lon or lat
        Temp(:,:) = squeeze(cut_data(:,min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1))));
        Yi(:,:)= squeeze(cut_depth(:,min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1))));
    else
        for k=1:1:gd.N
            lon_range=cut_lon_rho(min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1)));
            lat_range=cut_lat_rho(min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1)));
            data_range=squeeze(cut_data(k,min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1))));  %%data(zlevel, xmin : xmax, ymin : ymax)
            depth_range=squeeze(cut_depth(k,min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1)))); %%depth(zlevel, xmin : xmax, ymin : ymax)
            Temp(k,:)=griddata(lon_range,lat_range,data_range,lon_line,lat_line); %%get interpolated data section
            Yi(k,:)=griddata(lon_range,lat_range,depth_range,lon_line,lat_line); %%get interpolated depth section
        end
    end

    data=Temp;
    
    
    figure('position',[400 200 1500 550],'PaperUnits','inches','PaperPosition',[0 0 9.5 5]); %%   %%figure window, figure file 
    set(gca,'Position',[0.2 0.25 0.65 0.5]);  %% figure
    text_posi_x=(section(4)-section(3))/20+section(3);
    text_posi_y1=(section(6)-section(5))/20+section(5);
    text_posi_y2=2*(section(6)-section(5))/20+section(5);
    text_posi_y3=3*(section(6)-section(5))/20+section(5);
    %             switch section
    %                 case 1
            hold on
            pcolor(x,Yi,data)
            if (lon2-lon1) >= (lat2-lat1)        
                axis([section(1) section(2) section(5) section(6)]);
            else
                axis([section(3) section(4) section(5) section(6)]);
            end
            shading flat;
    %                 caxis(val_caxis)
            caxis(shadlev)
            set(gca,'box','on','linewidth',1.5,'fontsize',17)
            xlabel(x_label,'color','k','FontSize',17,'fontweight','bold')
            ylabel('Depth(m)','color','k','FontSize',17,'fontweight','bold')
    %                 title('Vertical Temperature','fontsize',17); % stlee
            title(titlename,'fontsize',17); % stlee
    %                 out_name_1=['YellowSea',out_name_1];
    %                 if (plot_contour)
              hold on
    %                   level_c =31:0.5:35;
%               [C,h2]=contour(x,Yi,data,conlev,'k','linewidth',1);
              [C,h2]=contour(x,Yi, data, conlev, m_contour_color, 'linewidth', m_contour_linewidth);
        clabel(C,h2,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
            'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);              
%               clabel(C,h,'FontSize',15,'Color','k','labelspacing',100000,'Rotation',0,'fontweight','bold');
    %                 end
    % set colorbar 
        h = colorbar;
        colormap(colormap_style);
        set(h,'fontsize',colorbar_fontsize);
        title(h,colorbar_title,'fontsize',colorbar_title_fontsize);
        caxis(shadlev);

        jpgname=strcat(outfile, '_', testname, '_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg
%         pause(2); %% prevent too short time between previous command and save command
        drawnow; %% prevent too short time between previous command and save command
        saveas(gcf,jpgname,'jpg');
        close all
end

status = 1;
end