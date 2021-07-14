function status=get_ROMS_section_data(testname, filedir, section, year, inputmonth, var, param_script)
% clear all;close all;
%==========================================================================
% % Updated 03-Apr-2018 by Y.Y.Kim

%==========================================================================
% vert_param;
romstools_param;
run(param_script);



% calendar{1} = ' Jan'; calendar{2} = ' Feb'; calendar{3} = ' Mar'; 
% calendar{4} = ' Apr'; calendar{5} = ' May'; calendar{6} = ' Jun';
% calendar{7} = ' Jul'; calendar{8} = ' Aug'; calendar{9} = ' Sep'; 
% calendar{10} = ' Oct'; calendar{11} = ' Nov'; calendar{12} = ' Dec';
tind=1;
for i=1:length(year)   
    tempyear=year(i);
    for monthij=1:length(inputmonth)
        tempmonth = inputmonth(monthij);
        filename = strcat(filedir, num2str(tempyear,'%04i'), '\', ...
                    testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');

    % %  slice data
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

        cut_lon_rho = lon_rho(lat_min:lat_max, lon_min:lon_max);
        cut_lat_rho = lat_rho(lat_min:lat_max, lon_min:lon_max);
        cut_depth = depth(:,lat_min:lat_max, lon_min:lon_max);

    % %  read data
        data_info = ncinfo(filename, varname); 
        data = ncread(filename,varname,[lon_min lat_min 1 1], [lon_max-lon_min+1 lat_max-lat_min+1 data_info.Size(3) 1]); %% read data[x y z t]
        data = permute(data, [3 2 1]); %% permute data[x y z t] -> [z y x]

    % %  slice section   
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

    % %     data[z hor] : variable data
    % %     x[z hor]    : lon or lat data
    % %     Yi[z hor]   : depth data

        sectionyear(tind)=tempyear;
        sectionmonth(tind)=tempmonth;
        sectionx=x;
        sectionz=Yi;
        sectiondata(tind,:,:)=data(:,:);
        tind=tind+1;
    end
end
sectionname = [num2str(section(1)),'_', num2str(section(2)), '_',num2str(section(3)), '_', num2str(section(4))];
fname = [filedir, 'isodepth\', 'ROMS_', testname, '_', var, '_section_', sectionname, '_', num2str(year(1)), '_', num2str(year(end)), '.mat'];
save(fname, 'sectionx', 'sectionz', 'var', 'section', 'sectionyear', 'sectionmonth', 'sectiondata'); 
status=1;
return