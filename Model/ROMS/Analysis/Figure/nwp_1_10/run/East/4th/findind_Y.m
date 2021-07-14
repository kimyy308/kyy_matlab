function [indw, inde, inds, indn]=findind(dl,section,lon_rho,lat_rho)
% clear all;close all;
%==========================================================================
% %  Updated 11-May-2018 by Yong-Yub Kim
% %  Updated 23-Aug-2018 by Yong-Yub Kim
% %  Updated 03-Jun-2019 by Yong-Yub Kim, consider cuvelinear grid (UM)
%==========================================================================
% % section --> (lon1, lon2, lat1, lat2) or (lon1, lat1)
% % lon_rho --> [y x]
% % lat_rho --> [y x]

    if length(section)==2
        section(3:4)=section(2);
        section(2)=section(1);
    end

    if(size(lon_rho,1)==1 || size(lon_rho,2)==1)
        lon_west = abs(lon_rho - (section(1)-dl));
        min_lon_west=min(lon_west(:));
        lon_east = abs(lon_rho - (section(2)+dl));
        min_lon_east=min(lon_east(:));
        lat_south = abs(lat_rho - (section(3)-dl));
        min_lat_south=min(lat_south(:));
        lat_north = abs(lat_rho - (section(4)+dl));
        min_lat_north=min(lat_north(:));

        indw = find(lon_west(:) == min_lon_west);
        inde = find(lon_east(:) == min_lon_east);
        inds = find(lat_south(:) == min_lat_south);
        indn = find(lat_north(:) == min_lat_north);
        
        indw = indw(1);
        inde = inde(end);
        inds = inds(1);
        indn = indn(end);
    else
        %% for meshed grid (UM, ~)
        
        max_row_lon=max(lon_rho,[],2);
        min_row_lon=min(lon_rho,[],2);
        max_col_lat=max(lat_rho,[],1);
        min_col_lat=min(lat_rho,[],1);
        
%         max_row_lon=max(lon_rho,[],1);
%         min_row_lon=min(lon_rho,[],1);
%         max_col_lat=max(lat_rho,[],2);
%         min_col_lat=min(lat_rho,[],2);
    
        lon_west = abs(max_row_lon - (section(1)-dl));
        min_lon_west=min(lon_west);
        lon_east = abs(min_row_lon - (section(2)+dl));
        max_lon_east=min(lon_east);
        lat_south = abs(max_col_lat - (section(3)-dl));
        min_lat_south=min(lat_south);
        lat_north = abs(min_col_lat - (section(4)+dl));
        max_lat_north=min(lat_north);

        indw = find(lon_west(:) == min_lon_west);
        inde = find(lon_east(:) == max_lon_east);
        inds = find(lat_south(:) == min_lat_south);
        indn = find(lat_north(:) == max_lat_north);
        
        indw = indw(1);
        inde = inde(end);
        inds = inds(1);
        indn = indn(end);
    end
end