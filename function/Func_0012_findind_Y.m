function [indw, inde, inds, indn]=Func_0012_findind_Y(dl,section,lon_rho,lat_rho, all)
% clear all;close all;
%==========================================================================
% %  Updated 11-May-2018 by Yong-Yub Kim
% %  Updated 23-Aug-2018 by Yong-Yub Kim
% %  Updated 03-Jun-2019 by Yong-Yub Kim, consider cuvelinear grid (UM)
% %  Updated 05-Jul-2021 by Yong-Yub Kim, modified to Func_0012
%==========================================================================
% % section --> (lon1, lon2, lat1, lat2) or (lon1, lat1)
% % lon_rho --> [y x]
% % lat_rho --> [y x]
% % all == 1 : judge valid grid from all of them (2-dimension)

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

%         lon_west = abs(lon_rho - (section(1)));
%         min_lon_west=min(lon_west(:));
%         lon_east = abs(lon_rho - (section(2)));
%         min_lon_east=min(lon_east(:));
%         lat_south = abs(lat_rho - (section(3)));
%         min_lat_south=min(lat_south(:));
%         lat_north = abs(lat_rho - (section(4)));
%         min_lat_north=min(lat_north(:));

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
%---------------- version 1
        max_row_lon=max(lon_rho,[],2);
        min_row_lon=min(lon_rho,[],2);
        max_col_lat=max(lat_rho,[],1);
        min_col_lat=min(lat_rho,[],1);
    
        lon_west = abs(max_row_lon - (section(1)-dl));
        min_lon_west=min(lon_west);
        lon_east = abs(min_row_lon - (section(2)+dl));
        max_lon_east=min(lon_east);
        lat_south = abs(max_col_lat - (section(3)-dl));
        min_lat_south=min(lat_south);
        lat_north = abs(min_col_lat - (section(4)+dl));
        max_lat_north=min(lat_north);
%         lon_west = abs(max_row_lon - (section(1)));
%         min_lon_west=min(lon_west);
%         lon_east = abs(min_row_lon - (section(2)));
%         max_lon_east=min(lon_east);
%         lat_south = abs(max_col_lat - (section(3)));
%         min_lat_south=min(lat_south);
%         lat_north = abs(min_col_lat - (section(4)));
%         max_lat_north=min(lat_north);
        

        indw = find(lon_west(:) == min_lon_west);
        inde = find(lon_east(:) == max_lon_east);
        inds = find(lat_south(:) == min_lat_south);
        indn = find(lat_north(:) == max_lat_north);
        
        indw = indw(1);
        inde = inde(end);
        inds = inds(1);
        indn = indn(end);

% --------------- version 2
% % % %         lon_rho(x,y), lat_rho(x,y)
        if (nargin ==5)
            lon_rho_temp=lon_rho;
            lat_rho_temp=lat_rho;
            lon_rho_temp(lon_rho<(section(1)-dl))=NaN;
            lon_rho_temp(lon_rho>(section(2)+dl))=NaN;
            lon_rho_temp(lat_rho<(section(3)-dl))=NaN;
            lon_rho_temp(lat_rho>(section(4)+dl))=NaN;
            lat_rho_temp(lon_rho<(section(1)-dl))=NaN;
            lat_rho_temp(lon_rho>(section(2)+dl))=NaN;
            lat_rho_temp(lat_rho<(section(3)-dl))=NaN;
            lat_rho_temp(lat_rho>(section(4)+dl))=NaN;
%             lon_rho_temp(lon_rho<(section(1)))=NaN;
%             lon_rho_temp(lon_rho>(section(2)))=NaN;
%             lon_rho_temp(lat_rho<(section(3)))=NaN;
%             lon_rho_temp(lat_rho>(section(4)))=NaN;
%             lat_rho_temp(lon_rho<(section(1)))=NaN;
%             lat_rho_temp(lon_rho>(section(2)))=NaN;
%             lat_rho_temp(lat_rho<(section(3)))=NaN;
%             lat_rho_temp(lat_rho>(section(4)))=NaN;
            
            
%             lon_rho(isnan(lat_rho))=NaN;
%             lat_rho(isnan(lon_rho))=NaN;
%             find(isnan(lon_rho)==0,2);
            lon_y=sum(lon_rho_temp,2,'omitnan');
            lon_y_min=min(find(lon_y~=0));
            lon_y_max=max(find(lon_y~=0));
            lat_x=sum(lat_rho_temp,1,'omitnan');
            lat_x_min=min(find(lat_x~=0));
            lat_x_max=max(find(lat_x~=0));
            
            if isempty(lon_y_min) && isempty(lon_y_max) && isempty(lat_x_min) && isempty(lat_x_max)
                indw = 1;
                inde = 1;
                inds = 1;
                indn = 1;
            else
                indw = lon_y_min(1);
                inde = lon_y_max(end);
                inds = lat_x_min(1);
                indn = lat_x_max(end);
            end
        end
    end
end