theta = lon_end - lon_beg;

dist = sin(deg2rad(lat_beg)) * sin(deg2rad(lat_end)) ...
        + cos(deg2rad(lat_beg)) * cos(deg2rad(lat_end)) * cos(deg2rad(theta)) ;
dist = acos(dist);
dist = rad2deg(dist);
dist = dist * 60 * 1.1515;
dist = dist * 1.609344;    %% change the unit as mile to km
dist = dist * 1000.0;  %% change the unit as km to m