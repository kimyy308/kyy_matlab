close all; clc; clear all;

ncload('D:\Data\Model\ROMS\Results_Jung_and_Cho_2020\roms_grid_EYECS_20190904.nc');

% 127.90, 39.52, -0.25, 172800, 101001




for loni=1:200
    for lati=1:500
        if mask_rho(loni,lati)==1
            fdata = [lon_rho(loni,lati), lat_rho(loni,lati), -0.25, 86400*2.5, 101001];
            if (exist(['D:\Data\Model\ROMS\Results_Jung_and_Cho_2020\loc.txt'])~=2)
                fid=fopen(['D:\Data\Model\ROMS\Results_Jung_and_Cho_2020\loc.txt'], 'w');
            else
                fid=fopen(['D:\Data\Model\ROMS\Results_Jung_and_Cho_2020\loc.txt'], 'a');
            end
            fprintf(fid,'%6.2f,%5.2f,%5.2f,%5d,%6d \r\n',fdata);
            fclose(fid);
        end
    end
end