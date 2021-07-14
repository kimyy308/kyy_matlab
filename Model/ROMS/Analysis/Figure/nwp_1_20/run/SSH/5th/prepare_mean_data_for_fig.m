function [cut_lon_rho, cut_lat_rho, mean_data, land] = prepare_mean_data_for_fig (testname, scenname, varname, inputyear);
    
    switch testname
        case {'test61', 'test62', 'test63', 'test64'}
            drivename='H:\';
        case {'test57', 'test58', 'test59', 'test60'}
            drivename='I:\';
        case {'test65', 'test66', 'test67', 'test68'}
            drivename='G:\';
        case {'ens08', 'ens09', 'ens10'}
            drivename='E:\';
    end

    filedir = strcat(drivename,'Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
    matdir = strcat('D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\mean\');
    matname = [matdir, testname, '_', regionname, '_', variable, '_mean_', num2str(min(inputyear),'%04i'), '-', num2str(max(inputyear),'%04i'), '.mat'];
    if (exist(matname , 'file') ~= 2 || fig_flag==2)
        for yearij=1:length(inputyear)
            tempyear=inputyear(yearij);
            yearstr=num2str(tempyear, '%04i');
            for monthij=1:length(inputmonth)
                tempmonth=inputmonth(monthij);
                monthstr=num2str(tempmonth, '%02i');
                filename=[filedir, yearstr, '\', testname, '_monthly_', yearstr, '_', monthstr, '.nc'];

                if (exist('lon_rho' , 'var') ~= 1)
                    lon_rho=ncread(filename, 'lon_rho');
                    lat_rho=ncread(filename, 'lat_rho');
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                    cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                    cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end
                data_info = ncinfo(filename, varname); 

                if (strcmp(variable,'SST')==1 || strcmp(variable,'SSS')==1)
                    data = ncread(filename,varname,[lon_min(1) lat_min(1) data_info.Size(3) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                elseif (strcmp(variable,'SSH')==1)
                    data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                elseif (strcmp(variable,'BT')==1)
                    data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
                end
                if (exist('mean_data' , 'var') ~= 1)
                    mean_data=zeros(size(data));
                end
                mean_data=mean_data + (data / (length(inputyear1) * length(inputmonth)));
            end
        end
        save(matname, 'mean_data', 'cut_lon_rho', 'cut_lat_rho');
    else
        load(matname);
    end










end