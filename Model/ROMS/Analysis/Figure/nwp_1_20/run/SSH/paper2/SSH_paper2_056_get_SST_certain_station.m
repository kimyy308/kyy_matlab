close all; clear all;  clc;

all_testname = {'ens10'};

% all_testname = {'test67', 'ens10'};

all_region ={'AKP4'};
% all_region ={'TEST'};
all_stationname = {'Donghae_Gas_Field', 'Wolseong', 'Guryongpo', 'Hyunpo', 'Deep_water_center'};
all_lonlat = [130.015138, 35.448027; ...
              129.544439, 35.695533; ...
              129.829711, 35.988450; ...
              130.823408, 37.549889; ...
              128.533238, 38.326033];
              
for testnameind=1:length(all_testname)
    system_name=computer;
    if (strcmp(system_name,'PCWIN64'))
        % % for windows
        dropboxpath='C:\Users\User\Dropbox';
        addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
        addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
        addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
        addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
        addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path
    elseif (strcmp(system_name,'GLNXA64'))
        dropboxpath='/home/kimyy/Dropbox';
        addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
        addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
        addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
        addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
    end

    shadlev = [-2 2];
    % rms_shadlev = [0 4];
    trend_shadlev = [0 4];
    % bias_shadlev = [-4 4];
    conlev  = 0:5:35;
    dl=1/10;

    % for snu_desktop
    testname=all_testname{testnameind} 
    inputyear = [2006:2100]; % % put year which you want to plot [year year ...]
    inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
    switch testname
        case {'test61', 'test62', 'test63', 'test64'}
            drivename='E:\';
        case {'test57', 'test58', 'test59', 'test60'}
            drivename='I:\';
        case {'test65', 'test66', 'test67', 'test68'}
            drivename='G:\';
        case {'ens08', 'ens09', 'ens10'}
            drivename='E:\';
    end
    system_name=computer;
    if (strcmp(system_name,'PCWIN64'))
        % % for windows
%             figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\4th_year\figure\nwp_1_20\',testname,'\',regionname,'\'); % % where figure files will be saved
        param_script ='C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m'
        filedir = strcat(drivename,'Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
        savedir = strcat('D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\');
        cmemsdir='Z:\내 드라이브\Data\Observation\CMEMS\';
    elseif (strcmp(system_name,'GLNXA64'))
    end
    
    variable = 'SST';
    varname ='temp';
% % %         get model data
    for ind_sta=1:length(all_stationname)
        sta_name = all_stationname{ind_sta};
        txtname = [savedir,testname,'_',sta_name,'_sst_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.txt'];
        tgtlon=all_lonlat(ind_sta,1);
        tgtlat=all_lonlat(ind_sta,2);
        filename = strcat(filedir, num2str(min(inputyear),'%04i'), '\', ...
                        testname, '_monthly_', num2str(min(inputyear),'%04i'), '_', num2str(min(inputmonth),'%02i'), '.nc');
        lon_rho=ncread(filename, 'lon_rho');
        lat_rho=ncread(filename, 'lat_rho');
        s_rho=ncread(filename, 's_rho');
        lonsize= size(lon_rho,1);
        latsize= size(lat_rho,2);
        zsize= length(s_rho);
        
        data_2d=ncread(filename, varname, [1 1 zsize 1], [lonsize, latsize, 1 1]);
        for i=1:lonsize
            for j=1:latsize
                if isfinite(data_2d(i,j))
                    dist(i,j)=m_lldist([lon_rho(i,j), tgtlon], [lat_rho(i,j), tgtlat]);
                else
                    dist(i,j)=NaN;
                end
            end
        end

        ind_sta_model(ind_sta)=find(dist(:)==min(dist(:)));
        ind_lon=mod(ind_sta_model(ind_sta),size(lon_rho,1));
        ind_lat=floor(ind_sta_model(ind_sta)/size(lon_rho,1))+1;
        lon_sta = lon_rho(ind_lon, ind_lat);
        lat_sta = lat_rho(ind_lon, ind_lat);
%         lon_rho(ind_lon,ind_lat)
%         lat_rho(ind_lon,ind_lat)
        ind=1;
        fid = fopen(txtname,'w+');
        fprintf(fid,[sta_name, ' \n']);
        fprintf(fid,'Year Mon Longitude Latitude    SST\n');
        fclose(fid);
        for yearij = 1:length(inputyear)
            for monthij = 1:length(inputmonth)
                disp([num2str(yearij), 'y_',num2str(monthij),'m'])
                tempyear = inputyear(yearij);
                tempmonth = inputmonth(monthij);
                yearstr = num2str(tempyear, '%04i');
                monthstr = num2str(tempmonth, '%02i');
                % ex : rootdir\test37\data\2001\test37_monthly_2001_01.nc
                filename = strcat(filedir, yearstr, '\', ...
                    testname, '_monthly_', yearstr, '_', monthstr, '.nc');
%                 data_info = ncinfo(filename, varname);  %% [lon lat depth time] -> [1601 1201 33 1]
                data = ncread(filename,varname,[ind_lon ind_lat zsize 1], [1 1 1 1]);
                fid = fopen(txtname,'a+');
                fprintf(fid,'%5.0f %3.0f %8.3f %8.3f %8.3f \n', tempyear, tempmonth, lon_sta, lat_sta, data);
                fclose(fid);
                comb_data(ind) = single(data);

                ind = ind + 1;
            end
        end      
    end
end
