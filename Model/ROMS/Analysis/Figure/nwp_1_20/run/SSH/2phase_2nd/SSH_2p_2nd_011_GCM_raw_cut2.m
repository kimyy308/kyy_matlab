close all; clear all;  clc;
warning off;

all_testname2 = {'CNRM-ESM2-1', 'EC-Earth3-Veg', 'ACCESS-CM2', 'CNRM-CM6-1-HR', 'CMCC-ESM2'};
% all_testname2 = {'CNRM-ESM2-1'};
% all_testname2 = {'EC-Earth3-Veg', 'ACCESS-CM2', 'CNRM-CM6-1-HR', 'CMCC-ESM2'};
% all_testname2 = {'CMCC-ESM2'};

all_region2 ={'NWP_GCM'};
all_var2 = {'zos', 'uo', 'vo', 'thetao', 'so'};

% all_region2 ={'NWP'}
for testnameind2=1:length(all_testname2)
    for regionind2=1:length(all_region2)
        close all;
        clearvars '*' -except regionind2 testnameind2 all_region2 all_testname2 all_var2
        % % % 
        % % % Read Model SST
        % % % interp
        % % % get RMS
        % % % get BIAS
        system_name=computer;
        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            dropboxpath='C:\Users\KYY\Dropbox';
            addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
            addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
            addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
            addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
            addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Roms_tools\Preprocessing_tools']));
            addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path
        elseif (strcmp(system_name,'GLNXA64'))
            dropboxpath='/home/kimyy/Dropbox';
            addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
            addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
            addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
            addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
        end


        % for snu_desktopd
        testname=all_testname2{testnameind2}    % % need to change
%         inputyear = [1985:2014]; % % put year which you want to plot [year year ...]
        inputyear = [2051:2100]; % % put year which you want to plot [year year ...]
        inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
        scenname ='ssp585';
%         variable ='zeta'
        run('nwp_polygon_point.m');
        regionname=all_region2{regionind2};
        switch(regionname)
            case('NWP_GCM') %% North western Pacific
                lonlat = [110, 170, 10, 60];  %% whole data area
                refpolygon(1,1)=lonlat(1);
                refpolygon(2,1)=lonlat(2);
                refpolygon(1,2)=lonlat(3);
                refpolygon(2,2)=lonlat(4);
            otherwise
                ('?')
        end
        lonlat(1)=min(refpolygon(:,1));
        lonlat(2)=max(refpolygon(:,1));
        lonlat(3)=min(refpolygon(:,2));
        lonlat(4)=max(refpolygon(:,2));

        % % % for EKB
        % regionname='EKB';
        % lonlat = [127, 129.5, 38, 40.5];

        if (strcmp(system_name,'PCWIN64'))
            % % for windows

        elseif (strcmp(system_name,'GLNXA64'))
            cmip6dir = ['/data2/CMIP6/public/', scenname, '/'];
            saverootdir = ['/data1/kimyy/Model/CMIP6/'];
        end
        
        

% start-------------------- later decadal SST, SSS plot
        for varind2=1:length(all_var2)
            variable=all_var2{varind2};
            savedir=[saverootdir,variable, '/', scenname, '/Omon/', testname, '/'];
            if (exist(strcat(savedir) , 'dir') ~= 7)
                mkdir(strcat(savedir));
            end 
            for yearij=1:length(inputyear)
                tempyear=inputyear(yearij);
                yearstr=num2str(tempyear, '%04i');
                
                switch testname
                    case {'CNRM-ESM2-1', 'CNRM-CM6-1-HR'}
                        ensname='r1i1p1f2';
                    case {'ACCESS-CM2', 'EC-Earth3-Veg', 'CMCC-ESM2'}
                        ensname='r1i1p1f1';
                end
                filedir = strcat(cmip6dir, variable, '/Omon/', testname, filesep, ensname, filesep, 'gn', filesep); % % where data files are
                flag_file_in = false;
                list = dir( [ filedir, filesep, variable, '*' ]); 
                for kk = 1 : length( list )
                    fname_in    = list(kk).name;
                    fname_split = strsplit( fname_in, {'_','.'} );
                    fyear_str   = strsplit( fname_split{end-1}, '-' );
                    fyear_start = str2num( fyear_str{1}(1:4) );
                    fyear_end   = str2num( fyear_str{2}(1:4) );
                    if( tempyear >= fyear_start && tempyear <= fyear_end &&     ...                 
                            strcmp( fname_split{2}, 'Omon' ) &&         ...
                            strcmp( fname_split{3}, testname ) &&      ...                 
                            strcmp( fname_split{4}, 'ssp585' ) )
                        flag_file_in = true;            break;
                    end         
                end         
                if( ~flag_file_in )
                    fprintf('Source File for %04i does not Exist. Continue...\n',tempyear);   continue;
                end
                filename=[filedir, filesep, fname_in];
                tind_min=(tempyear-fyear_start)*12+0;
                tind_max=(tempyear-fyear_start)*12+11;
                savefilename= [savedir, variable, '_Omon_',scenname, '_', testname, '_', yearstr,'.nc'];
                savelonname = [savedir, 'lon', '_Omon_',scenname, '_', testname,'.nc'];
                savelatname = [savedir, 'lat', '_Omon_',scenname, '_', testname,'.nc'];
                savedepthname = [savedir, 'depth', '_Omon_',scenname, '_', testname,'.nc'];
                savelev_boundsname = [savedir, 'lev_bounds', '_Omon_',scenname, '_', testname,'.nc'];
                switch(testname)
                    case{'CNRM-ESM2-1', 'CNRM-CM6-1-HR'}
                        lonname = 'lon';
                        latname = 'lat';
                        lonboundname = 'bounds_lon';
                        latboundname = 'bounds_lat';
                        xname = 'x';
                        yname = 'y';
                        zname = 'lev';
                        zboundname = 'lev_bounds';
                    case{'CMCC-ESM2'}
                        lonname = 'longitude';
                        latname = 'latitude';
                        lonboundname ='vertices_longitude';
                        latboundname ='vertices_latitude';
                        zname='lev';
                        xname = 'j';
                        yname = 'i';
                        zboundname='lev_bnds';
                    otherwise
                        lonname = 'longitude';
                        latname = 'latitude';
                        lonboundname ='vertices_longitude';
                        latboundname ='vertices_latitude';                        
                        zname='lev';
                        zboundname='lev_bnds';
                        xname = 'i';
                        yname = 'j';
                end
                
                if (exist('lon_rho' , 'var') ~= 1)
                    lon_rho=ncread(filename, lonname);
                    lat_rho=ncread(filename, latname);
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho, 1);
%                     cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
%                     cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                end
                lon_min_str=num2str(lon_min-1);
                lon_max_str=num2str(lon_max-1);
                lat_min_str=num2str(lat_min-1);
                lat_max_str=num2str(lat_max-1);
                
                
                
                    
%                 if (exist(savelonname, 'file') ~=2)
                    
                    system(['ncks --overwrite -C -v ', lonname, ',', lonboundname, ' ', ...
                   ' -d ', xname,',', lon_min_str, ',', lon_max_str, ' -d ', yname ',', lat_min_str, ',', lat_max_str, ' ', ...
                   filename, ' ', savelonname]);

                    system(['ncks --overwrite -C -v ', latname, ',', latboundname, ' ', ...
                   ' -d ', xname, ',', lon_min_str, ',', lon_max_str, ' -d ', yname, ',', lat_min_str, ',', lat_max_str, ' ', ...
                   filename, ' ', savelatname]);
%                 end
%                 if (exist(savedepthname, 'file') ~=2)
                    depthname ='lev';
                    if(strcmp(variable, 'zos')==1)
                    else
                        system(['ncks --overwrite -C -v ', zname, ',',zboundname, ' ',  filename, ' ', savedepthname]);
%                         system(['ncks --overwrite -C -v ', 'lev_bounds', ' ', filename, ' ', savelev_boundsname]);
                    end
%                 end
                
               system(['ncks --overwrite -C -v ', variable, ' -d time,',num2str(tind_min),',',num2str(tind_max), ...
                   ' -d ', xname, ',', lon_min_str, ',', lon_max_str, ' -d ', yname, ',', lat_min_str, ',', lat_max_str, ' ', ...
                   filename, ' ', savefilename]);
                disp([variable, ', ', yearstr])
            end
        end
    end
end