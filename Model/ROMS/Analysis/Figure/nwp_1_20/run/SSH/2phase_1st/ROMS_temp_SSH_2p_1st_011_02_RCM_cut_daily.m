close all; clear all;  clc;
warning off;

% all_testname2 = {'test2117', 'test2118', 'test2119', 'test2120', 'test2121'};
% all_testname2 = {'test2117', 'test2118'};
all_testname2 = {'test2121'};

all_region2 ={'ES'};

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
        inputyear = [1985:2014]; % % put year which you want to plot [year year ...]
%         inputyear = [1985]; % % put year which you want to plot [year year ...]
%         inputday = [1]; % % put month which you want to plot [month month ...]
        
%         scenname ='ssp585';
        scenname ='historical';
%         variable ='zeta'
        run('nwp_polygon_point.m');
        regionname=all_region2{regionind2};

        [refpolygon, lonlat, tmp.error_status] = Func_0007_get_polygon_data_from_regionname(regionname);


%         cmip6dir = ['/data2/RCM/CMIP6/Model/ROMS/nwp_1_20/output/monthly/'];
        cmip6dir = ['/data1/RCM/CMIP6/output_hdd3/nwp_1_20/output_daily/'];
        saverootdir = ['/data2/RCM/CMIP6/Model/ROMS/nwp_1_20/cut', '_',regionname, '_daily', filesep];



% start-------------------- later decadal SST, SSS plot
            
            for yearij=1:length(inputyear)
                tempyear=inputyear(yearij);
                
                inputday = [1:91, 334:eomday(tempyear,2)+337]; % % put month which you want to plot [month month ...]

                yearstr=num2str(tempyear, '%04i');
                savedir=[saverootdir, testname, '/', yearstr, filesep];
                    if (exist(strcat(savedir) , 'dir') ~= 7)
                        mkdir(strcat(savedir));
                    end 
%                 filedir = strcat(cmip6dir, testname, filesep, 'run', filesep, ...
%                     'packed_monthly', filesep, yearstr, filesep); % % where data files are
                filedir = strcat(cmip6dir, testname, filesep, 'run', filesep, ...
                    'packed_daily_avg', filesep, yearstr, filesep); % % where data files are
                for dayij = 1:length(inputday)
                    tempday=inputday(dayij);
                    daystr=num2str(tempday, '%04i');
                    flag_file_in = false;
                    list = dir( [ filedir,  '*' ]); 
                    for kk = 3 : length( list )
                        fname_in    = list(kk).name;
%                         fname_in    = list(3).name;

                        fname_split = strsplit( fname_in, {'_','.'} );
                        fyear_str   = strsplit( fname_split{end-4}, '-' );
                        fmonth_str   = strsplit( fname_split{end-1}, '-' );
                        fyear= str2num( fyear_str{1}(1:4) );
                        fday=str2num( fmonth_str{1}(1:4) );
                        if( tempyear == fyear  &&     ...    
                                tempday == fday && ...
                                strcmp( fname_split{2}, testname ) &&         ...
                                strcmp( fname_split{4}, 'daily' )   )              
                            flag_file_in = true;            break;
                        end         
                    end         
                    if( ~flag_file_in )
                        fprintf('Source File for %04i does not Exist. Continue...\n',tempyear);   continue;
                    end
                    filename=[filedir, fname_in];
                    
                    savefilename= [savedir, 'pck_', regionname, '_', testname, '_daily_', yearstr, '_', daystr, '.nc'];

                    lonname = 'lon_rho';
                    latname = 'lat_rho';
                    xrhoname = 'xi_rho';
                    yrhoname = 'eta_rho';
                    xuname = 'xi_u';
                    yuname = 'eta_u';
                    xvname = 'xi_v';
                    yvname = 'eta_v';

                    if (exist('lon_rho' , 'var') ~= 1)
                        lon_rho=ncread(filename, lonname);
                        lat_rho=ncread(filename, latname);
                        [lon_min, lon_max, lat_min, lat_max]=Func_0012_findind_Y(1/20,lonlat(1:4),lon_rho,lat_rho, 1);
                    end
                    lon_min_str=num2str(lon_min-1);
                    lon_max_str=num2str(lon_max-1);
                    lat_min_str=num2str(lat_min-1);
                    lat_max_str=num2str(lat_max-1);
                    
                    lonu_max_str=num2str(lon_max-2);
                    latv_max_str=num2str(lat_max-2);

%                    system(['ncks --overwrite -C -v ', variable, ' -d time,',num2str(tind_min),',',num2str(tind_max), ...
%                        ' -d ', xname, ',', lon_min_str, ',', lon_max_str, ' -d ', yname, ',', lat_min_str, ',', lat_max_str, ' ', ...
%                        filename, ' ', savefilename]);

                    system(['ncks --overwrite -C ', ...
                       ' -d ', xrhoname, ',', lon_min_str, ',', lon_max_str, ' -d ', yrhoname, ',', lat_min_str, ',', lat_max_str, ' ', ...
                       ' -d ', xuname, ',', lon_min_str, ',', lonu_max_str, ' -d ', yuname, ',', lat_min_str, ',', lat_max_str, ' ', ...
                       ' -d ', xvname, ',', lon_min_str, ',', lon_max_str, ' -d ', yvname, ',', lat_min_str, ',', latv_max_str, ' ', ...
                       filename, ' ', savefilename]);
                    disp([yearstr, ', ', daystr])
                end
            end
    end
end