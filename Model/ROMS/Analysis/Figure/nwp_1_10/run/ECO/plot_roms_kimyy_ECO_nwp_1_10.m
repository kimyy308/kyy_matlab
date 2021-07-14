% % This code based on MATLAB R2016b.
% % Updated 15-Jun-2018 by Yong-Yub Kim

clc;close all;clear all;
warning off;

linux=0; windows=1;

if (windows ==1)
    % % for windows
    dropboxpath='C:\Users\KYY\Dropbox'; %% SNU_desktop
    addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
    addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
    addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run']));
elseif (linux==1)
    % % for linux
    dropboxpath='/home01/kimyy/Dropbox'; %% ROMS server
%     dropboxpath='/home/kimyy/Dropbox'; %% DAMO server
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
    addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_10/run']));
end


% fid=fopen('E:\Data\Model\ROMS\nwp_1_10\test42\spinup\1989\files\modelinfo');
% modelinfo=textscan(fid,'%s');
% fclose(fid);

% % % for snu_desktop_auto
% % testname= modelinfo{1,1}{1,1};
% % year = str2num(modelinfo{1,1}{3,1}); 

testname= 'test06'
year = [1980:1985]; % % put year which you want to plot [year year ...]

month = [1:12]; % % put month which you want to plot [month month ...]
inputdir = strcat(['E:\Data\Model\ROMS\nwp_1_10\input\',testname,'\']);
% figdir =strcat([outputdir,'\figures\']); % % where figure files will be saved

figdir =strcat(['D:\OneDrive - 서울대학교\MEPL\project\etc\eco\2019\figure\nwp_1_10\',testname,'\']); % % where figure files will be saved

if (exist(strcat(figdir) , 'dir') ~= 7)
    mkdir(strcat(figdir));
end 


% % % % % % % % % % plot topography
% % % % 
% % % % % etopodir ='D:\MEPL\project\SSH\2nd_year\figure\etopo1\';
% % % % % plot_roms_kimyy_etopo
% % % % 
% % % % if (exist(strcat(figdir,'\Bathy\bathy_nwp_1_10.jpg'),'file') ~= 2)
% % % %     if (exist(strcat(figdir,'Bathy') , 'dir') ~= 7)
% % % %         mkdir(strcat(figdir, 'Bathy'));
% % % %     end 
% % % %     plot_roms_kimyy_topo_nwp_1_10;
% % % % end
% % % % 
% % % % 
% % % % %  plot ECO
lonlat_ECO = [126.1 129.1 33.1 35.1]; % [lon_start lon_end lat_start lat_end]
% for i=1:length(year)   

%     tempyear=year(i);
    outputdir = strcat(['E:\Data\Model\ROMS\nwp_1_10\',testname,'\run\']); % % where data files are

    if (exist(strcat(figdir,'ECO'),'dir')~=7)
        mkdir(strcat(figdir,'ECO'));
    end
    
    % %  plot clim SST (ECO)
    if (exist(strcat(figdir,'ECO\','SST') , 'dir') ~= 7)
        mkdir(strcat(figdir,'ECO\','SST'));
    end
    sstfile =[figdir, 'ECO\', 'SST\', 'ECO_SST']; % % ~/SST_ECO_testname_year_month.jpg
    status=plot_ROMS_monthly_surface_UV_data_clim(testname, sstfile, outputdir, lonlat_ECO, year, month, ...
                                    [0 30], 0:5:30, 'SST', 'fig_param_kyy_pollack', 1, 0);    
                                
%     plot TUV (pollack)                            
    if (exist(strcat(figdir,'ECO\','TUV\','TUV_0m') , 'dir') ~= 7)
        mkdir(strcat(figdir,'ECO\','TUV\','TUV_0m'));
    end
    suvfile =[figdir, 'ECO\','TUV\TUV_0m\', 'pollack_0m_TUV']; % % ~/SUV_ECO_testname_year_month.jpg
    status=plot_ROMS_monthly_surface_UV_data_clim(testname, suvfile, outputdir, lonlat_ECO, year, month, ...
                                    [0 30], 0:10:30, 'SST', 'fig_param_kyy_eco', 0, 1);

    % %  plot clim SSS (ECO)
    if (exist(strcat(figdir,'ECO\','SSS') , 'dir') ~= 7)
        mkdir(strcat(figdir,'ECO\','SSS'));
    end
    sstfile =[figdir, 'ECO\', 'SSS\', 'ECO_SSS']; % % ~/SST_ECO_testname_year_month.jpg
    status=plot_ROMS_monthly_surface_UV_data_clim(testname, sstfile, outputdir, lonlat_ECO, year, month, ...
                                    [30 35], 30:1:35, 'SSS', 'fig_param_kyy_pollack', 1, 0);    
    
    
% end