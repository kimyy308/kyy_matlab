

close all; clear all;  clc;
all_region ={'AKP2'}
regionind=1;
%for regionind=1:length(all_region)
    clearvars '*' -except regionind all_region
     system_name=computer;
    if (strcmp(system_name,'PCWIN64'))
        % % for windows
        dropboxpath='C:\Users\KYY\Dropbox';
        addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
        addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
        addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
        addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
    elseif (strcmp(system_name,'GLNXA64'))
        dropboxpath='/home/kimyy/Dropbox';
        addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
        addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
        addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
        addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
    end
    
    inputyear = [1993:2017]; % % put year which you want to plot [year year ...]
    inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
    
    
    varname ='thetao';
    run('nwp_polygon_point.m');
    regionname=all_region{regionind}
    switch(regionname) % 
        case('AKP2')
            refpolygon=akp2polygon;
        otherwise
            ('?')
    end
    
    % ¿øÇÏ´Â ¹üÀ§ÀÇ À§°æµµÀÇ ÃÖ´ëÃÖ¼Ò¹üÀ§ÁöÁ¤
    lonlat(1)=min(refpolygon(:,1));
    lonlat(2)=max(refpolygon(:,1));
    lonlat(3)=min(refpolygon(:,2));
    lonlat(4)=max(refpolygon(:,2));
    
    if (strcmp(system_name,'PCWIN64'))
        % % for windows
        figrawdir =strcat('D:\OneDrive - ????????\MEPL\project\MICT_pollack\3rd year\figure\nwp_1_20\',testname,'\',regionname,'\'); % % where figure files will be saved
        param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m'
        filedir = strcat('E:\Data\Observation\OIssh\monthly\'); % % where data files are
        MyOceandir='E:\Data\Observation\OIssh\monthly\';
    elseif (strcmp(system_name,'GLNXA64'))
        figrawdir =strcat('/data1/stlee/ext_hdd/MyOcean/',regionname,'/'); % % where figure files will be saved
        param_script ='/home/kimyy/Dropbox/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_20/run/fig_param/fig_param_kyy_EKB_RMS.m'
        filedir = strcat('/data1/stlee/ext_hdd/MyOcean/'); % % where data files are
        MyOceandir='/data1/stlee/ext_hdd/MyOcean/';
    end
    
    run(param_script); % ±×¸²±×¸®´Â ¼¼ÆÃÀÌ µé¾î°¡ ÀÖÀ½
    tt=[];
    tt1=[];
    tt2=[];
    ind=1;
    comb_MyOcean_data=NaN(306,240,300,50);
    for yearij = 1:length(inputyear)
        for monthij = 1:length(inputmonth)
            tempyear = inputyear(yearij);
            tempmonth = inputmonth(monthij);
            %1993/mercatorglorys12v1_gl12_mean_199301.nc
            % read OItemp DATA
            MyOceanfilename = strcat(MyOceandir,num2str(tempyear,'%04i'),'/','mercatorglorys12v1_gl12_mean_', num2str(tempyear,'%04i'),num2str(tempmonth,'%02i'),'.nc');
            if (exist('MyOcean_lon')==0)
                MyOceaninfo=ncinfo(MyOceanfilename);
                % ÀüÃ¼ À§°æµµ ÃßÃâ
                MyOcean_lon = ncread(MyOceanfilename,'longitude',[1],[MyOceaninfo.Dimensions(1).Length]);
                MyOcean_lat = ncread(MyOceanfilename,'latitude',[1],[MyOceaninfo.Dimensions(2).Length]);
                % ¹üÀ§ ÁöÁ¤
                MyOcean_lon_west = abs(MyOcean_lon - (lonlat(1)));
                min_MyOcean_lon_west=min(MyOcean_lon_west);
                MyOcean_lon_east = abs(MyOcean_lon - (lonlat(2)));
                min_MyOcean_lon_east=min(MyOcean_lon_east);
                MyOcean_lat_south = abs(MyOcean_lat - (lonlat(3)));
                min_MyOcean_lat_south=min(MyOcean_lat_south);
                MyOcean_lat_north = abs(MyOcean_lat - (lonlat(4)));
                min_MyOcean_lat_north=min(MyOcean_lat_north);
                
                
                MyOcean_lon_min = find(MyOcean_lon_west == min_MyOcean_lon_west);
                MyOcean_lon_max = find(MyOcean_lon_east == min_MyOcean_lon_east);
                MyOcean_lat_min = find(MyOcean_lat_south == min_MyOcean_lat_south);
                MyOcean_lat_max = find(MyOcean_lat_north == min_MyOcean_lat_north);
            
                % ¿øÇÏ´Â ºÎºÐÀÇ À§°æµµ ÃßÃâ
                MyOcean_lon = ncread(MyOceanfilename,'longitude', [MyOcean_lon_min(1)], [MyOcean_lon_max(1)-MyOcean_lon_min(1)]);
                MyOcean_lat = ncread(MyOceanfilename,'latitude', [MyOcean_lat_min(1)], [MyOcean_lat_max(1)-MyOcean_lat_min(1)]);
                MyOcean_st = ncread(MyOceanfilename,'depth', 1,50);
                
                %comb_spatial_meanrms=(zeros([length(MyOcean_lon),length(MyOcean_lat),12]));
                %comb_spatial_meanbias=(zeros([length(MyOcean_lon),length(MyOcean_lat),12]));
                comb_spatial_meanMyOcean=(zeros([length(MyOcean_lon),length(MyOcean_lat),length(MyOcean_st),12]));
                %comb_spatial_meanmodel=(zeros([length(MyOcean_lon),length(MyOcean_lat),12]));
                % À§°æµµ ¸Þ½¬±×¸®µå ¸¸µé±â
                [MyOcean_lat2 MyOcean_lon2 MyOcean_st2]=meshgrid(MyOcean_lat, MyOcean_lon, MyOcean_st);
                switch(regionname)
                    case('NWP') %% North western Pacific
                        mask_MyOcean(1:size(MyOcean_lon,1),1:size(MyOcean_lon,2))=1;
                    otherwise
                        mask_MyOcean = double(inpolygon(MyOcean_lon2,MyOcean_lat2,refpolygon(:,1),refpolygon(:,2)));
                        mask_MyOcean(mask_MyOcean==0)=NaN;
                end
                %                 % mask¸¦ 3Â÷¿øÀ¸·Î ¸¸µé±â // ÀÌ¹Ì 51x 40 x 50 À¸·Î µÇ¾îÀÖÀ½
                %                 for xy=1:length(MyOcean_st)
                %                     mask_MyOcean2(:,:,xy)=mask_MyOcean;
                %                 end
            end
            len_lon = length(MyOcean_lon(:));
            len_lat = length(MyOcean_lat(:));
            len_st = length(MyOcean_st(:));
            
            MyOcean_data = ncread(MyOceanfilename,varname,[MyOcean_lon_min(1) MyOcean_lat_min(1) 1 1], [MyOcean_lon_max(1)-MyOcean_lon_min(1) MyOcean_lat_max(1)-MyOcean_lat_min(1) 50 1]);
            MyOcean_data(MyOcean_data<-1000)=NaN;
            MyOcean_data(MyOcean_data>1000)=NaN;
            MyOcean_data=MyOcean_data.*mask_MyOcean; % Æ¯Á¤³âµµ,¿ùÀÇ ½ÇÁ¦ ¼ö¿Â°ªÀÌ ¼ö½Éº°·Î ±â·ÏµÊ
            comb_MyOcean_data(:,:,ind,:) = MyOcean_data;
            ind = ind+1;

            
        end
     end
%      set(gca,'ydir','reverse');

% save([filedir,regionname,'Myocean_temp_raw_data','.mat'],'MyOcean_lon','MyOcean_lat',...
%     'lonlat','varname','MyOcean_st','comb_MyOcean_data','all_region',...
%     'akp2polygon', 'inputyear', 'inputmonth', 'lonlat', 'filedir');
