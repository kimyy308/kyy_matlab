close all;clear all;clc

% all_region ={'NWP','ES', 'SS', 'YS', 'ECS'}
all_region ={'AKP2','ES','SS','YS'}


for regionind=1%:length(all_region)
    
    %%
    clearvars '*' -except regionind all_region
    
    % % %
    % % % Read Model ssh
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
    elseif (strcmp(system_name,'GLNXA64'))
        dropboxpath='/home/kimyy/Dropbox';
        addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
        addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
        addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
        addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
    end
    
    
    shadlev = [0 35];
    rms_shadlev = [0 4];
    bias_shadlev = [-4 4];
    conlev  = 0:5:35;
    dl=1/20;
    % for snu_desktop
    %     testname='test49'   % % need to change
    inputyear = [1993:2017]; % % put year which you want to plot [year year ...]
    inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
    
    
    
    varname ='thetao';  %% reference variable -> temperature
    run('nwp_polygon_point.m');
    regionname=all_region{regionind}
    switch(regionname)
        case('NWP') %% North western Pacific
            lonlat = [115, 164, 15, 52];  %% whole data area
            refpolygon(1,1)=lonlat(1);
            refpolygon(2,1)=lonlat(2);
            refpolygon(1,2)=lonlat(3);
            refpolygon(2,2)=lonlat(4);
        case('NWP2') %% North western Pacific
            lonlat = [115, 145, 25, 52];  %% whole data area
            refpolygon(1,1)=lonlat(1);
            refpolygon(2,1)=lonlat(2);
            refpolygon(1,2)=lonlat(3);
            refpolygon(2,2)=lonlat(4);
        case('ES') %% East Sea
            refpolygon=espolygon;
        case('SS') %% South Sea
            refpolygon=sspolygon;
        case('YS') %% Yellow Sea
            refpolygon=yspolygon;
        case('ECS') %% East China Sea
            refpolygon=ecspolygon;
        case('EKB') %% North western Pacific
            lonlat = [127, 131, 37, 42];  %% East Korea Bay
            refpolygon(1,1)=lonlat(1);
            refpolygon(2,1)=lonlat(2);
            refpolygon(1,2)=lonlat(3);
            refpolygon(2,2)=lonlat(4);
        case('AKP2')
            refpolygon=akp2polygon;
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
    
    
    run(param_script);
    tt1=[];
    tt2=[];
    ind=1;
    for yearij = 1:length(inputyear)
        for monthij = 1:length(inputmonth)
            %% raw ÀÚ·á »Ì±â
            disp([num2str(yearij), 'y_',num2str(monthij),'m'])
            tic;
            tempyear = inputyear(yearij);
            tempmonth = inputmonth(monthij);
            % ex : rootdir\test37\data\2001\test37_monthly_2001_01.nc
            
            % read OIssh DATA
            MyOceanfilename = strcat(MyOceandir,num2str(tempyear,'%04i'),'/','mercatorglorys12v1_gl12_mean_', ...
                num2str(tempyear,'%04i'),num2str(tempmonth,'%02i'),'.nc');
            if (exist('MyOcean_lon')==0)
                MyOceaninfo=ncinfo(MyOceanfilename);
                MyOcean_lon = ncread(MyOceanfilename, 'longitude',[1],[MyOceaninfo.Dimensions(1).Length]);
                MyOcean_lat = ncread(MyOceanfilename, 'latitude',[1],[MyOceaninfo.Dimensions(2).Length]);
                MyOcean_st = ncread(MyOceanfilename,'depth', 1,50);
                
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
                
                %         ncinfo('E:\Data\Observation\OIssh\monthly\MyOcean_monthly1983_11.nc');
                
                MyOcean_lon = ncread(MyOceanfilename,'longitude', [MyOcean_lon_min(1)], [MyOcean_lon_max(1)-MyOcean_lon_min(1)]);
                MyOcean_lat = ncread(MyOceanfilename,'latitude', [MyOcean_lat_min(1)], [MyOcean_lat_max(1)-MyOcean_lat_min(1)]);
                
                comb_spatial_meanrms=(zeros([length(MyOcean_lon),length(MyOcean_lat),12]));
                comb_spatial_meanbias=(zeros([length(MyOcean_lon),length(MyOcean_lat),12]));
                comb_spatial_meanMyOcean=(zeros([length(MyOcean_lon),length(MyOcean_lat),12]));
                comb_spatial_meanmodel=(zeros([length(MyOcean_lon),length(MyOcean_lat),12]));
                
                [MyOcean_lat2 MyOcean_lon2 MyOcean_st2]=meshgrid(MyOcean_lat, MyOcean_lon, MyOcean_st);
                switch(regionname)
                    case('NWP') %% North western Pacific
                        mask_MyOcean(1:size(MyOcean_lon,1),1:size(MyOcean_lon,2))=1;
                    otherwise
                        mask_MyOcean = double(inpolygon(MyOcean_lon2,MyOcean_lat2,refpolygon(:,1),refpolygon(:,2)));
                        mask_MyOcean(mask_MyOcean==0)=NaN;
                end
            end
            len_lon = length(MyOcean_lon(:));
            len_lat = length(MyOcean_lat(:));
            len_st = length(MyOcean_st(:));
            
            MyOcean_data = ncread(MyOceanfilename,varname,[MyOcean_lon_min(1) MyOcean_lat_min(1) 1 1], [MyOcean_lon_max(1)-MyOcean_lon_min(1) MyOcean_lat_max(1)-MyOcean_lat_min(1) 50 1]);
            MyOcean_data(MyOcean_data<-1000)=NaN;
            MyOcean_data(MyOcean_data>1000)=NaN;
            MyOcean_data=MyOcean_data.*mask_MyOcean;
            
            %% raw data profile
            comb_temp_raw_data(:,:,:,ind) = MyOcean_data;
            
            %% interpolation ÀÛ¾÷
            % ¼ö½É ÀÚ·á¿¡ Ç¥Ãþ(0m)Ãß°¡
            %         MyOcean_st3=[0;MyOcean_st];
            %         v5=zeros(length(MyOcean_lon),length(MyOcean_lat),300);% vertical mean of temp
            %         v6=zeros(300,1); % total spatial mean
            %         for aa=1:length(MyOcean_lon)
            %             for bb=1:length(MyOcean_lat)
            %                 num_v=zeros(length(MyOcean_lon),length(MyOcean_lat),300);
            %                 v4=zeros(length(MyOcean_lon),length(MyOcean_lat),300);
            %                 end_depth=zeros(length(MyOcean_lon),length(MyOcean_lat),300);
            %                 %Æ¯Á¤Áö¿ª ¼ö¿Â»Ì¾Æ¼­ Ç¥Ãþ(0m)¼ö¿Â Ãß°¡
            %                 point_temp = MyOcean_data(aa,bb,:);
            %                 point_temp2 = squeeze(point_temp);%1(Ç¥Ãþ)~50(½ÉÇØ)
            %                 if (~isnan(point_temp2(1))) % NaNÀÌ ¾Æ´Ï¶ó¸é
            %                     point_temp3 = [];
            %                     point_temp3(1) = point_temp2(1);
            %                     point_temp3(2:length(point_temp2)+1)=point_temp2(1:end);
            %                     %
            %                     a = find(isnan(point_temp3));
            %                     b = find(~isnan(point_temp3)); % NaNÀÌ ¾Æ´Ñ°÷ Ã£±â
            %                     std_dep = 0:1:fix(MyOcean_st3(b(end)));
            %                     point_temp4 = interp1(MyOcean_st3(1:b(end)),point_temp3(1:b(end)),std_dep,'spline');
            %                     v4(aa,bb,ind) = sum(point_temp4); % sum of vertical temp
            %                     num_v(aa,bb,ind) = length(point_temp4); % ÇÑ pointÀÇ ºÐÇÒ°³¼ö
            %                     end_depth(aa,bb,ind)=MyOcean_st3(b(end));
            %                     clear a b point_temp point_temp2 std_dep point_temp4
            %                 else % NaNÀÌ¶ó¸é
            %                     v4(aa,bb,ind)=NaN; % sum of vertical temp
            %                     num_v(aa,bb,ind)=NaN; % ÇÑ pointÀÇ ºÐÇÒ°³¼ö
            %                 end
            %             end
            %         end
            %         v5(:,:,ind) = v4(:,:,ind)./num_v(:,:,ind); % vertical mean of temp
            %         v6(ind) = sum(nansum(v4(:,:,ind)))/sum(nansum(num_v(:,:,ind))); % total spatial mean
            %         clear v4 num_v
                     ind=ind+1;
            %         tt1 = [tt1;tempyear];
            %         tt2 = [tt2;tempmonth];
                end
             end
            % comb_MyOcean_data = v5;
            % comb_MyOcean_temp_divided=reshape(comb_MyOcean_data,[len_lon, len_lat, 12, length(inputyear)]);
            
            %  save([filedir,regionname,'New_Myocean_temp','.mat'],'MyOcean_lon','MyOcean_lat',...
            %      'v5','v6','MyOcean_st','MyOcean_st3','comb_MyOcean_data','comb_MyOcean_temp_divided',...
            %      'len_lon','len_lat','len_st','lonlat','varname','tt1','tt2','regionname','end_depth')
            
%             save([filedir,regionname,'_New_Myocean_comb_raw','.mat'],'varname','regionname',...
%                 'MyOcean_lon','MyOcean_lat','MyOcean_st','comb_temp_raw_data')
%             
        end
        
        %% ÀúÀå
        % outpath='F:\¸ðµ¨_ÀÚ·á\MyOcean\';
        % outname = [outpath,'MyOcean_',poly_region,'_temp.nc'];
        % xlen = length(lon);ylen=length(lat);tlen=length(t);
        % nccreate(outname,'time','dimensions',{'time',tlen});
        % nccreate(outname,'longitude','dimensions',{'lon',xlen});
        % nccreate(outname,'latitude','dimensions',{'lat',ylen});
        % nccreate(outname,'vertical_mean','dimensions',{'lon',xlen,'lat',ylen,'time',tlen});
        % nccreate(outname,'total_mean','dimensions',{'time',tlen});
        % ncwrite(outname,'time',t);
        % ncwrite(outname,'longitude',lon);
        % ncwrite(outname,'latitude',lat);
        % ncwrite(outname,'vertical_mean',v5);
        % ncwrite(outname,'total_mean',v6);
