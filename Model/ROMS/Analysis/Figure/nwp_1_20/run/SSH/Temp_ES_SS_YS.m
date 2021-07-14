close all; clear all;  clc;
% all_region ={'NWP','ES', 'SS', 'YS', 'ECS'}
all_region ={'ES','SS','YS'}

for regionind=1:length(all_region)
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
    tt=[];
    tt1=[];
    tt2=[];
    ind=1;
    %% ÇØ´ç À§°æµµ¸¦ µû¿Í¼­ ¿øÇÏ´Â µ¥ÀÌÅÍ¸¦ ÃßÃâ, ÀüÃ¼ ±â°£(³â,¿ù)À» ÀÏ·Ä·Î ÀúÀå
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
                comb_spatial_meanMyOcean=(zeros([length(MyOcean_lon),length(MyOcean_lat),12]));
                %comb_spatial_meanmodel=(zeros([length(MyOcean_lon),length(MyOcean_lat),12]));
                % À§°æµµ ¸Þ½¬±×¸®µå ¸¸µé±â (µ¥ÀÌÅÍ Â÷¿ø¿¡ µû¶ó ´Ù¸§)
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

            
            MyOcean_data = ncread(MyOceanfilename,varname, ...
                [MyOcean_lon_min(1) MyOcean_lat_min(1) 1 1], ...
                [MyOcean_lon_max(1)-MyOcean_lon_min(1) MyOcean_lat_max(1)-MyOcean_lat_min(1) 50 1]);
            MyOcean_data(MyOcean_data<-1000)=NaN;
            MyOcean_data(MyOcean_data>1000)=NaN;
            MyOcean_data=MyOcean_data.*mask_MyOcean;
               %% interpolationÇØ¼­ Æò±Õ ±¸ÇØ¾ß
            % ¿©±â¼­ Æò±ÕÀ» ±¸ÇØ¾ß ÇÒµí
            % °¢ ¿ùº° °ÝÀÚ³» ÇØ¼ö¸é ¿ÂµµÀÇ ¼ö½É¿¡ µû¸¥ Æò±Õ
            point_temp3=[];
            for aa=1:length(MyOcean_lon)
                for bb=1:length(MyOcean_lat)
                    point_temp = MyOcean_data(aa,bb,:);
                    point_temp2 = squeeze(point_temp);%1(Ç¥Ãþ)~50(½ÉÇØ)
                    list_dep = find(isnan(point_temp2)==0);
                    % interpolation(Æ¯Á¤ ³âµµ, Æ¯Á¤ ¿ù, Æ¯Á¤ Àå¼Ò)
                    if length(list_dep)==0
                        end_depth=0;
                        mean_MyOcean=NaN;
                    else
                        i_dep = list_dep(end);
                        end_depth = round(MyOcean_st(i_dep),-1)+10;
                    % ÇÏ³ªÇÏ³ª interpolation
                    dep= 0:5:end_depth;
                    interp_MyOcean_data=interp1(MyOcean_st(1:i_dep),point_temp2(1:i_dep),dep,'spline');
                    MyOcean_data2 = interp_MyOcean_data;
                    mean_MyOcean = nanmean(MyOcean_data2);
                    end
                    mean_MyOcean_data(aa,bb)= mean_MyOcean;
                end
            end
                comb_MyOcean_data(:,:,ind) = mean_MyOcean_data;
                %
                %comb_spatial_meanMyOcean(:,:,monthij)=comb_spatial_meanMyOcean(:,:,monthij)+MyOcean_data/double(length(inputyear));???
                %%
                ind = ind + 1;
                %             toc;
                tt1 = [tt1;tempyear];
                tt2 = [tt2;tempmonth];
            end
    end
     comb_MyOcean_temp_divided=reshape(comb_MyOcean_data,[len_lon, len_lat, 12, length(inputyear)]);


%% unit ÁöÁ¤
    if strcmp(varname,'zos')==1
        unit = 'mm/year';
    elseif strcmp(varname,'thetao')==1
        unit = '\circC/year';
    elseif strcmp(varname,'temp')==1
        unit = '\circC/year';
    else
        disp('unit error')
        pause;
    end
    %
%%
        yy=tt1;
        mm=tt2;
   %% ¿ùº° ¼ö¿Â Æ®·»µå
        data3=[];
        month_MyOcean_temp_divided = (zeros([length(MyOcean_lon),length(MyOcean_lat),12]));
        for aa=1:length(MyOcean_lon)
            for bb=1:length(MyOcean_lat)
                for cc=1:12
                    i_mon = find(tt2==cc);
                    yy2 = yy(i_mon);
                    mm2 = mm(i_mon);
                    tt3 = datenum(yy2,mm2,15,0,0,0);
                    data2 = comb_MyOcean_temp_divided(aa,bb,cc,:); %¿Âµµ
                    data = squeeze(data2);                    
                    [r,m,d]=regression(double(tt3'),double(data'));%tt¿Ídata°¡ ¼¼·ÎÇàÇüÅÂ·Î µé¾î°¡¾ßÇÔ
                    %             period=(tt3(end)-tt3(1))/365;
                    %             height=m*tt3+d;
                    %             var=height(end)-height(1);
                    %             rr=(var/period);
                    %             rr=round(rr,3);
                    month_MyOcean_temp_divided(aa,bb,cc)=m*365; % Æ¯Á¤À§°æµµ¿ùÀÇ ¼ö¿ÂÆ®·»µå
                    clearvars rr data2 data
                end
            end
        end
         month_MyOcean_temp_divided(month_MyOcean_temp_divided==0)=NaN;

  save([filedir,'Temp_Myocean_',regionname,'.mat'],'MyOcean_lon','MyOcean_lat','month_MyOcean_temp_divided',...
    'comb_MyOcean_temp_divided','comb_MyOcean_data','len_lon','len_lat','len_st','lonlat','varname',...
    'tt1','tt2','unit','regionname');
       
end