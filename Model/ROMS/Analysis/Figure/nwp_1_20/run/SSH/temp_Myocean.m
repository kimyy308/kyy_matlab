%%



close all; clear all;  clc;
all_region ={'AKP2'}

 regionind=1%:length(all_region)
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
    
    
    varname ='depth';
    run('nwp_polygon_point.m');
    regionname=all_region{regionind}
    switch(regionname) % if¹®°ú ºñ½ÁÇÏ°Ô Á¶°Ç¿¡ ¸Â´Â ÇØ´ç ÄÉÀÌ½º¿¡ ´ëÇØ ½ÇÇà
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
        figrawdir =strcat('D:\OneDrive - ������б�\MEPL\project\MICT_pollack\3rd year\figure\nwp_1_20\',testname,'\',regionname,'\'); % % where figure files will be saved
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
                
                %         ncinfo('E:\Data\Observation\OIssh\monthly\MyOcean_monthly1983_11.nc');
                
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
                        if end_depth>=50
                            MyOcean_data2 = interp_MyOcean_data(11);
                            mean_MyOcean = nanmean(MyOcean_data2);
                        else
                            mean_MyOcean=NaN;
                        end
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
        % 51x40x12x38 size·Î ºÐ·ùµÊ
        %% ¿©±â´Â °øÅë
        
        yy=tt1;
        mm=tt2;
        %tt3 = datenum(yy,mm,15,0,0,0);

        
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

    
% for mon=1:12
%     figure
%     lonlat=[115 145 30 52];
%     val_min=-0.15; val_max=0.15; val_con=0.05; % color
%     m_proj('miller','lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
%     hold on;
%     m_gshhs_c('patch',[.7 .7 .7],'edgecolor','none');
%     m_grid('linestyle','none','Tickdir','out', 'xtick',lonlat(1):5:lonlat(2),...
%         'ytick',lonlat(3):3:lonlat(4), 'fontsize', 8, 'XaxisLocation', 'bottom');
%     colormap('jet'); caxis([val_min val_max]);
%     m_pcolor(MyOcean_lon,MyOcean_lat,month_MyOcean_temp_divided(:,:,mon)')
%     h=colorbar;  title(h,'\circC/year');
%     %shading flat
%     titlename =  ['Month Mean Temp Trend(5m) ',num2str(mon),'¿ù']  ;
%     title(titlename,'fontsize',15)
%     
%     savepath='C:\Users\½ÅÇÑ¼Ö\Desktop\Á¹¾÷³í¹®\figure\';
%     savename = strcat('month_temp_trend(5m interp)','_',num2str(mon),'¿ù');
%     print([savepath,savename],'-r200','-dtiffn')
%     close;
% end
% %     
%     save([filedir,regionname,'Myocean_temp(5m)_surface','.mat'],'MyOcean_lon','MyOcean_lat','month_MyOcean_temp_divided',...
%     'comb_MyOcean_temp_divided','comb_MyOcean_data','len_lon','len_lat','len_st','lonlat','MyOcean_data',...
%     'tt1','tt2','tt3','yy','yy2','mm','mm2');
    