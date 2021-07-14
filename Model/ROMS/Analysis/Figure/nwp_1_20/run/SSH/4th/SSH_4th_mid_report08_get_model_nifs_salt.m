close all; clear all;  clc;   

system_name=computer;
if (strcmp(system_name,'PCWIN64'))
    % % for windows
    dropboxpath='C:\Users\KYY\Dropbox';
    addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
    addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
%     addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
    addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path
elseif (strcmp(system_name,'GLNXA64'))
    dropboxpath='/home/kimyy/Dropbox';
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
    addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
end

shadlev = [0 35];
% rms_shadlev = [0 4];
%     trendlev = [-3 3];  %% trend lev
trendlev = [-5 5];  %% trend lev
% conlev  = 0:5:35;
meanplotlev =[-0.2 0.2];

% for snu_desktopd
% all_testname2 = {'ens01', 'ens02', 'ens03'};
all_testname2 = {'test53'};

for testnameind2=1:length(all_testname2)

testname=all_testname2{testnameind2}  

% testname='test49'   % % need to change
inputyear = [1976:2005]; % % put year which you want to plot [year year ...]
% inputyear = [2006:2005]; % % put year which you want to plot [year year ...]

inputmonth = [2 4 6 8 10 12]; % % put month which you want to plot [month month ...]
varname ='salt'
run('nwp_polygon_point.m');

load(['E:\Data\Observation\NIFS\xls\NIFS_salt_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);

NIFS6.lat(NIFS6.lat==0)=NaN;
NIFS6.lon(NIFS6.lon==0)=NaN;
maxlat=max(max(max(max(max(NIFS6.lat)))));
minlat=min(min(min(min(min(NIFS6.lat)))));
maxlon=max(max(max(max(max(NIFS6.lon)))));
minlon=min(min(min(min(min(NIFS6.lon)))));

refpolygon=[maxlon,maxlat; ...
            maxlon, minlat; ...
            minlon, minlat; ...
            minlon,maxlat];
lonlat(1)=min(refpolygon(:,1));
        lonlat(2)=max(refpolygon(:,1));
        lonlat(3)=min(refpolygon(:,2));
        lonlat(4)=max(refpolygon(:,2));
        
        
system_name=computer;
if (strcmp(system_name,'PCWIN64'))
% % for windows
    figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\4th_year\figure\nwp_1_20\',testname,'\'); % % where figure files will be saved
    param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m'
    filedir = strcat('I:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
elseif (strcmp(system_name,'GLNXA64'))
end 
for yearij = 1:length(inputyear)
    for monthij = 1:length(inputmonth)
        disp([num2str(yearij), 'y_',num2str(monthij),'m'])
        tic;
        tempyear = inputyear(yearij);
        tempmonth = inputmonth(monthij);
        % ex : rootdir\test37\data\2001\test37_monthly_2001_01.nc
        filename = strcat(filedir, num2str(tempyear,'%04i'), '\', ...
                testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');
        % read model data
        if (exist('lon')==0)
            modelinfo=ncinfo(filename);
            lon = ncread(filename,'lon_rho',[1 1],[modelinfo.Dimensions(5).Length,1]);
            lat = ncread(filename,'lat_rho',[1 1],[1,modelinfo.Dimensions(6).Length]);

            lon_west = abs(lon - (lonlat(1)-1));
            min_lon_west=min(lon_west);
            lon_east = abs(lon - (lonlat(2)+1));
            min_lon_east=min(lon_east);
            lat_south = abs(lat - (lonlat(3)-1));
            min_lat_south=min(lat_south);
            lat_north = abs(lat - (lonlat(4)+1));
            min_lat_north=min(lat_north);

            lon_min = find(lon_west == min_lon_west);
            lon_max = find(lon_east == min_lon_east);
            lat_min = find(lat_south == min_lat_south);
            lat_max = find(lat_north == min_lat_north);

            lon = ncread(filename,'lon_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
            lat = ncread(filename,'lat_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);

            mask_model = double(inpolygon(lon,lat,refpolygon(:,1),refpolygon(:,2)));
            mask_model(mask_model==0)=NaN;
            
            model_depth = ncread(filename,'h',[lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
            
            Vtransform= ncread(filename,'Vtransform');
            Vstretching= ncread(filename,'Vstretching');
            theta_s= ncread(filename,'theta_s');
            theta_b= ncread(filename,'theta_b');
            hc= ncread(filename,'hc');
            N = modelinfo.Dimensions(2).Length;
        end
        zeta = ncread(filename,'zeta',[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
        model_depth_z=zlevs(Vtransform,Vstretching,model_depth,zeta,theta_s,theta_b,hc,N,'r');
        model_depth_z=permute(model_depth_z,[2, 3, 1]);
        data_info = ncinfo(filename, varname);  %% [lon lat depth time] -> [1601 1201 33 1]

        data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 N 1]);
%         data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
        data=data.*mask_model;
        
        stddepth=unique(NIFS6.depth(~isnan(NIFS6.depth)));
        for lonij=1:size(data,1)
           for latij=1:size(data,2)
              temp_stddepth=stddepth(stddepth<=max(-model_depth_z(lonij,latij,end:-1:1)));
               if (isnan(data(lonij,latij,40))==1 || length(temp_stddepth)==1)
                   data_zdepth(lonij,latij,:)=NaN;
               else
                   %% depth --> 0 10 ...
                   model_depth_z(lonij,latij,end)=0;
                   data_zdepth(lonij,latij,1:length(temp_stddepth)) = ...
                       interp1( squeeze(-model_depth_z(lonij,latij,end:-1:1)), ...
                       squeeze(data(lonij,latij,end:-1:1)), ...
                       temp_stddepth);
               end
           end
        end
        data_zdepth(data_zdepth==0)=NaN;
        
        for numline=1:size(NIFS6.line,2)
            for numsta=1:size(NIFS6.line,3)
%                 for numdepth=1:size(NIFS6.line,5)
                    if (isnan(NIFS6.lon(yearij,numline,numsta,monthij,1))==1)
                       MODEL6.salt(yearij,numline,numsta,monthij,:)=NaN;
                    else
%                         MODEL6.salt(yearij,numline,numsta,monthij,numdepth) = griddata(lon,lat, ...
%                         squeeze(data_zdepth(:,:,numdepth)), ...
%                         squeeze(NIFS6.lon(yearij,numline,numsta,monthij,numdepth)), ...
%                         squeeze(NIFS6.lat(yearij,numline,numsta,monthij,numdepth)));

                        lon_sta_diff = abs(lon - NIFS6.lon(yearij,numline,numsta,monthij,1));
%                         min_lon_sta_diff=min(min(lon_sta_diff));
%                         lon_sta_min = find(lon_sta_diff == min_lon_sta_diff,1);
                        
                        lat_sta_diff = abs(lat - NIFS6.lat(yearij,numline,numsta,monthij,1));
%                         min_lat_sta_diff=min(lat_sta_diff);
%                         lat_sta_min = find(lat_sta_diff == min_lat_sta_diff,1);
                        
                        sta_diff=lon_sta_diff + lat_sta_diff;
                        sta_min = find(sta_diff == min(min(sta_diff)),1);
                        lon_sta_min=mod(sta_min,size(lon,1));
                        lat_sta_min=floor(sta_min/size(lon,1));
                        
                        MODEL6.salt(yearij,numline,numsta,monthij,1:length(stddepth))=squeeze(data_zdepth(lon_sta_min,lat_sta_min,:));
                    end
%                 end
            end
        end

        toc;
    end
end

for numyear=1:size(MODEL6.salt,1)
    for numline=1:size(MODEL6.salt,2)
        for nummon=1:size(MODEL6.salt,4)
            for numsta=1:size(MODEL6.salt,3)
                MODEL6.salt75_sta(numyear,numline,numsta,nummon)=mean(MODEL6.salt(numyear,numline,numsta,nummon,1:6),5,'omitnan');
                MODEL6.salt_surf_sta(numyear,numline,numsta,nummon)=MODEL6.salt(numyear,numline,numsta,nummon,1);
            end
            MODEL6.salt75_line(numyear,numline,nummon)=mean(MODEL6.salt75_sta(numyear,numline,:,nummon),3,'omitnan');
            MODEL6.salt_surf_line(numyear,numline,nummon)=mean(MODEL6.salt_surf_sta(numyear,numline,:,nummon),3,'omitnan');
        end
    end
end

MODEL6.salt75_sta(MODEL6.salt75_sta==0)=NaN;
MODEL6.salt_surf_sta(MODEL6.salt_surf_sta==0)=NaN;
MODEL6.salt75_line=squeeze(mean(MODEL6.salt75_sta,3,'omitnan'));
MODEL6.salt_surf_line=squeeze(mean(MODEL6.salt_surf_sta,3,'omitnan'));

numi=1;
for i=1:size(MODEL6.salt75_line,1)  %% year
    for j=1:size(MODEL6.salt75_line,3)  %% month
%         MODEL6.salt75_line_comb=reshape(MODEL6.salt75_line,size(MODEL6.salt75_line,1)*size(MODEL6.salt75_line,3),size(MODEL6.salt75_line,2));
        MODEL6.salt75_line_comb(numi,:)=MODEL6.salt75_line(i,:,j);
        MODEL6.salt_surf_line_comb(numi,:)=MODEL6.salt_surf_line(i,:,j);
        numi=numi+1;
    end
end
MODEL6.salt75_line_comb(MODEL6.salt75_line_comb==0)=NaN;

% plot(squeeze(MODEL6.salt75_line_comb(:,3)))  %%206 line

for numline=1:size(MODEL6.salt75_line,2)
    ttt=1/6:1/6:length(MODEL6.salt75_line_comb(:,numline))/6;
    idx=isnan(MODEL6.salt75_line_comb(:,numline)');
    p=polyfit(ttt(~idx),squeeze(MODEL6.salt75_line_comb(~idx,numline))',1);
    MODEL6.salt75_line_comb_trend(numline)=p(1);
    idx=isnan(MODEL6.salt_surf_line_comb(:,numline)');
    p=polyfit(ttt(~idx),squeeze(MODEL6.salt_surf_line_comb(~idx,numline))',1);
    MODEL6.salt_surf_line_comb_trend(numline)=p(1);
end
MODEL6.salt75_line_comb_trend(MODEL6.salt75_line_comb_trend==0)=NaN;
MODEL6.salt_surf_line_comb_trend(MODEL6.salt_surf_line_comb_trend==0)=NaN;

% % % [year line sta mon dep]
for i=1:size(MODEL6.salt,2)  %% line
    for j=1:size(MODEL6.salt,3)  %% sta
        numi=1;
        for k=1:size(MODEL6.salt,1) %% year
            for l=1:size(MODEL6.salt,4) %% month
        %         MODEL6.salt75_line_comb=reshape(MODEL6.salt75_line,size(MODEL6.salt75_line,1)*size(MODEL6.salt75_line,3),size(MODEL6.salt75_line,2));
                MODEL6.salt_surf_comb(numi,i,j)=MODEL6.salt(k,i,j,l,1);
                MODEL6.salt75_comb(numi,i,j)=MODEL6.salt75_sta(k,i,j,l);
                numi=numi+1;
            end
        end
    end
end
MODEL6.salt_surf_comb(MODEL6.salt_surf_comb==0)=NaN;

for i=1:size(MODEL6.salt,2)  %% line
    for j=1:size(MODEL6.salt,3)  %% sta
        idx=isnan(MODEL6.salt_surf_comb(:,i,j)');
        p=polyfit(ttt(~idx),squeeze(MODEL6.salt_surf_comb(~idx,i,j))',1);
        MODEL6.salt_surf_comb_trend(i,j)=p(1);
    end
end
MODEL6.salt_surf_comb_trend(MODEL6.salt_surf_comb_trend==0)=NaN;

for i=1:size(MODEL6.salt,2)  %% line
    for j=1:size(MODEL6.salt,3)  %% sta
        idx=isnan(MODEL6.salt75_comb(:,i,j)');
        p=polyfit(ttt(~idx),squeeze(MODEL6.salt75_comb(~idx,i,j))',1);
        MODEL6.salt75_comb_trend(i,j)=p(1);
    end
end
MODEL6.salt75_comb_trend(MODEL6.salt75_comb_trend==0)=NaN;


save(['E:\Data\Observation\NIFS\xls\',testname,'_nwp_1_20_salt_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat'], 'MODEL6', 'NIFS6');

end

% % --> dot figure

% plot(squeeze(MODEL6.salt_surf_line_comb(:,22)))
% hold on;
% plot(squeeze(NIFS6.salt_surf_line_comb(:,22)))
% hold off;
% ylim([25 35]);
% 
% 
% %% ss :15 , YS: 18 & 22, KS : 10, ES : 3
% plot(squeeze(MODEL6.salt75_line_comb(:,22)))
% hold on;
% plot(squeeze(NIFS6.salt75_line_comb(:,22)))
% hold off;
% ylim([25 35]);
