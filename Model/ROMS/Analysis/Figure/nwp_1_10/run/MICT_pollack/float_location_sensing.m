close all; clear all; clc;

system_name=computer;
if (strcmp(system_name,'PCWIN64'))
    % % for windows
    dropboxpath='C:\Users\KYY\Dropbox'; %% SNU_desktop, kyy_laptop
    addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
    addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
    addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run']));
elseif (strcmp(system_name,'GLNXA64'))
    dropboxpath='/HDD1/kimyy/Dropbox'; %% COCO
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
    addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_20/run']));
end

warning off;

run('nwp_polygon_point.m');

workdir='/HDD1/kimyy/Model/LTRANS/LTRANSv2b/LTRANSv2b/output/nwp_1_10_EnOI/';

totyear=1982:1991;
totmonth=1:3;
day=1;
totalday=90;

numfloat=96-32;


    for nyear=1:length(totyear)
        for nmonth =1:length(totmonth)
            year=totyear(nyear);
            month=totmonth(nmonth);
            if month==12 
                year=totyear(nyear)-1;
            end
            if (strcmp(system_name,'PCWIN64'))
                filename = ['E:\Data\Reanalysis\nwp_1_10_seo\avg_ens_10km_ens01\out_',num2str(year,'%04i'),'\output.nc'];
            elseif (strcmp(system_name,'GLNXA64'))
                filename = [workdir, ...
                    num2str(year,'%04i'),'/output_',num2str(year,'%04i'),'_',num2str(month,'%02i'),'_',num2str(day,'%02i'),'_',num2str(totalday,'%04i'),'.nc'];
            end
            lon = ncread(filename,'lon');
            lat = ncread(filename,'lat');
            if (strcmp(system_name,'PCWIN64'))
                gridname = ['E:\Data\Reanalysis\nwp_1_10_seo\avg_ens_10km_ens01\roms_grid02.nc'];
            elseif (strcmp(system_name,'GLNXA64'))
                gridname = ['/HDD1/kimyy/reanalysis_data/output/exp01/' ...
                    num2str(year,'%04i'),'/ocean_avg_0001.nc'];
            end
            mask_rho = ncread(gridname,'mask_rho');
            lon_rho = ncread(gridname,'lon_rho');
            lat_rho = ncread(gridname,'lat_rho');

    %         hold on;
            if (strcmp(system_name,'PCWIN64'))
                if (exist(['D:\OneDrive - ������б�\MEPL\project\MICT_pollack\3rd year\figure\ens_10km_ens01\drifter_LTRANS\',num2str(year,'%04i')] , 'dir') ~= 7)
                    mkdir(['D:\OneDrive - ������б�\MEPL\project\MICT_pollack\3rd year\figure\ens_10km_ens01\drifter_LTRANS\',num2str(year,'%04i')]);
                end 
                outfile = ['D:\OneDrive - ������б�\MEPL\project\MICT_pollack\3rd year\figure\ens_10km_ens01\drifter_LTRANS\',num2str(year,'%04i'),'\drifter'];
            elseif (strcmp(system_name,'GLNXA64'))
                if (exist([workdir,num2str(year,'%04i'),'/figures'] , 'dir') ~= 7)
                    mkdir([workdir,num2str(year,'%04i'),'/figures']);
                end 
                outfile = [workdir,num2str(year,'%04i'),'/figures','/drifter_',num2str(year,'%04i'),'_',num2str(month,'%02i'),'_',num2str(day,'%02i'),'_',num2str(totalday,'%04i')];
            end

            if (strcmp(system_name,'PCWIN64'))
                plot_dir = ['D:\OneDrive - ������б�\MEPL\project\MICT_pollack\3rd year\figure\ens_10km_ens01\drifter_LTRANS\',num2str(year,'%04i'),'\'];
            elseif (strcmp(system_name,'GLNXA64'))
                plot_dir = [workdir,num2str(year,'%04i'),'/figures','/'];
            end
            exp_name = 'pollack';
            k=90;
            jet96=jet(96);
            out_ind=0;
            t=1;
            while(t<totalday*24)
                wb_egg_in(nyear,nmonth,t)=0;
		nkc_egg_in(nyear,nmonth,t)=0;
		skc_egg_in(nyear,nmonth,t)=0;
		off_egg_in(nyear,nmonth,t)=0;
                for i=1:numfloat
                    wb_egg_location=double(inpolygon(lon(i,t),lat(i,t),wbpolygon(:,1),wbpolygon(:,2)));
                    wb_egg_in(nyear,nmonth,t)=wb_egg_in(nyear,nmonth,t)+wb_egg_location; 
                    nkc_egg_location=double(inpolygon(lon(i,t),lat(i,t),nkcpolygon(:,1),nkcpolygon(:,2)));
                    nkc_egg_in(nyear,nmonth,t)=nkc_egg_in(nyear,nmonth,t)+nkc_egg_location; 
                    skc_egg_location=double(inpolygon(lon(i,t),lat(i,t),skcpolygon(:,1),skcpolygon(:,2)));
                    skc_egg_in(nyear,nmonth,t)=skc_egg_in(nyear,nmonth,t)+skc_egg_location; 
                end
                off_egg_in(nyear,nmonth,t)=numfloat-wb_egg_in(nyear,nmonth,t)-nkc_egg_in(nyear,nmonth,t)-skc_egg_in(nyear,nmonth,t);
                t=t+1;
            end
	    disp(['year :', num2str(totyear(nyear)),' month :', num2str( totmonth(nmonth))]);
        end
    end
    save([workdir,'nwp_1_10_EnOI_comb_res_location','.mat'])
%end


