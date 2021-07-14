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

workdir='E:\Data\Model\LTRANS\output\test06_DA_3rd\';

% scalename={'e_folding', 'half', '2_3', 'tenth'};
scalename={'e_folding'};



totyear=1983:2019;
totmonth=1:3;
day=1;
totalday=90;

refpolygon(1,1)=127;
refpolygon(2,1)=128.3;
refpolygon(1,2)=38.4;
refpolygon(2,2)=40.0;
numfloat=96-32;

% refpolygon(1,1)=127;
% refpolygon(2,1)=128.6;
% refpolygon(1,2)=38.4;
% refpolygon(2,2)=40.0;
% numfloat=96;

e_fold_scale=numfloat * 1/exp(1);  %% numfloat * 1 / 2.7183 = 35.3164 (35)
half_scale=numfloat * 1/2;  %% numfloat * 1 / 2.7183 = 35.3164 (35)
two_over_three_scale=numfloat * 2/3;  %% numfloat * 1 / 2.7183 = 35.3164 (35)
tenth_scale=numfloat*1/10;
time_scale(1)=e_fold_scale;
time_scale(2)=half_scale;
time_scale(3)=two_over_three_scale;
time_scale(4)=tenth_scale;

for scaleind=1:length(scalename)
    tempscale=time_scale(scaleind);
    
    for nyear=1:length(totyear)
        for nmonth =1:length(totmonth)
            year=totyear(nyear);
            month=totmonth(nmonth);
            yearstr=num2str(year,'%04i');
            if month==12 
                year=totyear(nyear)-1;
            end
            if (strcmp(system_name,'PCWIN64'))
                filestr=[num2str(totalday,'%04i'),'d_',yearstr,num2str(nmonth,'%02i'),'01'];
                filename = ['E:\Data\Model\LTRANS\output\test06_DA_3rd\',yearstr,'\',...
                    filestr,'\output_',filestr,'.nc'];
            elseif (strcmp(system_name,'GLNXA64'))
                filename = [workdir, ...
                    num2str(year,'%04i'),'/output_',num2str(year,'%04i'),'_',num2str(month,'%02i'),'_',num2str(day,'%02i'),'_',num2str(totalday,'%04i'),'.nc'];
            end
            lon = ncread(filename,'lon');
            lat = ncread(filename,'lat');
            if (strcmp(system_name,'PCWIN64'))
                gridname = ['E:\Data\Model\ROMS\nwp_1_10\test06\DA\1983\test06_monthly_1983_01.nc'];
            elseif (strcmp(system_name,'GLNXA64'))
                gridname = ['/HDD1/kimyy/reanalysis_data/output/exp01/' ...
                    num2str(year,'%04i'),'/ocean_avg_0001.nc'];
            end
            mask_rho = ncread(gridname,'mask_rho');
            lon_rho = ncread(gridname,'lon_rho');
            lat_rho = ncread(gridname,'lat_rho');

    %         hold on;
%             if (strcmp(system_name,'PCWIN64'))
%                 if (exist(['D:\OneDrive - ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Ð±ï¿?\MEPL\project\MICT_pollack\3rd year\figure\ens_10km_ens01\drifter_LTRANS\',num2str(year,'%04i')] , 'dir') ~= 7)
%                     mkdir(['D:\OneDrive - ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Ð±ï¿?\MEPL\project\MICT_pollack\3rd year\figure\ens_10km_ens01\drifter_LTRANS\',num2str(year,'%04i')]);
%                 end 
%                 outfile = ['D:\OneDrive - ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Ð±ï¿?\MEPL\project\MICT_pollack\3rd year\figure\ens_10km_ens01\drifter_LTRANS\',num2str(year,'%04i'),'\drifter'];
%             elseif (strcmp(system_name,'GLNXA64'))
%                 if (exist([workdir,num2str(year,'%04i'),'/figures'] , 'dir') ~= 7)
%                     mkdir([workdir,num2str(year,'%04i'),'/figures']);
%                 end 
%                 outfile = [workdir,num2str(year,'%04i'),'/figures','/drifter_',num2str(year,'%04i'),'_',num2str(month,'%02i'),'_',num2str(day,'%02i'),'_',num2str(totalday,'%04i')];
%             end
            % pcolor
            lonlat = [128 131 37 41];

            if (strcmp(system_name,'PCWIN64'))
                plot_dir = ['D:\OneDrive - ¼­¿ï´ëÇÐ±³\MEPL\project\MICT_pollack\4th_year\Figure\nwp_1_10_EnOI\LTRANS\res_time\'];
            elseif (strcmp(system_name,'GLNXA64'))
                plot_dir = [workdir,num2str(year,'%04i'),'/figures','/'];
            end
            exp_name = 'pollack';
            k=90;
            jet96=jet(96);
            out_ind=0;
            t=25;
            while(t<(totalday+1)*24)
    %             hold on;
                egg_in=0;
                for i=1:numfloat
                    egg_location=double(inpolygon(lon(i,t),lat(i,t),refpolygon(:,1),refpolygon(:,2)));
                    egg_in=egg_in+egg_location;
                end
    %             if (egg_in<e_fold_scale)
                if (egg_in<tempscale)
                    out_ind=out_ind+1;
                else
                    out_ind=0;
                end
                if (out_ind>72)
                    disp(t/24)
                    comb_res_time(nyear,nmonth)=(t-72.0)/24.0;
    %                 return;
                    t=(totalday+1)*24;
                else
                   if (t==(totalday*24)-1)
                        disp((t+1)/24)
                        comb_res_time(nyear,nmonth)=(t+1)/24.0;
                   end
                   t=t+1;
                end

    %             disp(t)

    %          hold off;
            end
        end
    end
    save(['nwp_1_10_EnOI_comb_res_time_',scalename{scaleind},'_',num2str(min(totyear)),'_',num2str(max(totyear)),'.mat'])
end


