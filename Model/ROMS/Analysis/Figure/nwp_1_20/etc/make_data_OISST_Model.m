clear all; close all; clc;
windows=1;
if (windows ==1)
    % % for windows
    dropboxpath='C:\Users\KYY\Dropbox';
    addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
    addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
    addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
elseif (linux==1)
    dropboxpath='/home01/kimyy/Dropbox';
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
    addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
end



% % read model domain

name = 'monthly_spinup_';
filename_suffix = '.nc';
% ex : ~workdir\monthly_spinup_0001.nc
% filename = strcat('D:\MEPL\project\SSH\중간보고\smooth13_vtvs\',name,num2str(0001,'%04i'),filename_suffix);
filename = strcat('D:\MEPL\project\SSH\data\test24\',name,num2str(0049,'%04i'),filename_suffix); %% read grid

workdir = 'D:\MEPL\project\SSH\figure\avhrr';

lon = ncread(filename,'lon_rho');
lat = ncread(filename,'lat_rho');
% SST = ncread(filename,'temp',[1 1 40 1], [length(lon(:,1)) length(lat(1,:)) 1 1]);
% size(SST);
NWPdomain=[115 164 15 52];  %% nwp
ESdomain=[127 144 33 52];  %% east sea
  
    load D:\MEPL\project\SSH\중간보고\smooth13_vtvs\kyy_plot_subroutine\jet_mod
for month =1:12
modelname = ['D:\MEPL\project\SSH\data\test24\monthly_spinup_00',num2str(month+48,'%02i'), '.nc'];  %% read model data
oisstlon = ncread('D:\MEPL\project\SSH\data\OISST_monthly\avhrr-monthly_201312.nc','lon');  %% read oisst grid
oisstlat = ncread('D:\MEPL\project\SSH\data\OISST_monthly\avhrr-monthly_201312.nc','lat');  %% read oisst grid

% % for 2012. from cshwa
%     oiname = ['D:\MEPL\project\SSH\data\OISST_monthly\AVHRR_2012\2012', num2str(month,'%02i'), 'month'];
%     load (oiname)
%     oisst=monthly_mean;

% % for climate data, from MEPL Ocean server(stlee)
    oiname = ['D:\MEPL\project\SSH\data\OISST_monthly\OISST_MonthClim_1982-2011.nc'];
    oisst_scale_factor=0.01;
%     oisst(:,:) = ncread(oiname,'interpolated_sst', [1 1 month], [length(oisstlon), length(oisstlat) 1]) * oisst_scale_factor; 
    oisst(:,:) = ncread(oiname,'interpolated_sst', [1 1 month], [length(oisstlon), length(oisstlat) 1]); 
    for i=1:length(oisstlon)
        for j =1:length(oisstlat)
            if (oisst(i,j)<-9)
                oisst(i,j)=NaN;
            end
        end
    end
    interp_OISST = griddata(double(oisstlon), double(oisstlat), oisst',lon(:,1),lat(1,:));  %% 920 * 980 (transposed)
    for i=1:920
        for j =1:980
            if (interp_OISST(i,j)<-9)
                interp_OISST(i,j)=NaN;
            end
        end
    end
    comb_OISST(:,:,month)=interp_OISST(:,:)';
    comb_modelSST(:,:,month)=ncread(modelname,'temp',[1 1 40 1], [length(lon(:,1)) length(lat(1,:)) 1 1]);
    comb_modelSSS(:,:,month)=ncread(modelname,'salt',[1 1 40 1], [length(lon(:,1)) length(lat(1,:)) 1 1]);
    comb_modelSSH(:,:,month)=ncread(modelname,'zeta',[1 1 1], [length(lon(:,1)) length(lat(1,:)) 1]);
    month
end
save('D:\MEPL\project\SSH\data\OISST_monthly\avhrr_n_model_SST.mat','comb_modelSST','comb_OISST','comb_modelSSS','comb_modelSSH','lon','lat','NWPdomain','ESdomain');