

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
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run']));
elseif (linux==1)
    % % for linux
    dropboxpath='/home01/kimyy/Dropbox'; %% ROMS server
%     dropboxpath='/home/kimyy/Dropbox'; %% DAMO server
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
    addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_20/run']));
end


testname='test39';
inputdir = strcat('E:\Data\Model\ROMS\nwp_1_20\');
workdir =[inputdir, 'input\', testname, '\'];

load C:\Users\KYY\Dropbox\source\matlab\Common\Figure\jet_mod
% % name = 'roms_grid_combine2_';
name = 'roms_grid_nwp_1_20_';
filename_suffix = '.nc';
% ex : ~workdir\roms_grid_combine2_test37.nc
filename = strcat(workdir,name,testname,filename_suffix);

grd.lon_rho = ncread(filename,'lon_rho')';
grd.lat_rho = ncread(filename,'lat_rho')';
grd.mask_rho = ncread(filename,'mask_rho')';
grd.h = ncread(filename,'h')';


Vtransform = 2;
Vstretching = 4;
theta_s = 10;
theta_b =1;
hc = 5;
N = 40;

h = grd.h(495,281:381);  %   129E ~ 132E, 37N 

% zeta  
grd.zeta = zeros(size(h)); % default




% % %   grd.zeta = zeta;


% for theta_s = 10:10
%     for theta_b = 0:0.1:0.1
%         for hc = [5 100 250]
            grd.z_r=squeeze(zlevs(h,grd.zeta,theta_s,theta_b,hc,N,'r'));  
            grd.z_w=squeeze(zlevs(h,grd.zeta,theta_s,theta_b,hc,N+1,'w'));  

            close all;
            hold on
            for i=1:40
            %     contour(squeeze(grd.lon_rho(495,281:381)),0:-3000/39:-3000,-squeeze(grd.z_r(i,495,281:381)), 'k')
                plot(grd.lon_rho(495,281:381),squeeze(grd.z_w(i,:)), 'k')
            end
            ylim([-2500 0])
            a=area(grd.lon_rho(495,281:381),-grd.h(495,281:381),-2500, 'FaceColor',[0.8 0.8 0.8]);
            title(strcat('ts=',num2str(theta_s),' tb=',num2str(theta_b), ' hc=', num2str(hc)));
            hold off

% % % %             jpgname=strcat('D:\MEPL\project\SSH\3rd_year\figure\test37\coordinates\vertical\ES\', testname, '_','ts_',num2str(theta_s),'_tb_',num2str(theta_b),'_hc_',num2str(hc), '.jpg'); %% ~_year_month.jpg
           jpgname=strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\3rd_year\figure\test39\coordinates\vertical\ES\', testname, '_','ts_',num2str(theta_s),'_tb_',num2str(theta_b),'_hc_',num2str(hc), '.jpg'); %% ~_year_month.jpg

            saveas(gcf,jpgname,'jpg');
%         end
%     end
% end
