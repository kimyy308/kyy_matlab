clear all; clc; close all;

warning off;
testname='test23';
var='vert_temp';
linux=1; windows=0;
if (windows ==1)
    % % for windows
    dropboxpath='C:\Users\KYY\Dropbox'; %% SNU_desktop
    addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
    addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
    addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run']));
    rootdir='C:\Users\kyy\Desktop\conv\';
    filedir=[rootdir,testname,'\'];
    outfile=[filedir,'fig\'];
    run('C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\convection\fig_param_kyy_conv');
elseif (linux==1)
    % % for linux
    dropboxpath='/home01/kimyy/Dropbox'; %% ROMS server
%     dropboxpath='/home/kimyy/Dropbox'; %% DAMO server
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
    addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_20/run']));
    
    rootdir='/data/kimyy/todgist/conv/output/';
    filedir=[rootdir,testname,'/'];
    outfile=[filedir,'fig/'];
    if (exist(outfile,'dir')~=7)
        mkdir(outfile);
    end
    run('/home01/kimyy/Dropbox/source/matlab/Model/ROMS/Analysis/Figure/convection/fig_param_kyy_conv.m');
end


inputminute=120:120:720;
time=inputminute;
% shadlev=[-1 -0.05];
shadlev=[-0.5 -0.025];
% shadlev=[-0.2 -0.01];
% shadlev=[-0.1 -0.005];

conlev=10:10;

tempminute=inputminute(1);
filename = strcat(filedir,'ocean_his_', num2str(tempminute,'%04i'), '.nc');
Vstretching = ncread(filename,'Vstretching')';
Vtransform = ncread(filename,'Vtransform')';
theta_s = ncread(filename,'theta_s')';
theta_b = ncread(filename,'theta_b')';
s_rho = ncread(filename,'s_rho')';
N=length(s_rho);
hc = ncread(filename,'hc')';

jet_conv=jet;
close all;
fig1= figure;
fig2= figure;
for minuteij=1:length(inputminute)
    ll=minuteij;
    tempminute = inputminute(minuteij);
    t_seq = tempminute/10;
    fontsize = 15;
    cold_area_limit = -inf;
    delta_Q_limit = inf;
    plot_type = 'heat_gain';
     
    figName1=strcat(outfile, testname, '_',plot_type,'_Area', '.png');
     
    figName2=strcat(outfile, testname, '_',plot_type, '.png');
    t_air=shadlev(1);
    fname = strcat(filedir,'ocean_his_', num2str(tempminute,'%04i'), '.nc');
    
    grd_file=filename;
    if (exist('grd.lon_rho' , 'var') ~= 1)
        grd.lon_rho = ncread(grd_file,'x_rho')';
        grd.lat_rho = ncread(grd_file,'y_rho')';
        grd.mask_rho = ncread(grd_file,'mask_rho')';
        grd.lon_u = ncread(grd_file,'x_u')';
        grd.lat_v = ncread(grd_file,'y_v')';
        grd.h = ncread(grd_file,'h')';
        grd.mask_rho_nan = grd.mask_rho;
        land = find(grd.mask_rho_nan==0);
        grd.mask_rho_nan(land) = NaN;

        h = grd.h;
%         kgrid=0;
%         [sc_r,Cs_r]=stretching(Vstretching, theta_s, theta_b, hc, N, kgrid);
%         kgrid=1;
%         [sc_w,Cs_w]=stretching(Vstretching, theta_s, theta_b, hc, N, kgrid);

        % zeta  
          
        grd.N = N;
        lon_rho  = grd.lon_rho;
        lat_rho  = grd.lat_rho;
        lon_u  = grd.lon_u;
        lat_v  = grd.lat_v; 
        mask_rho = grd.mask_rho;
        h = grd.h;
        N = grd.N;
        pi=3.141592;
        R=6392;
        for j=2:511
            for i=2:511
                dA(j,i)= (lon_u(j,i)-lon_u(j,i-1))*(lat_v(j,i)-lat_v(j-1,i));
            end
        end
        dA(1,2:511)=dA(2,2:511);
        dA(512,2:511)=dA(511,2:511);
        dA(2:511,1)=dA(2:511,2);
        dA(2:511,512)=dA(2:511,511);
        dA(1,1)=dA(2,2); dA(1,512)=dA(2,511); dA(512,1)=dA(511,2); dA(512,512)=dA(511,511);
        mdA=mean(mean(dA));
    end
    pt_eps = -1e-6;
    prho=1027.306;
    c_p=3.99e3;
    zeta = zeros(size(grd.h)); % default
      if (strcmp(grd_file,'E:\Data\Reanalysis\nwp_1_10_seo\avg_ens_10km_mean\roms_grid_final.nc')==1)
        grd.zeta = zeta;
      else
%         grd.zeta = ncread(grd_file,'zeta')'
        grd.zeta=zeros(size(h));
      end
%       grd.z_r=zlevs(Vtransform,Vstretching,h,grd.zeta,theta_s,theta_b,hc,N,'r');  
%     depth=grd.z_r;
    
    lon_min=1; lon_max=512; lat_min=256; lat_max=256;
    grd.z_r=zlevs(Vtransform,Vstretching,h(lat_min(1):lat_max(1), lon_min(1):lon_max(1)),grd.zeta(lat_min(1):lat_max(1), lon_min(1):lon_max(1)),theta_s,theta_b,hc,N,'r');  
%     depth=grd.z_r;
    data_info = ncinfo(filename, varname); 
    if (length(data_info.Dimensions)==4)
        data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 data_info.Size(3) 1]);  %% cut horizontal area [x,y,z] (wider than target area)
    else
        data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 data_info.Size(3)]);  %% cut horizontal area [x,y,z] (wider than target area)    end
    end
    data=squeeze(data);
    cut_data = permute(data, [3 2 1]);  %% permute [x y z] -> [z y x]
    cut_data=reshape(cut_data,[208,1,512]);
    
    cut_mask_rho = mask_rho(lat_min(1):lat_max(1), lon_min(1):lon_max(1));
%     cut_depth = depth(:,lat_min(1):lat_max(1), lon_min(1):lon_max(1));
    cut_depth = grd.z_r;
    z=cut_depth(:,1,1);
    ylimit = [min(z) 0];
    
    plot_idx = length(t_seq)+ll;
    t_idx=ll;
    pt=ncread(fname,varname,[1,1,1,1], [512,512,208,1]);
    
    cold_area = squeeze (sum( sum( (pt < pt_eps)    , 2 ), 1 ) ) *mdA;
    delta_Q   = squeeze (sum( sum( (pt < pt_eps).*pt, 2 ), 1 ) ) *mdA*prho*c_p;
    cold_idx = find(cold_area, 1);
    fprintf('z = %7.1f, area = %7.1f, dQ = %9.1f\n', z(cold_idx), cold_area(cold_idx), delta_Q(cold_idx));
    figure(fig1); plot( cold_area, z, 'linewidth',2 ); ax1=gca; hold on;
    figure(fig2); plot( delta_Q,   z, 'linewidth',2 ); ax2=gca; hold on;
    legend_str{ll} = ['t=', num2str(time(t_idx)/60, '%5.0f'), 'hr.'];
    cold_area_limit = max( max(cold_area(1:end-2)), cold_area_limit );
    delta_Q_limit   = min( min(delta_Q(1:end-2)), delta_Q_limit );
end
figure(fig1);
set(ax1,'fontsize',fontsize);
text(0.7,0.8,['T_n=',num2str(t_air),'\circC'],'units','normalized','edgecolor','k','fontsize',fontsize);    % if length(exp_num)==1
legend( legend_str, 'location','east' );
ylim(ylimit); xlim([0 cold_area_limit]);
ylabel('z (m)'); xlabel(['area of T<T_0 (m^2)']);
saveas( fig1,[figName1,'.png'],'png');
figure(fig2);
set(ax2,'fontsize',fontsize);
text(0.7,0.8,['T_n=',num2str(t_air),'\circC'],'units','normalized','edgecolor','k','fontsize',fontsize);    % if length(exp_num)==1
legend( legend_str, 'location','east' );
ylim(ylimit); xlim([delta_Q_limit 0]); set(gca,'xdir','reverse');
ylabel('z (m)'); xlabel('heat gain, c\rhoA\DeltaT (J/m)');
saveas( fig2,[figName2,'.png'],'png');