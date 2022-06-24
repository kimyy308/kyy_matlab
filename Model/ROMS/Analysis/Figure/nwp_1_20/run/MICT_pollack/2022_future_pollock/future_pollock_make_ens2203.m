close all; clear all; clc;

% %  Updated 05-May-2022 by Yong-Yub Kim
tic;

% RCM_info.years=2081:2100;
RCM_info.years=1995:2013;
RCM_info.days=1:99;
RCM_info.ensname_raw='ens2202';
RCM_info.ensname_wind='ens2204';
RCM_info.testnames=testnames_from_ens(RCM_info.ensname_raw);
dirs.rootdir='/home/kimyy/Model/ROMS/nwp_1_20/output';
RCM_grid.dl=1/20;

for yearij = 1:length(RCM_info.years)
    tmp.tempyear = RCM_info.years(yearij);
    tmp.tempyear_str= num2str(tmp.tempyear, '%04i');
    dirs.ensdir_raw = [dirs.rootdir, '/', RCM_info.ensname_raw, '/run/', tmp.tempyear_str];
    dirs.ensdir_wind = [dirs.rootdir, '/', RCM_info.ensname_wind, '/run/', tmp.tempyear_str];
    if (exist(dirs.ensdir_wind, 'dir') ~=7) mkdir(dirs.ensdir_wind); end
    
    for dayij = 1:length(RCM_info.days)
        tmp.tempday = RCM_info.days(dayij);
        tmp.tempday_str = num2str(tmp.tempday, '%04i');
        tmp.filename_ens_raw=[dirs.ensdir_raw, '/', ...
            'pck_', RCM_info.ensname_raw, '_', tmp.tempyear_str, '_daily_avg_', tmp.tempday_str, '.nc'];
        tmp.filename_ens_wind=[dirs.ensdir_wind, '/', ...
            'pck_', RCM_info.ensname_wind, '_', tmp.tempyear_str, '_daily_avg_', tmp.tempday_str, '.nc'];
        system(['cp ', tmp.filename_ens_raw, ' ', tmp.filename_ens_wind]);
        
        for testnameind=1:length(RCM_info.testnames)
            tmp.testname=RCM_info.testnames{testnameind};
            dirs.testdir = [dirs.rootdir, '/', tmp.testname, '/run/', tmp.tempyear_str];
            tmp.filename_test=[dirs.testdir, '/', ...
            'pck_', tmp.testname, '_', tmp.tempyear_str, '_daily_avg_', tmp.tempday_str, '.nc'];
%             ncinfo(tmp.filename_test)
            if isfield(RCM_grid, 'lon_rho')~=1
                RCM_grid.lon_rho=ncread(tmp.filename_test, 'lon_rho');
                RCM_grid.xlen=size(RCM_grid.lon_rho,1);
                RCM_grid.ylen=size(RCM_grid.lon_rho,2);
                RCM_grid.zlen=length(ncread(tmp.filename_test, 's_rho'));
            end
            RCM_data.Uwind(1:RCM_grid.xlen, 1:RCM_grid.ylen, testnameind)=ncread(tmp.filename_test, 'Uwind');
            RCM_data.Vwind(1:RCM_grid.xlen, 1:RCM_grid.ylen, testnameind)=ncread(tmp.filename_test, 'Vwind');
%             RCM_data.Ubar(1:RCM_grid.xlen, 1:RCM_grid.ylen, testnameind)=ncread(tmp.filename_test, 'ubar');
%             RCM_data.Vbar(1:RCM_grid.xlen, 1:RCM_grid.ylen, testnameind)=ncread(tmp.filename_test, 'vbar');
%             RCM_data.u(1:RCM_grid.xlen-1, 1:RCM_grid.ylen, :, testnameind)= ncread(tmp.filename_test, 'u');
%             RCM_data.v(1:RCM_grid.xlen, 1:RCM_grid.ylen-1, :,testnameind)= ncread(tmp.filename_test, 'v');
            RCM_data.u_rho(1:RCM_grid.xlen, 1:RCM_grid.ylen, 1:RCM_grid.zlen, testnameind) = ...
                u2rho_3d(ncread(tmp.filename_test, 'u'));
            RCM_data.v_rho(1:RCM_grid.xlen, 1:RCM_grid.ylen, 1:RCM_grid.zlen, testnameind) = ...
                v2rho_3d(ncread(tmp.filename_test, 'v'));
%             for loni=1:RCM_grid.xlen
%                 for lati=1:RCM_grid.ylen
%                     if isfinite(RCM_data.Uwind(loni,lati,testnameind))
%                         RCM_data.vel_wind=sqrt(RCM_data.Uwind(loni,lati,testnameind).^2 + ...
%                             RCM_data.Vwind(loni,lati,testnameind).^2);
%                         RCM_data.angle_wind= ...
%                             mod(atan2d(RCM_data.Vwind(loni,lati,testnameind), RCM_data.Uwind(loni,lati,testnameind));
%                     end
%                 end
%             end
            RCM_data.vel_wind= sqrt(RCM_data.Uwind.^2 + RCM_data.Vwind.^2);
            RCM_data.angle_wind=mod(atan2d(RCM_data.Vwind, RCM_data.Uwind),360);
%             RCM_data.vel_bar= sqrt(RCM_data.Ubar.^2 + RCM_data.Vbar.^2);
%             RCM_data.angle_bar=mod(atan2d(RCM_data.Vbar, RCM_data.Ubar),360);
            RCM_data.vel_flow= sqrt(RCM_data.u_rho.^2 + RCM_data.v_rho.^2);
            RCM_data.angle_flow=mod(atan2d(RCM_data.v_rho, RCM_data.u_rho),360);
        end
        RCM_data.ens_Uwind=mean(RCM_data.Uwind,3);
        RCM_data.ens_Vwind=mean(RCM_data.Vwind,3);        
        RCM_data.ens_vel_wind=mean(RCM_data.vel_wind,3);
        RCM_data.ens_angle_wind=mod(atan2d(RCM_data.ens_Vwind, RCM_data.ens_Uwind),360);
        RCM_data.recon_Uwind=RCM_data.ens_vel_wind .* cosd(RCM_data.ens_angle_wind);
        RCM_data.recon_Vwind=RCM_data.ens_vel_wind .* sind(RCM_data.ens_angle_wind);
        
%         RCM_data.ens_Ubar=mean(RCM_data.Ubar,3);
%         RCM_data.ens_Vbar=mean(RCM_data.Vbar,3);        
%         RCM_data.ens_vel_bar=mean(RCM_data.vel_bar,3);
%         RCM_data.ens_angle_bar=mod(atan2d(RCM_data.ens_Vbar, RCM_data.ens_Ubar),360);
%         RCM_data.recon_Ubar=RCM_data.ens_vel_bar .* cosd(RCM_data.ens_angle_bar);
%         RCM_data.recon_Vbar=RCM_data.ens_vel_bar .* sind(RCM_data.ens_angle_bar);
        
        RCM_data.ens_u_rho=mean(RCM_data.u_rho,4);
        RCM_data.ens_v_rho=mean(RCM_data.v_rho,4);
        RCM_data.ens_vel_flow=mean(RCM_data.vel_flow,4);
        RCM_data.ens_angle_flow=mod(atan2d(RCM_data.ens_v_rho, RCM_data.ens_u_rho),360);
        RCM_data.recon_u_rho=RCM_data.ens_vel_flow .* cosd(RCM_data.ens_angle_flow);
        RCM_data.recon_v_rho=RCM_data.ens_vel_flow .* sind(RCM_data.ens_angle_flow);
        RCM_data.recon_u=rho2u_3d_new(RCM_data.recon_u_rho);
        RCM_data.recon_v=rho2v_3d_new(RCM_data.recon_v_rho);

        ncwrite(tmp.filename_ens_wind, 'Uwind', RCM_data.recon_Uwind);
        ncwrite(tmp.filename_ens_wind, 'Vwind', RCM_data.recon_Vwind);
%         ncwrite(tmp.filename_ens_bar, 'ubar', RCM_data.recon_Ubar);
%         ncwrite(tmp.filename_ens_bar, 'vbar', RCM_data.recon_Vbar);
        ncwrite(tmp.filename_ens_wind, 'u', RCM_data.recon_u);
        ncwrite(tmp.filename_ens_wind, 'v', RCM_data.recon_v);
        disp([tmp.tempyear_str, 'year ,', tmp.tempday_str, 'day'])
    end
    toc;
end



% quiver(RCM_data.ens_Uwind(1:10:end, 1:10:end)', ...
%     RCM_data.ens_Vwind(1:10:end, 1:10:end)', 'autoscale', 'off')
% figure; 
% quiver(RCM_data.recon_Uwind(1:10:end, 1:10:end)', ...
%     RCM_data.recon_Vwind(1:10:end, 1:10:end)', 'autoscale', 'off')
% pcolor(RCM_data.recon_Vwind'); shading flat; colorbar;
% figure; pcolor(RCM_data.ens_Vwind'); shading flat; colorbar;

% u=1; v=1; % 45d
% u=0.0001; v=1; % 90d
% u=-1; v=-1; %225d

% mod(180+180/pi*atan2(v,u),360)  % for wind
% mod(atan2d(v,u),360)  % for flow

% /home/kimyy/Model/ROMS/nwp_1_20/output/ens2201/run/2097/pck_ens2201_2097_daily_avg_0058.nc

% example for flow
% u = V cos(pi)
% v = V sin(pi)
% % atan2d(v,u)
% mod(atan2d(-1,1), 360)  % 315
% mod(atan2d(-1,-1), 360)  % 225
% mod(atan2d(1,-1), 360)  % 135
% mod(atan2d(1,1), 360)  % 45
% mod(atan2d(0,1), 360)  % 0


function testnames=testnames_from_ens(ensname)
    switch(ensname)
        case{'ens2201'}
            testnames={'test2127', 'test2128', 'test2129', 'test2130', 'test2131'};
        case{'ens2202'}
            testnames={'test2117', 'test2118', 'test2119', 'test2120', 'test2121'};
    end
end

% var_u=ncread(tmp.filename_test, 'u');

function var_rho=u2rho_3d(var_u)
    [L,Mp,vert]=size(var_u);
    Lp=L+1;
    Lm=L-1;
    var_rho=zeros(Lp,Mp,vert);
    var_rho(2:L,:,:)=0.5*(var_u(1:Lm,:,:)+var_u(2:L,:,:));
    var_rho(1,:,:)=var_rho(2,:,:);
    var_rho(Lp,:,:)=var_rho(L,:,:);
end

function var_rho=v2rho_3d(var_v)
    [Lp,M,vert]=size(var_v);
    Mp=M+1;
    Mm=M-1;
    var_rho=zeros(Lp,Mp,vert);
    var_rho(:,2:M,:)=0.5*(var_v(:,1:Mm,:)+var_v(:,2:M,:));
    var_rho(:,1,:)=var_rho(:,2,:);
    var_rho(:,Mp,:)=var_rho(:,M,:);
end


function var_u=rho2u_3d_new(var_rho)
    [Lp,Mp,N]=size(var_rho);
    L=Lp-1;
    var_u=0.5*(var_rho(1:L,:,:)+var_rho(2:Lp,:,:));
end

function var_v=rho2v_3d_new(var_rho)
    [N,Mp,Lp]=size(var_rho);
    M=Mp-1;
    var_v=0.5*(var_rho(:,1:M,:)+var_rho(:,2:Mp,:));
end