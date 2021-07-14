function grd = read_grid(grd_file,Vtransform,Vstretching,theta_s,theta_b,hc,N)

grd.lon_rho = ncread(grd_file,'lon_rho')';
grd.lat_rho = ncread(grd_file,'lat_rho')';
grd.mask_rho = ncread(grd_file,'mask_rho')';
grd.angle = ncread(grd_file,'angle')';
grd.h = ncread(grd_file,'h')';
grd.lon_psi = ncread(grd_file,'lon_psi')';
grd.lat_psi = ncread(grd_file,'lat_psi')';
grd.mask_psi = ncread(grd_file,'mask_psi')';
grd.lon_v = ncread(grd_file,'lon_v')';
grd.lat_v = ncread(grd_file,'lat_v')';
grd.mask_v = ncread(grd_file,'mask_v')';
grd.lon_u = ncread(grd_file,'lon_u')';
grd.lat_u = ncread(grd_file,'lat_u')';
grd.mask_u = ncread(grd_file,'mask_u')';
grd.pm = ncread(grd_file,'pm')';
grd.pn = ncread(grd_file,'pn')';

grd.mask_rho_nan = grd.mask_rho;
land = find(grd.mask_rho_nan==0);
grd.mask_rho_nan(land) = NaN;

h = grd.h;
kgrid=0;
[sc_r,Cs_r]=stretching(Vstretching, theta_s, theta_b, hc, N, kgrid);
kgrid=1;
[sc_w,Cs_w]=stretching(Vstretching, theta_s, theta_b, hc, N, kgrid);

% zeta  
  zeta = zeros(size(grd.h)); % default
  if (strcmp(grd_file,'E:\Data\Reanalysis\nwp_1_10_seo\avg_ens_10km_mean\roms_grid_final.nc')==1)
    grd.zeta = zeta;
  elseif (strcmp(grd_file,'E:\Data\Model\ROMS\es_1_40\input\etopo1_Eastsea_40.nc')==1)
    grd.zeta = zeta;
  elseif (strcmp(grd_file,'/data2/kimyy/Reanalysis/nwp_1_10_seo/ens_mean_monthly_mod/roms_grid_final.nc')==1)
    grd.zeta = zeta;
  elseif (strcmp(grd_file,'/data2/kimyy/Reanalysis/nwp_1_10_seo/ens_mean_monthly/roms_grid_final.nc')==1)
    grd.zeta = zeta;
  else
    grd.zeta = ncread(grd_file,'zeta')'
  end
  grd.z_r=zlevs(Vtransform,Vstretching,h,grd.zeta,theta_s,theta_b,hc,N,'r');  
  grd.z_w=zlevs(Vtransform,Vstretching,h,grd.zeta,theta_s,theta_b,hc,N+1,'w');  
  
  grd.theta_s = theta_s;
  grd.theta_b = theta_b;
  grd.Tcline = hc;
  grd.N = N;
  grd.hc = hc;
  grd.sc_w = sc_w;
  grd.Cs_w = Cs_w;
  grd.sc_r = sc_r;
  grd.Cs_r = Cs_r;
end