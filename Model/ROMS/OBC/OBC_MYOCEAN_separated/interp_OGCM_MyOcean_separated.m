function interp_OGCM_MyOcean_north(OGCM_dir,OGCM_prefix,year,month,Roa,interp_method,...
    rmdepth,tin,...
    nc_bry,lon,lon_south, lon_north, lon_west, lon_east, ...
    lat,         angle,h,tout,Vtransform,Vstretching, conserv)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Read the local OGCM files and perform the interpolations
%
% Ok, I am lazy and I did not do something special for the bry files...
%
%
%  Further Information:
%  http://www.brest.ird.fr/Roms_tools/
%
%  This file is part of ROMSTOOLS
%
%  ROMSTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  ROMSTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Copyright (c) 2005-2006 by Pierrick Penven
%  e-mail:Pierrick.Penven@ird.fr
%
%  Updated    6-Sep-2006 by Pierrick Penven : Nothing special for the bry file
%  Update    13-Sep-2009 by Gildas Cambon :   Begin treatments case  for the bry
%  file, no to be continued ...
%  Updated    5-Nov-2006 by Pierrick Penven : A bit of cleaning...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updated 20-Feb-2019 by Y.Y.Kim

grid_name = [OGCM_path, OGCM_prefix, num2str(Ymin), '.nc'];
nc = netcdf(grid_name);

% % MyOcean
LAT = nc{'latitude'}(:);
LON = nc{'longitude'}(:);
Z = -nc{'depth'}(:);

lonT = LON; latT = LAT;
lonU = LON; latU = LAT;
lonV = LON; latV = LAT;

NZ = length(Z);
NZ = NZ - rmdepth;
Z = Z(1:NZ);
close(nc)


%
disp(['  Horizontal interpolation: ',...
    'MyOcean ','Y',num2str(year),'M',num2str(month),'.nc'])
%
%
% ROMS grid angle
%
cosa = cos(angle);
sina = sin(angle);
%
% Open the OGCM file
%
nc = netcdf([OGCM_dir,OGCM_prefix,num2str(year),'.nc']);
missvalue = -32767;
%
% Interpole data on the OGCM Z grid and ROMS horizontal grid
%
% Read and extrapole the 2D variables
%
zeta = ext_data_OGCM_MyOcean(nc,lonT,latT,'ssh',tin,lon,lat,1,missvalue,Roa,interp_method);
%
% Read and extrapole the 3D variables
%
NZ = length(Z);
[M,L] = size(lon);
dz = gradient(Z);
temp = zeros(NZ,M,L);
salt = zeros(NZ,M,L);
u = zeros(NZ,M,L-1);
v = zeros(NZ,M-1,L);

% % % set ROMS grid
lon_u=rho2u_2d(lon);
lon_v=rho2u_2d(lon);
lat_u=rho2u_2d(lat);
lat_v=rho2u_2d(lat);

lon_south=squeeze(lon(1,:));
lon_east=squeeze(lon(:,end));
lon_north=squeeze(lon(end,:));
lon_west=squeeze(lon(:,1));

lat_south=squeeze(lat(1,:));
lat_east=squeeze(lat(:,end));
lat_north=squeeze(lat(end,:));
lat_west=squeeze(lat(:,1));

lon_u_south=squeeze(lon_u(1,:));
lon_u_east=squeeze(lon_u(:,end));
lon_u_north=squeeze(lon_u(end,:));
lon_u_west=squeeze(lon_u(:,1));

lat_u_south=squeeze(lat_u(1,:));
lat_u_east=squeeze(lat_u(:,end));
lat_u_north=squeeze(lat_u(end,:));
lat_u_west=squeeze(lat_u(:,1));

lon_v_south=squeeze(lon_v(1,:));
lon_v_east=squeeze(lon_v(:,end));
lon_v_north=squeeze(lon_v(end,:));
lon_v_west=squeeze(lon_v(:,1));

lat_v_south=squeeze(lat_v(1,:));
lat_v_east=squeeze(lat_v(:,end));
lat_v_north=squeeze(lat_v(end,:));
lat_v_west=squeeze(lat_v(:,1));

dzu_3d_north = repmat(-dz, [1, L-1]);
dzu_3d_south = dzu_3d_north;
dzu_3d_east = repmat(-dz, [1, M]);
dzu_3d_west = dzu_3d_east;

dzv_3d_north = repmat(-dz, [1, L]);
dzv_3d_south = dzv_3d_north;
dzv_3d_east = repmat(-dz, [1, M-1]);
dzv_3d_west = dzv_3d_east;


for k = 1:NZ
    if rem(k,10) == 0
        disp(['  Level ',num2str(k),' of ',num2str(NZ)])
    end
    u2d = ext_data_OGCM_MyOcean(nc,lonU,latU,'u',tin,lon,lat,...
        k,missvalue,Roa,interp_method);
    v2d = ext_data_OGCM_MyOcean(nc,lonV,latV,'v',tin,lon,lat,...
        k,missvalue,Roa,interp_method);
    u(k,:,:) = rho2u_2d(u2d.*cosa+v2d.*sina);
    v(k,:,:) = rho2v_2d(v2d.*cosa-u2d.*sina);
    temp(k,:,:) = ext_data_OGCM_MyOcean(nc,lonT,latT,'temp',tin,lon,lat,...
        k,missvalue,Roa,interp_method);
    salt(k,:,:) = ext_data_OGCM_MyOcean(nc,lonT,latT,'salt',tin,lon,lat,...
        k,missvalue,Roa,interp_method); 
end

%
% Close the OGCM file
%
close(nc)
%
% Calculate Ubar and Vbar
%
dzu_3d = repmat(-dz, [1, M, L-1]);
dzv_3d = repmat(-dz, [1, M-1, L]);

ubar = squeeze(sum(u.*dzu_3d)/sum(-dz));
vbar = squeeze(sum(v.*dzv_3d)/sum(-dz));
% %
% %Initialisation in case of bry files
% %
if ~isempty(nc_bry)
    if obc(1)==1
        zeta_south = squeeze(ext_data_OGCM_MyOcean(nc,lonT,latT,'ssh',tin,lon_south,lat_south,1,missvalue,Roa,interp_method,1));
        u_south_temp = squeeze(ext_data_OGCM_MyOcean(nc,lonU,latU,'u',tin,lon_u_south,lat_u_south,1,missvalue,Roa,interp_method,1));
        v_south_temp = squeeze(ext_data_OGCM_MyOcean(nc,lonV,latV,'v',tin,lon_v_south,lat_v_south,1,missvalue,Roa,interp_method,1));
        u_south = u_south_temp.*cosa + v_south_temp.*sina;
        v_south = v_south_temp.*cosa - u_south_temp.*sina;
        ubar_south = squeeze(sum(u_south.*dzu_3d_south)/sum(-dz));
        vbar_south = squeeze(sum(v_south.*dzv_3d_south)/sum(-dz));
        temp_south = squeeze(ext_data_OGCM_MyOcean(nc,lonT,latT,'temp',tin,lon_south,lat_south,1,missvalue,Roa,interp_method,1));
        salt_south = squeeze(ext_data_OGCM_MyOcean(nc,lonT,latT,'salt',tin,lon_south,lat_south,1,missvalue,Roa,interp_method,1));
    end
    if obc(2)==1
        zeta_east = squeeze(ext_data_OGCM_MyOcean(nc,lonT,latT,'ssh',tin,lon_east,lat_east,1,missvalue,Roa,interp_method,2));
        u_east_temp = squeeze(ext_data_OGCM_MyOcean(nc,lonU,latU,'u',tin,lon_u_east,lat_u_east,1,missvalue,Roa,interp_method,2));
        v_east_temp = squeeze(ext_data_OGCM_MyOcean(nc,lonV,latV,'v',tin,lon_v_east,lat_v_east,1,missvalue,Roa,interp_method,2));
        u_east = u_east_temp.*cosa + v_east_temp.*sina;
        v_east = v_east_temp.*cosa - u_east_temp.*sina;
        ubar_east = squeeze(sum(u_east.*dzu_3d_east)/sum(-dz));
        vbar_east = squeeze(sum(v_east.*dzv_3d_east)/sum(-dz));
        temp_east = squeeze(ext_data_OGCM_MyOcean(nc,lonT,latT,'temp',tin,lon_east,lat_east,1,missvalue,Roa,interp_method,2));
        salt_east = squeeze(ext_data_OGCM_MyOcean(nc,lonT,latT,'salt',tin,lon_east,lat_east,1,missvalue,Roa,interp_method,2));
    end
    if obc(3)==1
        zeta_north = squeeze(ext_data_OGCM_MyOcean(nc,lonT,latT,'ssh',tin,lon_north,lat_north,1,missvalue,Roa,interp_method,3));
        u_north_temp = squeeze(ext_data_OGCM_MyOcean(nc,lonU,latU,'u',tin,lon_u_north,lat_u_north,1,missvalue,Roa,interp_method,3));
        v_north_temp = squeeze(ext_data_OGCM_MyOcean(nc,lonV,latV,'v',tin,lon_v_north,lat_v_north,1,missvalue,Roa,interp_method,3));
        u_north = u_north_temp.*cosa + v_north_temp.*sina;
        v_north = v_north_temp.*cosa - u_north_temp.*sina;
        ubar_north = squeeze(sum(u_north.*dzu_3d_north)/sum(-dz));
        vbar_north = squeeze(sum(v_north.*dzv_3d_north)/sum(-dz));
        temp_north = squeeze(ext_data_OGCM_MyOcean(nc,lonT,latT,'temp',tin,lon_north,lat_north,1,missvalue,Roa,interp_method,3));
        salt_north = squeeze(ext_data_OGCM_MyOcean(nc,lonT,latT,'salt',tin,lon_north,lat_north,1,missvalue,Roa,interp_method,3));
    end
    if obc(4)==1
        zeta_west = squeeze(ext_data_OGCM_MyOcean(nc,lonT,latT,'ssh',tin,lon_west,lat_west,1,missvalue,Roa,interp_method,4));
        u_west_temp = squeeze(ext_data_OGCM_MyOcean(nc,lonU,latU,'u',tin,lon_u_west,lat_u_west,1,missvalue,Roa,interp_method,4));
        v_west_temp = squeeze(ext_data_OGCM_MyOcean(nc,lonV,latV,'v',tin,lon_v_west,lat_v_west,1,missvalue,Roa,interp_method,4));
        u_west = u_west_temp.*cosa + v_west_temp.*sina;
        v_west = v_west_temp.*cosa - u_west_temp.*sina;
        ubar_west = squeeze(sum(u_west.*dzu_3d_west)/sum(-dz));
        vbar_west = squeeze(sum(v_west.*dzv_3d_west)/sum(-dz));
        temp_west = squeeze(ext_data_OGCM_MyOcean(nc,lonT,latT,'temp',tin,lon_west,lat_west,1,missvalue,Roa,interp_method,4));
        salt_west = squeeze(ext_data_OGCM_MyOcean(nc,lonT,latT,'salt',tin,lon_west,lat_west,1,missvalue,Roa,interp_method,4));
        
    end
end

%
% Get the ROMS vertical grid
%
disp('  Vertical interpolations')
if ~isempty(nc_clm)
    theta_s = nc_clm{'theta_s'}(:);
    theta_b = nc_clm{'theta_b'}(:);
    hc = nc_clm{'hc'}(:);
    N = length(nc_clm('s_rho'));
end
if ~isempty(nc_bry)
    theta_s = nc_bry{'theta_s'}(:);
    theta_b = nc_bry{'theta_b'}(:);
    hc = nc_bry{'hc'}(:);
    N = length(nc_bry('s_rho'));
end
%
% Add an extra bottom layer (-100000m) and an extra surface layer (+100m)
% to prevent vertical extrapolations
%
Z = [100;Z;-100000];
%
% ROMS vertical grid
%
zr = zlevs(Vtransform,Vstretching,h,zeta,theta_s,theta_b,hc,N,'r');
zu = rho2u_3d(zr);
zv = rho2v_3d(zr);
zw = zlevs(Vtransform,Vstretching,h,zeta,theta_s,theta_b,hc,N,'w');
dzr = zw(2:end,:,:)-zw(1:end-1,:,:);
dzu = rho2u_3d(dzr);
dzv = rho2v_3d(dzr);

%
%
% Vertical interpolation in case of bry files
%
%
if ~isempty(nc_bry)
    %
    %South
    %
    if obc(1) == 1
        [u_south,v_south,...
            ubar_south,vbar_south,...
            temp_south,salt_south]=vinterp_OGCM_MyOcean(zr(:,1,:),zu(:,1,:),zv(:,1,:),...
            dzr(:,1,:),dzu(:,1,:),dzv(:,1,:),...
            u_south,v_south,...
            ubar_south,vbar_south,...
            temp_south,salt_south,...
            N,Z,conserv);
    end
    if obc(2) == 1
        [u_east,v_east,...
            ubar_east,vbar_east,...
            temp_east,salt_east]=vinterp_OGCM_MyOcean(zr(:,:,end),zu(:,:,end),zv(:,:,end),...
            dzr(:,:,end),dzu(:,:,end),dzv(:,:,end),...
            u_east,v_east,...
            ubar_east,vbar_east,...
            temp_east,salt_east,...
            N,Z,conserv);
    end
    if obc(3) == 1
        [u_north,v_north,...
            ubar_north,vbar_north,...
            temp_north,salt_north]=vinterp_OGCM_MyOcean(zr(:,end,:),zu(:,end,:),zv(:,end,:),...
            dzr(:,end,:),dzu(:,end,:),dzv(:,end,:),...
            u_north,v_north,...
            ubar_north,vbar_north,...
            temp_north,salt_north,...
            N,Z,conserv);
    end
    if obc(4) == 1
        [u_west,v_west,...
            ubar_west,vbar_west,...
            temp_west,salt_west]=vinterp_OGCM_MyOcean(zr(:,:,1),zu(:,:,1),zv(:,:,1),...
            dzr(:,:,1),dzu(:,:,1),dzv(:,:,1),...
            u_west,v_west,...
            ubar_west,vbar_west,...
            temp_west,salt_west,...
            N,Z,conserv);
    end
end   %~isempty(nc_bry)
%--------------------------------------------------------------

%
%  fill the %  Boundary files
%
if ~isempty(nc_bry)
    if obc(1) == 1
        nc_bry{'zeta_south'}(tout,:) = zeta_south;
        nc_bry{'temp_south'}(tout,:,:) = temp_south;
        nc_bry{'salt_south'}(tout,:,:) = salt_south;
        nc_bry{'u_south'}(tout,:,:) = u_south;
        nc_bry{'v_south'}(tout,:,:) = v_south;
        nc_bry{'ubar_south'}(tout,:,:) = ubar_south;
        nc_bry{'vbar_south'}(tout,:,:) = vbar_south;
    end
    if obc(2) == 1
        nc_bry{'zeta_east'}(tout,:) = zeta_east;
        nc_bry{'temp_east'}(tout,:,:) = temp_east;
        nc_bry{'salt_east'}(tout,:,:) = salt_east;
        nc_bry{'u_east'}(tout,:,:) = u_east;
        nc_bry{'v_east'}(tout,:,:) = v_east;
        nc_bry{'ubar_east'}(tout,:,:) = ubar_east;
        nc_bry{'vbar_east'}(tout,:,:) = vbar_east;
    end
    if obc(3) == 1
        nc_bry{'zeta_north'}(tout,:) = zeta_north;
        nc_bry{'temp_north'}(tout,:,:) = temp_north;
        nc_bry{'salt_north'}(tout,:,:) = salt_north;
        nc_bry{'u_north'}(tout,:,:) = u_north;
        nc_bry{'v_north'}(tout,:,:) = v_north;
        nc_bry{'ubar_north'}(tout,:,:) = ubar_north;
        nc_bry{'vbar_north'}(tout,:,:) = vbar_north;
    end
    if obc(4) == 1
        nc_bry{'zeta_west'}(tout,:) = zeta_west;
        nc_bry{'temp_west'}(tout,:,:) = temp_west;
        nc_bry{'salt_west'}(tout,:,:) = salt_west;
        nc_bry{'u_west'}(tout,:,:) = u_west;
        nc_bry{'v_west'}(tout,:,:) = v_west;
        nc_bry{'ubar_west'}(tout,:,:) = ubar_west;
        nc_bry{'vbar_west'}(tout,:,:) = vbar_west;
    end
end