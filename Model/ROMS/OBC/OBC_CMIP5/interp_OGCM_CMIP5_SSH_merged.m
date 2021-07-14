function interp_OGCM_SODA3_SSH_merged(OGCM_dir,OGCM_prefix,year,month,Roa,interp_method,...
    lonU,latU,lonV,latV,lonT,latT,tin,...
    nc_clm,nc_bry,lon,lat,angle,h,tout,obc)
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
% Updated 14-Jun-2018 by Y.Y.Kim


conserv = 0; % same barotropic velocities as the OGCM,
% % if switch is on, barotropic velocity is obtained by divided volume transport by maximum depth --> error



%
disp(['  Horizontal interpolation: ',...
    'CMIP5_CMEMS_SSH - ','Y',num2str(year),' - M',num2str(month),'.nc'])
%
%
% ROMS grid angle
%
cosa = cos(angle);
sina = sin(angle);
%
% Open the OGCM file
%
nc = netcdf([OGCM_dir,OGCM_prefix,'.nc']);
% missvalue = nc{'ssh'}.missing_value(:);
missvalue = NaN;

%
% Interpole data on the OGCM Z grid and ROMS horizontal grid
%
%
% Read and extrapole the 2D variables
%
zeta = ext_data_OGCM_CMIP5_total(nc,lonT,latT,'cmems_sla',tin,lon,lat,1,missvalue,Roa,interp_method);  %(cm -> m)

%
% Close the OGCM file (CSEOF reconstructed SSH file)
%
close(nc)

%
%Initialisation in case of bry files
%

zeta_south_soda =  nc_bry{'zeta_south'}(tout,:);
zeta_north_soda =  nc_bry{'zeta_north'}(tout,:);
zeta_west_soda =  nc_bry{'zeta_west'}(tout,:);
zeta_east_soda =  nc_bry{'zeta_east'}(tout,:);


if ~isempty(nc_bry)
    if obc(1)==1
        zeta_south = squeeze(zeta(1,:));
    end
    if obc(2)==1
        zeta_east = squeeze(zeta(:,end));
    end
    if obc(3)==1
        zeta_north = squeeze(zeta(end,:));
    end
    if obc(4)==1
        zeta_west = squeeze(zeta(:,1));
    end
end





%
% Boundary file
%
if ~isempty(nc_bry)
    if obc(1) == 1
        nc_bry{'zeta_south'}(tout,:) = zeta_south + zeta_south_soda;
    end
    if obc(2) == 1
        nc_bry{'zeta_east'}(tout,:) = zeta_east + zeta_east_soda';
    end
    if obc(3) == 1
        nc_bry{'zeta_north'}(tout,:) = zeta_north + zeta_north_soda;
    end
    if obc(4) == 1
        nc_bry{'zeta_west'}(tout,:) = zeta_west + zeta_west_soda';
    end
end