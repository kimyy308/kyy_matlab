function  create_grid_Y(L,M,grdname,title)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Create an empty netcdf gridfile
%       L: total number of psi points in x direction  
%       M: total number of psi points in y direction  
%       grdname: name of the grid file
%       title: title in the netcdf file  
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
%  Copyright (c) 2001-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%  Updated 11-May-2018 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lp=L+1;
Mp=M+1;

ncid = netcdf.create(grdname, 'CLOBBER');
% ncid = netcdf.create(grdname, 'NETCDF4');

%
%  Create dimensions
%

xi_u_dimid = netcdf.defDim(ncid, 'xi_u', L);
xi_v_dimid = netcdf.defDim(ncid, 'xi_v', Lp);
xi_rho_dimid = netcdf.defDim(ncid, 'xi_rho', Lp);
xi_psi_dimid = netcdf.defDim(ncid, 'xi_psi', L);
eta_u_dimid = netcdf.defDim(ncid, 'eta_u', Mp);
eta_v_dimid = netcdf.defDim(ncid, 'eta_v', M);
eta_rho_dimid = netcdf.defDim(ncid, 'eta_rho', Mp);
eta_psi_dimid = netcdf.defDim(ncid, 'eta_psi', M);
one_dimid = netcdf.defDim(ncid, 'one', 1);
two_dimid = netcdf.defDim(ncid, 'two', 2);
four_dimid = netcdf.defDim(ncid, 'four', 4);
bath_dimid = netcdf.defDim(ncid, 'bath', 1);


%
%  Create variables and attributes
%
xl_varid = netcdf.defVar(ncid, 'xl', 'double', [one_dimid]);
netcdf.putAtt(ncid, xl_varid, 'long_name', 'domain length in the XL-direction');
netcdf.putAtt(ncid, xl_varid, 'units', 'meter');

el_varid = netcdf.defVar(ncid, 'el', 'double', [one_dimid]);
netcdf.putAtt(ncid, el_varid, 'long_name', 'domain length in the ETA-direction');
netcdf.putAtt(ncid, el_varid, 'units', 'meter');

depthmin_varid = netcdf.defVar(ncid, 'depthmin', 'double', [one_dimid]);
netcdf.putAtt(ncid, depthmin_varid, 'long_name', 'Shallow bathymetry clipping depth');
netcdf.putAtt(ncid, depthmin_varid, 'units', 'meter');

depthmax_varid = netcdf.defVar(ncid, 'depthmax', 'double', [one_dimid]);
netcdf.putAtt(ncid, depthmax_varid, 'long_name', 'Deep bathymetry clipping depth');
netcdf.putAtt(ncid, depthmax_varid, 'units', 'meter');

spherical_varid = netcdf.defVar(ncid, 'spherical', 'char', [one_dimid]);
netcdf.putAtt(ncid, spherical_varid, 'long_name', 'Grid type logical switch');
netcdf.putAtt(ncid, spherical_varid, 'units', 'spherical');

angle_varid = netcdf.defVar(ncid, 'angle', 'double', [xi_rho_dimid, eta_rho_dimid]);
netcdf.putAtt(ncid, angle_varid, 'long_name', 'angle between xi axis and east');
netcdf.putAtt(ncid, angle_varid, 'units', 'degree');

h_varid = netcdf.defVar(ncid, 'h', 'double', [xi_rho_dimid, eta_rho_dimid]);
netcdf.putAtt(ncid, h_varid, 'long_name', 'Final bathymetry at RHO-points');
netcdf.putAtt(ncid, h_varid, 'units', 'meter');

hraw_varid = netcdf.defVar(ncid, 'hraw', 'double', [xi_rho_dimid, eta_rho_dimid bath_dimid]);
netcdf.putAtt(ncid, hraw_varid, 'long_name', 'Working bathymetry at RHO-points');
netcdf.putAtt(ncid, hraw_varid, 'units', 'meter');

alpha_varid = netcdf.defVar(ncid, 'alpha', 'double', [xi_rho_dimid, eta_rho_dimid]);
netcdf.putAtt(ncid, alpha_varid, 'long_name', 'Weights between coarse and fine grids at RHO-points');

f_varid = netcdf.defVar(ncid, 'f', 'double', [xi_rho_dimid, eta_rho_dimid]);
netcdf.putAtt(ncid, f_varid, 'long_name', 'Coriolis parameter at RHO-points');
netcdf.putAtt(ncid, f_varid, 'units', 'second-1');

pm_varid = netcdf.defVar(ncid, 'pm', 'double', [xi_rho_dimid, eta_rho_dimid]);
netcdf.putAtt(ncid, pm_varid, 'long_name', 'curvilinear coordinate metric in XI');
netcdf.putAtt(ncid, pm_varid, 'units', 'meter-1');

pn_varid = netcdf.defVar(ncid, 'pn', 'double', [xi_rho_dimid, eta_rho_dimid]);
netcdf.putAtt(ncid, pn_varid, 'long_name', 'curvilinear coordinate metric in ETA');
netcdf.putAtt(ncid, pn_varid, 'units', 'meter-1');

dndx_varid = netcdf.defVar(ncid, 'dndx', 'double', [xi_rho_dimid, eta_rho_dimid]);
netcdf.putAtt(ncid, dndx_varid, 'long_name', 'xi derivative of inverse metric factor pn');
netcdf.putAtt(ncid, dndx_varid, 'units', 'meter');

dmde_varid = netcdf.defVar(ncid, 'dmde', 'double', [xi_rho_dimid, eta_rho_dimid]);
netcdf.putAtt(ncid, dmde_varid, 'long_name', 'eta derivative of inverse metric factor pm');
netcdf.putAtt(ncid, dmde_varid, 'units', 'meter');

x_rho_varid = netcdf.defVar(ncid, 'x_rho', 'double', [xi_rho_dimid, eta_rho_dimid]);
netcdf.putAtt(ncid, x_rho_varid, 'long_name', 'x location of RHO-points');
netcdf.putAtt(ncid, x_rho_varid, 'units', 'meter');

x_u_varid = netcdf.defVar(ncid, 'x_u', 'double', [xi_u_dimid, eta_u_dimid]);
netcdf.putAtt(ncid, x_u_varid, 'long_name', 'x location of U-points');
netcdf.putAtt(ncid, x_u_varid, 'units', 'meter');

x_v_varid = netcdf.defVar(ncid, 'x_v', 'double', [xi_v_dimid, eta_v_dimid]);
netcdf.putAtt(ncid, x_v_varid, 'long_name', 'x location of V-points');
netcdf.putAtt(ncid, x_v_varid, 'units', 'meter');

x_psi_varid = netcdf.defVar(ncid, 'x_psi', 'double', [xi_psi_dimid, eta_psi_dimid]);
netcdf.putAtt(ncid, x_psi_varid, 'long_name', 'x location of PSI-points');
netcdf.putAtt(ncid, x_psi_varid, 'units', 'meter');

y_rho_varid = netcdf.defVar(ncid, 'y_rho', 'double', [xi_rho_dimid, eta_rho_dimid]);
netcdf.putAtt(ncid, y_rho_varid, 'long_name', 'y location of RHO-points');
netcdf.putAtt(ncid, y_rho_varid, 'units', 'meter');

y_u_varid = netcdf.defVar(ncid, 'y_u', 'double', [xi_u_dimid, eta_u_dimid]);
netcdf.putAtt(ncid, y_u_varid, 'long_name', 'y location of U-points');
netcdf.putAtt(ncid, y_u_varid, 'units', 'meter');

y_v_varid = netcdf.defVar(ncid, 'y_v', 'double', [xi_v_dimid, eta_v_dimid]);
netcdf.putAtt(ncid, y_v_varid, 'long_name', 'y location of V-points');
netcdf.putAtt(ncid, y_v_varid, 'units', 'meter');

y_psi_varid = netcdf.defVar(ncid, 'y_psi', 'double', [xi_psi_dimid, eta_psi_dimid]);
netcdf.putAtt(ncid, y_psi_varid, 'long_name', 'y location of PSI-points');
netcdf.putAtt(ncid, y_psi_varid, 'units', 'meter');

lon_rho_varid = netcdf.defVar(ncid, 'lon_rho', 'double', [xi_rho_dimid, eta_rho_dimid]);
netcdf.putAtt(ncid, lon_rho_varid, 'long_name', 'longitude of RHO-points');
netcdf.putAtt(ncid, lon_rho_varid, 'units', 'degree_east');

lon_u_varid = netcdf.defVar(ncid, 'lon_u', 'double', [xi_u_dimid, eta_u_dimid]);
netcdf.putAtt(ncid, lon_u_varid, 'long_name', 'longitude of U-points');
netcdf.putAtt(ncid, lon_u_varid, 'units', 'degree_east');

lon_v_varid = netcdf.defVar(ncid, 'lon_v', 'double', [xi_v_dimid, eta_v_dimid]);
netcdf.putAtt(ncid, lon_v_varid, 'long_name', 'longitude of V-points');
netcdf.putAtt(ncid, lon_v_varid, 'units', 'degree_east');

lon_psi_varid = netcdf.defVar(ncid, 'lon_psi', 'double', [xi_psi_dimid, eta_psi_dimid]);
netcdf.putAtt(ncid, lon_psi_varid, 'long_name', 'longitude of PSI-points');
netcdf.putAtt(ncid, lon_psi_varid, 'units', 'degree_east');

lat_rho_varid = netcdf.defVar(ncid, 'lat_rho', 'double', [xi_rho_dimid, eta_rho_dimid]);
netcdf.putAtt(ncid, lat_rho_varid, 'long_name', 'latitude of RHO-points');
netcdf.putAtt(ncid, lat_rho_varid, 'units', 'degree_north');

lat_u_varid = netcdf.defVar(ncid, 'lat_u', 'double', [xi_u_dimid, eta_u_dimid]);
netcdf.putAtt(ncid, lat_u_varid, 'long_name', 'latitude of U-points');
netcdf.putAtt(ncid, lat_u_varid, 'units', 'degree_north');

lat_v_varid = netcdf.defVar(ncid, 'lat_v', 'double', [xi_v_dimid, eta_v_dimid]);
netcdf.putAtt(ncid, lat_v_varid, 'long_name', 'latitude of V-points');
netcdf.putAtt(ncid, lat_v_varid, 'units', 'degree_north');

lat_psi_varid = netcdf.defVar(ncid, 'lat_psi', 'double', [xi_psi_dimid, eta_psi_dimid]);
netcdf.putAtt(ncid, lat_psi_varid, 'long_name', 'latitude of PSI-points');
netcdf.putAtt(ncid, lat_psi_varid, 'units', 'degree_north');

mask_rho_varid = netcdf.defVar(ncid, 'mask_rho', 'double', [xi_rho_dimid, eta_rho_dimid]);
netcdf.putAtt(ncid, mask_rho_varid, 'long_name', 'mask on RHO-points');
netcdf.putAtt(ncid, mask_rho_varid, 'option_0', 'land');
netcdf.putAtt(ncid, mask_rho_varid, 'option_1', 'water');

mask_u_varid = netcdf.defVar(ncid, 'mask_u', 'double', [xi_u_dimid, eta_u_dimid]);
netcdf.putAtt(ncid, mask_u_varid, 'long_name', 'mask on U-points');
netcdf.putAtt(ncid, mask_u_varid, 'option_0', 'land');
netcdf.putAtt(ncid, mask_u_varid, 'option_1', 'water');

mask_v_varid = netcdf.defVar(ncid, 'mask_v', 'double', [xi_v_dimid, eta_v_dimid]);
netcdf.putAtt(ncid, mask_v_varid, 'long_name', 'mask on V-points');
netcdf.putAtt(ncid, mask_v_varid, 'option_0', 'land');
netcdf.putAtt(ncid, mask_v_varid, 'option_1', 'water');

mask_psi_varid = netcdf.defVar(ncid, 'mask_psi', 'double', [xi_psi_dimid, eta_psi_dimid]);
netcdf.putAtt(ncid, mask_psi_varid, 'long_name', 'mask on PSI-points');
netcdf.putAtt(ncid, mask_psi_varid, 'option_0', 'land');
netcdf.putAtt(ncid, mask_psi_varid, 'option_1', 'water');

%
% Create global attributes
%
type = 'GRID file' ; 
history = 'ROMS' ;

varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid,varid,'title', title);
netcdf.putAtt(ncid,varid,'date',datestr(date, 'yyyymmdd'));
netcdf.putAtt(ncid,varid,'grd_file', grdname);
netcdf.putAtt(ncid,varid,'type', type);
netcdf.putAtt(ncid,varid,'history', history);
netcdf.putAtt(ncid,varid,'author', 'Created by Yong-Yub Kim');

%
% Leave define mode
%
netcdf.endDef(ncid);
netcdf.close(ncid);
result = 1;
