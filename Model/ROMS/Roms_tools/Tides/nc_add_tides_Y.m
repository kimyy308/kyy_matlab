function nc_add_tides(fname,Ntides,start_tide_mjd,components)
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
%  Copyright (c) 2001-2006 by Patrick Marchesiello
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ncid= netcdf.open(fname, 'WRITE');
netcdf.reDef(ncid);

%
%  Add dimension
%
tide_period_dimid=netcdf.defDim(ncid, 'tide_period', Ntides);
xi_rho_dimid=netcdf.inqDimID(ncid,'xi_rho');
eta_rho_dimid=netcdf.inqDimID(ncid,'eta_rho');
%
%  Add variables and attributes
%
tide_period_varid=netcdf.defVar(ncid, 'tide_period', 'double', tide_period_dimid);
netcdf.putAtt(ncid, tide_period_varid, 'long_name', 'Tide angular period');
netcdf.putAtt(ncid, tide_period_varid, 'units', 'Hours');

tide_Ephase_varid=netcdf.defVar(ncid, 'tide_Ephase', 'double', [xi_rho_dimid, eta_rho_dimid, tide_period_dimid]);
netcdf.putAtt(ncid, tide_Ephase_varid, 'long_name', 'Tidal elevation phase angle');
netcdf.putAtt(ncid, tide_Ephase_varid, 'units', 'Degrees');

tide_Eamp_varid=netcdf.defVar(ncid, 'tide_Eamp', 'double', [xi_rho_dimid, eta_rho_dimid, tide_period_dimid]);
netcdf.putAtt(ncid, tide_Eamp_varid, 'long_name', 'Tidal elevation amplitude');
netcdf.putAtt(ncid, tide_Eamp_varid, 'units', 'Meter');

tide_Cmin_varid=netcdf.defVar(ncid, 'tide_Cmin', 'double', [xi_rho_dimid, eta_rho_dimid, tide_period_dimid]);
netcdf.putAtt(ncid, tide_Cmin_varid, 'long_name', 'Tidal current ellipse semi-minor axis');
netcdf.putAtt(ncid, tide_Cmin_varid, 'units', 'Meter second-1');

tide_Cmax_varid=netcdf.defVar(ncid, 'tide_Cmax', 'double', [xi_rho_dimid, eta_rho_dimid, tide_period_dimid]);
netcdf.putAtt(ncid, tide_Cmax_varid, 'long_name', 'Tidal current, ellipse semi-major axis');
netcdf.putAtt(ncid, tide_Cmax_varid, 'units', 'Meter second-1');

tide_Cangle_varid=netcdf.defVar(ncid, 'tide_Cangle', 'double', [xi_rho_dimid, eta_rho_dimid, tide_period_dimid]);
netcdf.putAtt(ncid, tide_Cangle_varid, 'long_name', 'Tidal current inclination angle');
netcdf.putAtt(ncid, tide_Cangle_varid, 'units', 'Degrees between semi-major axis and East');

tide_Cphase_varid=netcdf.defVar(ncid, 'tide_Cphase', 'double', [xi_rho_dimid, eta_rho_dimid, tide_period_dimid]);
netcdf.putAtt(ncid, tide_Cphase_varid, 'long_name', 'Tidal current phase angle');
netcdf.putAtt(ncid, tide_Cphase_varid, 'units', 'Degrees');

netcdf.endDef(ncid);
netcdf.close(ncid);

nc=netcdf(fname,'write');
nc.date = ncchar(date);
nc.date = date;
nc.start_tide_mjd=start_tide_mjd;
nc.components = ncchar(components);
nc.components = components;
close(nc)
