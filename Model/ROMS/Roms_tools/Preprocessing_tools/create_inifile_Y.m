function create_inifile_Y(inifile,gridfile,title,...
                         Vtransform, Vstretching, theta_s,theta_b,hc,N,time,clobber)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  function nc=create_inifile(inifile,gridfile,theta_s,...
%                  theta_b,hc,N,ttime,stime,utime,... 
%                  cycle,clobber)
%
%   This function create the header of a Netcdf climatology 
%   file.
%
%   Input: 
% 
%   inifile      Netcdf initial file name (character string).
%   gridfile     Netcdf grid file name (character string).
%   theta_s      S-coordinate surface cncid ontrol parameter.(Real)
%   theta_b      S-coordinate bottom control parameter.(Real)
%   hc           Width (m) of surface or bottom boundary layer
%                where higher vertical resolution is required 
%                during stretching.(Real)
%   N            Number of vertical levels.(Integer)  
%   time         Initial time.(Real) 
%   clobber      Switch to allow or not writing over an existing
%                file.(character string) 
%
%   Output
%
%   nc       Output netcdf object.
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% romstools_param;

disp(' ')
disp([' Creating the file : ',inifile])
%
%  Read the grid file
%
nc=netcdf(gridfile);
h=nc{'h'}(:);  
mask=nc{'mask_rho'}(:);
close(nc);
hmin=min(min(h(mask==1)));
if hc > hmin
%   error([' hc (',num2str(hc),' m) > hmin (',num2str(hmin),' m)'])
[' hc (',num2str(hc),' m) > hmin (',num2str(hmin),' m)']
end
[Mp,Lp]=size(h);
L=Lp-1;
M=Mp-1;
Np=N+1;
%
%  Create the initial file
%
type = 'INITIAL file' ; 
history = 'ROMS' ;
% ncid = netcdf.create(inifile, 'NETCDF4');
ncid = netcdf.create(inifile, 'CLOBBER');

%
%  Create dimensions
%
xi_u_dimid = netcdf.defDim(ncid, 'xi_u', L);
xi_v_dimid = netcdf.defDim(ncid, 'xi_v', Lp);
xi_rho_dimid = netcdf.defDim(ncid, 'xi_rho', Lp);
eta_u_dimid = netcdf.defDim(ncid, 'eta_u', Mp);
eta_v_dimid = netcdf.defDim(ncid, 'eta_v', M);
eta_rho_dimid = netcdf.defDim(ncid, 'eta_rho', Mp);
s_rho_dimid = netcdf.defDim(ncid, 's_rho', N);
s_w_dimid = netcdf.defDim(ncid, 's_w', Np);
tracer_dimid = netcdf.defDim(ncid, 'tracer', 2);
time_dimid = netcdf.defDim(ncid, 'time', 0);
one_dimid = netcdf.defDim(ncid, 'one', 1);

%
%  Create variables
%

tstart_varid = netcdf.defVar(ncid, 'tstart', 'double', [one_dimid]);
netcdf.putAtt(ncid, tstart_varid, 'long_name', 'start processing day');
netcdf.putAtt(ncid, tstart_varid, 'units', 'day');  

tend_varid = netcdf.defVar(ncid, 'tend', 'double', [one_dimid]);
netcdf.putAtt(ncid, tend_varid, 'long_name', 'end processing day');
netcdf.putAtt(ncid, tend_varid, 'units', 'day');

theta_s_varid = netcdf.defVar(ncid, 'theta_s', 'double', [one_dimid]);
netcdf.putAtt(ncid, theta_s_varid, 'long_name', 'S-coordinate surface control parameter');
netcdf.putAtt(ncid, theta_s_varid, 'units', 'nondimensional');

theta_b_varid = netcdf.defVar(ncid, 'theta_b', 'double', [one_dimid]);
netcdf.putAtt(ncid, theta_b_varid, 'long_name', 'S-coordinate bottom control parameter');
netcdf.putAtt(ncid, theta_b_varid, 'units', 'nondimensional');

Tcline_varid = netcdf.defVar(ncid, 'Tcline', 'double', [one_dimid]);
netcdf.putAtt(ncid, Tcline_varid, 'long_name', 'S-coordinate surface/bottom layer width');
netcdf.putAtt(ncid, Tcline_varid, 'units', 'meter');

hc_varid = netcdf.defVar(ncid, 'hc', 'double', [one_dimid]);
netcdf.putAtt(ncid, hc_varid, 'long_name', 'S-coordinate parameter, critical depth');
netcdf.putAtt(ncid, hc_varid, 'units', 'meter');

sc_r_varid = netcdf.defVar(ncid, 'sc_r', 'double', [s_rho_dimid]);
netcdf.putAtt(ncid, sc_r_varid, 'long_name', 'S-coordinate at RHO-points');
netcdf.putAtt(ncid, sc_r_varid, 'valid_min', -1);
netcdf.putAtt(ncid, sc_r_varid, 'valid_max', 0);

Cs_r_varid = netcdf.defVar(ncid, 'Cs_r', 'double', [s_rho_dimid]);
netcdf.putAtt(ncid, Cs_r_varid, 'long_name', 'S-coordinate stretching curves at RHO-points');
netcdf.putAtt(ncid, Cs_r_varid, 'valid_min', -1);
netcdf.putAtt(ncid, Cs_r_varid, 'valid_max', 0);

ocean_time_varid = netcdf.defVar(ncid, 'ocean_time', 'double', [time_dimid]);
netcdf.putAtt(ncid, ocean_time_varid, 'long_name', 'time since initialization');
netcdf.putAtt(ncid, ocean_time_varid, 'units', 'second');

scrum_time_varid = netcdf.defVar(ncid, 'scrum_time', 'double', [time_dimid]);
netcdf.putAtt(ncid, scrum_time_varid, 'long_name', 'time since initialization');
netcdf.putAtt(ncid, scrum_time_varid, 'units', 'second');

u_varid = netcdf.defVar(ncid, 'u', 'double', [xi_u_dimid, eta_u_dimid, s_rho_dimid, time_dimid]);
netcdf.putAtt(ncid, u_varid, 'long_name', 'u-momentum component');
netcdf.putAtt(ncid, u_varid, 'units', 'meter second-1');

v_varid = netcdf.defVar(ncid, 'v', 'double', [xi_v_dimid, eta_v_dimid, s_rho_dimid, time_dimid]);
netcdf.putAtt(ncid, v_varid, 'long_name', 'v-momentum component');
netcdf.putAtt(ncid, v_varid, 'units', 'meter second-1');

ubar_varid = netcdf.defVar(ncid, 'ubar', 'double', [xi_u_dimid, eta_u_dimid, time_dimid]);
netcdf.putAtt(ncid, ubar_varid, 'long_name', 'vertically integrated u-momentum component');
netcdf.putAtt(ncid, ubar_varid, 'units', 'meter second-1');

vbar_varid = netcdf.defVar(ncid, 'vbar', 'double', [xi_v_dimid, eta_v_dimid, time_dimid]);
netcdf.putAtt(ncid, vbar_varid, 'long_name', 'vertically integrated v-momentum component');
netcdf.putAtt(ncid, vbar_varid, 'units', 'meter second-1');

zeta_varid = netcdf.defVar(ncid, 'zeta', 'double', [xi_rho_dimid, eta_rho_dimid, time_dimid]);
netcdf.putAtt(ncid, zeta_varid, 'long_name', 'free-surface');
netcdf.putAtt(ncid, zeta_varid, 'units', 'meter');

temp_varid = netcdf.defVar(ncid, 'temp', 'double', [xi_rho_dimid, eta_rho_dimid, s_rho_dimid, time_dimid]);
netcdf.putAtt(ncid, temp_varid, 'long_name', 'potential temperature');
netcdf.putAtt(ncid, temp_varid, 'units', 'Celsius');

salt_varid = netcdf.defVar(ncid, 'salt', 'double', [xi_rho_dimid, eta_rho_dimid, s_rho_dimid, time_dimid]);
netcdf.putAtt(ncid, salt_varid, 'long_name', 'salinity');
netcdf.putAtt(ncid, salt_varid, 'units', ' ');

Vtransform_varid = netcdf.defVar(ncid, 'Vtransform', 'double', [one_dimid]);
netcdf.putAtt(ncid, Vtransform_varid, 'long_name', 'vertical terrain-following transformation equation');

Vstretching_varid = netcdf.defVar(ncid, 'Vstretching', 'double', [one_dimid]);
netcdf.putAtt(ncid, Vstretching_varid, 'long_name', 'vertical terrain-following stretching function');

%
% Create global attributes
%
varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid,varid,'title', title);
netcdf.putAtt(ncid,varid,'date',datestr(date, 'yyyymmdd'));
netcdf.putAtt(ncid,varid,'clim_file', inifile);
netcdf.putAtt(ncid,varid,'grd_file', gridfile);
netcdf.putAtt(ncid,varid,'type', type);
netcdf.putAtt(ncid,varid,'history', history);
netcdf.putAtt(ncid,varid,'author', 'Created by Yong-Yub Kim');


%
% Leave define mode
%
netcdf.endDef(ncid);
netcdf.close(ncid);

% Write variables
%
nc = netcdf(inifile, 'write');

nc{'Vtransform'}(:)=Vtransform;
nc{'Vstretching'}(:)=Vstretching;
kgrid=0;
[sc,Cs]=stretching(Vstretching, theta_s, theta_b, hc, N, kgrid);

nc{'tstart'}(:) =  time; 
nc{'tend'}(:) =  time; 
nc{'theta_s'}(:) =  theta_s; 
nc{'theta_b'}(:) =  theta_b; 
nc{'Tcline'}(:) =  hc; 
nc{'hc'}(:) =  hc; 
nc{'sc_r'}(:) =  sc; 
nc{'Cs_r'}(:) =  Cs; 
nc{'scrum_time'}(1) =  time*24*3600; 
nc{'ocean_time'}(1) =  time*24*3600; 
nc{'u'}(:) =  0; 
nc{'v'}(:) =  0; 
nc{'zeta'}(:) =  0; 
nc{'ubar'}(:) =  0; 
nc{'vbar'}(:) =  0; 
nc{'temp'}(:) =  0; 
nc{'salt'}(:) =  0; 
%
% Synchronize on disk
%
close(nc);
return


