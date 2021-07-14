function z=get_depths_ltrans(fname,gname,type, lonlat_ind);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Get the depths of the sigma levels
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
%  Copyright (c) 2002-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin <= 3                           

    % nc = netcdf(gname);
    % h = nc{'h'}(:);
    % close(nc);
    % nc=netcdf(fname);
    % zeta=nc{'zeta'}(:);
    % theta_s=nc{'theta_s'}(:);
    % s_rho=nc{'s_rho'}(:);
    % N=length(s_rho);
    % theta_b=nc{'theta_b'}(:);
    % Tcline=nc{'Tcline'}(:);
    % Vtransform=nc{'Vtransform'}(:);
    % Vstretching=nc{'Vstretching'}(:);
    % hc=nc{'hc'}(:);
    % close(nc);

    h=ncread(gname, 'h')';
    zeta=ncread(fname, 'zeta')';
    theta_s=ncread(fname, 'theta_s');
    s_rho =ncread(fname, 's_rho');
    N=length(s_rho);
    theta_b = ncread(fname, 'theta_b');
    % Tcline=ncread(fname, 'Tcline');
    Vtransform=ncread(fname, 'Vtransform');
    Vstretching=ncread(fname, 'Vstretching');
    hc=ncread(fname, 'hc');


    if isempty(zeta)
      zeta=0.*h;
    end

    vtype=type;
    if (type=='u')|(type=='v')
      vtype='r';
    end
    z = zlevs(Vtransform, Vstretching, h,zeta,theta_s,theta_b,hc,N,vtype);
    if type=='u'
      z=rho2u_3d(z);
    end
    if type=='v'
      z=rho2v_3d(z);
    end

else
    loncount=lonlat_ind(2)-lonlat_ind(1)+1;
    latcount=lonlat_ind(4)-lonlat_ind(3)+1;
    h=ncread(gname, 'h', [lonlat_ind(1), lonlat_ind(3)], [loncount, latcount])';
    zeta=ncread(fname, 'zeta', [lonlat_ind(1), lonlat_ind(3) 1], [loncount, latcount 1])';
    theta_s=ncread(fname, 'theta_s');
    s_rho =ncread(fname, 's_rho');
    N=length(s_rho);
    theta_b = ncread(fname, 'theta_b');
    % Tcline=ncread('Tcline');
    Vtransform=ncread(fname, 'Vtransform');
    Vstretching=ncread(fname, 'Vstretching');
    hc=ncread(fname, 'hc');


    if isempty(zeta)
      zeta=0.*h;
    end

    vtype=type;
    if (type=='u')|(type=='v')
      vtype='r';
    end
    z = zlevs(Vtransform, Vstretching, h,zeta,theta_s,theta_b,hc,N,vtype);
    if type=='u'
      z=rho2u_3d(z);
    end
    if type=='v'
      z=rho2v_3d(z);
    end

end
return
