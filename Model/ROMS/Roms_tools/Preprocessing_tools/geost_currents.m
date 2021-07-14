function geost_currents(clmname,grdname,oaname,frcname,zref,obc,tini)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% compute SSH and the geostrophic currents from Hydrology data
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
%  Adapted from a previous program of Patrick Marchesiello (IRD).
%
%  Copyright (c) 2001-2006 by Patrick Marchesiello and Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%  Updated 04-Apr-2018 by Y.Y.Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clearvars -except clmname grdname oaname frcname zref obc tini Vtransform Vstretching

rho0=1025; % Bousinesq background density [kg.m-3]
g=9.8;     % Gravity acceleration [m.s-2]
De=40;     % Ekman depth [m]
conserv=1;
%
%  grid parameters
%
%disp(' Read grid parameters ...');
nc=netcdf(grdname);
pm=nc{'pm'}(:);
pn=nc{'pn'}(:);
h=nc{'h'}(:);
f=nc{'f'}(:);
lon=nc{'lon_rho'}(:);
lat=nc{'lat_rho'}(:);
rmask=nc{'mask_rho'}(:);
umask=nc{'mask_u'}(:);
vmask=nc{'mask_v'}(:);
[M,L]=size(rmask);
close(nc)
%
%  Levitus vertical levels
%
noa=netcdf(oaname);
Z=-noa{'Z'}(:);
NL=length(Z);
time=noa{'tclm_time'}(:);
tlen=length(time);
%
%  Model grid vertical levels
%
nc=netcdf(clmname,'write');
Vtransform = nc{'Vtransform'}(:);
Vstretching = nc{'Vstretching'}(:);
theta_s = nc{'theta_s'}(:);
theta_b = nc{'theta_b'}(:);
hc  =  nc{'hc'}(:);
N =  length(nc('s_rho'));
%
% get the reference level
%
kref=min(find(Z<=zref));
if isempty(kref);
  kref=length(Z);
  disp(['Warning zref not found. Taking :',num2str(Z(kref))])
end
Z=Z(1:kref);
z=reshape(Z,kref,1,1);
z=repmat(z,[1 M L]);
% % % % % %
% % % % % % Open the forcing file
% % % % % %
% % % % % if ~isempty(frcname)
% % % % %   nfrc=netcdf(frcname);
% % % % % end
%%%%%%%%%%%%%%%%%%%
% START MAIN LOOP
%%%%%%%%%%%%%%%%%%%
%
% loop on time
%
if tini ~= 0  % initial file 
  tlen=1;
end
for l=1:tlen
%for l=1:1
  if tini == 0
    disp(['time index: ',num2str(l),' of total: ',num2str(tlen)])
  else
    l=find(time==tini);
    disp(['time index: ',num2str(l)])
  end
%
% read T and S
%
%   save 'E:\Data\Model\ROMS\nwp_1_20\input\test39\geost_currents_var.mat' %%Y.Y.Kim
  T3d=squeeze(noa{'temp'}(l,1:kref,:,:));
  S3d=squeeze(noa{'salt'}(l,1:kref,:,:));
  Ts=squeeze(T3d(1,:,:));
  Ss=squeeze(S3d(1,:,:));
  rhos=rho_eos(Ts,Ss,0);
  rho=rho_eos(T3d,S3d,z);
  rho_w=.5*(rho(1:kref-1,:,:)+rho(2:kref,:,:));
  z_w=.5*(z(1:kref-1,:,:)+z(2:kref,:,:));
  dz_w=z(1:kref-1,:,:)-z(2:kref,:,:);
%
%  COMPUTE PRESSURE
%
%  disp('Pressure field and sea level ...')
  pres=0*T3d;
%  initialize pressure at kref in Pascal
  pres(kref,:,:)=-zref*1.e4;
  for k=kref-1:-1:1;
    pres(k,:,:)=pres(k+1,:,:)-rho_w(k,:,:).*g.*dz_w(k,:,:);
  end
%
%  compute SSH and remove area averaged
%
  ssh=squeeze(pres(1,:,:))./(rhos*g);
  avgssh=sum(sum(rmask.*ssh./(pm.*pn)))/sum(sum(rmask./(pm.*pn)));
  ssh=ssh-avgssh;
  avgp=tridim(avgssh.*rhos*g,kref);
  pres=pres-avgp;
%
%  COMPUTE GEOSTROPHIC BAROCLINIC VELOCITIES
%
%  disp('Baroclinic geostrophic component ...')

  p_u=rho2u_3d(pres);
  p_v=rho2v_3d(pres);
  
  px(:,:,2:L-1)=p_u(:,:,2:L-1)-p_u(:,:,1:L-2);
  px(:,:,1)=2.*px(:,:,2)-px(:,:,3);
  px(:,:,L)=2.*px(:,:,L-1)-px(:,:,L-2);
  
  py(:,2:M-1,:)=p_v(:,2:M-1,:)-p_v(:,1:M-2,:);
  py(:,1,:)=2.*py(:,2,:)-py(:,3,:);
  py(:,M,:)=2.*py(:,M-1,:)-py(:,M-2,:);

  pn3d=tridim(pn,kref); 
  pm3d=tridim(pm,kref); 
  f3d=tridim(f,kref); 
  m3d=tridim(rmask,kref);   
  u_r=-m3d.*pn3d.*py./(rho0*f3d);
  v_r=m3d.*pm3d.*px./(rho0*f3d);
  
  
%
% Ekman transport  %% not calculated Y.Y. Kim.
%
% % % %   if ~isempty(frcname)
% % % % %    disp('Add the Ekman transport')
% % % %     sustr(:,:)=nfrc{'sustr'}(l,:,:);
% % % %     svstr(:,:)=nfrc{'svstr'}(l,:,:);
% % % % %    rhoA=1.3; Cd=1.4e-3;
% % % % %    W=sqrt(u2rho_2d(sustr).^2+v2rho_2d(svstr).^2);
% % % % %    De=mean(mean(  sqrt(W).* ...
% % % % %       min(250,4.3./sqrt(rhoA.*Cd.*sin(abs(lat.*pi/180)))) ));
% % % % save 'E:\Data\Model\ROMS\nwp_1_20\input\test39\geost_currents_var.mat' %%Y.Y.Kim
% % % %     k_ekman=min(find(Z<=-De));
% % % %     u_r(1:k_ekman,:,:)=u_r(1:k_ekman,:,:)+tridim(rmask.*...
% % % %                        v2rho_2d(svstr)./(rho0*De*f),k_ekman);
% % % %     v_r(1:k_ekman,:,:)=v_r(1:k_ekman,:,:)+tridim(-rmask.*...
% % % %                        u2rho_2d(sustr)./(rho0*De*f),k_ekman);
% % % %   end
  

%
% Replace/interpolate Equatorial values 
%
  equatlat=(lat>=-2 & lat<=2);
  if sum(sum(equatlat))~=0
%    disp('Extrapole values outside the Equator')
    D=find(~equatlat);
    if length(D)~=0
      for k=1:kref
        u_r(k,:,:)=griddata(lon(D),lat(D),u_r(k,D),lon,lat);
        v_r(k,:,:)=griddata(lon(D),lat(D),v_r(k,D),lon,lat);
      end
    else
%      disp('No values outside the Equator to extrapole')
    end
  end
%
%  Masking
%  
  umask3d=tridim(umask,kref); 
  vmask3d=tridim(vmask,kref); 
  u=umask3d.*rho2u_3d(u_r);
  v=vmask3d.*rho2v_3d(v_r);
  ssh=ssh.*rmask; 
%
% Vertical interpolation of baroclinic fields
%
%  disp('Vertical interpolation ...')
  zroms=zlevs(Vtransform,Vstretching,h,ssh,theta_s,theta_b,hc,N,'r');
  zu=0.5*(zroms(:,:,1:end-1)+zroms(:,:,2:end));
  zv=0.5*(zroms(:,1:end-1,:)+zroms(:,2:end,:));
%
% Add non gradient velocities at the top and nul velocities 
% at -10000m for vertical extrapolation.
%
  u=cat(1,u(1,:,:),u);
  v=cat(1,v(1,:,:),v);
  u=flipdim(cat(1,u,0*u(1,:,:)),1);
  v=flipdim(cat(1,v,0*v(1,:,:)),1);
  z1=flipud([100;Z;-10000]);
%
% Do the interpolation
%
  u=ztosigma(u,zu,z1);
  v=ztosigma(v,zv,z1);
%
%  Barotropic velocities
%
%  disp('Barotropic component ...')
  zw=zlevs(Vtransform,Vstretching,h,ssh,theta_s,theta_b,hc,N,'w');
  dz=zw(2:end,:,:)-zw(1:end-1,:,:);
  dzu=0.5*(dz(:,:,1:end-1)+dz(:,:,2:end));
  dzv=0.5*(dz(:,1:end-1,:)+dz(:,2:end,:));
  hu(:,:)=sum(dzu.*u);
  hv(:,:)=sum(dzv.*v);
  D_u(:,:)=sum(dzu);
  D_v(:,:)=sum(dzv);
  ubar(:,:)=hu./D_u;
  vbar(:,:)=hv./D_v;
  u=u-tridim(ubar,N);
  v=v-tridim(vbar,N);
%
% Mass conservation
%
  if conserv==1
%    disp('Volume conservation enforcement ...')
    [hu,hv]=get_obcvolcons(hu,hv,pm,pn,rmask,obc);
%
% Get the stream function
%
    psi=get_psi(hu,hv,pm,pn,rmask); 
    hu(2:end-1,1:end)=-0.5*umask(2:end-1,1:end).*...
                      (psi(2:end,1:end)-psi(1:end-1,1:end)).*...
                      (pn(2:end-1,2:end)+pn(2:end-1,1:end-1));
    hv(1:end,2:end-1)=0.5*vmask(1:end,2:end-1).*...
                     (psi(1:end,2:end)-psi(1:end,1:end-1)).*...
                     (pm(2:end,2:end-1)+pm(1:end-1,2:end-1));
    [hu,hv]=get_obcvolcons(hu,hv,pm,pn,rmask,obc);
    ubar(:,:)=hu./D_u;
    vbar(:,:)=hv./D_v;
  end
  u=u+tridim(ubar,N);
  v=v+tridim(vbar,N);
%
% corners
%
  ubar(1,1)=0.5*(ubar(1,2)+ubar(2,1)); 
  ubar(end,1)=0.5*(ubar(end,2)+ubar(end-1,1)); 
  ubar(1,end)=0.5*(ubar(1,end-1)+ubar(2,end)); 
  ubar(end,end)=0.5*(ubar(end,end-1)+ubar(end-1,end)); 
  vbar(1,1)=0.5*(vbar(1,2)+vbar(2,1)); 
  vbar(end,1)=0.5*(vbar(end,2)+vbar(end-1,1)); 
  vbar(1,end)=0.5*(vbar(1,end-1)+vbar(2,end)); 
  vbar(end,end)=0.5*(vbar(end,end-1)+vbar(end-1,end)); 
  u(:,1,1)=0.5*(u(:,1,2)+u(:,2,1)); 
  u(:,end,1)=0.5*(u(:,end,2)+u(:,end-1,1)); 
  u(:,1,end)=0.5*(u(:,1,end-1)+u(:,2,end)); 
  u(:,end,end)=0.5*(u(:,end,end-1)+u(:,end-1,end)); 
  v(:,1,1)=0.5*(v(:,1,2)+v(:,2,1)); 
  v(:,end,1)=0.5*(v(:,end,2)+v(:,end-1,1)); 
  v(:,1,end)=0.5*(v(:,1,end-1)+v(:,2,end)); 
  v(:,end,end)=0.5*(v(:,end,end-1)+v(:,end-1,end)); 
%
%  Write into file
%
%  disp('Writes into climatology file ...')
  nc{'u'}(l,:,:,:)=u;
  nc{'v'}(l,:,:,:)=v;
  nc{'ubar'}(l,:,:)=ubar;
  nc{'vbar'}(l,:,:)=vbar;
  if tini == 0
    nc{'SSH'}(l,:,:)=ssh;
    nc{'zeta'}(l,:,:)=ssh;
  else
    nc{'zeta'}(l,:,:)=ssh;
  end
end
close(nc)
% % % % % %
% % % % % % Close the forcing file
% % % % % %
% % % % % if ~isempty(frcname)
% % % % %   close(nfrc);
% % % % % end
if nargin==0
  close all
  test_clim(clmname,grdname,'temp',1)
end
