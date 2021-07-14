%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Build a ROMS initial file from Levitus Data
%
%  Extrapole and interpole temperature and salinity from a
%  Climatology to get initial conditions for
%  ROMS (initial netcdf files) .
%  Get the velocities and sea surface elevation via a 
%  geostrophic computation.
%
%  Data input format (netcdf):
%     temperature(T, Z, Y, X)
%     T : time [Months]
%     Z : Depth [m]
%     Y : Latitude [degree north]
%     X : Longitude [degree east]
%
%  Data source : IRI/LDEO Climate Data Library (World Ocean Atlas 1998)
%    http://ingrid.ldgo.columbia.edu/
%    http://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NODC/.WOA98/
%
%  P. Marchesiello & P. Penven - IRD 2005
%
%  Version of 21-Sep-2005
%  made    10-Jul-2018 by Yong-Yub Kim 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
%%%%%%%%%%%%%%%%%%%%% USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%
% Common parameters
%
romstools_param

%
%  Title 
%
title='ROMS initial file from restart data(different topography)';

%
%
%%%%%%%%%%%%%%%%%%% END USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%
%
% Title
%
disp(' ')
disp([' Making initial file: ',ininame])
disp(' ')
disp([' Title: ',title])
%
% Initial file
%
% create_inifile(ininame,grdname,title,...
%                theta_s,theta_b,hc,N,...
%                tini,'clobber');

create_inifile_Y(ininame,grdname,title,...
               theta_s,theta_b,hc,N,...
               tini,'clobber');
%
% Horizontal and vertical interp/extrapolations 
%
disp(' ')
disp(' Interpolations / extrapolations')
disp(' ')


% Open initial file
%
ncnew=netcdf(ininame,'write');
new_Vtransform = ncnew{'Vtransform'}(:);
new_Vstretching =  ncnew{'Vstretching'}(:);
new_theta_s = ncnew{'theta_s'}(:);
new_theta_b =  ncnew{'theta_b'}(:);
new_hc  =  ncnew{'hc'}(:);
new_zeta  =  ncnew{'zeta'}(:);
new_N =  length(ncnew('s_rho'));
%
% Open and Read grid file  
% 
ng=netcdf(grdname);
lon=ng{'lon_rho'}(:);
lat=ng{'lat_rho'}(:);
lon_u=ng{'lon_u'}(:);
lat_u=ng{'lat_u'}(:);
lon_v=ng{'lon_v'}(:);
lat_v=ng{'lat_v'}(:);
new_h=ng{'h'}(:);
close(ng);
[M,L]=size(lon);  %% M -> size of lat(y), L -> size of lon(x)

disp(' Read data from old rst file...')
tic;
ncold=netcdf(rstname);
old_lon = ncold{'lon_rho'}(:);
old_lat = ncold{'lat_rho'}(:);
old_lon_u = ncold{'lon_u'}(:);
old_lat_u = ncold{'lat_u'}(:);
old_lon_v = ncold{'lon_v'}(:);
old_lat_v= ncold{'lat_v'}(:);
old_Vtransform = ncold{'Vtransform'}(:);
old_Vstretching =  ncold{'Vstretching'}(:);
old_theta_s = ncold{'theta_s'}(:);
old_theta_b =  ncold{'theta_b'}(:);
old_hc  =  ncold{'hc'}(:);
old_N =  length(ncold('s_rho'));
old_h=ncold{'h'}(:);
old_zeta=ncold{'zeta'}(:);
old_temp=ncold{'temp'}(:);
old_salt=ncold{'salt'}(:);
old_u=ncold{'u'}(:);
old_v=ncold{'v'}(:);
old_ubar=ncold{'ubar'}(:);
old_vbar=ncold{'vbar'}(:);
close(ncold)

[old_M,old_L]=size(old_h);  %% M -> size of lat(y), L -> size of lon(x)


old_temp(old_temp<-3)=NaN;
old_temp(old_temp>40)=NaN;
old_salt(old_salt<-5)=NaN;
old_salt(old_salt>40)=NaN;
old_u(old_u<-10)=NaN;
old_u(old_u>10)=NaN;
old_v(old_v<-10)=NaN;
old_v(old_v>10)=NaN;


disp(' horizontal interpolation from old data to new data (zeta)...')
new_zeta=griddata(old_lon,old_lat,old_zeta,lon,lat);
new_zeta(new_zeta>5)=0;
new_zeta(new_zeta<-5)=0;
new_z_r=zlevs(new_Vtransform,new_Vstretching,new_h,new_zeta,new_theta_s,new_theta_b,new_hc,new_N,'r');
new_z_w=zlevs(new_Vtransform,new_Vstretching,new_h,new_zeta,new_theta_s,new_theta_b,new_hc,new_N,'w');
new_z_r(new_z_r>0)=0;
new_z_w(new_z_w>0)=0;


% % new_h_v=0.5*(new_h(1:M-1,:)+new_h(2:M,:));
% % new_h_u=0.5*(new_h(:,1:L-1)+new_h(:,2:L));
% % for level=1:N   
% %     new_thickness(level,:,:)=new_z_w(level+1,:,:)-new_z_w(level,:,:); 
% % end
% % new_thickness_u=0.5*(new_thickness(:,:,1:L-1)+new_thickness(:,:,2:L)); % each layer thickness
% % new_z_u(1,:,:)=-new_h_u(:,:);             % z @ bottom of each layer
% % for k=2:+1:N
% %     new_z_u(k,:,:)=new_z_u(k-1,:,:)+new_thickness_u(k-1,:,:);
% % end
% %                 
% % new_thickness_v=0.5*(new_thickness(:,1:M-1,:)+new_thickness(:,2:M,:)); % each layer thickness
% % new_z_v(1,:,:)=-new_h_v(:,:);             % z @ bottom of each layer
% % for k=2:+1:N
% %     new_z_v(k,:,:)=new_z_v(k-1,:,:)+new_thickness_v(k-1,:,:); 
% % end
new_z_u=rho2u_3d(new_z_r);
new_z_v=rho2v_3d(new_z_r);
new_z_u(new_z_u>0)=0;
new_z_v(new_z_v>0)=0;


old_z_r=zlevs(old_Vtransform,old_Vstretching,old_h,old_zeta,old_theta_s,old_theta_b,old_hc,old_N,'r');                
old_z_w=zlevs(old_Vtransform,old_Vstretching,old_h,old_zeta,old_theta_s,old_theta_b,old_hc,old_N,'w');

old_h_v=0.5*(old_h(1:old_M-1,:)+old_h(2:old_M,:));
old_h_u=0.5*(old_h(:,1:old_L-1)+old_h(:,2:old_L));

for level=1:N   
    old_thickness(level,:,:)=old_z_w(level+1,:,:)-old_z_w(level,:,:); 
end
old_thickness_u=0.5*(old_thickness(:,:,1:old_L-1)+old_thickness(:,:,2:old_L)); % each layer thickness
old_z_u(1,:,:)=-old_h_u(:,:);             % z @ bottom of each layer
for k=2:+1:N
    old_z_u(k,:,:)=old_z_u(k-1,:,:)+old_thickness_u(k-1,:,:);
end
                
old_thickness_v=0.5*(old_thickness(:,1:old_M-1,:)+old_thickness(:,2:old_M,:)); % each layer thickness
old_z_v(1,:,:)=-old_h_v(:,:);             % z @ bottom of each layer
for k=2:+1:N
    old_z_v(k,:,:)=old_z_v(k-1,:,:)+old_thickness_v(k-1,:,:); 
end            
toc;

zcorlev=[-10000 -5000 -4500 -4000 -3500 -3000 -2500 -2000 -1700 -1400 -1150 -1000 ...
         -950 -900 -850 -800 -750 -700 -650 -600 -550 -500 ...
         -480 -460 -440 -420 -400 -380 -360 -340 -320 -300 ...
         -280 -260 -240 -220 -200 -190 -180 -170 -160 -150 ...
         -140 -130 -120 -110 -100 -90 -80 -70 -60 -50 ...
         -45 -40 -35 -30 -25 -20 -15 -10 -5 -0];
zcor_z_r(1:length(zcorlev),1:old_M,1:old_L)=NaN;
zcor_z_u(1:length(zcorlev),1:old_M,1:old_L-1)=NaN;     
zcor_z_v(1:length(zcorlev),1:old_M-1,1:old_L)=NaN;     

for k=1:length(zcorlev)
    zcor_z_r(k,1:size(old_z_r,2),1:size(old_z_r,3))=zcorlev(k);
    zcor_z_u(k,1:size(old_z_u,2),1:size(old_z_u,3))=zcorlev(k);
    zcor_z_v(k,1:size(old_z_v,2),1:size(old_z_v,3))=zcorlev(k);
end

disp(' vertical interpolation from old coordinate to z coordinate (T, S)...')
tic;
z_temp(1:length(zcorlev),1:old_M,1:old_L)=NaN;
z_salt(1:length(zcorlev),1:old_M,1:old_L)=NaN;
for i= 1:old_M
    for j=1:old_L
        z_temp(:,i,j)=interp1(old_z_r(:,i,j),old_temp(:,i,j),zcor_z_r(:,i,j),'linear');
        z_salt(:,i,j)=interp1(old_z_r(:,i,j),old_salt(:,i,j),zcor_z_r(:,i,j),'linear');
    end
end
toc;
num=1;
while(num<length(zcorlev))  
    nanind=mean(mean(isnan(z_temp(num,:,:))));
    if (nanind==1)
        num=num+1;
    else
        valid_bot_num=num;
        num=length(zcorlev);
    end
end
% while(num>1)
%     nanind=mean(mean(isnan(z_temp(num,:,:))));
%     if (nanind==1)
%         num=num-1;
%     else
%         valid_surf_num=num;
%         num=1;
%     end
% end
if (valid_bot_num>1)
    for i=1:valid_bot_num-1
        z_temp(i,:,:)=z_temp(valid_bot_num,:,:);
        z_salt(i,:,:)=z_salt(valid_bot_num,:,:);
    end
end

valid_surf_num=length(zcorlev)-2;
if (valid_surf_num<length(zcorlev))
    for i=valid_surf_num+1:length(zcorlev)
        z_temp(i,:,:)=z_temp(valid_surf_num,:,:);
        z_salt(i,:,:)=z_salt(valid_surf_num,:,:);
    end
end

disp(' horizontal interpolation to fill from the nearest value (T, S)...')
tic;
for k=1:length(zcorlev)
    squeezed_var=squeeze(z_temp(k,:,:));
    ismask=isnan(squeezed_var);
    squeezed_var(ismask)=griddata(old_lon(~ismask),old_lat(~ismask),squeezed_var(~ismask),old_lon(ismask),old_lat(ismask),'nearest'); %#ok<GRIDD>
    z_temp(k,:,:)=squeezed_var;
    
    squeezed_var=squeeze(z_salt(k,:,:));
    ismask=isnan(squeezed_var);
    squeezed_var(ismask)=griddata(old_lon(~ismask),old_lat(~ismask),squeezed_var(~ismask),old_lon(ismask),old_lat(ismask),'nearest'); %#ok<GRIDD>
    z_salt(k,:,:)=squeezed_var;
end

toc;
disp(' horizontal interpolation from old grid to new grid (T, S)...')
tic;
z_temp2(1:length(zcorlev),1:M,1:L)=NaN;
z_salt2(1:length(zcorlev),1:M,1:L)=NaN;
for k=1:length(zcorlev)
    z_temp2(k,1:M,1:L)=griddata(old_lon,old_lat,squeeze(z_temp(k,:,:)),lon,lat); %#ok<GRIDD>
    z_salt2(k,1:M,1:L)=griddata(old_lon,old_lat,squeeze(z_salt(k,:,:)),lon,lat); %#ok<GRIDD>
end
toc;

disp(' vertical interpolation from z coordinate to new coordinate (T, S)...')
tic;
% for i= 1:M
%     for j=1:L
%         new_temp(:,i,j)=interp1(zcor_z_r(:,i,j),z_temp(:,i,j),new_z_r(:,i,j),'pchip','extrap');
%         new_salt(:,i,j)=interp1(zcor_z_r(:,i,j),z_salt(:,i,j),new_z_r(:,i,j),'linear',mean(z_salt(:,i,j),'omitnan'));
%     end
% end
% new_temp(new_temp>36)=36;
% new_temp(new_temp<-2)=-2;
% new_salt(new_salt>36)=36;
% new_salt(new_salt<-1)=-1;
new_temp=ztosigma(z_temp2,new_z_r,zcorlev);
new_salt=ztosigma(z_salt2,new_z_r,zcorlev);

toc;

disp(' vertical interpolation from old coordinate to z coordinate (U)...')
tic;
z_u(1:length(zcorlev),1:old_M,1:old_L-1)=NaN;
for i= 1:old_M
    for j=1:old_L-1
        z_u(:,i,j)=interp1(old_z_u(:,i,j),old_u(:,i,j),zcor_z_u(:,i,j),'nearest');
    end
end
toc;
if (valid_bot_num>1)
    for i=1:valid_bot_num-1
        z_u(i,:,:)=z_u(valid_bot_num,:,:);
    end
end
valid_surf_num=length(zcorlev)-2;
if (valid_surf_num<length(zcorlev))
    for i=valid_surf_num+1:length(zcorlev)
        z_u(i,:,:)=z_u(valid_surf_num,:,:);
    end
end

disp(' horizontal interpolation to fill from the nearest value (U)...')
tic;
for k= 1:length(zcorlev)
    squeezed_var_u=squeeze(z_u(k,:,:));
    ismask=isnan(squeezed_var_u);
    squeezed_var_u(ismask)=griddata(old_lon_u(~ismask),old_lat_u(~ismask),squeezed_var_u(~ismask),old_lon_u(ismask),old_lat_u(ismask),'nearest'); %#ok<GRIDD>
    z_u(k,:,:)=squeezed_var_u;
end
toc;

disp(' horizontal interpolation from old grid to new grid (U)...')
tic;
z_u2(1:length(zcorlev),1:M,1:L-1)=NaN;
for k=1:length(zcorlev)
    z_u2(k,1:M,1:L-1)=griddata(old_lon_u,old_lat_u,squeeze(z_u(k,:,:)),lon_u,lat_u); %#ok<GRIDD>
end

clear new_ubar
new_ubar(1:M,1:L-1)=0;
for i= 1:M
    for j=1:L-1
        new_ubar(i,j)= sum(z_u2(1:length(zcorlev)-1,i,j) .* diff(zcorlev)')/5000.0;
    end
end
toc;


disp(' vertical interpolation from z coordinate to new coordinate (U)...')
tic;
% for i= 1:M
%     for j=1:L-1
%         new_u(:,i,j)=interp1(zcor_z_u(:,i,j),z_u(:,i,j),new_z_u(:,i,j),'pchip',0);
%     end
% end
% new_u(new_u>5)=5;
% new_u(new_u<-5)=-5;
new_u=ztosigma(z_u2,new_z_u,zcorlev);
toc;

disp(' vertical interpolation from old coordinate to z coordinate (V)...')
tic;
z_v(1:length(zcorlev),1:old_M-1,1:old_L)=NaN;

for i= 1:old_M-1
    for j=1:old_L
        z_v(:,i,j)=interp1(old_z_v(:,i,j),old_v(:,i,j),zcor_z_v(:,i,j),'nearest');
    end
end
toc;
if (valid_bot_num>1)
    for i=1:valid_bot_num-1
        z_v(i,:,:)=z_v(valid_bot_num,:,:);
    end
end
valid_surf_num=length(zcorlev)-2;
if (valid_surf_num<length(zcorlev))
    for i=valid_surf_num+1:length(zcorlev)
        z_v(i,:,:)=z_v(valid_surf_num,:,:);
    end
end

disp(' horizontal interpolation to fill from the nearest value (V)...')
tic;
for k=1:length(zcorlev)
    squeezed_var_v=squeeze(z_v(k,:,:));
    ismask=isnan(squeezed_var_v);
    squeezed_var_v(ismask)=griddata(old_lon_u(~ismask),old_lat_u(~ismask),squeezed_var_v(~ismask),old_lon_u(ismask),old_lat_u(ismask),'nearest'); %#ok<GRIDD>
    z_v(k,:,:)=squeezed_var_v;
end
toc;

disp(' horizontal interpolation from old grid to new grid (V)...')
tic;
z_v2(1:length(zcorlev),1:M-1,1:L)=NaN;
for k=1:length(zcorlev)
    z_v2(k,1:M-1,1:L)=griddata(old_lon_v,old_lat_v,squeeze(z_v(k,:,:)),lon_v,lat_v); %#ok<GRIDD>
end

clear new_vbar
new_vbar(1:M-1,1:L)=0;
for i= 1:M-1
    for j=1:L
        new_vbar(i,j)= sum(z_v2(1:length(zcorlev)-1,i,j) .* diff(zcorlev)' )/5000.0;
    end
end
toc;

disp(' vertical interpolation from z coordinate to new coordinate (V)...')
tic;
% for i= 1:M-1
%     for j=1:L
%         new_v(:,i,j)=interp1(zcor_z_v(:,i,j),z_v(:,i,j),new_z_v(:,i,j),'pchip',0);
%     end
% end
% new_v(new_v>5)=5;
% new_v(new_v<-5)=-5;
new_v=ztosigma(z_v2,new_z_v,zcorlev);
toc;


disp(' write data to new initial file...')
tic;
ncnew{'temp'}(1,:,:,:)=new_temp;
ncnew{'salt'}(1,:,:,:)=new_salt;
ncnew{'u'}(1,:,:,:)=new_u;
ncnew{'v'}(1,:,:,:)=new_v;
ncnew{'ubar'}(1,:,:)=new_ubar;
ncnew{'vbar'}(1,:,:)=new_vbar;
ncnew{'zeta'}(1,:,:)=new_zeta;
close(ncnew);
toc;
%
% End
%0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
