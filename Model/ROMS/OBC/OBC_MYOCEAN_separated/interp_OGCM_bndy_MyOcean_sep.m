function interp_OGCM_bndy_MyOcean_sep(OGCM_path_day, Roa, interp_method, ...
        OGCM, year, month, ROMS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read the separated OGCM files(East, West, South, North) and perform the interpolations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updated 15-Oct-2020 by Y.Y.Kim

% 
% present time setting for filename
% 

% today.all=OGCM_path_day(end-8:end-1);
% today.year=today.all(1:4);
% today.month=today.all(5:6);
% today.day=today.all(7:8);
% today.num=datenum([num2str(today.month),'-',num2str(today.day),'-',num2str(today.year)]); % mm-dd-yyyy
% 
% presentday.vec=datevec(datestr(today.num+tndx_OGCM-1));
% presentday.year=presentday.vec(1);
% presentday.month=presentday.vec(2);
% presentday.day=presentday.vec(3);
% presentday.all=[num2str(presentday.year,'%04i'), num2str(presentday.month,'%02i'), num2str(presentday.day,'%02i')];


OGCM_path = [OGCM_path_day];
% if (exist(OGCM_path, 'dir')~=7)
%     OGCM_path = [OGCM_path_day(1:end-9), presentday.all, '_forecast/'];
% end
%

disp(['  Horizontal interpolation: ' ])

%
% ROMS grid angle
%

ROMS.cosa = cos(ROMS.angle);
ROMS.sina = sin(ROMS.angle);

%
% Open the OGCM file and read data(value)
% var# in OGCM info
% var1 = zeta, var2 = time, var3 = salt, var4 = temp, 
% var5 = lon, var6 = lat, var7 = u, var8 = depth, var9= v
% you don't have to consider scale_factor and add_offset and missing value
% when using ncread function (in recent version, R2018b)
%  --> wrong


yearstr=num2str(year,'%04i');
monthstr=num2str(month, '%02i');

for diri=1:4
    OGCM.data(diri).filename = [OGCM_path, 'myocean_', ROMS.data(diri).direction, '_', yearstr, monthstr,'.nc'];
    OGCM.data(diri).info = ncinfo(OGCM.data(diri).filename);
    OGCM.num_var = length(OGCM.data(1).info.Variables(:));
    OGCM.data(diri).missvalue = OGCM.data(diri).info.Variables(1).Attributes(8).Value;
    for vari = 1:OGCM.num_var
        OGCM.data(diri).var(vari).value = ...
            ncread(OGCM.data(diri).filename, OGCM.data(diri).info.Variables(vari).Name);
%         ROMS.var(vari).name=OGCM.data(diri).info.Variables(vari).Name;
    end
end

for vari=1:OGCM.num_var
%    ROMS.var(vari).name=OGCM.data(diri).info.Variables(vari).Name;
%    if strcmp(ROMS.var(vari).name, 'zeta')
%        ind_zeta=vari;
%    elseif strcmp(ROMS.var(vari).name, 'temp')
%        ind_temp=vari;
%    elseif strcmp(ROMS.var(vari).name, 'salt')
%        ind_salt=vari;
%    elseif strcmp(ROMS.var(vari).name, 'u')
%        ind_u=vari;
%    elseif strcmp(ROMS.var(vari).name, 'v')
%        ind_v=vari;
%    end
    ROMS.var(vari).name=OGCM.data(diri).info.Variables(vari).Name;
    if strcmp(OGCM.data(diri).info.Variables(vari).Name, 'zos')
        ind_zeta=vari;
    elseif strcmp(OGCM.data(diri).info.Variables(vari).Name, 'thetao')
        ind_temp=vari;
    elseif strcmp(OGCM.data(diri).info.Variables(vari).Name, 'so')
        ind_salt=vari;
    elseif strcmp(OGCM.data(diri).info.Variables(vari).Name, 'uo')
        ind_u=vari;
    elseif strcmp(OGCM.data(diri).info.Variables(vari).Name, 'vo')
        ind_v=vari;
    end
end

ROMS.var(ind_zeta).name='zeta';
ROMS.var(ind_temp).name='temp';
ROMS.var(ind_salt).name='salt';
ROMS.var(ind_u).name='u';
ROMS.var(ind_v).name='v';

% ROMS.var(1).name = 'zeta';
% ROMS.var(2).name = 'time';
% ROMS.var(3).name = 'salt';
% ROMS.var(4).name = 'temp';
% ROMS.var(5).name = 'lon';
% ROMS.var(6).name = 'lat';
% ROMS.var(7).name = 'u';
% ROMS.var(8).name = 'depth';
% ROMS.var(9).name = 'v';
ROMS.var(10).name = 'ubar';
ROMS.var(11).name = 'vbar';


for diri=1
    for vari=[ind_zeta,ind_salt,ind_temp] % 1 : zeta, 3 : salt, 4 : temp
        ROMS.data(diri).var(vari).validgrid = ROMS.lat.rho.all(:,end);
    end
    ROMS.data(diri).var(ind_u).validgrid = ROMS.lat.u.all(:,end);  % 7 : u
    ROMS.data(diri).var(10).validgrid = ROMS.lat.u.all(:,end);  % 10 : ubar
    ROMS.data(diri).var(ind_v).validgrid = ROMS.lat.v.all(:,end);  % 9 : v
    ROMS.data(diri).var(11).validgrid = ROMS.lat.v.all(:,end);  % 9 : vbar
    ROMS.data(diri).h= ROMS.h(:,end);
    ROMS.data(diri).hu= ROMS.hu(:,end);
    ROMS.data(diri).hv= ROMS.hv(:,end);
end
for diri=2
    for vari=[ind_zeta,ind_salt,ind_temp] % 1 : zeta, 3 : salt, 4 : temp
        ROMS.data(diri).var(vari).validgrid = ROMS.lat.rho.all(:,1);
    end
    ROMS.data(diri).var(ind_u).validgrid = ROMS.lat.u.all(:,1);  % 7 : u
    ROMS.data(diri).var(10).validgrid = ROMS.lat.u.all(:,1);  % 10 : ubar
    ROMS.data(diri).var(ind_v).validgrid = ROMS.lat.v.all(:,1);  % 9 : v
    ROMS.data(diri).var(11).validgrid = ROMS.lat.v.all(:,1);  % 9 : vbar
    ROMS.data(diri).h= ROMS.h(:,1);
    ROMS.data(diri).hu= ROMS.hu(:,1);
    ROMS.data(diri).hv= ROMS.hv(:,1);
end
for diri=3
    for vari=[ind_zeta,ind_salt,ind_temp] % 1 : zeta, 3 : salt, 4 : temp
        ROMS.data(diri).var(vari).validgrid = ROMS.lon.rho.all(1,:);
    end
    ROMS.data(diri).var(ind_u).validgrid = ROMS.lon.u.all(1,:);  % 7 : u
    ROMS.data(diri).var(10).validgrid = ROMS.lon.u.all(1,:);  % 10 : ubar
    ROMS.data(diri).var(ind_v).validgrid = ROMS.lon.v.all(1,:);  % 9 : v
    ROMS.data(diri).var(11).validgrid = ROMS.lon.v.all(1,:);  % 11 : vbar
    ROMS.data(diri).h= ROMS.h(1,:);
    ROMS.data(diri).hu= ROMS.hu(1,:);
    ROMS.data(diri).hv= ROMS.hv(1,:);
end
for diri=4
    for vari=[ind_zeta,ind_salt,ind_temp] % 1 : zeta, 3 : salt, 4 : temp
        ROMS.data(diri).var(vari).validgrid = ROMS.lon.rho.all(end,:);
    end
    ROMS.data(diri).var(ind_u).validgrid = ROMS.lon.u.all(end,:);  % 7 : u
    ROMS.data(diri).var(10).validgrid = ROMS.lon.u.all(end,:);  % 10 : ubar
    ROMS.data(diri).var(ind_v).validgrid = ROMS.lon.v.all(end,:);  % 9 : v
    ROMS.data(diri).var(11).validgrid = ROMS.lon.v.all(end,:);  % 11 : vbar
    ROMS.data(diri).h= ROMS.h(end,:);
    ROMS.data(diri).hu= ROMS.hu(end,:);
    ROMS.data(diri).hv= ROMS.hv(end,:);
end

%
% Read and extrapole the 2D variable (zeta) 
% and 3D variables (temp, salt, u, v) of OGCM data
%
for diri = 1:4
    OGCM.data(diri).var(ind_zeta).value = ...
        ext_line_auto(OGCM.data(diri).var(ind_zeta).value, OGCM.data(diri).var(ind_zeta).value, OGCM.data(diri).validgrid); 
    for depthi = 1:OGCM.NZ
        for vari = [ind_salt, ind_temp, ind_u, ind_v] % 3 : salt, 4 : temp, 7 : u, 9 : v
            if (depthi == 1)
                OGCM.data(diri).var(vari).value(:,:,depthi) = ...
                    ext_line_auto(OGCM.data(diri).var(vari).value(:,:,depthi), ...
                    OGCM.data(diri).var(vari).value(:,:,depthi), OGCM.data(diri).validgrid);
            else
                OGCM.data(diri).var(vari).value(:,:,depthi) = ...
                    ext_line_auto(OGCM.data(diri).var(vari).value(:,:,depthi), ...
                    OGCM.data(diri).var(vari).value(:,:,depthi-1), OGCM.data(diri).validgrid);
            end
        end
    end
end

% get dz array for u and v at each direction
for diri=1:4
    OGCM.data(diri).dz=gradient(OGCM.depth);
    OGCM.data(diri).dzu_3d = repmat(reshape(-OGCM.data(diri).dz, ...
        [1,1,OGCM.NZ]), [size(OGCM.data(diri).var(ind_u).value(:,:,1)), 1]);
    OGCM.data(diri).dzv_3d = repmat(reshape(-OGCM.data(diri).dz, ...
        [1,1,OGCM.NZ]), [size(OGCM.data(diri).var(ind_v).value(:,:,1)), 1]);
end

% caculate barotropic u and v
for diri=1:4
    OGCM.data(diri).var(10).value=sum(OGCM.data(diri).var(ind_u).value .* OGCM.data(diri).dzu_3d, 3) ./ sum(OGCM.data(diri).dzu_3d, 3);
    OGCM.data(diri).var(11).value=sum(OGCM.data(diri).var(ind_v).value .* OGCM.data(diri).dzv_3d, 3) ./ sum(OGCM.data(diri).dzv_3d, 3);
end

% 
% % Curved u and v considering angle will not be calculated.
% 


%
% Interpolation of OGCM data on ROMS model grid at each direction
%
for diri = 1:4
    for vari = [ind_zeta, 10, 11]  % 1: zeta, 10: ubar, 11: vbar
        ROMS.data(diri).var(vari).value = ...
            interp1(OGCM.data(diri).validgrid, OGCM.data(diri).var(vari).value, ...
            ROMS.data(diri).var(vari).validgrid);
    end
    for depthi = 1:OGCM.NZ
        for vari = [ind_salt, ind_temp, ind_u, ind_v] % 3: salt, 4: temp, 7: u, 9: v
            ROMS.data(diri).var(vari).value(:,:,depthi) = ...
            interp1(OGCM.data(diri).validgrid, OGCM.data(diri).var(vari).value(:,:,depthi), ...
                ROMS.data(diri).var(vari).validgrid);
        end
    end
end

%
%Initialisation in case of bry files
%

%
% Get the ROMS vertical grid
%
disp('  Vertical interpolations')

%
% Add an extra bottom layer (-100000m) and an extra surface layer (+100m)
% to prevent vertical extrapolations
%
OGCM.depth = [100; -OGCM.depth;-100000];

%
% ROMS vertical grid
%
for diri=1:4
    ROMS.data(diri).zr = zlevs(ROMS.Vtransform, ROMS.Vstretching, ...
        ROMS.data(diri).h, ROMS.data(diri).var(ind_zeta).value, ...
        ROMS.theta_s, ROMS.theta_b, ROMS.hc, ROMS.N, 'r');
    ROMS.data(diri).zw = zlevs(ROMS.Vtransform, ROMS.Vstretching, ...
        ROMS.data(diri).h, ROMS.data(diri).var(ind_zeta).value, ...
        ROMS.theta_s, ROMS.theta_b, ROMS.hc, ROMS.N, 'w');
end
for diri=1:2
    ROMS.data(diri).zu = zlevs(ROMS.Vtransform, ROMS.Vstretching, ...
        ROMS.data(diri).hu, ROMS.data(diri).var(ind_zeta).value, ...
        ROMS.theta_s, ROMS.theta_b, ROMS.hc, ROMS.N, 'r');
    ROMS.data(diri).zv = zlevs(ROMS.Vtransform, ROMS.Vstretching, ...
        ROMS.data(diri).hv, ROMS.data(diri).var(ind_zeta).value(1:end-1), ...
        ROMS.theta_s, ROMS.theta_b, ROMS.hc, ROMS.N, 'r');
    ROMS.data(diri).zwu = zlevs(ROMS.Vtransform, ROMS.Vstretching, ...
        ROMS.data(diri).hu, ROMS.data(diri).var(ind_zeta).value, ...
        ROMS.theta_s, ROMS.theta_b, ROMS.hc, ROMS.N, 'w');
    ROMS.data(diri).zwv = zlevs(ROMS.Vtransform, ROMS.Vstretching, ...
        ROMS.data(diri).hv, ROMS.data(diri).var(ind_zeta).value(1:end-1), ...
        ROMS.theta_s, ROMS.theta_b, ROMS.hc, ROMS.N, 'w');
end
for diri=3:4
    ROMS.data(diri).zu = zlevs(ROMS.Vtransform, ROMS.Vstretching, ...
        ROMS.data(diri).hu, ROMS.data(diri).var(ind_zeta).value(1:end-1), ...
        ROMS.theta_s, ROMS.theta_b, ROMS.hc, ROMS.N, 'r');
    ROMS.data(diri).zv = zlevs(ROMS.Vtransform, ROMS.Vstretching, ...
        ROMS.data(diri).hv, ROMS.data(diri).var(ind_zeta).value, ...
        ROMS.theta_s, ROMS.theta_b, ROMS.hc, ROMS.N, 'r');
    ROMS.data(diri).zwu = zlevs(ROMS.Vtransform, ROMS.Vstretching, ...
        ROMS.data(diri).hu, ROMS.data(diri).var(ind_zeta).value(1:end-1), ...
        ROMS.theta_s, ROMS.theta_b, ROMS.hc, ROMS.N, 'w');
    ROMS.data(diri).zwv = zlevs(ROMS.Vtransform, ROMS.Vstretching, ...
        ROMS.data(diri).hv, ROMS.data(diri).var(ind_zeta).value, ...
        ROMS.theta_s, ROMS.theta_b, ROMS.hc, ROMS.N, 'w');
end
for diri=1:4
    ROMS.data(diri).dzr = ROMS.data(diri).zw(2:end,:,:) - ROMS.data(diri).zw(1:end-1,:,:);
    ROMS.data(diri).dzu = ROMS.data(diri).zwu(2:end,:,:) - ROMS.data(diri).zwu(1:end-1,:,:);
    ROMS.data(diri).dzv = ROMS.data(diri).zwv(2:end,:,:) - ROMS.data(diri).zwv(1:end-1,:,:);
end
% zr = zlevs(Vtransform,Vstretching,h,zeta,theta_s,theta_b,hc,N,'r');
% zu = rho2u_3d(zr);
% zv = rho2v_3d(zr);
% zw = zlevs(Vtransform,Vstretching,h,zeta,theta_s,theta_b,hc,N,'w');
% dzr = zw(2:end,:,:)-zw(1:end-1,:,:);
% dzu = rho2u_3d(dzr);
% dzv = rho2v_3d(dzr);

%
%
% Vertical interpolation in case of bry files
%
%
% pcolor(squeeze(ROMS.data(diri).var(vari).value))
% pcolor(squeeze(OGCM.data(diri).var(vari).value))
ROMS = vinterp_OGCM_auto(ROMS,  OGCM.depth, ind_temp, ind_salt, ind_u, ind_v);


%--------------------------------------------------------------

%
%  fill the files
%
%
% Boundary file
%
for diri=1:4
    if ~isempty(ROMS.brync)
        for vari=[ind_zeta,10,11]
            ROMS.brync{[ROMS.var(vari).name, '_', ROMS.data(diri).direction]}(1,:) = ...
                squeeze(ROMS.data(diri).var(vari).value);
        end
    end
end
for diri=1:4
    if ~isempty(ROMS.brync)
        for vari=[ind_salt, ind_temp, ind_u, ind_v]
            ROMS.brync{[ROMS.var(vari).name, '_', ROMS.data(diri).direction]}(1,:,:) = ...
                squeeze(ROMS.data(diri).var(vari).value);
        end
    end
end
% pcolor(squeeze(ROMS.data(diri).var(vari).value))
% pcolor(squeeze(OGCM.data(diri).var(vari).value))
% abc = ROMS.brync{[ROMS.var(vari).name, '_', ROMS.data(diri).direction]}(tndx_OGCM,:)
