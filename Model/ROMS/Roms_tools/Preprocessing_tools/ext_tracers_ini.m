function ext_tracers_ini(ininame,grdname,seas_datafile,ann_datafile,...
                         dataname,vname,type,tini);
 warning off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% P. Marchesiello - 2005. Adapted from P. Penven's ext_tracers.m 
%
%  Ext tracers in a ROMS initial file
%  take seasonal data for the upper levels and annual data for the
%  lower levels
%
%  input:
%    ininame       : ROMS initial file name
%    grdname       : ROMS grid file name    
%    seas_datafile : regular longitude - latitude - z seasonal data 
%                    file used for the upper levels  (netcdf)
%    ann_datafile  : regular longitude - latitude - z annual data 
%                    file used for the lower levels  (netcdf)
%    dataname      : variable name in data file
%    vname         : variable name in ROMS file
%    type          : position on C-grid ('r', 'u', 'v', 'p')
%    tini          : initialisation time [days]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
%
% open the vertical coordinate parameter
% 
% vert_param;
romstools_param;

%
% set the value of ro (oa decorrelation scale [m]) 
% and default (value if no data)
%
ro=0;
default=NaN;
disp([' Ext tracers: ro = ',num2str(ro/1000),...
      ' km - default value = ',num2str(default)])

% Open initial file
%
nc=netcdf(ininame,'write');
theta_s = nc{'theta_s'}(:);
theta_b =  nc{'theta_b'}(:);
hc  =  nc{'hc'}(:);
N =  length(nc('s_rho'));
%
% Open and Read grid file  
% 
ng=netcdf(grdname);
lon=ng{'lon_rho'}(:);
lat=ng{'lat_rho'}(:);
h=ng{'h'}(:);
close(ng);
[M,L]=size(lon);
%
% Read seasonal datafile 
%
if (WOA_switch==0)
    % for WOA 1998
    ncseas=netcdf(seas_datafile);
    X=ncseas{'X'}(:);
    Y=ncseas{'Y'}(:);
    Zseas=-ncseas{'Z'}(:);
    T=ncseas{'T'}(:).*30;
else 
    % % for WOA 2013
    ncseas=netcdf(seas_datafile);
    X=ncseas{'LON'}(:);
    Y=ncseas{'LAT'}(:);
    Zseas=-ncseas{'DEPTH'}(:);
    T=ncseas{'TIME'}(:).*30;
end

tlen=length(T);
Nzseas=length(Zseas);

if (WOA_switch==0)
% % % Y.Y.Kim doesn't use annual data -> comment
%
% Read annual datafile
%
ncann=netcdf(ann_datafile);
Zann=-ncann{'Z'}(:);
Nz=length(Zann);
% % Y.Y.Kim doesn't use annual data -> comment
end
%
% Determine time index to process
%
ll=find(T<=tini);
if (size(ll,1) ~= 0)
 l=ll(size(ll,1));
else
 l=1;
end
disp(['   ext_tracers_ini: time index: ',num2str(l),' of total: ',num2str(tlen)])
%
% get a subgrid
%
dl=1;
lonmin=min(min(lon))-dl;
lonmax=max(max(lon))+dl;
latmin=min(min(lat))-dl;
latmax=max(max(lat))+dl;
%
j=find(Y>=latmin & Y<=latmax);
if (WOA_switch==0)
    i1=find(X-360>=lonmin & X-360<=lonmax);  %% for WOA 1998
    i2=find(X>=lonmin & X<=lonmax);
    i3=find(X+360>=lonmin & X+360<=lonmax);  %% for WOA 1998
    x=cat(1,X(i1)-360,X(i2),X(i3)+360);  %% for WOA 1998
else
    i1=find(X-1440>=lonmin & X-1440<=lonmax);  %% for WOA 2013
    i2=find(X>=lonmin & X<=lonmax);
    i3=find(X+1440>=lonmin & X+1440<=lonmax);  %% for WOA 2013
    x=cat(1,X(i1)-1440,X(i2),X(i3)+1440); %% for WOA 2013
end
y=Y(j);
%
%------------------------------------------------------------
% Horizontal interpolation
%------------------------------------------------------------
%
% % %
% % % Interpole annual dataset on horizontal ROMS grid
% % %

if (WOA_switch==0)
% %     Y.Y.Kim doesn't use annual data -> comment
    if Nz > Nzseas
      if Zseas~=Zann(1:length(Zseas)) 
        error('vertical levels dont match')
      end
      datazgrid=zeros(Nz,M,L);
      missval=ncann{dataname}.missing_value(:);
      for k=Nzseas+1:Nz
        if ~isempty(i2)
          data=squeeze(ncann{dataname}(k,j,i2));
        else
          data=[];
        end
        if ~isempty(i1)
          data=cat(2,squeeze(ncann{dataname}(k,j,i1)),data);
        end
        if ~isempty(i3)
                    disp('tt');

          data=cat(2,data,squeeze(ncann{dataname}(k,j,i3)));
        end
        data=get_missing_val(x,y,data,missval,ro,default);
        datazgrid(k,:,:)=interp2(x,y,data,lon,lat,'cubic');
      end
    end
    close(ncann);
% %     Y.Y.Kim doesn't use annual data -> comment
end


%
% interpole seasonal dataset on horizontal roms grid
%
disp(['   ext_tracers_ini: horizontal interpolation of seasonal data'])
if (WOA_switch==0)
    missval=ncseas{dataname}.missing_value(:);  %%for WOA 1998  
else
    % missval=ncseas{dataname}.FillValue(:);  %% for WOA 2013
    missval=single(9.969209968386869e+36);  %% for WOA 2013
end
if (WOA_switch==0)
    if Nz <= Nzseas
      datazgrid=zeros(Nz,M,L);
    end
else
    datazgrid=zeros(Nzseas,M,L);
end
if (WOA_switch==0)
    Nzseas = min([Nz Nzseas]);
end
for k=1:Nzseas
  if ~isempty(i2)
    data=squeeze(ncseas{dataname}(l,k,j,i2));
  else
    data=[];
  end
  if ~isempty(i1)
    data=cat(2,squeeze(ncseas{dataname}(l,k,j,i1)),data);
  end
  if ~isempty(i3)
    data=cat(2,data,squeeze(ncseas{dataname}(l,k,j,i3)));
  end
%   'kyy'
  data=get_missing_val(x,y,data,missval,ro,default);
  datazgrid(k,:,:)=interp2(x,y,data,lon,lat,'cubic');
end
close(ncseas);
%
%----------------------------------------------------
%  Vertical interpolation
%-----------------------------------------------------
%
disp('   ext_tracers_ini: vertical interpolation')
%
% Get the sigma depths
%
zroms=zlevs(Vtransform,Vstretching,h,0.*h,theta_s,theta_b,hc,N,'r');
if type=='u'
  zroms=rho2u_3d(zroms);
end
if type=='v'
  zroms=rho2v_3d(zroms);
end
zmin=min(min(min(zroms)));
zmax=max(max(max(zroms)));
%
% Check if the min z level is below the min sigma level
%    (if not add a deep layer)
%
if (WOA_switch==0)
    z=Zann;  %% for WOA 1998;
else
    z=Zseas; %% for WOA 2013
end
addsurf=max(z)<zmax;
addbot=min(z)>zmin;
if addsurf
 z=[100;z];
end
if addbot
 z=[z;-100000];
end
Nz=min(find(z<zmin));
z=z(1:Nz);
var=datazgrid; clear datazgrid;
if addsurf
  var=cat(1,var(1,:,:),var);
end
if addbot
  var=cat(1,var,var(end,:,:));
end
var=var(1:Nz,:,:);
%
% Do the vertical interpolation and write in inifile
%
nc{vname}(1,:,:,:)=ztosigma(flipdim(var,1),zroms,flipud(z));
close(nc);

return
