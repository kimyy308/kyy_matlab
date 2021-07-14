function ext_datas_ini_MyOcean(ininame,grdname,seas_datafile, ...
                         dataname,vname,type,tini);
 warning off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Ext tracers in a ROMS initial file
% 
%    input:
%    ininame       : ROMS initial file name
%    grdname       : ROMS grid file name    
%    seas_datafile : regular longitude - latitude - z 5-day averaged data
%                    file used for the all levels  (netcdf)
%    dataname      : variable name in data file
%    vname         : variable name in ROMS file
%    type          : position on C-grid ('r', 'u', 'v')
%    tini          : initialisation time [days]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
%
% open the vertical coordinate parameter
% 
% vert_param;
romstools_param;

addpath(genpath([matlabroot,'/toolbox/matlab/imagesci/']));

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
% Read SODA3 5-day averaged data file 
%
% if (type=='r')
    ncseas=netcdf(seas_datafile);
    X=ncseas{'longitude'}(:);
    Y=ncseas{'latitude'}(:);
    Zseas=-ncseas{'depth'}(:);
    T=ncseas{'time'}(:);
% elseif (type=='u' || type=='v')
%     % % for WOA 2013
%     ncseas=netcdf(seas_datafile);
%     X=ncseas{'xu_ocean'}(:);
%     Y=ncseas{'yu_ocean'}(:);
%     Zseas=-ncseas{'st_ocean'}(:);
%     T=ncseas{'time'}(:);
% end

tlen=length(T);
Nzseas=length(Zseas);

%
% Determine time index to process
%
l=1;
disp(['   ext_tracers_ini: time index: ',num2str(l),' of total: ',num2str(tlen)])
%
% get a subgrid
%
% dl=1;
% lonmin=min(min(lon))-dl;
% lonmax=max(max(lon))+dl;
% latmin=min(min(lat))-dl;
% latmax=max(max(lat))+dl;
%
% j=find(Y>=latmin & Y<=latmax);
% % if (WOA_switch==0)
%     i1=find(X-360>=lonmin & X-360<=lonmax);  
%     i2=find(X>=lonmin & X<=lonmax);
%     i3=find(X+360>=lonmin & X+360<=lonmax); 
%     x=cat(1,X(i1)-360,X(i2),X(i3)+360);  
% y=Y(j);
x=X;
y=Y;

%
%------------------------------------------------------------
% Horizontal interpolation
%------------------------------------------------------------
%
% % %
% % % Interpole annual dataset on horizontal ROMS grid
% % %

%
% interpole seasonal dataset on horizontal roms grid
%
disp(['   ext_datas_ini_MyOcean: horizontal interpolation of seasonal data'])
info=ncinfo(seas_datafile);
% missval=ncseas{dataname}.FillValue(:);  
missval=NaN;
% missval=info.Variables(4).Attributes(5).Value;
% scale_factor=info.Variables(4).Attributes(7).Value;
% add_offset=info.Variables(4).Attributes(6).Value;

% data=permute(ncread(seas_datafile,'thetao'),[3,2,1]);
switch(dataname)
    case('zos')
        datazgrid=zeros(M,L);
%         if ~isempty(i2)
%             data=squeeze(ncseas{dataname}(l,j,i2));
%         else
%             data=[];
%         end
%         if ~isempty(i1)
%             data=cat(2,squeeze(ncseas{dataname}(l,j,i1)),data);
%         end
%         if ~isempty(i3)
%             data=cat(2,data,squeeze(ncseas{dataname}(l,j,i3)));
%         end
        data=squeeze(permute(ncread(seas_datafile,dataname),[2,1]));
        scale_factor=ncreadatt(seas_datafile,dataname,'scale_factor');
        add_offset=ncreadatt(seas_datafile,dataname,'add_offset');
%         data=squeeze(data.*scale_factor + add_offset);
        data=get_missing_val(x,y,data,missval,ro,default);
        datazgrid(:,:)=interp2(x,y,data,lon,lat,'cubic');
    
    otherwise
        datazgrid=zeros(Nzseas,M,L);
        for k=1:Nzseas
%           if ~isempty(i2)
%             data=squeeze(ncseas{dataname}(l,k,j,i2));
%           else
%             data=[];
%           end
%           if ~isempty(i1)
%             data=cat(2,squeeze(ncseas{dataname}(l,k,j,i1)),data);
%           end
%           if ~isempty(i3)
%             data=cat(2,data,squeeze(ncseas{dataname}(l,k,j,i3)));
%           end
          data=squeeze(permute(ncread(seas_datafile,dataname,[1, 1, k, 1], [721, 601, 1, 1]),[3, 2, 1]));
            scale_factor=ncreadatt(seas_datafile,dataname,'scale_factor');
            add_offset=ncreadatt(seas_datafile,dataname,'add_offset');
%             data=squeeze(data*scale_factor + add_offset);
          data=get_missing_val(x,y,data,missval,ro,default);
          if (sum(sum(isnan(data)))~=0)
            datazgrid(k,:,:)=datazgrid(k-1,:,:);
          else
            datazgrid(k,:,:)=interp2(x,y,data,lon,lat,'cubic');
          end
        end
end
close(ncseas);

switch(dataname)
    case('zos')
        nc{vname}(1,:,:)=datazgrid(:,:);
        close(nc);
    otherwise
        %
        %----------------------------------------------------
        %  Vertical interpolation
        %-----------------------------------------------------
        %
        disp('   ext_datas_ini_MyOcean : vertical interpolation')
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
        z=Zseas;
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
end

return
