clear all
temp=ncread('temp.nc','TEMP');
x=ncread('temp.nc','XT_OCEAN');
y=ncread('temp.nc','YT_OCEAN');
time=ncread('vort2.nc','TIME1');
rel_vort=ncread('vort2.nc','SURF_VORT');
for i=1:120
    for j=1:110
        for l=1:360
            k=1;
            while(k>0)
                if (temp(i,j,k,l)>19.5);
                    k=k+1;
                elseif (isnan(temp(i,j,1,l))==1)
                    isodepth(i,j,l)=NaN;
                    k=0;
                else
                    isodepth(i,j,l)=(k-1)*10;
                    k=0;
                end
            end
        end
    end
    i
end

pi = 3.14159265358979323846;
radius = 6371.0e3;
radian = 180./pi;
omega  = 7.292e-5;
cp_ocean = 3989.24495292815;  %(J/kg/K), J =joule, kg= kilogram, K=kelvin temperature, McDougall(2002)
% 1J = 1N*m = 1kg*m^2/s^2 = 10^7 erg
% 1N = 1kg * m /s^2 
% In dimensional analysis -> F=ML/T^2
% It is equal to the energy transferred (or work done) to an object when a force of one newton 
% acts on that object in the direction of its motion through a distance of one metre (1 newton metre or N¡¤m). 
% It is also the energy dissipated as heat when an electric current of one ampere passes through a resistance of one ohm for one second. 
% It is named after the English physicist James Prescott Joule (1818?1889).
deg2m=radius/radian;
for i=1:120
    for j=1:110
sin1(i,j)  = sin(y(j)*pi/180.0);
cos1(i,j)  = cos(y(j)*pi/180.0);
% = (yu(i,j)-lat)*deg2m
beta(i,j) = 2.0*omega*cos(y(j)*pi/180.0)/radius;
f(i,j)    = 2.0*omega*sin1(i,j) + beta(i,j);
%fstar(i,j) = 2.0*omega*cos1  + beta(i,j)
    end
end
rel_vort2(:,:,:)=rel_vort(:,:,1,:);

% nccreate('isodepth.nc',{'isodepth','XT_OCEAN','YT_OCEAN','f','sfc_rel_vort'}, 'Dimensions',{'XT_OCEAN',120,'YT_OCEAN',110,'time',inf}, ...
%           'Format','classic')
% ncwrite('isodepth.nc','XT_OCEAN',x);
% ncwrite('isodepth.nc','YT_OCEAN',y);
% ncwrite('isodepth.nc','f',f);
% ncwrite('isodepth.nc','isodepth',isodepth);
% ncwrite('isodepth.nc','sfc_rel_vort',rel_vort2);

ncid= netcdf.create('isodepth.nc','NOCLOBBER');
xdimid = netcdf.defDim(ncid,'XT_OCEAN',120);
ydimid = netcdf.defDim(ncid,'YT_OCEAN',110);
tdimid = netcdf.defDim(ncid,'time',netcdf.getConstant('NC_UNLIMITED'));
xvarid = netcdf.defVar(ncid,'XT_OCEAN','double',xdimid);
yvarid = netcdf.defVar(ncid,'YT_OCEAN','double',ydimid);
tvarid = netcdf.defVar(ncid,'time','double',tdimid);
fvarid = netcdf.defVar(ncid,'f','double',[xdimid,ydimid]);
isovarid = netcdf.defVar(ncid,'isodepth','double',[xdimid,ydimid,tdimid]);
relvarid = netcdf.defVar(ncid,'sfc_rel_vort','double',[xdimid,ydimid,tdimid]);
netcdf.endDef(ncid);
netcdf.putVar(ncid,xvarid,x);
netcdf.putVar(ncid,yvarid,y);
netcdf.putVar(ncid,tvarid,0,360,time);
netcdf.putVar(ncid,fvarid,f);
netcdf.putVar(ncid,isovarid,isodepth);
netcdf.putVar(ncid,relvarid,rel_vort2);
netcdf.close(ncid);