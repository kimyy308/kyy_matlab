clear all; close all; clc;


testname='avg_ens_10km_mean' 
filedir = strcat('E:\Data\Reanalysis\nwp_1_10_seo\', testname, '\'); % % where data files are
varname='temp';
inputyear=2006:2010;
for month=5:3:11
    for i=1:length(inputyear)
    %     for tempmonth=1:12
            tempyear=inputyear(i);
            filename01=strcat(filedir, testname, '_monthly_', num2str(tempyear,'%04i'), '_',num2str(month,'%02i'), '.nc');
            data_info = ncinfo(filename01, varname); 
            if (length(data_info.Dimensions)==3)
                temp(:,:,:,1)=squeeze(ncread(filename01,'temp',[1 1 1], [334 332 20]));
                salt(:,:,:,1)=squeeze(ncread(filename01,'salt',[1 1 1], [334 332 20]));
                u(:,:,:,1)=squeeze(ncread(filename01,'u',[1 1 1], [333 332 20]));
                v(:,:,:,1)=squeeze(ncread(filename01,'v',[1 1 1], [334 331 20]));
                zeta(:,:,1)=squeeze(ncread(filename01,'zeta',[1 1], [334 331]));
            else
                temp(:,:,:,1)=squeeze(ncread(filename01,'temp',[1 1 1 1], [334 332 20 1]));
                salt(:,:,:,1)=squeeze(ncread(filename01,'salt',[1 1 1 1], [334 332 20 1]));
                u(:,:,:,1)=squeeze(ncread(filename01,'u',[1 1 1 1], [333 332 20 1]));
                v(:,:,:,1)=squeeze(ncread(filename01,'v',[1 1 1 1], [334 331 20 1]));
                zeta(:,:,1)=squeeze(ncread(filename01,'zeta',[1 1 1], [334 331 1]));
            end
            combtemp(:,:,:,i)=temp;
            combsalt(:,:,:,i)=salt;
            combu(:,:,:,i)=u;
            combv(:,:,:,i)=v;
            combzeta(:,:,i)=zeta;
      %     end
        disp(['year : ',num2str(tempyear)])
        disp(['month : ',num2str(month)])
    end
    meantemp=mean(temp,4);
    meansalt=mean(salt,4);
    meanu=mean(u,4);
    meanv=mean(v,4);
    meanzeta=mean(zeta,3);

filename=strcat(filedir, testname, '_monthly_', num2str(min(inputyear),'%04i'), '_',num2str(max(inputyear),'%04i'), '_',num2str(month,'%02i'), '.nc')
ncid = netcdf.create(filename,'CLOBBER');

tdimid=netcdf.defDim(ncid, 'ocean_time', 0);
xdimid=netcdf.defDim(ncid, 'xi_rho', 334);
ydimid=netcdf.defDim(ncid, 'eta_rho', 332);
sdimid=netcdf.defDim(ncid, 's_rho', 20);

tvarid=netcdf.defVar(ncid, 'ocean_time','double',tdimid);
netcdf.putAtt(ncid,tvarid,'long_name','averaged time since initialization');
netcdf.putAtt(ncid,tvarid,'units','seconds since 1980-01-01 00:00:00');
netcdf.putAtt(ncid,tvarid,'calendar','gregorian');
netcdf.putAtt(ncid,tvarid,'field','time, scalar, series');

tempvarid=netcdf.defVar(ncid, 'temp','nc_float',[xdimid ydimid sdimid tdimid]);
saltvarid=netcdf.defVar(ncid, 'salt','nc_float',[xdimid ydimid sdimid tdimid]);
uvarid=netcdf.defVar(ncid, 'u','nc_float',[xdimid ydimid sdimid tdimid]);
vvarid=netcdf.defVar(ncid, 'v','nc_float',[xdimid ydimid sdimid tdimid]);
zetavarid=netcdf.defVar(ncid, 'zeta','nc_float',[xdimid ydimid tdimid]);

netcdf.endDef(ncid);
temptime=datenum(round((min(inputyear)+max(inputyear))/2),month,15) - datenum(1980,1,1);
netcdf.putVar(ncid,tvarid,0,1,temptime);
netcdf.putVar(ncid,tempvarid,[0 0 0 0],[334 332 20 1], meantemp);
netcdf.putVar(ncid,saltvarid,[0 0 0 0],[334 332 20 1], meansalt);
netcdf.putVar(ncid,uvarid,[0 0 0 0],[333 332 20 1], meanu);
netcdf.putVar(ncid,vvarid,[0 0 0 0],[334 331 20 1], meanv);
netcdf.putVar(ncid,zetavarid,[0 0 0],[334 331 1], meanzeta);
netcdf.close(ncid);
end