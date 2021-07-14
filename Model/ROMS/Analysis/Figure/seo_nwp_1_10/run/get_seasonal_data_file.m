clear all; close all; clc;


testname='avg_ens_10km_mean' 
filedir = strcat('E:\Data\Reanalysis\nwp_1_10_seo\', testname, '\'); % % where data files are
varname='temp';
for tempyear=1981:2010
%     for tempmonth=1:12
        
        filename12=strcat(filedir, testname, '_monthly_', num2str(tempyear-1,'%04i'), '_12', '.nc');
        filename01=strcat(filedir, testname, '_monthly_', num2str(tempyear,'%04i'), '_01', '.nc');
        filename02=strcat(filedir, testname, '_monthly_', num2str(tempyear,'%04i'), '_02', '.nc');
        data_info = ncinfo(filename01, varname); 
        if (length(data_info.Dimensions)==3)
            temp(:,:,:,1)=squeeze(ncread(filename01,'temp',[1 1 1], [334 332 20]));
            u(:,:,:,1)=squeeze(ncread(filename01,'u',[1 1 1], [333 332 20]));
            v(:,:,:,1)=squeeze(ncread(filename01,'v',[1 1 1], [334 331 20]));
        else
            temp(:,:,:,1)=squeeze(ncread(filename01,'temp',[1 1 1 1], [334 332 20 1]));
            u(:,:,:,1)=squeeze(ncread(filename01,'u',[1 1 1 1], [333 332 20 1]));
            v(:,:,:,1)=squeeze(ncread(filename01,'v',[1 1 1 1], [334 331 20 1]));
        end
        data_info = ncinfo(filename02, varname); 
        if (length(data_info.Dimensions)==3)
            temp(:,:,:,2)=squeeze(ncread(filename02,'temp',[1 1 1], [334 332 20]));
            u(:,:,:,2)=squeeze(ncread(filename02,'u',[1 1 1], [333 332 20]));
            v(:,:,:,2)=squeeze(ncread(filename02,'v',[1 1 1], [334 331 20]));
        else
            temp(:,:,:,2)=squeeze(ncread(filename02,'temp',[1 1 1 1], [334 332 20 1]));
            u(:,:,:,2)=squeeze(ncread(filename02,'u',[1 1 1 1], [333 332 20 1]));
            v(:,:,:,2)=squeeze(ncread(filename02,'v',[1 1 1 1], [334 331 20 1]));
        end
        data_info = ncinfo(filename12, varname); 
        if (length(data_info.Dimensions)==3)
            temp(:,:,:,3)=squeeze(ncread(filename12,'temp',[1 1 1], [334 332 20]));
            u(:,:,:,3)=squeeze(ncread(filename12,'u',[1 1 1], [333 332 20]));
            v(:,:,:,3)=squeeze(ncread(filename12,'v',[1 1 1], [334 331 20]));
        else
            temp(:,:,:,3)=squeeze(ncread(filename12,'temp',[1 1 1 1], [334 332 20 1]));
            u(:,:,:,3)=squeeze(ncread(filename12,'u',[1 1 1 1], [333 332 20 1]));
            v(:,:,:,3)=squeeze(ncread(filename12,'v',[1 1 1 1], [334 331 20 1]));
        end
        meantemp=mean(temp,4);
        meanu=mean(u,4);
        meanv=mean(v,4);
        
        
        filename=strcat(filedir, testname, '_monthly_', num2str(tempyear,'%04i'), '_13', '.nc')
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
        uvarid=netcdf.defVar(ncid, 'u','nc_float',[xdimid ydimid sdimid tdimid]);
        vvarid=netcdf.defVar(ncid, 'v','nc_float',[xdimid ydimid sdimid tdimid]);
        
        netcdf.endDef(ncid);
        temptime=datenum(tempyear,1,15) - datenum(1980,1,1);
        netcdf.putVar(ncid,tvarid,0,1,temptime);
        netcdf.putVar(ncid,tempvarid,[0 0 0 0],[334 332 20 1], meantemp);
        netcdf.putVar(ncid,uvarid,[0 0 0 0],[333 332 20 1], meanu);
        netcdf.putVar(ncid,vvarid,[0 0 0 0],[334 331 20 1], meanv);
        netcdf.close(ncid);
%     end
end

for tempyear=1980:1980
%     for tempmonth=1:12
        
%         filename12=strcat(filedir, testname, '_monthly_', num2str(tempyear-1,'%04i'), '_12', '.nc');
        filename01=strcat(filedir, testname, '_monthly_', num2str(tempyear,'%04i'), '_01', '.nc');
        filename02=strcat(filedir, testname, '_monthly_', num2str(tempyear,'%04i'), '_02', '.nc');
        data_info = ncinfo(filename01, varname); 
        if (length(data_info.Dimensions)==3)
            temp(:,:,:,1)=squeeze(ncread(filename01,'temp',[1 1 1], [334 332 20]));
            u(:,:,:,1)=squeeze(ncread(filename01,'u',[1 1 1], [333 332 20]));
            v(:,:,:,1)=squeeze(ncread(filename01,'v',[1 1 1], [334 331 20]));
        else
            temp(:,:,:,1)=squeeze(ncread(filename01,'temp',[1 1 1 1], [334 332 20 1]));
            u(:,:,:,1)=squeeze(ncread(filename01,'u',[1 1 1 1], [333 332 20 1]));
            v(:,:,:,1)=squeeze(ncread(filename01,'v',[1 1 1 1], [334 331 20 1]));
        end
        data_info = ncinfo(filename02, varname); 
        if (length(data_info.Dimensions)==3)
            temp(:,:,:,2)=squeeze(ncread(filename02,'temp',[1 1 1], [334 332 20]));
            u(:,:,:,2)=squeeze(ncread(filename02,'u',[1 1 1], [333 332 20]));
            v(:,:,:,2)=squeeze(ncread(filename02,'v',[1 1 1], [334 331 20]));
        else
            temp(:,:,:,2)=squeeze(ncread(filename02,'temp',[1 1 1 1], [334 332 20 1]));
            u(:,:,:,2)=squeeze(ncread(filename02,'u',[1 1 1 1], [333 332 20 1]));
            v(:,:,:,2)=squeeze(ncread(filename02,'v',[1 1 1 1], [334 331 20 1]));
        end
%         data_info = ncinfo(filename12, varname); 
%         if (length(data_info.Dimensions)==3)
%             temp(:,:,:,3)=squeeze(ncread(filename12,'temp',[1 1 1], [334 332 20]));
%             u(:,:,:,3)=squeeze(ncread(filename12,'u',[1 1 1], [334 332 20]));
%             v(:,:,:,3)=squeeze(ncread(filename12,'v',[1 1 1], [334 332 20]));
%         else
%             temp(:,:,:,3)=squeeze(ncread(filename12,'temp',[1 1 1 1], [334 332 20 1]));
%             u(:,:,:,3)=squeeze(ncread(filename12,'u',[1 1 1 1], [334 332 20 1]));
%             v(:,:,:,3)=squeeze(ncread(filename12,'v',[1 1 1 1], [334 332 20 1]));
%         end
        meantemp=mean(temp,4);
        meanu=mean(u,4);
        meanv=mean(v,4);
        
        
        filename=strcat(filedir, testname, '_monthly_', num2str(tempyear,'%04i'), '_13', '.nc')
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
        uvarid=netcdf.defVar(ncid, 'u','nc_float',[xdimid ydimid sdimid tdimid]);
        vvarid=netcdf.defVar(ncid, 'v','nc_float',[xdimid ydimid sdimid tdimid]);
        
        netcdf.endDef(ncid);
        temptime=datenum(tempyear,1,15) - datenum(1980,1,1);
        netcdf.putVar(ncid,tvarid,0,1,temptime);
        netcdf.putVar(ncid,tempvarid,[0 0 0 0],[334 332 20 1], meantemp);
        netcdf.putVar(ncid,uvarid,[0 0 0 0],[333 332 20 1], meanu);
        netcdf.putVar(ncid,vvarid,[0 0 0 0],[334 331 20 1], meanv);
        netcdf.close(ncid);
%     end
end
