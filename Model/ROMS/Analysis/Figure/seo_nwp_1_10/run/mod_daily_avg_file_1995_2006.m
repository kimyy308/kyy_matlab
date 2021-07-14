clear all; close all; clc;

testname='avg_ens_10km_mean_'   %% object file
filedir2 = strcat('/data2/kimyy/Reanalysis/nwp_1_10_seo/packed_tak_daily/'); % % where data files are

for tempyear=1995:2006
    for tempday=1:370
        filedir = strcat('/data2/kimyy/Reanalysis/nwp_1_10_seo/ens_mean_daily/', num2str(tempyear,'%04i'), '/'); % % where data files are
        filedir3 = strcat('/data2/kimyy/Reanalysis/nwp_1_10_seo/ens_mean_daily/', num2str(1990,'%04i'), '/'); % % where data files are

        filename=strcat(filedir,'avg_',num2str(tempyear,'%04i'),'_ens_mean_',num2str(tempday,'%04i'), '.nc')
        if (tempday>366)
            filename2=strcat(filedir2, testname, num2str(tempyear+1,'%04i'), '_', num2str(tempday-366,'%04i'), '.nc');
        else
            filename2=strcat(filedir2, testname, num2str(tempyear,'%04i'), '_', num2str(tempday,'%04i'), '.nc');
        end
        filename3=strcat(filedir3,'avg_',num2str(1990,'%04i'),'_ens_mean_',num2str(tempday,'%04i'), '.nc');

        reftime=squeeze(ncread(filename3,'ocean_time'));
        
        temp=squeeze(ncread(filename2,'temp'));
        salt=squeeze(ncread(filename2,'salt'));
        zeta=squeeze(ncread(filename2,'zeta'));
        ubar=squeeze(ncread(filename2,'ubar'));
        vbar=squeeze(ncread(filename2,'vbar'));
        u=squeeze(ncread(filename2,'u'));
        v=squeeze(ncread(filename2,'v'));

        ncid = netcdf.open(filename,'WRITE');
        netcdf.reDef(ncid);

        ncinfo(filename);
        tvarid=netcdf.inqVarID(ncid,'ocean_time');
        tempvarid=netcdf.inqVarID(ncid,'temp');
        saltvarid=netcdf.inqVarID(ncid,'salt');
        zetavarid=netcdf.inqVarID(ncid,'zeta');
        ubarvarid=netcdf.inqVarID(ncid,'ubar');
        vbarvarid=netcdf.inqVarID(ncid,'vbar');
        uvarid=netcdf.inqVarID(ncid,'u');
        vvarid=netcdf.inqVarID(ncid,'v');
        
        netcdf.putAtt(ncid,tvarid,'units',['seconds since ',num2str(tempyear,'%04i'), '-01-01 00:00:00']);
        netcdf.endDef(ncid);
        
        netcdf.putVar(ncid,tvarid,0,1,reftime);
        netcdf.putVar(ncid,tempvarid,[0 0 0 0],[334 332 20 1], temp);
        netcdf.putVar(ncid,saltvarid,[0 0 0 0],[334 332 20 1], salt);
        netcdf.putVar(ncid,zetavarid,[0 0 0],[334 332 1], zeta);
        netcdf.putVar(ncid,ubarvarid,[0 0 0],[333 332 1], ubar);
        netcdf.putVar(ncid,vbarvarid,[0 0 0],[334 331 1], vbar);
        netcdf.putVar(ncid,uvarid,[0 0 0 0],[333 332 20 1], u);
        netcdf.putVar(ncid,vvarid,[0 0 0 0],[334 331 20 1], v);
        
        netcdf.close(ncid);
    end
end


% for tempyear=1991:2010
%     for tempmonth=1:12
%         filename=strcat(filedir, testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');
%         ncid = netcdf.open(filename,'WRITE');
%         netcdf.reDef(ncid);
% %         if tempyear==1980
% %             tdimid=netcdf.defDim(ncid, 'ocean_time', 0);
% %         else
%             tdimid=netcdf.inqDimID(ncid,'ocean_time');
% %         end
%         netcdf.renameDim(ncid,tdimid,'ocean_time2');
%         tvarid=netcdf.inqVarID(ncid,'ocean_time');
%         netcdf.renameVar(ncid,tvarid,'ocean_time2');
%         ncinfo(filename)
% %         tdimid=netcdf.defDim(ncid, 'ocean_time', 0);
%         tvarid=netcdf.defVar(ncid, 'ocean_time','double',tdimid);
%         netcdf.putAtt(ncid,tvarid,'long_name','averaged time since initialization');
%         netcdf.putAtt(ncid,tvarid,'units','seconds since 1980-01-01 00:00:00');
%         netcdf.putAtt(ncid,tvarid,'calendar','gregorian');
%         netcdf.putAtt(ncid,tvarid,'field','time, scalar, series');
%         netcdf.endDef(ncid);
%         temptime=datenum(tempyear,tempmonth,15) - datenum(1980,1,1);
%         netcdf.putVar(ncid,tvarid,0,1,temptime);
%         netcdf.close(ncid);
%     end
% end

