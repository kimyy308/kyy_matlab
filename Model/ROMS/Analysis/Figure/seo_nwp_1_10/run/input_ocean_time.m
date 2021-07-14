clear all; close all; clc;


testname='avg_ens_10km_mean' 
filedir = strcat('E:\Data\Reanalysis\nwp_1_10_seo\', testname, '\'); % % where data files are

for tempyear=1982:1984
    for tempmonth=1:12
        filename=strcat(filedir, testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');
        ncid = netcdf.open(filename,'WRITE');
        netcdf.reDef(ncid);
%         if tempyear==1980
%             tdimid=netcdf.defDim(ncid, 'ocean_time', 0);
%         else
%             tdimid=netcdf.inqDimID(ncid,'ocean_time');
%         end
%         netcdf.renameDim(ncid,tdimid,'ocean_time2');
        ncinfo(filename)
        tdimid=netcdf.defDim(ncid, 'ocean_time', 0);
        tvarid=netcdf.inqVarID(ncid,'ocean_time');
        netcdf.renameVar(ncid,tvarid,'ocean_time2');
        tvarid=netcdf.defVar(ncid, 'ocean_time','double',tdimid);
        netcdf.putAtt(ncid,tvarid,'long_name','averaged time since initialization');
        netcdf.putAtt(ncid,tvarid,'units','seconds since 1980-01-01 00:00:00');
        netcdf.putAtt(ncid,tvarid,'calendar','gregorian');
        netcdf.putAtt(ncid,tvarid,'field','time, scalar, series');
        netcdf.endDef(ncid);
        temptime=datenum(tempyear,tempmonth,15) - datenum(1980,1,1);
        netcdf.putVar(ncid,tvarid,0,1,temptime);
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

