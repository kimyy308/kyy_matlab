clear all; close all; clc;
% % % 
% % % Read the KODC data[line_st yyyymmdd lon lat depth temp salt], 1961 - 2015 
% % % 

% filename = 'C:\Users\kyy\Desktop\KODC data\KODC1961-2015.txt';
% startRow = 2;
% 
% formatSpec = '%5f%9f%11f%11f%4f%11f%f%[^\n\r]';
% 
% fileID = fopen(filename,'r');
% 
% dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
% 
% fclose(fileID);
% 
% line_st = dataArray{:, 1};
% yyyymmdd = dataArray{:, 2};
% lon = dataArray{:, 3};
% lat = dataArray{:, 4};
% depth = dataArray{:, 5};
% temp = dataArray{:, 6};
% salt = dataArray{:, 7};
% 
% clearvars filename startRow formatSpec fileID dataArray ans;
% 
% tempstr=num2str(yyyymmdd);
% year=str2num(tempstr(:,1:4));
% month=str2num(tempstr(:,5:6));
% day=str2num(tempstr(:,7:8));
% 
% clearvars tempstr
% 
% save('KODC_rawvar.mat','year','month','day','lon','lat','line_st','depth','temp','salt')

load E:\Data\Observation\NIFS\KODC_rawvar.mat

% time = juliandate(year,month,day);


% % % stddepth info
% % %  0 m = 0 m ~ 5 m
% % %  10 m = 6 m ~ 15 m  
% % %  20 m = 16 m ~ 25 m  
% % %  30 m = 26 m ~ 40 m  
% % %  50 m = 41 m ~ 62 m  
% % %  75 m = 63 m ~ 87 m  
% % %  100 m = 88 m ~ 112 m  
% % %  125 m = 113 m ~ 137 m  
% % %  150 m = 138 m ~ 175 m  
% % %  200 m = 176 m ~ 225 m  
% % %  250 m = 226 m ~ 275 m
% % %  300 m = 276 m ~ 350 m
% % %  400 m = 351 m ~ 450 m
% % %  500 m = 451 m ~ 551 m
latsize=length(unique(lat))
lonsize=length(unique(lon))
stddepth=[0 10 20 30 50 75 100 125 150 200 250 300 400 500];
tempstddepth=[-5 0 10 20 30 50 75 100 125 150 200 250 300 400 500 600];


i=12;
% % % error in the 1993year
for tempyear=1961:2015
ind_year=find(year==tempyear);
y_month=month(ind_year);
y_temp=temp(ind_year);
y_salt=salt(ind_year);
y_lat=lat(ind_year);
y_lon=lon(ind_year);
y_depth=depth(ind_year);
    for tempmonth=1:12
        ind_month=find(y_month==tempmonth);
        if (isfinite(ind_month))
            if (i<12)
                i=i+1;
            else
                i=1;
            end
            m_time(i) = mjuliandate(tempyear,tempmonth,15);
            m_temp=y_temp(ind_month);
            m_salt=y_salt(ind_month);
            m_lat=y_lat(ind_month);
            m_lon=y_lon(ind_month);
            m_depth=y_depth(ind_month);            
            for d=1:14
                ind_depth=find(m_depth>=(tempstddepth(d+1)+tempstddepth(d))/2.0+1.0);
                d_depth1=m_depth(ind_depth);
                d_temp1=m_temp(ind_depth);
                d_salt1=m_salt(ind_depth);
                d_lat1=m_lat(ind_depth);
                d_lon1=m_lon(ind_depth);
                ind_depth2=find(d_depth1<=(tempstddepth(d+2)+tempstddepth(d+1))/2.0);
                d_depth2=m_depth(ind_depth2);
                d_temp2=m_temp(ind_depth2);
                d_salt2=m_salt(ind_depth2);
                d_lat2=m_lat(ind_depth2);
                d_lon2=m_lon(ind_depth2);
%                 if (isfinite(griddata(d_lon2,d_lat2,d_temp2,unique(lon),unique(lat)')))
                if(length(unique(d_lat2))>3 && length(unique(d_lon2))>3)
                    m_d_temp(1:179,1:374,d,i)=griddata(d_lon2,d_lat2,d_temp2,unique(lon),unique(lat)');
                    m_d_salt(1:179,1:374,d,i)=griddata(d_lon2,d_lat2,d_salt2,unique(lon),unique(lat)');
                else
                    m_d_temp(1:179,1:374,d,i)=NaN;
                    m_d_salt(1:179,1:374,d,i)=NaN;
                end
                clearvars d_depth1 d_temp1 d_salt1 d_lat1 d_lon1
                clearvars d_depth2 d_temp2 d_salt2 d_lat2 d_lon2
            end
            clearvars m_depth m_temp m_salt m_lat m_lon
        else
            if (i<12)
                i=i+1;
            else
                i=1;
            end
            m_time(i) = mjuliandate(tempyear,tempmonth,1);
            for d=1:14
                m_d_temp(1:179,1:374,d,i)=NaN;
                m_d_salt(1:179,1:374,d,i)=NaN;
            end
        end
    end
%     ncfilename=strcat('KODC_',num2str(tempyear),'_',num2str(tempmonth,'%02i'),'ts.nc');
    ncfilename=strcat('KODC_',num2str(tempyear),'_','ts.nc');
    ncid = netcdf.create(ncfilename,'CLOBBER');
    londimid = netcdf.defDim(ncid,'lon',lonsize);
    latdimid = netcdf.defDim(ncid,'lat',latsize);
    depthdimid = netcdf.defDim(ncid,'depth',length(stddepth));
    timedimid =netcdf.defDim(ncid,'Time',netcdf.getConstant('NC_UNLIMITED'));
    
    lonvarid = netcdf.defVar(ncid,'lon','NC_FLOAT',londimid);
    netcdf.putAtt(ncid,lonvarid,'long_name','longitude');
    netcdf.putAtt(ncid,lonvarid,'units','degrees_E');
    netcdf.putAtt(ncid,lonvarid,'cartesian_axis','X');
    latvarid = netcdf.defVar(ncid,'lat','NC_FLOAT',latdimid);
    netcdf.putAtt(ncid,latvarid,'long_name','latitude');
    netcdf.putAtt(ncid,latvarid,'units','degrees_N');
    netcdf.putAtt(ncid,latvarid,'cartesian_axis','Y');
    depthvarid = netcdf.defVar(ncid,'depth','NC_FLOAT',depthdimid);
    netcdf.putAtt(ncid,depthvarid,'long_name','depth');
    netcdf.putAtt(ncid,depthvarid,'units','meters');
    netcdf.putAtt(ncid,depthvarid,'cartesian_axis','Z');
    netcdf.putAtt(ncid,depthvarid,'positive','down');
    
    timevarid = netcdf.defVar(ncid,'Time','NC_DOUBLE',timedimid);
    netcdf.putAtt(ncid,timevarid,'long_name','time');
    netcdf.putAtt(ncid,timevarid,'units','days since 1858-11-17 00:00:00');
    netcdf.putAtt(ncid,timevarid,'axis','T');
    netcdf.putAtt(ncid,timevarid,'calendar_type','JULIAN');
    netcdf.putAtt(ncid,timevarid,'calendar','JULIAN');
    
    tempvarid = netcdf.defVar(ncid,'temp','NC_FLOAT',[londimid latdimid depthdimid timedimid]);
    netcdf.putAtt(ncid,tempvarid,'long_name','KODC in-situ temperature');
    saltvarid = netcdf.defVar(ncid,'salt','NC_FLOAT',[londimid latdimid depthdimid timedimid]);
    netcdf.putAtt(ncid,saltvarid,'long_name','KODC in-situ salinity');
    glvarid = netcdf.getConstant('GLOBAL');
    netcdf.putAtt(ncid,glvarid,'made by','Y.Y. Kim');
    netcdf.endDef(ncid)
    
    netcdf.putVar(ncid,lonvarid,unique(lon));
    netcdf.putVar(ncid,latvarid,unique(lat));
    netcdf.putVar(ncid,depthvarid,stddepth);
%     You must use the start and count form, 
%     becuase stupid matlab doesn't recognize dimension of unlimited variable
    netcdf.putVar(ncid,timevarid,0,12,m_time(1:12));  
    netcdf.putVar(ncid,tempvarid,[0 0 0 0],size(permute(m_d_temp,[2 1 3 4])),permute(m_d_temp,[2 1 3 4]));
    netcdf.putVar(ncid,saltvarid,[0 0 0 0],size(permute(m_d_salt,[2 1 3 4])),permute(m_d_salt,[2 1 3 4]));
    netcdf.close(ncid);
    clearvars y_year y_month y_temp y_salt y_lat y_lon y_depth m_d_temp m_d_salt
end

