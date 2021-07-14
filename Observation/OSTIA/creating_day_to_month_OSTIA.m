clear all;clc;close all;

%read reference setting file
% filepath1 = '/home01/kimyy/1992/'
% filepath2 = '/home01/kimyy/1992/monthly'
filepath = 'E:\Data\Observation\ostia\'
filepath2 = 'E:\Data\Observation\ostia\monthly\'

% filepath1 = '/data1/kimyy/Observation/avhrr/L4_OISST/monthly/podaac.jpl.nasa.gov/allData/avhrr/L4/reynolds_oi/v2/daily/1992/';%daily folder
% filepath2 ='/data1/kimyy/Observation/avhrr/L4_OISST/monthly/podaac.jpl.nasa.gov/allData/avhrr/L4/reynolds_oi/v2/daily/1992/monthly/';   %monthly folder
name1 = '-UKMO-L4HRfnd-GLOB-v01-fv02-OSTIA';
outname = 'OSTIA_monthly';
file1 = strcat(filepath,num2str(2007),'/', '20070101-UKMO-L4HRfnd-GLOB-v01-fv02-OSTIA.nc');

% % read longitude, latitude
row_lat = ncread(file1,'lat');
row_lon = ncread(file1,'lon');
% row_temp = ncread(file1,'sst');

% nc = netcdf(file1);
% row_lat = nc{'lat'}(:); %unit: degrees_E  
% row_lon = nc{'lon'}(:);%unit: degrees_N 
% %row_temp= nc{'temp'}(:);
% clear nc;

%creating nc file : yearly temp
ind=1;
for year = 2007 : 2009
    filepath1 = [filepath, num2str(year,'%04i'),'\'];
    for i = 1:12  % month
        disp([num2str(year),num2str(i,'%02i')])
%         initialization
        temp_sum = zeros(length(row_lon),length(row_lat));
        err_sum = zeros(length(row_lon),length(row_lat));
        time_sum = 0;
        
        %define name1(month name)
        name2 = num2str(i,'%02i');
        
        %define mx
        if (i == 2)
            mx = 28;
            if (mod(year,4)==0)
                mx =29;
                if (mod(year,100)==0)
                    mx = 28;
                end
            end
        elseif i==4 || 6 || 9 || 11;  mx=30;
        elseif i==1 || 3 || 5 || 7 || 8 ||10||12; mx=31;
        end
        for j = 1:mx
            name3 = num2str(j,'%02i');
            file_name = strcat(filepath1,num2str(year,'%04i'),num2str(i,'%02i'),num2str(mx,'%02i'),name1,'.nc');
    %         %make sst data
            sst = ncread(file_name,'analysed_sst');
            err = ncread(file_name,'analysis_error');
            time = ncread(file_name, 'time');
    %         sst = nc{'sst'}(1,1,:,:);
%             sst = 0.01 * sst + 0 ; %scale factor, add offset 
%             err = 0.01 * err + 0 ; %scale factor, add offset 
            sst(sst<=-30) = nan;  %missing process

            err(err<=-30) = nan;  %missing process
            err(err>=100) = nan;  %missing process
            %%%sst = sst + 0 ;  
            temp_sum = temp_sum + sst;
            err_sum = err_sum + err;
            time_sum = time_sum + time;
            clear nc;
        end
        temp_monthly_mean = temp_sum / mx;
        err_monthly_mean = err_sum / mx;
        time_monthly_mean = time_sum /mx;
        comb_temp(:,:,ind)=temp_monthly_mean(:,:) - 273.15;
        comb_temp(comb_temp>=100) = nan;  %missing process
        comb_err(:,:,ind)=err_monthly_mean(:,:);
        comb_err(comb_err>=100) = nan;  %missing process
        comb_time(ind)=time_monthly_mean;
        
        ind=1;
        %--------------------creating .nc code-------------------------
        nx = 7200 ; ny = 3600;
        filenc=strcat(filepath2, outname,'_', num2str(year), num2str(i,'%02i'), '.nc');
        ncid = netcdf.create(filenc,'NC_WRITE');
         % Crear dimensiones
        dimid_eta_rho = netcdf.defDim(ncid,'lat',ny);
        dimid_xi_rho = netcdf.defDim(ncid,'lon',nx); 
        dimid_time= netcdf.defDim(ncid,'time',0); 
         % Crear variables atributos
        varid_lon = netcdf.defVar(ncid,'lon','double',[dimid_xi_rho]);
        netcdf.putAtt(ncid,varid_lon,'long_name','longitude')
        netcdf.putAtt(ncid,varid_lon,'units','degrees_east')
         % 
        varid_lat = netcdf.defVar(ncid,'lat','double',[dimid_eta_rho]);
        netcdf.putAtt(ncid,varid_lat,'long_name','latitude')
        netcdf.putAtt(ncid,varid_lat,'units','degrees_north')
         % 
        varid_time = netcdf.defVar(ncid,'time','double',[dimid_time]);
        netcdf.putAtt(ncid,varid_time,'long_name','reference time of sst field')
        netcdf.putAtt(ncid,varid_time,'units','seconds since 1981-01-01 00:00:00')
        % 
        varid_temp = netcdf.defVar(ncid,'temp','float',[dimid_xi_rho, dimid_eta_rho, dimid_time]);
        netcdf.putAtt(ncid,varid_temp,'long_name','Monthly Temperature')
        netcdf.putAtt(ncid,varid_temp,'units','Celsius')
        % 
        varid_err = netcdf.defVar(ncid,'err','double',[dimid_xi_rho, dimid_eta_rho, dimid_time]);
        netcdf.putAtt(ncid,varid_err,'long_name','Estimated error standard deviation of analysed_sst')
        netcdf.putAtt(ncid,varid_err,'units','Celsius')
        netcdf.endDef(ncid)

         % % % Agregar datos de coordenadas
        netcdf.putVar(ncid,varid_lon,row_lon');
        netcdf.putVar(ncid,varid_lat,row_lat');
        netcdf.putVar(ncid,varid_time,0,1,comb_time);
        netcdf.putVar(ncid,varid_temp,comb_temp);
        netcdf.putVar(ncid,varid_err,comb_err);
        netcdf.close(ncid);
    end
end



