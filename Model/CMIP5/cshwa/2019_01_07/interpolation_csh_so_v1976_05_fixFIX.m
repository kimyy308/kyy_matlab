%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NorESM1-M
close all; clear all;  clc;

Modelname='NorESM1-M'
varname='so'
rootdir=['/data1/cshwa/tmp/history/',Modelname,'/ocean/'];
dyear=4;


for year=1976:1976
    if(mod(year,dyear)==2)
        t_start = year - round(dyear/2);  %% 1974, 1978, ...
        t_end = num2str(t_start + dyear-1, '%04i');
    else
        t_start = year-dyear+mod(year+2,dyear);
        t_end = num2str(t_start + dyear-1, '%04i');
    end
    
% %     get CMIP5 filename, grid
    t_start = num2str(t_start,'%04i');
    mo_name = [rootdir, varname,'_Omon_',Modelname,'_historical_r1i1p1_',t_start,'01-',t_end,'12.nc'];
    lon = double(ncread(mo_name,'lon'));
    lat = double(ncread(mo_name,'lat'));
    var = ncread(mo_name,varname);
    lev = ncread(mo_name,'lev');

% %   get ECMWF grid
    tic;
    airvar = 'airT'
    if(year>=1979)
        airyear = year;
    else
        airyear = 1979;
    end
    ECMWF_dir = ['/data1/temp/ECMWF_interim/',airvar,'/'];
    ECMWF_filename = [ECMWF_dir,'ECMWF_Interim_',airvar,'_',num2str(airyear,'%04i'),'.nc'];
    lon_ec = double(ncread(ECMWF_filename,'longitude'));
    lat_ec = double(ncread(ECMWF_filename,'latitude'));
    load('standard_dep.txt');
    len_stddep=length(standard_dep);
    toc;
    
    var(var >  10000)=NaN;
    var(var < -10000)=NaN;
    %%% make surface value
    if lev(1) > 0
        h(2:length(lev)+1) = lev;
        h(1) = 0;
        var_1 = var(:,:,1,:); %%[x y z t]
        var_temp(:,:,2:length(lev)+1,:) = var;
        var_temp(:,:,1,:) = var_1;
    else
        h = lev;
        var_temp = var;
    end
    gridnum=size(var_temp);
    gridnum_ec(1)=length(lon_ec);
    gridnum_ec(2)=length(lat_ec);
    yearind=1+12*(year-str2num(t_start)); %#ok<ST2NM>

    disp('put valid value on the NaN grid')
    tic;
    var_temp2=NaN(gridnum(1),gridnum(2),gridnum(3),12);
    for k=1:gridnum(3)
        for t= yearind : yearind+12-1
            sq_var_temp=double(squeeze(var_temp(:,:,k,t-yearind+1)));
            nanind=isnan(sq_var_temp);
            sq_var_temp(nanind)=griddata(lon(~nanind),lat(~nanind),sq_var_temp(~nanind),lon(nanind),lat(nanind),'nearest'); %#ok<GRIDD>
            var_temp2(:,:,k,t-yearind+1)=sq_var_temp;
        end
    end
    toc;
    disp('vertical interpolation')
    tic;
    var_temp3=NaN(gridnum(1),gridnum(2),len_stddep,12);
    for t= yearind : yearind+12-1
        for i=1:gridnum(1) %% x
            for j=1:gridnum(2) %%y
                var_temp3(i,j,:,t-yearind+1)=interp1(h,squeeze(var_temp2(i,j,:,t-yearind+1)),standard_dep,'linear');
            end
        end
    end
    toc;
    clear sq_var_temp var_temp2
    disp('horizontal linear interpolation')
    tic;
    var_temp4=NaN(gridnum_ec(1),gridnum_ec(2),len_stddep,12);
    for k= 1:len_stddep
        for t= yearind : yearind+12-1
            sq_var_temp=squeeze(var_temp3(:,:,k,t-yearind+1));
            var_temp4(:,:,k,t-yearind+1)=(griddata(lon,lat,sq_var_temp,lon_ec,lat_ec'))';
        end
    end
    toc;
    
    outputdir=[rootdir,'../../',num2str(year,'%04i'),'/'];
    if(exist(outputdir,'dir')~=7) %#ok<EXIST>
        mkdir(outputdir);
    end
    ncid = netcdf.create(strcat(outputdir,Modelname,'_',varname,'_',num2str(year,'%04i'),'.nc'),'CLOBBER');
    lon_dimid = netcdf.defDim(ncid,'lon',gridnum_ec(1));
    lat_dimid = netcdf.defDim(ncid, 'lat', gridnum_ec(2));
    dep_dimid = netcdf.defDim(ncid, 'depth', len_stddep);
    time_dimid = netcdf.defDim(ncid, 'time', 0);  

    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'type', ' NorESM1-M interpolated data');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'title', 'sea_water_salinity');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'source', 'CMIP5 data');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'author', 'Created by S.H Chae');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'date', date);

    lonvarid = netcdf.defVar(ncid, 'lon', 'NC_float', lon_dimid);
    netcdf.putAtt(ncid,lonvarid,'standard_name','lon');
    netcdf.putAtt(ncid,lonvarid,'long_name','lon');
    netcdf.putAtt(ncid,lonvarid,'units','degrees_east');

    latvarid = netcdf.defVar(ncid, 'lat', 'NC_float', lat_dimid);
    netcdf.putAtt(ncid,latvarid,'standard_name','lat');
    netcdf.putAtt(ncid,latvarid,'long_name','lat');
    netcdf.putAtt(ncid,latvarid,'units','degrees_north');

    depvarid = netcdf.defVar(ncid, 'depth', 'NC_float', dep_dimid);
    netcdf.putAtt(ncid,depvarid,'standard_name','ocean depth');
    netcdf.putAtt(ncid,depvarid,'long_name','depth_below_ocean_surface');
    netcdf.putAtt(ncid,depvarid,'units','meter');

    timevarid = netcdf.defVar(ncid, 'time', 'NC_float', time_dimid);
    netcdf.putAtt(ncid,timevarid,'standard_name','time');
    netcdf.putAtt(ncid,timevarid,'long_name','time');
    netcdf.putAtt(ncid,timevarid,'units','months from 1974-jan');

    varid = netcdf.defVar(ncid, varname, 'NC_float', [lon_dimid lat_dimid dep_dimid time_dimid]);
    
    time=1+(year-1974)*12 : 12 +(year-1974)*12;
    switch(varname)
        case 'so'
            netcdf.putAtt(ncid,varid,'standard_name','salt');
            netcdf.putAtt(ncid,varid,'long_name','sea_water_salinity');
            netcdf.putAtt(ncid,varid,'units','PSU');
        otherwise
            disp('?')
            return
    end
            

    netcdf.endDef(ncid);

    netcdf.putVar(ncid, lonvarid, 0, gridnum_ec(1), lon_ec);
    netcdf.putVar(ncid, latvarid, 0, gridnum_ec(2), lat_ec);
    netcdf.putVar(ncid, depvarid, 0, len_stddep, standard_dep);
    netcdf.putVar(ncid, timevarid, 0, length(time), time);
    netcdf.putVar(ncid, varid, [0 0 0 0], [gridnum_ec(1) gridnum_ec(2) len_stddep length(time)], var_temp4);

    netcdf.close(ncid);
end
