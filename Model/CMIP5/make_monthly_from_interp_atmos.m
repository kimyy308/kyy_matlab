%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; %clc;
warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
check_file_out = false;
test_plot = false;

Experiment = 'historical';
% Variables = {'hur', 'rsds', 'ua', 'va', 'tas', 'psl'};
Variables = {'tas'};
% Models    = {'NorESM1-M'};
Models    = {'bcc-csm1-1-m', 'CanESM2', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'GFDL-ESM2G', 'IPSL-CM5A-MR', 'MIROC5', 'MPI-ESM-MR', 'MRI-CGCM3', 'NorESM1-M', 'IPSL-CM5A-LR', 'MIROC-ESM-CHEM', 'MPI-ESM-LR'};
year_o = 1976;      % the first year to count months, see definition of "time"
year_range = (1976:2005);   % for year loop, historical era
% year_range = (1976);   % for year loop, historical era
yearday = 365;      % depend on model, some models have 365 days even in a leap year.
monthday = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
sample_month = [1,8];

file_SODA = '/home/auto/IPCC/example/soda3.4.1_mn_ocean_reg_2009.nc';
var_SODA = ncread( file_SODA, 'temp', [1 1 1 1], [inf inf inf 1] );

for year = year_range
    if( leapyear( year ) ), yearday = 366;  end
for ii = 1 : length(Variables)
    Varname = Variables{ii};
var_tgt = zeros( size(var_SODA,1), size(var_SODA,2), 12 );
for jj = 1 : length(Models)
    Modelname = Models{jj};
% %     select CMIP5 file
    path_in  = [ '/home/auto/ext_hdd/CMIP5/',Varname,'/',Experiment,'/interp/',Modelname];
    fname_in = [ Varname, '_interp_', Modelname, '_', Experiment, '_', ...
                 'r1i1p1', '_', num2str(year,'%04i'), '.nc' ];
    file_in = fullfile( path_in, fname_in );
    if( ~exist( file_in, 'file' ) )
        fprintf('Source File for %s does not Exist. Continue...\n',Modelname);   continue;
    end
% %     Check whether the output file exists
    path_out = [ path_in, '/../../mon/',Modelname ];
    fname_out = [ Varname, '_mon_', Modelname, '_', Experiment, '_', ...
                  'r1i1p1', '_', num2str(year,'%04i'), '.nc' ];
    file_out = fullfile( path_out, fname_out );
    if( check_file_out && exist( file_out, 'file' ) )
        fprintf('Target File %s already exists. Continue...\n',file_out);    continue;
    end
    if( ~exist( path_out, 'dir' ) ),	mkdir(path_out);   end

    fprintf( ' %-6s %-15s %04i ', Varname, Modelname, year );
    lap_time = tic;

% %     get CMIP5 file
    var_src = ncread( file_in, Varname );
    num_ind = ~isnan( var_src );
    for month = 1:12
        monthday_range = 1+sum(monthday(1:month-1)) : sum(monthday(1:month));
        var_tgt(:,:,month) = nanmean( var_src(:,:,monthday_range),3 );
    end     % month

    lon_tgt  = ncread( file_in, 'lon' );
    lat_tgt  = ncread( file_in, 'lat' );
    time_tgt = (year - year_o)*12 + (1:12);

    clear var_src

if( test_plot )
    [lon,lat]=meshgrid(lon_tgt,lat_tgt);
    for month = sample_month
        var = var_tgt(:,:,month);
        figure('name',[Varname,'_',Modelname,'_',num2str(month)]);
        pcolor(lon, lat, var'); shading flat; colorbar; caxis([240 310]);
        drawnow
    end
    continue
end

    %%% prepare output information
    source_str = ['CMIP5 data: ',Models{jj}];
    %%% open output file
%     fprintf('%-40s','write target file')
%     tic
    ncid = netcdf.create(fullfile(path_out,fname_out),'CLOBBER');
    %%% define matrix dimensions
    dimid_lon  = netcdf.defDim( ncid, 'lon',   length( lon_tgt ) );
    dimid_lat  = netcdf.defDim( ncid, 'lat',   length( lat_tgt ) );
    dimid_time = netcdf.defDim( ncid, 'time',  0);  
    %%% write global attributes
    netcdf.putAtt( ncid,netcdf.getConstant('NC_GLOBAL'), 'type', 'interpolated monthly data');
    netcdf.putAtt( ncid,netcdf.getConstant('NC_GLOBAL'), 'title', '3D Ocean Data into the grid system of SODA3.4.1');
    netcdf.putAtt( ncid,netcdf.getConstant('NC_GLOBAL'), 'source', source_str);
    netcdf.putAtt( ncid,netcdf.getConstant('NC_GLOBAL'), 'author', 'Created by BG Kim');
    netcdf.putAtt( ncid,netcdf.getConstant('NC_GLOBAL'), 'date', date);
    %%% write variable attributes
    varid_lon = netcdf.defVar( ncid, 'lon', 'NC_float', dimid_lon);
    netcdf.putAtt( ncid, varid_lon, 'standard_name', ncreadatt(file_in,'lon','standard_name'));
    netcdf.putAtt( ncid, varid_lon, 'long_name',     ncreadatt(file_in,'lon','long_name') );
    netcdf.putAtt( ncid, varid_lon, 'units',         ncreadatt(file_in,'lon','units') );
    varid_lat = netcdf.defVar( ncid, 'lat', 'NC_float', dimid_lat);
    netcdf.putAtt( ncid, varid_lat, 'standard_name', ncreadatt(file_in,'lat','standard_name'));
    netcdf.putAtt( ncid, varid_lat, 'long_name',     ncreadatt(file_in,'lat','long_name') );
    netcdf.putAtt( ncid, varid_lat, 'units',         ncreadatt(file_in,'lat','units') );
    varid_time = netcdf.defVar( ncid, 'time', 'NC_float', dimid_time);
    netcdf.putAtt( ncid, varid_time, 'standard_name', ncreadatt(file_in,'time','standard_name'));
    netcdf.putAtt( ncid, varid_time, 'long_name',     ncreadatt(file_in,'time','long_name') );
    netcdf.putAtt( ncid, varid_time, 'units',         ncreadatt(file_in,'time','units') );
    varid_var = netcdf.defVar( ncid, Varname, 'NC_float', [dimid_lon dimid_lat dimid_time]);
    netcdf.putAtt( ncid, varid_var, 'standard_name', ncreadatt(file_in,Varname,'standard_name') );
    netcdf.putAtt( ncid, varid_var, 'long_name',     ncreadatt(file_in,Varname,'long_name') );
    netcdf.putAtt( ncid, varid_var, 'units',         ncreadatt(file_in,Varname,'units') );
    netcdf.endDef(ncid);    % end of definition
    %%% write variables
    netcdf.putVar(ncid, varid_lon,  0, length( lon_tgt ),  lon_tgt);
    netcdf.putVar(ncid, varid_lat,  0, length( lat_tgt ),  lat_tgt);
    netcdf.putVar(ncid, varid_time, 0, length( time_tgt ), time_tgt);
    netcdf.putVar(ncid, varid_var, [0 0 0], [length( lon_tgt ) length( lat_tgt ) length(time_tgt)], var_tgt);
    %%% close output file
    netcdf.close(ncid);
%     toc

    fprintf('%7.1f sec\n', toc(lap_time) );
end     % jj for Models
end     % ii for Variables
end     % year
