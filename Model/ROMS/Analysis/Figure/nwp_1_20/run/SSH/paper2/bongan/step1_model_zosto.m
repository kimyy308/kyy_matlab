%-- Calculate local themosteric sea level change
%-      using CMIP5 data of thetao and so.
%-                          created by bongan
%-                          last edit 2021-04-26
%-------------------------------------------------
close all; clearvars; %clc;
addpath(genpath('/home/bongan/Documents/MATLAB/sea_water'));
check_file_out = true;
dry_run = false;
flag_akp  = false;

% Tables = ["Omon","interp"];
% Experiments = ["historical","rcp85","rcp45"];
% Models = ["IPSL-CM5A-LR","IPSL-CM5A-MR","MPI-ESM-LR","NorESM1-M"];
Tables = ["Omon"];
Experiments = ["historical"];
Models = ["IPSL-CM5A-LR"];
TitleTgt  = 'Thermosteric Sea Level Change since January 1976';
VarTgt    = 'zosto';       % 2-D, not a global average
VarStdName = 'local thermosteric sea level change';
VarLongName= replace(VarStdName,{' ','-'},'_');
VarUnit    = 'm';

for Table = Tables
for Experiment = Experiments
for Model = Models
    if( contains(Experiment,'rcp') )
            base_dir = '/data1/CMIP/cmip5/rcp_extHDD/CMIP5';
            yr_range = (2006:2100);
    else,   base_dir = '/data1/CMIP/cmip5/historical_extHDD/CMIP5';
            yr_range = (1976:2005);
    end
    fprintf( '== %s %s %s %s %04i-%04i ==\n', VarTgt, Table, Experiment, Model, yr_range(1), yr_range(end) );

    path_out = fullfile(VarTgt,Experiment);
    fname_out = sprintf('%s_%s_%s_%s_%04i-%04i.nc', ...
                        VarTgt,Table,Model,Experiment,yr_range(1),yr_range(end));
    %------ Check whether the output file exists
    file_out = fullfile( path_out, fname_out );
    if( check_file_out && exist( file_out, 'file' ) && ~dry_run )
        try
            var_tmp = ncread(file_out,VarTgt,[1 1 length(yr_range)*12],[inf inf 1]);
            if( sum(~isnan(var_tmp),'all')>0 )
                fprintf('Target File %s already exists. Continue...\n',file_out);    continue;
            end
        catch
            fprintf('Target File %s is broken. Rewrite...\n',file_out);
        end
    end

%--- Define input grid
fprintf('  Define grid information..'); lap_time = tic;

path_in_S  = fullfile(base_dir,'so',Experiment,Table,Model);
path_in_T  = fullfile(base_dir,'thetao',Experiment,Table,Model);
file_ref  = select_file( dir(fullfile(path_in_S,['so','*'])), yr_range(1) );
fname_ref = fullfile( path_in_S, file_ref.fname );

%- cut lon, lat for regional domain
[lon,lonname] = ncread_var( fname_ref, {'longitude','lon','nav_lon'} );
[lat,latname] = ncread_var( fname_ref, {'latitude','lat','nav_lat'} );
if( flag_akp && strcmp(Table,'interp') )
    refpolygon = akp4polygon;
    [~,i_min_lon] = min(abs( lon - (min(refpolygon(:,1))-1) ));
    [~,i_max_lon] = min(abs( lon - (max(refpolygon(:,1))+1) ));
    [~,i_min_lat] = min(abs( lat - (min(refpolygon(:,2))-1) ));
    [~,i_max_lat] = min(abs( lat - (max(refpolygon(:,2))+1) ));
    lon = ncread(fname_ref,lonname, i_min_lon, i_max_lon - i_min_lon + 1 );
    lat = ncread(fname_ref,latname, i_min_lat, i_max_lat - i_min_lat + 1 );
else
    i_min_lon = 1;
    i_min_lat = 1;
end

%- make lon, lat 2-D
if( length( lon ) == numel( lon ) )     % if lon is a 1-D array
    [lon, lat] = ndgrid( lon, lat );
end
[n_lon, n_lat] = size( lon );

lev = ncread_var( fname_ref, {'olevel','lev','depth','deptht'} );
if( strcmp(Table,'Omon') )
    lev_bnds = ncread_var( fname_ref, {'olevel_bnds','lev_bnds','depth_bnds','deptht_bnds'} )';
elseif( strcmp(Table,'interp') )
    lev_bnds = [0; 0.5*(lev(1:end-1)+lev(2:end)); 5500];
    lev_bnds = [lev_bnds(1:end-1), lev_bnds(2:end)];
end
n_lev = length( lev );

fprintf('%7.1f sec\n', toc(lap_time) );

%--- Define reference rho
fprintf('  Calculate reference state..'); lap_time = tic;

yr = 1976;
%------ read S
file_in = select_file( dir(fullfile('/data1/CMIP/cmip5/historical_extHDD/CMIP5','so','historical',Table,Model,['so','*'])), yr );
if( isempty(file_in) ), error('Source File for %04i does not Exist.',yr); end
fname_in = fullfile( file_in.folder, file_in.fname );
t_offset = 12 * ( yr - file_in.fyr_start );
time_src = ncread( fname_in, 'time', 1+t_offset, 12 );
S_src    = ncread( fname_in, 'so', [i_min_lon,i_min_lat,1,1+t_offset], [n_lon,n_lat,inf,12] );
S_src    = mean( S_src, 4 );
%------ read T
file_in = select_file( dir(fullfile('/data1/CMIP/cmip5/historical_extHDD/CMIP5','thetao','historical',Table,Model,['thetao','*'])), yr );
if( isempty(file_in) ), error('Source File for %04i does not Exist.',yr); end
fname_in = fullfile( file_in.folder, file_in.fname );
t_offset = 12 * ( yr - file_in.fyr_start );
PT_src   = ncread( fname_in, 'thetao', [i_min_lon,i_min_lat,1,1+t_offset], [n_lon,n_lat,inf,12] );
PT_src   = PT_src - 273.15;     % K -> degC
PT_src   = mean( PT_src, 4 );
%------ read Area
if( strcmp(Table,'Omon') )
    file_A = dir(fullfile('/data1/CMIP/cmip5/historical_extHDD/CMIP5/areacello/piControl/fx',Model,'areacello*'));
    dA     = ncread(fullfile(file_A.folder,file_A.name),'areacello');
elseif( strcmp(Table,'interp') )
    dA = (deg2km(0.5)*1e3)^2 * cosd(lat);
end

%------ Calculate rho_0
dz = diff(lev_bnds,1,2);
dz = reshape(repmat(dz', [n_lon*n_lat,1]),[n_lon,n_lat,n_lev]);
P  = reshape(repmat(lev',[n_lon*n_lat,1]),[n_lon,n_lat,n_lev]);
T_src   = sw_temp( S_src, PT_src, P, zeros(size(P)) );
rho_src = sw_dens( S_src, T_src, P );
is_ocean = ~isnan( rho_src );
dep_0    = sum( is_ocean.*dz, 3 );                  % 2-D
vol_0    = nansum( is_ocean.*dz.*dA, 'all' );       % global sum, 0-D
area_0   = nansum( is_ocean(:,:,1).*dA, 'all' );    % global sum, 0-D
rho_0_va = nansum( rho_src.*dz, 3 ) ./ dep_0;       % vertical(local) average, 2-D
rho_0_ga = nansum( rho_src.*dz.*dA, 'all' ) ./ vol_0;   % global average

fprintf('%7.1f sec\n', toc(lap_time) );

%====== Open output file
if( ~dry_run )
if( ~exist( path_out, 'dir' ) ),	mkdir(path_out);   end
ncid = netcdf.create(file_out,'CLOBBER');
%%% define matrix dimensions
dimid_lon  = netcdf.defDim( ncid, 'lon',   n_lon );
dimid_lat  = netcdf.defDim( ncid, 'lat',   n_lat );
% dimid_dep  = netcdf.defDim( ncid, 'dep',   dim_tgt(3) );
dimid_time = netcdf.defDim( ncid, 'time',  netcdf.getConstant('NC_UNLIMITED'));  
%%% write global attributes
netcdf.putAtt( ncid,netcdf.getConstant('NC_GLOBAL'), 'type', strcat(Model,' derived data'));
netcdf.putAtt( ncid,netcdf.getConstant('NC_GLOBAL'), 'title', TitleTgt);
netcdf.putAtt( ncid,netcdf.getConstant('NC_GLOBAL'), 'source', 'CMIP5 data');
netcdf.putAtt( ncid,netcdf.getConstant('NC_GLOBAL'), 'author', 'Created by BG Kim');
netcdf.putAtt( ncid,netcdf.getConstant('NC_GLOBAL'), 'date', date);
%%% write variable attributes
varid_lon = netcdf.defVar( ncid, 'lon', 'NC_float', [dimid_lon dimid_lat]);
    netcdf.putAtt( ncid, varid_lon, 'standard_name', ncreadatt(fname_ref,lonname,'standard_name') );
    netcdf.putAtt( ncid, varid_lon, 'long_name',     ncreadatt(fname_ref,lonname,'long_name') );
    netcdf.putAtt( ncid, varid_lon, 'units',         ncreadatt(fname_ref,lonname,'units') );
try,netcdf.putAtt( ncid, varid_lon, 'axis',          ncreadatt(fname_ref,lonname,'axis') ); end
varid_lat = netcdf.defVar( ncid, 'lat', 'NC_float', [dimid_lon dimid_lat]);
    netcdf.putAtt( ncid, varid_lat, 'standard_name', ncreadatt(fname_ref,latname,'standard_name') );
    netcdf.putAtt( ncid, varid_lat, 'long_name',     ncreadatt(fname_ref,latname,'long_name') );
    netcdf.putAtt( ncid, varid_lat, 'units',         ncreadatt(fname_ref,latname,'units') );
try,netcdf.putAtt( ncid, varid_lat, 'axis',          ncreadatt(fname_ref,latname,'axis') ); end
% varid_dep = netcdf.defVar( ncid, 'depth', 'NC_float', dimid_dep);
%     netcdf.putAtt( ncid, varid_dep, 'standard_name', ncreadatt(file_WOD,'depth','standard_name') );
%     netcdf.putAtt( ncid, varid_dep, 'units',         ncreadatt(file_WOD,'depth','units') );
%     netcdf.putAtt( ncid, varid_dep, 'axis',          ncreadatt(file_WOD,'depth','axis') );
varid_time = netcdf.defVar( ncid, 'time', 'NC_float', dimid_time);
    netcdf.putAtt( ncid, varid_time, 'standard_name', ncreadatt(fname_ref,'time','standard_name') );
try,netcdf.putAtt( ncid, varid_time, 'long_name',     ncreadatt(fname_ref,'time','long_name') ); end
    netcdf.putAtt( ncid, varid_time, 'units',         ncreadatt(fname_ref,'time','units') );
try,netcdf.putAtt( ncid, varid_time, 'calendar',      ncreadatt(fname_ref,'time','calendar') ); end
try,netcdf.putAtt( ncid, varid_time, 'axis',          ncreadatt(fname_ref,'time','axis') ); end
varid_var = netcdf.defVar( ncid, VarTgt, 'NC_float', [dimid_lon dimid_lat dimid_time]);
    netcdf.putAtt( ncid, varid_var, 'standard_name', VarStdName );
try,netcdf.putAtt( ncid, varid_var, 'long_name', VarLongName ); end
    netcdf.putAtt( ncid, varid_var, 'units',         VarUnit );
netcdf.endDef(ncid);    % end of definition
%%% write dimensional variables
netcdf.putVar(ncid, varid_lon, [0 0], [n_lon n_lat], lon);
netcdf.putVar(ncid, varid_lat, [0 0], [n_lon n_lat], lat);
% netcdf.putVar(ncid, varid_dep, 0, dim_tgt(3), lev_tgt);
end

fprintf( '  Calculate n-th state..' ); lap_time = tic; lap_time_l = tic;
zosto = NaN( [n_lon, n_lat, length(yr_range)*12] );
zostoga = NaN( [1, length(yr_range)*12] );
time    = NaN( [1, length(yr_range)*12] );
for i_tot = 1:length(yr_range)*12
    i_yr = ceil(i_tot/12);
    i_mn = i_tot - 12*floor((i_tot-1)/12);
    yr = yr_range(i_yr);
    nchar(1) = fprintf('  %d-%02d..',yr,i_mn); elapsed = toc(lap_time_l);
    nchar(2) = fprintf(' %.0f sec. (%.0f sec. left)', elapsed, elapsed*(length(yr_range)*12-i_tot+1)/i_tot );
%------ read S
    file_in = select_file( dir(fullfile(path_in_S,['so','*'])), yr );
    if( isempty(file_in) ), error('Source File for %04i does not Exist.',yr); end
    fname_in = fullfile( path_in_S, file_in.fname );
    t_offset = 12 * ( yr - file_in.fyr_start );
    S_src    = ncread( fname_in, 'so', [1 1 1 i_mn+t_offset], [inf inf inf 1] );
    time_src = ncread( fname_in, 'time', i_mn+t_offset, 1 );
%------ read T
    file_in = select_file( dir(fullfile(path_in_T,['thetao','*'])), yr );
    if( isempty(file_in) ), error('Source File for %04i does not Exist.',yr); end
    fname_in = fullfile( path_in_T, file_in.fname );
    t_offset = 12 * ( yr - file_in.fyr_start );
    PT_src    = ncread( fname_in, 'thetao', [1 1 1 i_mn+t_offset], [inf inf inf 1] );
    PT_src   = PT_src - 273.15;     % K -> degC

%------ Define output time : days from 1850-01-01
    time_tgt = time_src;

%------ Calculate rho
    dz = diff(lev_bnds,1,2);
    dz = reshape(repmat(dz', [n_lon*n_lat,1]),[n_lon,n_lat,n_lev]);
    P  = reshape(repmat(lev',[n_lon*n_lat,1]),[n_lon,n_lat,n_lev]);
    T_src   = sw_temp( S_src, PT_src, P, zeros(size(P)) );
    rho_src = sw_dens( S_src, T_src, P );
    rho_n_va = nansum( rho_src.*dz, 3 ) ./ dep_0;      % vertical(local) average
    rho_n_ga = nansum( rho_src.*dz.*dA, 'all' ) ./ vol_0;  % global average
    zosto(:,:,i_tot) = dep_0       .* ( 1 - rho_n_va ./ rho_0_va );
    zostoga(i_tot)   = vol_0/area_0 * ( 1 - rho_n_ga  / rho_0_ga );
    time(i_tot)      = time_tgt;
    var_tgt = zosto(:,:,i_tot);

%--- Test Plot
    if( dry_run )
    lon_lim = [115,164]; %lon : 115 ~164
    lat_lim = [15,52];   %lat : 15 ~ 52
    tmp_var = zosto(:,:,i_tot);     c_lim = [-0.2 0.2];
%     tmp_var = rho_src(:,:,1);       c_lim = [1020 1027];
%     tmp_var = rho_src(:,:,21);      c_lim = [1030 1031];
%     tmp_var = S_src(:,:,1);         c_lim = [30 36];
%     tmp_var = T_src(:,:,1);         c_lim = [0 28];
    tmp_var( lon<lon_lim(1) | lon>lon_lim(2) | lat<lat_lim(1) | lat>lat_lim(2) ) = NaN;
    pcolor(lon,lat,tmp_var); shading flat; xlim(lon_lim); ylim(lat_lim); colorbar; caxis(c_lim);
    text(0.02,0.98,sprintf('%s\n%d-%02d',Model,yr,i_mn),'verticalalignment','top','horizontalalignment','left','units','normalized','fontsize',14);
    drawnow;
    end

%====== Write output file
    if( ~dry_run )
    netcdf.putVar(ncid, varid_time, i_tot-1, 1, time_tgt);
    netcdf.putVar(ncid, varid_var, [0 0 i_tot-1], [n_lon n_lat 1], var_tgt);
    end
    fprintf(repmat('\b',1,sum(nchar)));
end     % i_tot

fprintf('%7.1f sec\n', toc(lap_time) );

%%% close output file
if( ~dry_run ), netcdf.close(ncid); end

%------ save calculated zostoga
save(fullfile('zostoga',Experiment,sprintf('%s_%s_%s_%s_%04i-%04i.mat','zostoga',Table,Model,Experiment,yr_range(1),yr_range(end))), ...
    'zostoga','time');

end     % for Models
end     % for Experiments
end     % for Tables

%% function definition
function finfo = select_file( list, yr )

%     flag_exist = false;
    finfo = struct([]);
%     fyr_start = -1;
    for ll = 1 : length( list )
        fname = list(ll).name;
        fname_split = strsplit( fname, {'_','.'} );
        fyr_str   = strsplit( fname_split{end-1}, '-' );
        fyr_start = str2num( fyr_str{1}(1:4) );
        if( length( fyr_str ) > 1 )
            fyr_end   = str2num( fyr_str{2}(1:4) );
        else
            fyr_end   = fyr_start;
        end
        if( yr >= fyr_start && yr <= fyr_end )
            finfo = struct;
            finfo.folder     = list(ll).folder;
            finfo.fname      = list(ll).name;
            finfo.Table      = fname_split{2};
            finfo.Model      = fname_split{3};
            finfo.Experiment = fname_split{4};
            finfo.Variant    = fname_split{5};
            finfo.Grid       = fname_split{6};
            finfo.fyr_start  = fyr_start;
            finfo.fyr_end    = fyr_end;
            return;
        end
    end

end

function [var,varname] = ncread_var( fname, varnames, varargin)
    finfo   = ncinfo( fname );
    kk = find(ismember({finfo.Variables.Name},varnames));
    if( ~isempty(kk) ),	varname = finfo.Variables(kk).Name;
    else,               error('Can''t find variable name');
    end
    if( nargin > 2 ),	var = ncread( fname, varname, varargin );
    else,               var = ncread( fname, varname );
    end
end

function polygon = akp4polygon
% %  around Korean Peninsula polygon
    polygon = ...
           [117.0, 30.0;
            117.0, 52.0;
            142.5, 52;
            142.3, 47;
            142, 46.5;
            142, 45;
            142, 43;
            141, 43;
            141, 42.8;
            140.2, 42.6;
            140.2, 42.2;
            140.4, 41.8;
            140.5, 41;
            140.5, 38;
            137, 36;
            136, 35;
            133, 35;
            132, 34;
            131, 34;
            131, 33;
            131.0, 32.0;
            128.4, 30.0]; 
end

function polygon = nwppolygon
% % North Western Pacific polygon
    polygon = ...
        [115.0, 15.0;
         115.0, 52.0;
         164.0, 52.0;
         164.0, 15.0];
end
