function status=ROMS_surface_forcing_Qair(grdfile, year, tinterval, outfiledir, testname)
% % This function is based on MATLAB 2017a

%     outfiledir = '/data1/kimyy/Model/ROMS/roms_nwp/nwp_1_20/forcing_matlab/';
%     testname = 'nwp_1_20';
        
    varname_time = strcat('Qair','_time');
    pfile = strcat(outfiledir,testname,'_',num2str(year,'%04i'),'_Pair.nc');
    tfile = strcat(outfiledir,testname,'_',num2str(year,'%04i'),'_Tair.nc');
    dfile = strcat(outfiledir,testname,'_',num2str(year,'%04i'),'_Dair.nc');
    outfile = strcat(outfiledir,testname,'_',num2str(year,'%04i'),'_Qair.nc');
    
    disp(['forcing file is the ',outfile])
    
%     totalday  = yeardays(year); 
    [stat,totalday]=system(['date +%j -d "', num2str(year+1),'0101 -1 day"']);
    totalday=str2num(totalday);
    
    lon_rho = ncread(grdfile,'lon_rho');
    lat_rho = ncread(grdfile,'lat_rho');
    data_info = ncinfo(grdfile, 'lon_rho'); 

    time=[0.5:1:totalday-0.5];
 
    ncid = netcdf.create(outfile,'CLOBBER');

    eta_rho_dimid = netcdf.defDim(ncid,'eta_rho',data_info.Size(2));
    xi_rho_dimid = netcdf.defDim(ncid, 'xi_rho', data_info.Size(1));
    eta_u_dimid = netcdf.defDim(ncid,'eta_u',data_info.Size(2));
    xi_u_dimid = netcdf.defDim(ncid, 'xi_u', data_info.Size(1)-1);
    eta_v_dimid = netcdf.defDim(ncid,'eta_v',data_info.Size(2)-1);
    xi_v_dimid = netcdf.defDim(ncid, 'xi_v', data_info.Size(1));
    eta_psi_dimid = netcdf.defDim(ncid,'eta_psi',data_info.Size(2)-1);
    xi_psi_dimid = netcdf.defDim(ncid, 'xi_psi', data_info.Size(1)-1);
    time_dimid = netcdf.defDim(ncid, varname_time, 0);

    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'type', ' ROMS Surface forcing file ');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'title', ' Bulk Formular Forcing file ');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'source', ' ECMWF ERA Interim ');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'author', 'Created by Y.Y.Kim');
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
        'date', date);

    timevarid=netcdf.defVar(ncid, varname_time, 'NC_DOUBLE', time_dimid);
    netcdf.putAtt(ncid,timevarid,'long_name','cyclic day');
    netcdf.putAtt(ncid,timevarid,'units','DAYS');
    netcdf.putAtt(ncid,timevarid,'cycle_length', totalday);

    dvarid=netcdf.defVar(ncid,'Qair', 'NC_FLOAT', [xi_rho_dimid eta_rho_dimid time_dimid]);  %% [x y t]
    netcdf.putAtt(ncid,dvarid,'long_name','Relative humidity calculated from ECMWF P, T, D');
    netcdf.putAtt(ncid,dvarid,'units','Percentage');
    netcdf.putAtt(ncid,dvarid,'time',varname_time);

    netcdf.endDef(ncid);
    
    for i=1:totalday
        % %     index to the ncread starts from 1, 
        % %     but to the netcdf.putVar starts from 0.   
        p_data = ncread(pfile, 'Pair', [1 1 i], [data_info.Size(1) data_info.Size(2) 1]);  %% [x y t]
        t_data = ncread(tfile, 'Tair', [1 1 i], [data_info.Size(1) data_info.Size(2) 1]);  %% [x y t]
        d_data = ncread(dfile, 'Dair', [1 1 i], [data_info.Size(1) data_info.Size(2) 1]);  %% [x y t]
        
        q_data = (qsat(d_data,p_data)./qsat(t_data,p_data)).*100;
        
        netcdf.putVar(ncid, timevarid, i-1, 1, time(i));
        netcdf.putVar(ncid, dvarid, [0 0 i-1], [data_info.Size(1) data_info.Size(2) 1], q_data); %%[x y t], Z'(x,y)
    end
    netcdf.close(ncid);
    status=1;
end

function q=qsat(Ta,Pa)
% QSAT: computes specific humidity at saturation. 
% q=QSAT(Ta) computes the specific humidity (kg/kg) at satuation at
% air temperature Ta (deg C). Dependence on air pressure, Pa, is small,
% but is included as an optional input.
%
%    INPUT:   Ta - air temperature  [C]
%             Pa - (optional) pressure [mb]
%
%    OUTPUT:  q  - saturation specific humidity  [kg/kg]

% Version 1.0 used Tetens' formula for saturation vapor pressure 
% from Buck (1981), J. App. Meteor., 1527-1532.  This version 
% follows the saturation specific humidity computation in the COARE
% Fortran code v2.5b.  This results in an increase of ~5% in 
% latent heat flux compared to the calculation with version 1.0.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3/8/97: version 1.0
% 4/7/99: version 1.2 (revised as above by AA)
% 8/5/99: version 2.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==1,
  as_consts;
  Pa=P_default; % pressure in mb
end;

% original code
% a=(1.004.*6.112*0.6220)./Pa;
% q=a.*exp((17.502.*Ta)./(240.97+Ta))

% as in Fortran code v2.5b for COARE
ew = 6.1121*(1.0007+3.46e-6*Pa).*exp((17.502*Ta)./(240.97+Ta)); % in mb
q  = 0.62197*(ew./(Pa-0.378*ew));                         % mb -> kg/kg
end