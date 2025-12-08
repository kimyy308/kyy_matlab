clear all; clc; close all;

% file=dir('*.hdf')
% fn=file.name;
% info=hdfinfo(fn)
% info.SDS
% npp=hdfread(fn,'npp');

% For 1080 by 2160 data, the grid spacing is 1/6 of a degree in both latitude and longitude.
% 1080 rows * 1/6 degree per row = 180 degrees of latitude (+90 to -90).
% 2160 columns * 1/6 degree per column = 360 degrees of longitude (-180 to +180).

% The units of the npp files are: mg C / m**2 / day




% % fn='/Volumes/kyy_raid/kimyy/Observation/VGPM/cbpm/cbpm.2005001.hdf';
% % info=hdfinfo(fn);
% % npp=hdfread(fn,'npp');
% % npp(npp<0)=NaN;
% % pcolor(flip(npp)); shading flat; colorbar;
% % 
% % 
% % fn='/Volumes/kyy_raid/kimyy/Observation/VGPM/e_vgpm/eppley.2005001.hdf';
% % info=hdfinfo(fn);
% % npp=hdfread(fn,'npp');
% % npp(npp<0)=NaN;
% % pcolor(flip(npp)); shading flat; colorbar;
% % 
% % 
% % fn='/Volumes/kyy_raid/kimyy/Observation/VGPM/s_vgpm/vgpm.2005001.hdf';
% % info=hdfinfo(fn);
% % npp=hdfread(fn,'npp');
% % npp(npp<0)=NaN;
% % pcolor(flip(npp)); shading flat; colorbar;


dataroot='/Volumes/kyy_raid/kimyy/Observation/VGPM';
data_list= {'cbpm', 'e_vgpm', 's_vgpm'};
fname_start={'cbpm', 'eppley', 'vgpm'};

lon=-180+1/12:1/6:180-1/12;
lat=-90+1/12:1/6:90-1/12;

for di=1:length(data_list)
    for yi=2003:2020
        flist=dir([dataroot,'/', data_list{di},'/*.',num2str(yi),'*.hdf']);
        ncroot=[dataroot, filesep, 'netcdf',filesep, data_list{di}];
        mkdir(ncroot);
        for mi=1:length(flist)
            fname=[dataroot, filesep, data_list{di}, filesep, flist(mi).name];
            npp=hdfread(fname,'npp');
            npp(npp<0)=NaN;
            npp=flip(npp);
            npp=npp';
            
            savename=[ncroot, filesep, data_list{di}, '_', num2str(yi),num2str(mi,'%02i'), '.nc'];

            ncid=netcdf.create(savename, 'NETCDF4');

            lon_dimid = netcdf.defDim(ncid, 'lon', length(lon));
            lat_dimid = netcdf.defDim(ncid, 'lat', length(lat));
            time_dimid = netcdf.defDim(ncid, 'time', 0);

            netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                'author', 'Created by Y.Y. Kim');
            netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                'date', date);

            timevarid=netcdf.defVar(ncid, 'time', 'NC_DOUBLE', time_dimid);
                netcdf.putAtt(ncid,timevarid,'long_name','time');
                netcdf.putAtt(ncid,timevarid,'units','days since 0000-01-01 00:00:00');
    %             netcdf.putAtt(ncid,timevarid,'bounds','time_bound');           
                netcdf.putAtt(ncid,timevarid,'calendar','noleap');

            lonvarid=netcdf.defVar(ncid, 'lon', 'NC_DOUBLE', lon_dimid);
                netcdf.putAtt(ncid,lonvarid,'long_name','longitude');
                netcdf.putAtt(ncid,lonvarid,'units','degree_east');

            latvarid=netcdf.defVar(ncid, 'lat', 'NC_DOUBLE', lat_dimid);
                netcdf.putAtt(ncid,latvarid,'long_name','latitude');
                netcdf.putAtt(ncid,latvarid,'units','degree_north');

            nppvarid=netcdf.defVar(ncid, 'npp', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                netcdf.putAtt(ncid,nppvarid,'long_name', ['net primary productivity from ',data_list{di}]);
                netcdf.defVarChunking(ncid, nppvarid, 'CHUNKED', [length(lon)/10, length(lat)/10 1]);
                netcdf.defVarDeflate(ncid, nppvarid, true, true, 1);

            netcdf.endDef(ncid);
            netcdf.putVar(ncid, timevarid, 0, 1, datenum(yi-1,mi,15)-120);
            netcdf.putVar(ncid, lonvarid, 0, length(lon), lon);
            netcdf.putVar(ncid, latvarid, 0, length(lat), lat);
            netcdf.putVar(ncid, nppvarid, [0 0 0], [length(lon), length(lat) 1], npp);
            
            netcdf.close(ncid);
            
            disp(['saved file is ', savename]);
        end
    end
end



fname='/Volumes/kyy_raid/kimyy/Observation/VGPM/netcdf_regrid/corr_cbpm.nc';
fname='/Volumes/kyy_raid/kimyy/Observation/VGPM/netcdf_regrid/corr_s_vgpm.nc';
fname='/Volumes/kyy_raid/kimyy/Observation/VGPM/netcdf_regrid/corr_e_vgpm.nc';
fname='/Volumes/kyy_raid/kimyy/Observation/VGPM/netcdf_regrid/corr_ens_vgpm.nc';

gname='/Volumes/kyy_raid/kimyy/Observation/OC_CCI/monthly_reg_pop/ocn_cesm2_grid.nc';
nppcorr=ncread(fname, 'photoC_TOT_zint_100m');
tlong=ncread(gname, 'TLONG');
tlat=ncread(gname, 'TLAT');
Func_0011_get_area_weighted_mean(nppcorr, tlong, tlat)
