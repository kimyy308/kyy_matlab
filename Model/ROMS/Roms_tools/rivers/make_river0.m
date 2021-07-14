!copy roms_eas10_river.nc roms_eas10_river0_kcs.nc
nc=netcdf('roms_eas10_river0_kcs.nc','write');
tr=nc{'river_transport'}(:)
tr2=tr(:,1)*0;
tr=[tr2,tr(:,2)]
nc{'river_transport'}(:)=tr;
close(nc)
