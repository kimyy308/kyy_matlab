close all; clear all; clc;

testnames = {'test53', 'test54', 'test55', 'test56'};
for testnameind=1:4
    testname=testnames{testnameind};
    % 'J:\Data\Model\ROMS\nwp_1_20\test56\run\test56_NWPcmems_interped_sst_trend_1976_2005.nc'
    filename = ['D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\', testname, '_NWPcmems_interped_sst_trend_1976_2005.nc'];
    ncinfo(filename)
    lon_cmems=ncread(filename,'lon_cmems');
    lat_cmems=ncread(filename,'lat_cmems');
    interped_trend=ncread(filename, 'interped_trend');
    dA_cmems = (m_lldist([0 0.5], [0 0])*1e3)^2 * cosd(lat_cmems);
    [mask_cmems_ocean, mask_cmems_land, er_status] = Func_0010_get_mask_from_data(interped_trend);
    dA_cmems = dA_cmems .* mask_cmems_ocean;
    dA_cmems_sum = sum(dA_cmems(isfinite(dA_cmems)));
    hist_mean_trend_sst(testnameind,1)=sum(interped_trend(:,:).*dA_cmems, 'all', 'omitnan')./dA_cmems_sum;
end

testnames = {'test57', 'test58', 'test59', 'test60'};
for testnameind=1:4
    testname=testnames{testnameind};
    % 'J:\Data\Model\ROMS\nwp_1_20\test56\run\test56_NWPcmems_interped_sst_trend_1976_2005.nc'
    filename = ['D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\', testname, '_NWPcmems_interped_sst_trend_2006_2100.nc'];
    ncinfo(filename)
    lon_cmems=ncread(filename,'lon_cmems');
    lat_cmems=ncread(filename,'lat_cmems');
    interped_trend=ncread(filename, 'interped_trend');
    dA_cmems = (m_lldist([0 0.5], [0 0])*1e3)^2 * cosd(lat_cmems);
    [mask_cmems_ocean, mask_cmems_land, er_status] = Func_0010_get_mask_from_data(interped_trend);
    dA_cmems = dA_cmems .* mask_cmems_ocean;
    dA_cmems_sum = sum(dA_cmems(isfinite(dA_cmems)));
    rcp45_mean_trend_sst(testnameind,1)=sum(interped_trend(:,:).*dA_cmems, 'all', 'omitnan')./dA_cmems_sum;
end

testnames = {'test65', 'test66', 'test67', 'test68'};
for testnameind=1:4
    testname=testnames{testnameind};
    % 'J:\Data\Model\ROMS\nwp_1_20\test56\run\test56_NWPcmems_interped_sst_trend_1976_2005.nc'
    filename = ['D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\', testname, '_NWPcmems_interped_sst_trend_2006_2100.nc'];
    ncinfo(filename)
    lon_cmems=ncread(filename,'lon_cmems');
    lat_cmems=ncread(filename,'lat_cmems');
    interped_trend=ncread(filename, 'interped_trend');
    dA_cmems = (m_lldist([0 0.5], [0 0])*1e3)^2 * cosd(lat_cmems);
    [mask_cmems_ocean, mask_cmems_land, er_status] = Func_0010_get_mask_from_data(interped_trend);
    dA_cmems = dA_cmems .* mask_cmems_ocean;
    dA_cmems_sum = sum(dA_cmems(isfinite(dA_cmems)));
    rcp85_mean_trend_sst(testnameind,1)=sum(interped_trend(:,:).*dA_cmems, 'all', 'omitnan')./dA_cmems_sum;
end


testnames = {'test53', 'test54', 'test55', 'test56'};
for testnameind=1:4
    testname=testnames{testnameind};
%     test56_NWP_ssh_recon_analysis_1976_2005.nc
    filename = ['D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\', testname, '_NWP_ssh_recon_analysis_1976_2005.nc'];
    ncinfo(filename)
    lon_recon=ncread(filename,'lon');
    lat_recon=ncread(filename,'lat');
    interped_trend=ncread(filename, 'trend');
    [lat_recon2, lon_recon2]=meshgrid(lat_recon, lon_recon);
    dA_recon = (m_lldist([0 0.5], [0 0])*1e3)^2 * cosd(lat_recon2);
    [mask_recon_ocean, mask_recon_land, er_status] = Func_0010_get_mask_from_data(interped_trend);
    dA_recon = dA_recon .* mask_recon_ocean;
    dA_recon_sum = sum(dA_recon(isfinite(dA_recon)));
    hist_mean_trend_ssh(testnameind,1)=sum(interped_trend(:,:).*dA_recon, 'all', 'omitnan')./dA_recon_sum;
end

testnames = {'test57', 'test58', 'test59', 'test60'};
for testnameind=1:4
    testname=testnames{testnameind};
    filename = ['D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\', testname, '_NWPcmems_interped_ssh_trend_2006_2100.nc'];
     ncinfo(filename)
    lon_cmems=ncread(filename,'lon_cmems');
    lat_cmems=ncread(filename,'lat_cmems');
    interped_trend=ncread(filename, 'interped_trend');
    dA_cmems = (m_lldist([0 0.5], [0 0])*1e3)^2 * cosd(lat_cmems);
    [mask_cmems_ocean, mask_cmems_land, er_status] = Func_0010_get_mask_from_data(interped_trend);
    dA_cmems = dA_cmems .* mask_cmems_ocean;
    dA_cmems_sum = sum(dA_cmems(isfinite(dA_cmems)));
    rcp45_mean_trend_ssh(testnameind,1)=sum(interped_trend(:,:).*dA_cmems, 'all', 'omitnan')./dA_cmems_sum;
end

testnames = {'test65', 'test66', 'test67', 'test68'};
for testnameind=1:4
    testname=testnames{testnameind};
    filename = ['D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\', testname, '_NWPcmems_interped_ssh_trend_2006_2100.nc'];
     ncinfo(filename)
    lon_cmems=ncread(filename,'lon_cmems');
    lat_cmems=ncread(filename,'lat_cmems');
    interped_trend=ncread(filename, 'interped_trend');
    dA_cmems = (m_lldist([0 0.5], [0 0])*1e3)^2 * cosd(lat_cmems);
    [mask_cmems_ocean, mask_cmems_land, er_status] = Func_0010_get_mask_from_data(interped_trend);
    dA_cmems = dA_cmems .* mask_cmems_ocean;
    dA_cmems_sum = sum(dA_cmems(isfinite(dA_cmems)));
    rcp85_mean_trend_ssh(testnameind,1)=sum(interped_trend(:,:).*dA_cmems, 'all', 'omitnan')./dA_cmems_sum;
end