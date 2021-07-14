clc; close all; clear all;


cmip_model='IPSL-CM5A-LR';
ssh_cmip(1,:,:,:)=ncread(['H:\zos_interp_',cmip_model,'_rcp85_r1i1p1_2006.nc'], 'zos');
cmip_model='IPSL-CM5A-MR';
ssh_cmip(2,:,:,:)=ncread(['H:\zos_interp_',cmip_model,'_rcp85_r1i1p1_2006.nc'], 'zos');
cmip_model='NorESM1-M';
ssh_cmip(3,:,:,:)=ncread(['H:\zos_interp_',cmip_model,'_rcp85_r1i1p1_2006.nc'], 'zos');
cmip_model='MPI-ESM-LR';
ssh_cmip(4,:,:,:)=ncread(['H:\zos_interp_',cmip_model,'_rcp85_r1i1p1_2006.nc'], 'zos');
lon_cmip=ncread(['H:\zos_interp_',cmip_model,'_rcp85_r1i1p1_2006.nc'], 'lon');
lat_cmip=ncread(['H:\zos_interp_',cmip_model,'_rcp85_r1i1p1_2006.nc'], 'lat');

lon_west=find(abs(lon_cmip-115)==min(abs(lon_cmip-115)));
lon_east=find(abs(lon_cmip-164)==min(abs(lon_cmip-164)));
lat_south=find(abs(lat_cmip-15)==min(abs(lat_cmip-15)));
lat_north=find(abs(lat_cmip-52)==min(abs(lat_cmip-52)));

testname='test53';
ssh_model(1,:,:)=ncread(['H:\Data\Model\ROMS\nwp_1_20\',testname,'\run\2005\ocean_rst2.nc'], 'zeta');
testname='test54';
ssh_model(2,:,:)=ncread(['H:\Data\Model\ROMS\nwp_1_20\',testname,'\run\2005\ocean_rst2.nc'], 'zeta');
testname='test55';
ssh_model(3,:,:)=ncread(['H:\Data\Model\ROMS\nwp_1_20\',testname,'\run\2005\ocean_rst2.nc'], 'zeta');
testname='test56';
ssh_model(4,:,:)=ncread(['H:\Data\Model\ROMS\nwp_1_20\',testname,'\run\2005\ocean_rst2.nc'], 'zeta');

lon_rho=ncread(['H:\Data\Model\ROMS\nwp_1_20\',testname,'\run\2005\ocean_rst2.nc'], 'lon_rho');

for i=1:4
    ssh_east(i,1)=mean(ssh_model(i,980,:),'omitnan');
    ssh_west(i,1)=mean(ssh_model(i,1,:),'omitnan');
    ssh_south(i,1)=mean(ssh_model(i,:,1),'omitnan');
    ssh_north(i,1)=mean(ssh_model(i,:,920),'omitnan');
    
    for k=1:12
        ssh_east(i,k+1)=mean(ssh_cmip(i,lon_east(1),lat_south(1):lat_north(1),k),'omitnan');
        ssh_west(i,k+1)=mean(ssh_cmip(i,lon_west(1),lat_south(1):lat_north(1),k),'omitnan');
        ssh_south(i,k+1)=mean(ssh_cmip(i,lon_west(1):lon_east(1),lat_south(1),k),'omitnan');
        ssh_north(i,k+1)=mean(ssh_cmip(i,lon_west(1):lon_east(1),lat_north(1),k),'omitnan');
    end
end

for i=1:4
    jpgname=['D:\OneDrive - 서울대학교\MEPL\project\SSH\5th_year\figure\bndy\','bndy_',num2str(i),'_east.jpg'];
    mslplot=plot(0:12,ssh_east(i,:))
    ylabel('ssh')
    xlabel('time(month)')
    set(mslplot,'LineWidth',2);
    title(['ssh east ', num2str(i)])
    set(gca,'FontSize',20);
    saveas(gcf,jpgname,'jpg');
    
    jpgname=['D:\OneDrive - 서울대학교\MEPL\project\SSH\5th_year\figure\bndy\','bndy_',num2str(i),'_west.jpg'];
    mslplot=plot(0:12,ssh_west(i,:))
    ylabel('ssh')
    xlabel('time(month)')
    set(mslplot,'LineWidth',2);
    title(['ssh west ', num2str(i)])
    set(gca,'FontSize',20);
    saveas(gcf,jpgname,'jpg');
    
    jpgname=['D:\OneDrive - 서울대학교\MEPL\project\SSH\5th_year\figure\bndy\','bndy_',num2str(i),'_south.jpg'];
    mslplot=plot(0:12,ssh_south(i,:))
    ylabel('ssh')
    xlabel('time(month)')
    set(mslplot,'LineWidth',2);
    title(['ssh south ', num2str(i)])
    set(gca,'FontSize',20);
    saveas(gcf,jpgname,'jpg');
    
    jpgname=['D:\OneDrive - 서울대학교\MEPL\project\SSH\5th_year\figure\bndy\','bndy_',num2str(i),'_north.jpg'];
    mslplot=plot(0:12,ssh_north(i,:))
    ylabel('ssh')
    xlabel('time(month)')
    set(mslplot,'LineWidth',2);
    title(['ssh north ', num2str(i)])
    set(gca,'FontSize',20);
    saveas(gcf,jpgname,'jpg');
end