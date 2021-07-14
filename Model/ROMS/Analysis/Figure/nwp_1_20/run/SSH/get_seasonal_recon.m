close all; clc; clear all;

a=ncread('E:\Data\Reanalysis\CSEOF_reconstructed_SSHA\CCAR_recon_sea_level_monthly_merged.nc','ssha');
lat=ncread('E:\Data\Reanalysis\CSEOF_reconstructed_SSHA\CCAR_recon_sea_level_monthly_merged.nc','lat');
% ssha=ans(:,:,361:708);
% seasonal_ssha=reshape(ssha,[720,361,29,12]);
for i=1:361
    ssha(:,i,:)=a(:,i,:)*squeeze(sqrt(cos(lat(i))^2));
end
seasonal_ssha=reshape(ssha(:,:,1:708),[720,361,59,12]);
clim_ssha=mean(mean(mean(seasonal_ssha,1,'omitnan'),2,'omitnan'),3,'omitnan');
plot(squeeze(clim_ssha));
% plot(squeeze(mean(mean(ssha,1,'omitnan'),2,'omitnan'))*10.0);

% 
% ncread('E:\Data\Reanalysis\CSEOF_reconstructed_SSHA\weekly\CCAR_recon_sea_level_19900103_19991227_v1.nc','ssha');
% ssha=ans(:,:,:);
% seasonal_ssha=reshape(ssha,[720,361,10,52]);
% % ssha=ans(:,:,1:708);
% % seasonal_ssha=reshape(ssha,[720,361,59,12]);
% clim_ssha=mean(mean(mean(seasonal_ssha,1,'omitnan'),2,'omitnan'),3,'omitnan');
% plot(squeeze(clim_ssha));
% plot(squeeze(mean(mean(ssha,1,'omitnan'),2,'omitnan')));