% %  Created 12-Apr-2023 by Yong-Yub Kim
clc; clear all; close all;
warning off;

%% set path
[error_status, tmp.hostname] = system('hostname');
tmp.hostname=tmp.hostname(1:end-1);
switch tmp.hostname
    case 'Yong-Yubs-iMac-Pro.local'
        tmp.dropboxpath = '/Volumes/kyy_raid/kimyy/Dropbox';
    case {'da1', 'da2', 'da3', 'da4'}
        tmp.dropboxpath = '/mnt/lustre/proj/kimyy/Dropbox';
end
tmp.fs=filesep;
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);

%% model configuration

% dirs.hcstroot=['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/HCST_EXP/archive/atm/', cfg.var];
% dirs.obsroot=['/Volumes/kyy_raid/kimyy/Observation/GPCP/monthly_reg_cam'];
% dirs.figroot=['/Volumes/kyy_raid/kimyy/Figure/CESM2/ESP/HCST_EXP/archive/atm/', cfg.var];


% cfg.len_t_y = length(cfg.iyears);
cfg.var='photoC_TOT_zint_100m';
% cfg.var='SSH';

% grid.regions=[45 80 -30 -20];
% grid.regions=[135 200 10 20];
grid.regions=[240 290 -40 -20];
grid.regions=[210 240 30 60];
grid.regions=[160 160 20 20];

cfg.obs=f_obs_name(cfg.var);


str_regions = ['_lon',num2str(grid.regions(1), '%02i'), '_', num2str(grid.regions(2), '%02i'), ...
    '__lat',num2str(grid.regions(3),'%02i'),'_',num2str(grid.regions(4),'%02i')];

load(['/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/ts_data/', ...
    'ts_data_all_', cfg.var,str_regions,'_v01_v01_obs_', cfg.obs, '.mat'])
hold on
for mi=1:50
    plot(m_data_lens2.([cfg.var, '_ym'])(mi,:), 'g');
end
plot(mean(m_data_lens2.([cfg.var, '_ym']),1), 'g', 'linewidth', 3);
for mi=1:20
    plot(m_data_hcst.([cfg.var, '_ym']).ly1(mi,:), 'b');
    plot(m_data_assm.([cfg.var, '_ym'])(mi,:), 'r');
end
plot(mean(m_data_hcst.([cfg.var, '_ym']).ly1,1), 'b', 'linewidth', 3);
plot(mean(m_data_assm.([cfg.var, '_ym']),1), 'r', 'linewidth', 3);

plot(m_data_obs.([cfg.var, '_ym']), 'k', 'linewidth', 2);
hold off

abc(:,:,1)=rand(3,3);
abc(:,:,2)=rand(3,3);
abc(:,:,3)=rand(3,3);
abc(:,:,4)=rand(3,3);
plot(squeeze(mean(mean(abc,1),2)))
tmmean(1,1,1:4)=1;
abcd=abc-tmmean;
plot(squeeze(mean(mean(abcd,1),2)))




ly=4;
close all;
hold on
for mi=1:50
    plot(1960:2020, m_data_lens2.([cfg.var, '_ym'])(mi,:), 'g')
end
plot(1960:2020, mean( m_data_lens2.([cfg.var, '_ym']),1), 'g', 'linewidth', 3)
for mi=1:20
    plot(1960+ly-1:2020+ly-1, m_data_hcst.([cfg.var, '_ym']).(['ly',num2str(ly)])(mi,:), 'b')
    plot(1960:2020, m_data_assm.([cfg.var, '_ym'])(mi,:), 'r')
end
plot(1960+ly-1:2020+ly-1, mean(m_data_hcst.([cfg.var, '_ym']).(['ly',num2str(ly)]),1), 'b', 'linewidth', 3)
plot(1960:2020, mean(m_data_assm.([cfg.var, '_ym']),1), 'r', 'linewidth', 3)

plot(m_data_obs.([cfg.var, '_ym']), 'k', 'linewidth', 2)
hold off





load ('/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/statistics/clim_indices/clim_indices_all_AMOobs_ERSST.mat')

plot(data_AMO.time, data_AMO.obs_dseason_lp, 'color', 'k', 'linewidth', 2)
hold on
% for mi=1:20
    plot(squeeze(mean(data_AMO.assm_dseason_lp(:,:),1)), 'color','r', 'linewidth', 2)
    plot(squeeze(mean(data_AMO.hcst_dseason_lp(1,:,:),2)), 'color','b', 'linewidth', 2)
    plot(mean(data_AMO.lens2_dseason_lp(:,:),1), 'color','g', 'linewidth', 2)
hold off
% end

for mi=1:20
    hold on
    plot(data_AMO.time,squeeze(data_AMO.hcst_dseason_lp(1,mi,:)), 'color','b', 'linewidth', 2)
end



corrcoef(squeeze(mean(data_AMO.assm_dseason_lp(:,:),1)), data_AMO.obs_dseason_lp)

corrcoef(squeeze(median(data_AMO.hcst_dseason_lp(1,:,:),2)), data_AMO.obs_dseason_lp)
corrcoef(squeeze(median(data_AMO.hcst_dseason_lp(2,:,1:end-1),2)), data_AMO.obs_dseason_lp(2:end))
corrcoef(squeeze(median(data_AMO.hcst_dseason_lp(3,:,1:end-2),2)), data_AMO.obs_dseason_lp(3:end))
corrcoef(squeeze(median(data_AMO.hcst_dseason_lp(4,:,1:end-3),2)), data_AMO.obs_dseason_lp(4:end))
corrcoef(squeeze(median(data_AMO.hcst_dseason_lp(5,:,1:end-4),2)), data_AMO.obs_dseason_lp(5:end))




function obsname_simple = f_obs_name(varn)
    switch varn
        case 'SST'
            obsname_simple='ERSST';
        case 'PRECT'
            obsname_simple='GPCC';
        case 'RAIN'
            obsname_simple='GPCC';
        case 'PSL'
            obsname_simple='ERA5';
        case 'SOILWATER_10CM'
%             obsname_simple='CMEMS';
            obsname_simple='GLEAM';
        case 'TWS'
            obsname_simple='NOAA';
        case 'SSH'
            obsname_simple='CMEMS';
        case 'TS'
%             obsname_simple='HadCRUT5';
            obsname_simple='ERA5';
        case 'sumChl'
            obsname_simple='OC_CCI';
        case 'TLAI'
            obsname_simple='NOAA'; % NOAA LAI
        case 'FAREA_BURNED'
            obsname_simple='MODIS'; % MODIS Fire_cci v5.1
            obsname_simple='AVHRR'; % AVHRR-LTDR
        case 'COL_FIRE_CLOSS'
            obsname_simple='GFED'; % GFED
        case 'photoC_TOT_zint'
%             obsname_simple='VGPM'; % VGPM
            obsname_simple='CMEMS'; %Globcolour            
        case 'photoC_TOT_zint_100m'
%             obsname_simple='VGPM'; % VGPM
            obsname_simple='CMEMS'; %Globcolour
        case 'GPP'
%             obsname_simple='ORNL_DAAC';
            obsname_simple='VODCA2GPP';
        case 'TEMP'
            obsname_simple='EN4';
%             obsname_simple='projdv7.3';
        case 'SALT'
            obsname_simple='EN4';
%             obsname_simple='projdv7.3';
        otherwise
            obsname_simple='nan';
    end
end
