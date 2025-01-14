close all; clc; clear all;

warning off;

%% set path
[error_status, tmp.hostname] = system('hostname');
tmp.hostname=tmp.hostname(1:end-1);
switch tmp.hostname
    case 'Yong-Yubs-iMac-Pro.local'
        tmp.dropboxpath = '/Volumes/kyy_raid/kimyy/Dropbox';
    case {'da1', 'da2', 'da3', 'da4'}
        tmp.dropboxpath = '/mnt/lustre/proj/kimyy/Dropbox';
    otherwise
        tmp.dropboxpath = '/mnt/lustre/proj/kimyy/Dropbox';
end
tmp.fs=filesep;
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Common', tmp.fs, 'mca']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Common', tmp.fs, 'order']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);

% cfg.varnames={'zeta','zosto','zosto_thermo','zosto_halo','sustr','svstr','ubar','vbar','Uwind','Vwind'};
cfg.varnames={'zeta','zosto','zosto_thermo','zosto_halo','ubar','vbar','Uwind','Vwind', 'wcurl', 'sustr', 'svstr', 'wscurl', 'vec'};
% cfg.varnames={'Uwind', 'Vwind', 'wcurl', 'sustr', 'svstr', 'wscurl'};
cfg.varnames={'vec', 'wvec', 'zeta','zosto','zosto_thermo','zosto_halo', 'wsvec', 'ubar','vbar','Uwind','Vwind', 'sustr', 'svstr', 'wcurl', 'wscurl'};


dir.figdir='/Users/kimyy/Desktop/backup/Research/Postdoc/03_IBS/2024_TS_KS_connectivity/Figure';



%% data load
load('/Volumes/kyy_raid/kimyy/Model/ROMS/drifter_ROMS/mat/stlee_drifter_model_depth_1982_2019_01_12.mat');
load('/Volumes/kyy_raid/kimyy/Model/ROMS/drifter_ROMS/mat/stlee_drifter_model_steric_ssh_1983_2019_01_12.mat');

xn=size(lon,1);
yn=size(lat,2);
tn=size(comb_data.zeta,3);

% h2=h;
% mask_100=h2; 
% mask_50(mask_50>50)=NaN; mask_50=mask_50./mask_5J0;
% mask_100(mask_100>100)=NaN; mask_100=mask_100./mask_100;
% mask_200(mask_200>200)=NaN; mask_200=mask_200./mask_200;

depths=[100, 200, 300];
% depths=[150];
depths=[300];

for depthi=1:length(depths)
    dd=depths(depthi);
    str_d=num2str(dd);
    masks.(['d', str_d])=h;
    masks.(['d', str_d])(masks.(['d', str_d])>dd)=NaN;
%     masks.(['d', str_d])(masks.(['d', str_d])<50)=NaN;
    masks.(['d', str_d])=masks.(['d', str_d])./masks.(['d', str_d]);
end

[comb_data.wcurl] = Func_0037_get_curl(lon,lat,comb_data.Uwind,comb_data.Vwind);
[comb_data.wscurl] = Func_0037_get_curl(lon,lat,comb_data.sustr,comb_data.svstr);


% %% annual mean
% for vi=1:length(cfg.varnames)
%     tmp.varn=cfg.varnames{vi};
%     tmp.var=reshape(comb_data.(tmp.varn), [xn yn 12 tn/12]);
%     comb_data_an.(tmp.varn)=squeeze(mean(tmp.var,3));
% end




% monrange=4:10;
% monrange=3:8;
% monrange=6:8;


% %% seasonal mean
% for lm=1:12
%     monrange=lm:lm;
%     for vi=1:length(cfg.varnames)
%         tmp.varn=cfg.varnames{vi};
%         tmp.var=reshape(comb_data.(tmp.varn), [xn yn 12 tn/12]);
%         comb_data_an.(tmp.varn)(lm,:,:,:)=squeeze(mean(tmp.var(:,:,monrange,:),3));
%     end
% end
% 
% for lm=1:10
%     monrange=lm:lm+2;
%     for vi=1:length(cfg.varnames)
%         tmp.varn=cfg.varnames{vi};
%         tmp.var=reshape(comb_data.(tmp.varn), [xn yn 12 tn/12]);
%         comb_data_an_mov.(tmp.varn)(lm,:,:,:)=squeeze(mean(tmp.var(:,:,monrange,:),3));
%     end
% end
% lm=11;
% for vi=1:length(cfg.varnames)
%     tmp.varn=cfg.varnames{vi};
%     tmp.var=reshape(comb_data.(tmp.varn), [xn yn 12 tn/12]);
%     tmp.var1=squeeze(tmp.var(:,:,lm,:));
%     tmp.var2=squeeze(tmp.var(:,:,lm+1,:));
%     tmp.var3=squeeze(tmp.var(:,:,1,2:end));
%     tmp.var3(:,:,end+1)=NaN;
%     comb_data_an_mov.(tmp.varn)(lm,:,:,:)=(tmp.var1+tmp.var2+tmp.var3)./3;
% end
% lm=12;
% for vi=1:length(cfg.varnames)
%     tmp.varn=cfg.varnames{vi};
%     tmp.var=reshape(comb_data.(tmp.varn), [xn yn 12 tn/12]);
%     tmp.var1=squeeze(tmp.var(:,:,lm,:));
%     tmp.var2=squeeze(tmp.var(:,:,1,2:end));
%     tmp.var3=squeeze(tmp.var(:,:,2,2:end));
%     tmp.var2(:,:,end+1)=NaN;
%     tmp.var3(:,:,end+1)=NaN;
%     comb_data_an_mov.(tmp.varn)(lm,:,:,:)=(tmp.var1+tmp.var2+tmp.var3)./3;
% end


% seasonal mean

for vi=1:length(cfg.varnames)
    tmp.varn=cfg.varnames{vi};
    if strcmp(tmp.varn, 'vec') || strcmp(tmp.varn, 'wvec') || strcmp(tmp.varn, 'wsvec')
        disp('vec = no preprocessing here');
    else
        for lm=1:12
            monrange=lm:lm;
            tmp.var=reshape(comb_data.(tmp.varn), [xn yn 12 tn/12]);
            comb_data_an.(tmp.varn)(lm,:,:,:)=squeeze(mean(tmp.var(:,:,monrange,:),3));
        end
        
        for lm=1:10
            monrange=lm:lm+2;
            tmp.var=reshape(comb_data.(tmp.varn), [xn yn 12 tn/12]);
            comb_data_an_mov.(tmp.varn)(lm,:,:,:)=squeeze(mean(tmp.var(:,:,monrange,:),3));
        end
        
        lm=11;
        tmp.var=reshape(comb_data.(tmp.varn), [xn yn 12 tn/12]);
        tmp.var1=squeeze(tmp.var(:,:,lm,:));
        tmp.var2=squeeze(tmp.var(:,:,lm+1,:));
        tmp.var3=squeeze(tmp.var(:,:,1,2:end));
        tmp.var3(:,:,end+1)=NaN;
        comb_data_an_mov.(tmp.varn)(lm,:,:,:)=(tmp.var1+tmp.var2+tmp.var3)./3;
        
        lm=12;
        tmp.var=reshape(comb_data.(tmp.varn), [xn yn 12 tn/12]);
        tmp.var1=squeeze(tmp.var(:,:,lm,:));
        tmp.var2=squeeze(tmp.var(:,:,1,2:end));
        tmp.var3=squeeze(tmp.var(:,:,2,2:end));
        tmp.var2(:,:,end+1)=NaN;
        tmp.var3(:,:,end+1)=NaN;
        comb_data_an_mov.(tmp.varn)(lm,:,:,:)=(tmp.var1+tmp.var2+tmp.var3)./3;
    end
end



% part=[1800, 1750, 2050, 2030, 2000, 2400, 2200, 2010, 2500, 1730, 2150, 1930, 2700, 2030, 1510, ...
%     2230, 2400, 2500, 2600, 2700, 1980, 2700, 2150, 2430, 3150, ...
%     1930, 2400, 2250, 2400, 2380, 2370, 2000, 1800, 1500, 1900, 1950, 1880, 1490];

part=load('/Users/kimyy/Desktop/backup/Research/Postdoc/03_IBS/2024_TS_KS_connectivity/layer_particle_month_2.mat');

% part.years=1982:2020;
% part.months=1:12;

part_names= {'tcage', 'wcage', 'ecage', 'tcdepth', 'wcdepth', 'ecdepth', 'tccount', 'wccount', 'eccount'};

for ni=1:length(part_names)
    tmp.data=squeeze(sum(part.(part_names{ni}), 1, 'omitnan'));
    part2.(part_names{ni})=tmp.data(2:end-1,:); % subsample 1983~2019
    for lm=1:12
        monrange=lm:lm;
        part3.(part_names{ni})(:,lm)=squeeze(mean(part2.(part_names{ni})(:,monrange),2));
    end
    for lm=1:10
        monrange=lm:lm+2;
        part3_mov.(part_names{ni})(:,lm)=squeeze(mean(part2.(part_names{ni})(:,monrange),2));
    end
    lm=11;
    tmp.var1=squeeze(part2.(part_names{ni})(:,lm));
    tmp.var2=squeeze(part2.(part_names{ni})(:,lm+1));
    tmp.var3=squeeze(part2.(part_names{ni})(2:end,1));
    tmp.var3(end+1)=NaN;
    part3_mov.(part_names{ni})(:,lm)=(tmp.var1+tmp.var2+tmp.var3)./3;
    lm=12;
    tmp.var1=squeeze(part2.(part_names{ni})(:,lm));
    tmp.var2=squeeze(part2.(part_names{ni})(2:end,1));
    tmp.var3=squeeze(part2.(part_names{ni})(2:end,2));
    tmp.var2(end+1)=NaN;
    tmp.var3(end+1)=NaN;
    part3_mov.(part_names{ni})(:,lm)=(tmp.var1+tmp.var2+tmp.var3)./3;
end


regionname='stlee_drifter';
%% polygon set
%% stleedrifter

stlee_drifter_polygon = ...
    [120, 25;
    121, 25;
    124, 26;
    126, 26.5;
    128, 28;
    128, 30;
    130, 32;
    130, 35;
    120, 35;
    120, 32;
    122, 32;
    122, 29;
    120, 28;];  %% with onshore

switch(regionname)
    case('stlee_drifter') %% for debugging
        refpolygon=stlee_drifter_polygon;
    otherwise
        ('?');
end

%% masking
switch(regionname)
    case('NWP') %% North western Pacific
        mask_model(1:size(lon,1),1:size(lon,2))=1;
    otherwise
        mask_model = double(inpolygon(lon,lat,refpolygon(:,1),refpolygon(:,2)));
        mask_model(mask_model==0)=NaN;
end

stlee_orig_polygon = ...
    [121, 24;
    122, 24.5;
    124, 25;
    126, 26;
    129, 28;
    131, 32;
    122, 32;
    122, 30;
    120, 28;
    120, 27];

lonlat(1)=min(stlee_orig_polygon(:,1));
lonlat(2)=max(stlee_orig_polygon(:,1));
lonlat(3)=min(stlee_orig_polygon(:,2));
lonlat(4)=max(stlee_orig_polygon(:,2));

lonlat(1)=min(min(lon));
lonlat(2)=max(max(lon));
lonlat(3)=min(min(lat));
lonlat(4)=max(max(lat));

cfg.years=1983:2019;
[fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwr_20', tmp.dropboxpath);


%% depth masked
% tmp.varn='Vwind'

for lm=1:12
    for depthi=1:length(depths)
        dd=depths(depthi);
        str_d=num2str(dd);
        masks.(['d', str_d]);
    
        tmp.svd_modes=4;
        for vi=1:length(cfg.varnames)
            tmp.varn=cfg.varnames{vi};
%             'vec', 'wvec', 'wsvec'
            if strcmp(tmp.varn, 'vec') || strcmp(tmp.varn, 'wvec') || strcmp(tmp.varn, 'wsvec')
                switch tmp.varn
                    case 'vec'
                        part_u='ubar';
                        part_v='vbar';
                    case 'wvec'
                        part_u='Uwind';
                        part_v='Vwind';
                    case 'wsvec'
                        part_u='sustr';
                        part_v='svstr';
                end

                tmp.size_uv=size(comb_data_an.ubar);
                comb_data_an.(tmp.varn)(lm,1:tmp.size_uv(2),:,:)=squeeze(comb_data_an.(part_u)(lm,:,:,:)).*mask_model.*masks.(['d', str_d]);
                comb_data_an.(tmp.varn)(lm,tmp.size_uv(2)+1:tmp.size_uv(2)*2, :,:)=squeeze(comb_data_an.(part_v)(lm,:,:,:)).*mask_model.*masks.(['d', str_d]);
                comb_data_an_mov.(tmp.varn)(lm,1:tmp.size_uv(2),:,:)=squeeze(comb_data_an_mov.(part_u)(lm,:,:,:)).*mask_model.*masks.(['d', str_d]);
                comb_data_an_mov.(tmp.varn)(lm,tmp.size_uv(2)+1:tmp.size_uv(2)*2, :,:)=squeeze(comb_data_an_mov.(part_v)(lm,:,:,:)).*mask_model.*masks.(['d', str_d]);

                tmp.lat_vec(1:tmp.size_uv(2), :)=lat;
                tmp.lat_vec(tmp.size_uv(2)+1:tmp.size_uv(2)*2, :)=lat;
                tmp.lon_vec(1:tmp.size_uv(2), :)=lon;
                tmp.lon_vec(tmp.size_uv(2)+1:tmp.size_uv(2)*2, :)=lon;

                [EOF_result.(tmp.varn).lv(lm,:,:,:), ...
                EOF_result.(tmp.varn).pcs(lm,:,:), ...
                EOF_result.(tmp.varn).var_exp(lm,:)] = ...
                Func_0024_EOF_3d( squeeze(comb_data_an.(tmp.varn)(lm,:,:,:)), tmp.svd_modes, tmp.lat_vec);

                if lm<=10
                    [EOF_result_mov.(tmp.varn).lv(lm,:,:,:), ...
                    EOF_result_mov.(tmp.varn).pcs(lm,:,:), ...
                    EOF_result_mov.(tmp.varn).var_exp(lm,:)] = ...
                    Func_0024_EOF_3d( squeeze(comb_data_an_mov.(tmp.varn)(lm,:,:,:)), tmp.svd_modes, tmp.lat_vec);
                else
                    [EOF_result_mov.(tmp.varn).lv(lm,:,:,:), ...
                    EOF_result_mov.(tmp.varn).pcs(lm,1:end-1,:), ...
                    EOF_result_mov.(tmp.varn).var_exp(lm,:)] = ...
                    Func_0024_EOF_3d( squeeze(comb_data_an_mov.(tmp.varn)(lm,:,:,1:end-1)), tmp.svd_modes, tmp.lat_vec);
                    EOF_result_mov.(tmp.varn).pcs(lm,end,:)=NaN;
                end
            else
                    [EOF_result.(tmp.varn).lv(lm,:,:,:), ...
                EOF_result.(tmp.varn).pcs(lm,:,:), ...
                EOF_result.(tmp.varn).var_exp(lm,:)] = ...
                Func_0024_EOF_3d( squeeze(comb_data_an.(tmp.varn)(lm,:,:,:)) ...
                    .*mask_model.*masks.(['d', str_d]), tmp.svd_modes, lat);

                if lm<=10
                [EOF_result_mov.(tmp.varn).lv(lm,:,:,:), ...
                EOF_result_mov.(tmp.varn).pcs(lm,:,:), ...
                EOF_result_mov.(tmp.varn).var_exp(lm,:)] = ...
                Func_0024_EOF_3d( squeeze(comb_data_an_mov.(tmp.varn)(lm,:,:,:)) ...
                    .*mask_model.*masks.(['d', str_d]), tmp.svd_modes, lat);
                else
                    [EOF_result_mov.(tmp.varn).lv(lm,:,:,:), ...
                    EOF_result_mov.(tmp.varn).pcs(lm,1:end-1,:), ...
                    EOF_result_mov.(tmp.varn).var_exp(lm,:)] = ...
                    Func_0024_EOF_3d( squeeze(comb_data_an_mov.(tmp.varn)(lm,:,:,1:end-1)) ...
                        .*mask_model.*masks.(['d', str_d]), tmp.svd_modes, lat);
                    EOF_result_mov.(tmp.varn).pcs(lm,end,:)=NaN;
                end
            end
        end
    
        %% sign check
        for vi=1:length(cfg.varnames)
%         for vi=1:1
            tmp.varn=cfg.varnames{vi};
            for amodi=1:tmp.svd_modes
%                 if mean(EOF_result.(tmp.varn).pcs(lm,1:19,amodi))<0
                if mean(mean(squeeze(EOF_result.(tmp.varn).lv(lm,:,:,amodi)),'omitnan'),'omitnan')<0
                    EOF_result.(tmp.varn).lv(lm,:,:,amodi)=-EOF_result.(tmp.varn).lv(lm,:,:,amodi);
                    EOF_result.(tmp.varn).pcs(lm,:,amodi)=-EOF_result.(tmp.varn).pcs(lm,:,amodi);
                end
                if mean(mean(squeeze(EOF_result_mov.(tmp.varn).lv(lm,:,:,amodi)),'omitnan'),'omitnan')<0
                    EOF_result_mov.(tmp.varn).lv(lm,:,:,amodi)=-EOF_result_mov.(tmp.varn).lv(lm,:,:,amodi);
                    EOF_result_mov.(tmp.varn).pcs(lm,:,amodi)=-EOF_result_mov.(tmp.varn).pcs(lm,:,amodi);
                end
                
                if amodi==1
                    tmp.data_raw=squeeze(comb_data_an.(tmp.varn)(lm,:,:,:));
                    tmp.data_raw(abs(tmp.data_raw)>10e6)=NaN;
                    tmp.data_mov=squeeze(comb_data_an_mov.(tmp.varn)(lm,:,:,:));
                    tmp.data_mov(abs(tmp.data_mov)>10e6)=NaN;
                    if strcmp(tmp.varn, 'vec') || strcmp(tmp.varn, 'wvec') || strcmp(tmp.varn, 'wsvec')
                        tmp.data_m=Func_0011_get_area_weighted_mean( ...
                            tmp.data_raw, tmp.lon_vec,tmp.lat_vec);
                        tmp.data_mov_m=Func_0011_get_area_weighted_mean( ...
                            tmp.data_mov, tmp.lon_vec,tmp.lat_vec);
                    else
                        tmp.data_m=Func_0011_get_area_weighted_mean( ...
                            tmp.data_raw, lon, lat);
                        tmp.data_mov_m=Func_0011_get_area_weighted_mean( ...
                            tmp.data_mov, lon, lat);
                    end
                    
                end
                for lt=0:12
                    for ni=1:length(part_names)
                        if (lm+lt)<=12
                            tmp.part=part3.(part_names{ni})(:,lm+lt);
                        elseif (lm+lt)>12
                            tmp.part=part3.(part_names{ni})(2:end,lm+lt-12);
                            tmp.part=[tmp.part; NaN];
                        end
                        tmp.corr=corrcoef(EOF_result.(tmp.varn).pcs(lm,:,amodi), tmp.part, 'Rows', 'complete');
                        corr_EOF(lt+1,vi,ni,lm,amodi)=tmp.corr(1,2);
                        if amodi==1
                            tmp.corr=corrcoef(tmp.data_m, tmp.part, 'Rows', 'complete');
                            corr_raw(lt+1,vi,ni,lm)=tmp.corr(1,2);
                        end

                        if (lm+lt)<=12
                            tmp.part_mov=part3_mov.(part_names{ni})(:,lm+lt);
                        elseif (lm+lt)>12
                            tmp.part_mov=part3_mov.(part_names{ni})(2:end,lm+lt-12);
                            tmp.part_mov=[tmp.part_mov; NaN];
                        end

                        tmp.corr=corrcoef(EOF_result_mov.(tmp.varn).pcs(lm,:,amodi), tmp.part_mov, 'Rows', 'complete');
                        corr_EOF_mov(lt+1,vi,ni,lm,amodi)=tmp.corr(1,2);
                        if amodi==1
                            tmp.corr=corrcoef(tmp.data_mov_m, tmp.part, 'Rows', 'complete');
                            corr_raw_mov(lt+1,vi,ni,lm)=tmp.corr(1,2);
                        end
                    end
                    disp([tmp.varn, ', ', part_names{ni}, ', ', num2str(lt), 'lt'])
                end
            end
        end
    end
end
for vi=1:length(cfg.varnames)
    tmp.varn=cfg.varnames{vi};
    if strcmp(tmp.varn, 'vec') || strcmp(tmp.varn, 'wvec') || strcmp(tmp.varn, 'wsvec')
        EOF_result.([tmp.varn, '_u']).lv= EOF_result.(tmp.varn).lv(:,1:tmp.size_uv(2),:,:);
        EOF_result.([tmp.varn, '_v']).lv= EOF_result.(tmp.varn).lv(:,tmp.size_uv(2)+1:tmp.size_uv(2)*2,:,:);
        EOF_result.([tmp.varn, '_u']).pcs=EOF_result.(tmp.varn).pcs;
        EOF_result.([tmp.varn, '_v']).pcs=EOF_result.(tmp.varn).pcs;
        EOF_result.([tmp.varn, '_u']).var_exp=EOF_result.(tmp.varn).var_exp;
        EOF_result.([tmp.varn, '_v']).var_exp=EOF_result.(tmp.varn).var_exp;

        EOF_result_mov.([tmp.varn, '_u']).lv= EOF_result.(tmp.varn).lv(:,1:tmp.size_uv(2),:,:);
        EOF_result_mov.([tmp.varn, '_v']).lv= EOF_result.(tmp.varn).lv(:,tmp.size_uv(2)+1:tmp.size_uv(2)*2,:,:);
        EOF_result_mov.([tmp.varn, '_u']).pcs=EOF_result.(tmp.varn).pcs;
        EOF_result_mov.([tmp.varn, '_v']).pcs=EOF_result.(tmp.varn).pcs;
        EOF_result_mov.([tmp.varn, '_u']).var_exp=EOF_result.(tmp.varn).var_exp;
        EOF_result_mov.([tmp.varn, '_v']).var_exp=EOF_result.(tmp.varn).var_exp;    
    end
end


%% 95% significance low & high range
    conf = 0.90;
%     n=cfg.t_ini_max-cfg.t_ini_min+1;
    n=38;
    alpha = 1 - conf;
    pLo = alpha/2;
    pUp = 1 - alpha/2;
    crit = tinv([pLo pUp], n-1);
    % xbar = mean(r); % = 0
    xbar = 0; % = 0
    r_crit=sqrt((crit.^2)./(n-2+(crit).^2));



% % %% figure raw
% % for vi=1:length(cfg.varnames)
% %     for ni=7:9
% %         tmp.data=squeeze(corr_raw(1:13,vi,ni,:));
% % %         tmp.data2=[tmp.data, NaN(12,1)];
% % %         tmp.data2(13,:)=NaN;
% %         tmp.data(end+1,:)=NaN;
% %         tmp.data(:,end+1)=NaN;
% %         
% %         tmp.lt=0:13;
% %         tmp.lt=repmat(tmp.lt', [1 13]);
% % 
% %         tmp.lm=1:13;
% %         tmp.lm=repmat(tmp.lm, [14, 1]);
% %         
% % %         tmp.lm_str={'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'NaN'};
% %         tmp.lm_str={'Feb', 'Apr', 'Jun', 'Aug', 'Oct', 'Dec'};
% % 
% %         fig_h=figure('visible', 'off');
% %         pcolor(tmp.lt', tmp.lm', tmp.data', 'parent', fig_h); shading flat; colorbar;
% %         hold on;
% %         caxis([-1 1]); colormap(fig_cfg.c_map);
% %         tmp.data3=tmp.data;
% %         tmp.data3(abs(tmp.data)>r_crit(1))=NaN;
% %         if sum(isfinite(tmp.data3(:)))>0
% %             pp2=pcolor(tmp.lt', tmp.lm', tmp.data3', 'parent', fig_h);
% %             set(pp2,'linestyle','none','Tag','HatchingRegion');
% %             hp = findobj(pp2,'Tag','HatchingRegion');
% %             hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',300,'HatchColor','k','HatchLineWidth',0.5);
% %         end
% %         ylabel('month');
% %         xlabel('lead time');
% %         yticklabels(tmp.lm_str);
% % 
% %         set(gca, 'fontsize', 20);
% %         title([cfg.varnames{vi}, ' vs ', part_names{ni}]);
% % 
% %         dir.figtgdir=[dir.figdir, '/', 'corrs_raw'];
% %             mkdir(dir.figtgdir);
% %         cfg.figname=[dir.figtgdir, '/', 'corrs_raw_', ...
% %             str_d, 'm', '_', cfg.varnames{vi}, '_', part_names{ni} '.tif'];
% %             print(fig_h, cfg.figname, '-dpng');
% %         axis tight;
% %         RemoveWhiteSpace([], 'file', cfg.figname);
% %         hold off;
% %         close all;
% %     end
% % end
% % 
% % %% figure raw (movm, 3)
% % for vi=1:length(cfg.varnames)
% %     for ni=7:9
% %         tmp.data=squeeze(corr_raw_mov(1:13,vi,ni,:));
% % %         tmp.data2=[tmp.data, NaN(12,1)];
% % %         tmp.data2(13,:)=NaN;
% %         tmp.data(end+1,:)=NaN;
% %         tmp.data(:,end+1)=NaN;
% %         
% %         tmp.lt=0:13;
% %         tmp.lt=repmat(tmp.lt', [1 13]);
% % 
% %         tmp.lm=1:13;
% %         tmp.lm=repmat(tmp.lm, [14, 1]);
% %         
% % %         tmp.lm_str={'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'NaN'};
% %         tmp.lm_str={'Feb-Apr', 'Apr-Jun', 'Jun-Aug', 'Aug-Oct', 'Oct-Dec', 'Dec-Feb'};
% % 
% % 
% %         fig_h=figure('visible', 'off');
% %         pcolor(tmp.lt', tmp.lm', tmp.data', 'parent', fig_h); shading flat; colorbar;
% %         hold on;
% %         caxis([-1 1]); colormap(fig_cfg.c_map);
% %         tmp.data3=tmp.data;
% %         tmp.data3(abs(tmp.data)>r_crit(1))=NaN;
% %         if sum(isfinite(tmp.data3(:)))>0
% %             pp2=pcolor(tmp.lt', tmp.lm', tmp.data3', 'parent', fig_h);
% %             set(pp2,'linestyle','none','Tag','HatchingRegion');
% %             hp = findobj(pp2,'Tag','HatchingRegion');
% %             hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',300,'HatchColor','k','HatchLineWidth',0.5);
% %         end
% %         ylabel('season');
% %         xlabel('lead time');
% %         yticklabels(tmp.lm_str);
% %         set(gca, 'fontsize', 20);
% %         title([cfg.varnames{vi}, ' vs ', part_names{ni}]);
% % 
% %         dir.figtgdir=[dir.figdir, '/', 'corrs_raw_mov'];
% %             mkdir(dir.figtgdir);
% %         cfg.figname=[dir.figtgdir, '/', 'corrs_raw_mov_', ...
% %             str_d, 'm', '_', cfg.varnames{vi}, '_', part_names{ni} '.tif'];
% %             print(fig_h, cfg.figname, '-dpng');
% %         axis tight;
% %         RemoveWhiteSpace([], 'file', cfg.figname);
% %         hold off;
% %         close all;
% %     end
% % end




% Spring: 3~5;
% Summer: 6~8;
% Autumn: 9~11;
% Winter: 12~2;

%% EOF LV & PCT (1-4 modes) - seasonal
for pmon= [7, 6, 8, 9, 10, 11]
    for si=[3, 6, 9, 12]
    % for si=[9, 12]
        for vi=1:length(cfg.varnames)
            tmp.varn=cfg.varnames{vi};
%             pmon=9;
            if si<=pmon
                part_tc = part3.tccount(:,pmon); 
                part_wc = part3.wccount(:,pmon); 
                part_ec = part3.eccount(:,pmon); 
            else
                part_tc = part3.tccount(2:end,pmon); 
                part_tc = [part_tc; NaN];
                part_wc = part3.wccount(2:end,pmon); 
                part_wc = [part_wc; NaN];
                part_ec = part3.eccount(2:end,pmon); 
                part_ec = [part_ec; NaN];
            end
        
            if strcmp(tmp.varn, 'vec') || strcmp(tmp.varn, 'wvec') || strcmp(tmp.varn, 'wsvec')
                %% vector
                tmp.varn=cfg.varnames{vi};
                
                tmp.intv=5;
                tmp.amp_size=15;

                close all;
                fig_cfg.fig_size=[0 0 14 16];
                fig_h = figure('name', 'EOF','PaperUnits','inches', ...
                        'PaperPosition',(fig_cfg.fig_size),'position', fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
                amod=1;
                subplot(4,2,1);
                lv_u=squeeze(EOF_result_mov.([tmp.varn, '_u']).lv(si,:,:,amod));
                lv_v=squeeze(EOF_result_mov.([tmp.varn, '_v']).lv(si,:,:,amod));

                m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                m_quiver(lon(1:tmp.intv:end, 1:tmp.intv:end)', ...
                    lat(1:tmp.intv:end, 1:tmp.intv:end)', ...
                    lv_u(1:tmp.intv:end, 1:tmp.intv:end)' * tmp.amp_size, ...
                    lv_v(1:tmp.intv:end, 1:tmp.intv:end)' * tmp.amp_size, ...
                    'AutoScale','off','LineWidth', 1); 
                m_grid;  
                m_gshhs_i('color',[1 1 1]);
                m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
                title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
                
                subplot(4,2,2);
                plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod));
                tmp.corr1=corrcoef(part_tc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr2=corrcoef(part_wc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr3=corrcoef(part_ec, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                title([tmp.varn, ' pcs, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%, ', ...
                    'r(tc):', num2str(round(tmp.corr1(1,2),2)), ...
                    ', r(wc):', num2str(round(tmp.corr2(1,2),2)), ...
                    ', r(ec):', num2str(round(tmp.corr3(1,2),2))]);
                grid on;
                
                
                amod=2;
                subplot(4,2,3);
                lv_u=squeeze(EOF_result_mov.([tmp.varn, '_u']).lv(si,:,:,amod));
                lv_v=squeeze(EOF_result_mov.([tmp.varn, '_v']).lv(si,:,:,amod));

                m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                m_quiver(lon(1:tmp.intv:end, 1:tmp.intv:end)', ...
                    lat(1:tmp.intv:end, 1:tmp.intv:end)', ...
                    lv_u(1:tmp.intv:end, 1:tmp.intv:end)' * tmp.amp_size, ...
                    lv_v(1:tmp.intv:end, 1:tmp.intv:end)' * tmp.amp_size, ...
                    'AutoScale','off','LineWidth', 1); 
                m_grid;    
                m_gshhs_i('color',[1 1 1]);
                m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
                title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
                
                subplot(4,2,4);
                plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod));
                tmp.corr1=corrcoef(part_tc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr2=corrcoef(part_wc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr3=corrcoef(part_ec, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                title([tmp.varn, ' pcs, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%, ', ...
                    'r(tc):', num2str(round(tmp.corr1(1,2),2)), ...
                    ', r(wc):', num2str(round(tmp.corr2(1,2),2)), ...
                    ', r(ec):', num2str(round(tmp.corr3(1,2),2))]);
                grid on;
                
                amod=3;
                subplot(4,2,5)
                lv_u=squeeze(EOF_result_mov.([tmp.varn, '_u']).lv(si,:,:,amod));
                lv_v=squeeze(EOF_result_mov.([tmp.varn, '_v']).lv(si,:,:,amod));

                m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                m_quiver(lon(1:tmp.intv:end, 1:tmp.intv:end)', ...
                    lat(1:tmp.intv:end, 1:tmp.intv:end)', ...
                    lv_u(1:tmp.intv:end, 1:tmp.intv:end)' * tmp.amp_size, ...
                    lv_v(1:tmp.intv:end, 1:tmp.intv:end)' * tmp.amp_size, ...
                    'AutoScale','off','LineWidth', 1); 
                m_grid;  
                m_gshhs_i('color',[1 1 1]);
                m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
                title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
                
                subplot(4,2,6);
                plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod));
                tmp.corr1=corrcoef(part_tc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr2=corrcoef(part_wc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr3=corrcoef(part_ec, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                title([tmp.varn, ' pcs, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%, ', ...
                    'r(tc):', num2str(round(tmp.corr1(1,2),2)), ...
                    ', r(wc):', num2str(round(tmp.corr2(1,2),2)), ...
                    ', r(ec):', num2str(round(tmp.corr3(1,2),2))]);
                grid on;
                
                amod=4;
                subplot(4,2,7)
                lv_u=squeeze(EOF_result_mov.([tmp.varn, '_u']).lv(si,:,:,amod));
                lv_v=squeeze(EOF_result_mov.([tmp.varn, '_v']).lv(si,:,:,amod));

                m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                m_quiver(lon(1:tmp.intv:end, 1:tmp.intv:end)', ...
                    lat(1:tmp.intv:end, 1:tmp.intv:end)', ...
                    lv_u(1:tmp.intv:end, 1:tmp.intv:end)' * tmp.amp_size, ...
                    lv_v(1:tmp.intv:end, 1:tmp.intv:end)' * tmp.amp_size, ...
                    'AutoScale','off','LineWidth', 1); 
                m_grid;  
                m_gshhs_i('color',[1 1 1]);
                m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
                title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
                
                subplot(4,2,8);
                plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod));
                tmp.corr1=corrcoef(part_tc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr2=corrcoef(part_wc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr3=corrcoef(part_ec, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                title([tmp.varn, ' pcs, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%, ', ...
                    'r(tc):', num2str(round(tmp.corr1(1,2),2)), ...
                    ', r(wc):', num2str(round(tmp.corr2(1,2),2)), ...
                    ', r(ec):', num2str(round(tmp.corr3(1,2),2))]);
                grid on;
                
                
                dir.figtgdir=[dir.figdir, '/', 'EOF_mov', '/', 'par', num2str(pmon)];
                mkdir(dir.figtgdir);
                
                cfg.figname=[dir.figtgdir, '/', 'EOF_mov_',tmp.varn,'_season_',num2str(si), '_', str_d, '_p_', num2str(pmon) '.tif'];
                print(fig_h, cfg.figname, '-dpng');
                RemoveWhiteSpace([], 'file', cfg.figname);

                %% basic figure_v
                tmp.varn=[cfg.varnames{vi}, '_v'];

                close all;
                fig_cfg.fig_size=[0 0 14 16];
                fig_h = figure('name', 'EOF','PaperUnits','inches', ...
                        'PaperPosition',(fig_cfg.fig_size),'position', fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
                amod=1;
                subplot(4,2,1);
                lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
                m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
                % m_grid('fontsize', m_grid_fontsize-5, 'tickdir', m_grid_tickdir_type, 'xticklabels', [], 'xtick',(120:20:160), 'parent', ax1_1);  
                m_grid;  
                m_gshhs_i('color',[1 1 1]);
                m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
                title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
                mcax=max(abs(lv_l(:)));
                caxis([-mcax mcax]);
                colormap(fig_cfg.c_map);
                
                subplot(4,2,2);
                plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod));
                tmp.corr1=corrcoef(part_tc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr2=corrcoef(part_wc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr3=corrcoef(part_ec, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                title([tmp.varn, ' pcs, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%, ', ...
                    'r(tc):', num2str(round(tmp.corr1(1,2),2)), ...
                    ', r(wc):', num2str(round(tmp.corr2(1,2),2)), ...
                    ', r(ec):', num2str(round(tmp.corr3(1,2),2))]);
                grid on;
                
                
                amod=2;
                subplot(4,2,3);
                lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
                m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
                m_grid;  
                m_gshhs_i('color',[1 1 1]);
                m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
                title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
                mcax=max(abs(lv_l(:)));
                caxis([-mcax mcax]);
                colormap(fig_cfg.c_map);
                
                subplot(4,2,4);
                plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod));
                tmp.corr1=corrcoef(part_tc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr2=corrcoef(part_wc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr3=corrcoef(part_ec, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                title([tmp.varn, ' pcs, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%, ', ...
                    'r(tc):', num2str(round(tmp.corr1(1,2),2)), ...
                    ', r(wc):', num2str(round(tmp.corr2(1,2),2)), ...
                    ', r(ec):', num2str(round(tmp.corr3(1,2),2))]);
                grid on;
                
                amod=3;
                subplot(4,2,5)
                lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
                m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
                m_grid;  
                m_gshhs_i('color',[1 1 1]);
                m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
                title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
                mcax=max(abs(lv_l(:)));
                caxis([-mcax mcax]);
                colormap(fig_cfg.c_map);
                
                subplot(4,2,6);
                plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod));
                tmp.corr1=corrcoef(part_tc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr2=corrcoef(part_wc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr3=corrcoef(part_ec, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                title([tmp.varn, ' pcs, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%, ', ...
                    'r(tc):', num2str(round(tmp.corr1(1,2),2)), ...
                    ', r(wc):', num2str(round(tmp.corr2(1,2),2)), ...
                    ', r(ec):', num2str(round(tmp.corr3(1,2),2))]);
                grid on;
                
                amod=4;
                subplot(4,2,7)
                lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
                m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
                m_grid;  
                m_gshhs_i('color',[1 1 1]);
                m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
                title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
                mcax=max(abs(lv_l(:)));
                caxis([-mcax mcax]);
                colormap(fig_cfg.c_map);
                
                subplot(4,2,8);
                plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod));
                tmp.corr1=corrcoef(part_tc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr2=corrcoef(part_wc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr3=corrcoef(part_ec, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                title([tmp.varn, ' pcs, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%, ', ...
                    'r(tc):', num2str(round(tmp.corr1(1,2),2)), ...
                    ', r(wc):', num2str(round(tmp.corr2(1,2),2)), ...
                    ', r(ec):', num2str(round(tmp.corr3(1,2),2))]);
                grid on;
                
                
                dir.figtgdir=[dir.figdir, '/', 'EOF_mov', '/', 'par', num2str(pmon)];
                mkdir(dir.figtgdir);
                
                cfg.figname=[dir.figtgdir, '/', 'EOF_mov_',tmp.varn,'_season_',num2str(si), '_', str_d, '_p_', num2str(pmon) '.tif'];
                print(fig_h, cfg.figname, '-dpng');
                RemoveWhiteSpace([], 'file', cfg.figname);

                %% basic figure_u
                tmp.varn=[cfg.varnames{vi}, '_u'];

                close all;
                fig_cfg.fig_size=[0 0 14 16];
                fig_h = figure('name', 'EOF','PaperUnits','inches', ...
                        'PaperPosition',(fig_cfg.fig_size),'position', fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
                amod=1;
                subplot(4,2,1);
                lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
                m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
                % m_grid('fontsize', m_grid_fontsize-5, 'tickdir', m_grid_tickdir_type, 'xticklabels', [], 'xtick',(120:20:160), 'parent', ax1_1);  
                m_grid;  
                m_gshhs_i('color',[1 1 1]);
                m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
                title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
                mcax=max(abs(lv_l(:)));
                caxis([-mcax mcax]);
                colormap(fig_cfg.c_map);
                
                subplot(4,2,2);
                plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod));
                tmp.corr1=corrcoef(part_tc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr2=corrcoef(part_wc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr3=corrcoef(part_ec, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                title([tmp.varn, ' pcs, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%, ', ...
                    'r(tc):', num2str(round(tmp.corr1(1,2),2)), ...
                    ', r(wc):', num2str(round(tmp.corr2(1,2),2)), ...
                    ', r(ec):', num2str(round(tmp.corr3(1,2),2))]);
                grid on;
                
                
                amod=2;
                subplot(4,2,3);
                lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
                m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
                m_grid;  
                m_gshhs_i('color',[1 1 1]);
                m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
                title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
                mcax=max(abs(lv_l(:)));
                caxis([-mcax mcax]);
                colormap(fig_cfg.c_map);
                
                subplot(4,2,4);
                plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod));
                tmp.corr1=corrcoef(part_tc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr2=corrcoef(part_wc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr3=corrcoef(part_ec, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                title([tmp.varn, ' pcs, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%, ', ...
                    'r(tc):', num2str(round(tmp.corr1(1,2),2)), ...
                    ', r(wc):', num2str(round(tmp.corr2(1,2),2)), ...
                    ', r(ec):', num2str(round(tmp.corr3(1,2),2))]);
                grid on;
                
                amod=3;
                subplot(4,2,5)
                lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
                m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
                m_grid;  
                m_gshhs_i('color',[1 1 1]);
                m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
                title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
                mcax=max(abs(lv_l(:)));
                caxis([-mcax mcax]);
                colormap(fig_cfg.c_map);
                
                subplot(4,2,6);
                plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod));
                tmp.corr1=corrcoef(part_tc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr2=corrcoef(part_wc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr3=corrcoef(part_ec, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                title([tmp.varn, ' pcs, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%, ', ...
                    'r(tc):', num2str(round(tmp.corr1(1,2),2)), ...
                    ', r(wc):', num2str(round(tmp.corr2(1,2),2)), ...
                    ', r(ec):', num2str(round(tmp.corr3(1,2),2))]);
                grid on;
                
                amod=4;
                subplot(4,2,7)
                lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
                m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
                m_grid;  
                m_gshhs_i('color',[1 1 1]);
                m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
                title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
                mcax=max(abs(lv_l(:)));
                caxis([-mcax mcax]);
                colormap(fig_cfg.c_map);
                
                subplot(4,2,8);
                plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod));
                tmp.corr1=corrcoef(part_tc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr2=corrcoef(part_wc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr3=corrcoef(part_ec, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                title([tmp.varn, ' pcs, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%, ', ...
                    'r(tc):', num2str(round(tmp.corr1(1,2),2)), ...
                    ', r(wc):', num2str(round(tmp.corr2(1,2),2)), ...
                    ', r(ec):', num2str(round(tmp.corr3(1,2),2))]);
                grid on;
                
                
                dir.figtgdir=[dir.figdir, '/', 'EOF_mov', '/', 'par', num2str(pmon)];
                mkdir(dir.figtgdir);
                
                cfg.figname=[dir.figtgdir, '/', 'EOF_mov_',tmp.varn,'_season_',num2str(si), '_', str_d, '_p_', num2str(pmon) '.tif'];
                print(fig_h, cfg.figname, '-dpng');
                RemoveWhiteSpace([], 'file', cfg.figname);

                %% basic figure_v
                tmp.varn=[cfg.varnames{vi}, '_v'];

                close all;
                fig_cfg.fig_size=[0 0 14 16];
                fig_h = figure('name', 'EOF','PaperUnits','inches', ...
                        'PaperPosition',(fig_cfg.fig_size),'position', fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
                amod=1;
                subplot(4,2,1);
                lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
                m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
                % m_grid('fontsize', m_grid_fontsize-5, 'tickdir', m_grid_tickdir_type, 'xticklabels', [], 'xtick',(120:20:160), 'parent', ax1_1);  
                m_grid;  
                m_gshhs_i('color',[1 1 1]);
                m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
                title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
                mcax=max(abs(lv_l(:)));
                caxis([-mcax mcax]);
                colormap(fig_cfg.c_map);
                
                subplot(4,2,2);
                plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod));
                tmp.corr1=corrcoef(part_tc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr2=corrcoef(part_wc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr3=corrcoef(part_ec, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                title([tmp.varn, ' pcs, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%, ', ...
                    'r(tc):', num2str(round(tmp.corr1(1,2),2)), ...
                    ', r(wc):', num2str(round(tmp.corr2(1,2),2)), ...
                    ', r(ec):', num2str(round(tmp.corr3(1,2),2))]);
                grid on;
                
                
                amod=2;
                subplot(4,2,3);
                lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
                m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
                m_grid;  
                m_gshhs_i('color',[1 1 1]);
                m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
                title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
                mcax=max(abs(lv_l(:)));
                caxis([-mcax mcax]);
                colormap(fig_cfg.c_map);
                
                subplot(4,2,4);
                plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod));
                tmp.corr1=corrcoef(part_tc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr2=corrcoef(part_wc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr3=corrcoef(part_ec, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                title([tmp.varn, ' pcs, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%, ', ...
                    'r(tc):', num2str(round(tmp.corr1(1,2),2)), ...
                    ', r(wc):', num2str(round(tmp.corr2(1,2),2)), ...
                    ', r(ec):', num2str(round(tmp.corr3(1,2),2))]);
                grid on;
                
                amod=3;
                subplot(4,2,5)
                lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
                m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
                m_grid;  
                m_gshhs_i('color',[1 1 1]);
                m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
                title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
                mcax=max(abs(lv_l(:)));
                caxis([-mcax mcax]);
                colormap(fig_cfg.c_map);
                
                subplot(4,2,6);
                plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod));
                tmp.corr1=corrcoef(part_tc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr2=corrcoef(part_wc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr3=corrcoef(part_ec, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                title([tmp.varn, ' pcs, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%, ', ...
                    'r(tc):', num2str(round(tmp.corr1(1,2),2)), ...
                    ', r(wc):', num2str(round(tmp.corr2(1,2),2)), ...
                    ', r(ec):', num2str(round(tmp.corr3(1,2),2))]);
                grid on;
                
                amod=4;
                subplot(4,2,7)
                lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
                m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
                m_grid;  
                m_gshhs_i('color',[1 1 1]);
                m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
                title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
                mcax=max(abs(lv_l(:)));
                caxis([-mcax mcax]);
                colormap(fig_cfg.c_map);
                
                subplot(4,2,8);
                plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod));
                tmp.corr1=corrcoef(part_tc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr2=corrcoef(part_wc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr3=corrcoef(part_ec, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                title([tmp.varn, ' pcs, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%, ', ...
                    'r(tc):', num2str(round(tmp.corr1(1,2),2)), ...
                    ', r(wc):', num2str(round(tmp.corr2(1,2),2)), ...
                    ', r(ec):', num2str(round(tmp.corr3(1,2),2))]);
                grid on;
                
                
                dir.figtgdir=[dir.figdir, '/', 'EOF_mov', '/', 'par', num2str(pmon)];
                mkdir(dir.figtgdir);
                
                cfg.figname=[dir.figtgdir, '/', 'EOF_mov_',tmp.varn,'_season_',num2str(si), '_', str_d, '_p_', num2str(pmon) '.tif'];
                print(fig_h, cfg.figname, '-dpng');
                RemoveWhiteSpace([], 'file', cfg.figname);
            else
                %% basic figure
                close all;
                fig_cfg.fig_size=[0 0 14 16];
                fig_h = figure('name', 'EOF','PaperUnits','inches', ...
                        'PaperPosition',(fig_cfg.fig_size),'position', fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
                amod=1;
                subplot(4,2,1);
                lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
                m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
                % m_grid('fontsize', m_grid_fontsize-5, 'tickdir', m_grid_tickdir_type, 'xticklabels', [], 'xtick',(120:20:160), 'parent', ax1_1);  
                m_grid;  
                m_gshhs_i('color',[1 1 1]);
                m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
                title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
                mcax=max(abs(lv_l(:)));
                caxis([-mcax mcax]);
                colormap(fig_cfg.c_map);
                
                subplot(4,2,2);
                plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod));
                tmp.corr1=corrcoef(part_tc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr2=corrcoef(part_wc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr3=corrcoef(part_ec, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                title([tmp.varn, ' pcs, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%, ', ...
                    'r(tc):', num2str(round(tmp.corr1(1,2),2)), ...
                    ', r(wc):', num2str(round(tmp.corr2(1,2),2)), ...
                    ', r(ec):', num2str(round(tmp.corr3(1,2),2))]);
                grid on;
                
                
                amod=2;
                subplot(4,2,3);
                lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
                m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
                m_grid;  
                m_gshhs_i('color',[1 1 1]);
                m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
                title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
                mcax=max(abs(lv_l(:)));
                caxis([-mcax mcax]);
                colormap(fig_cfg.c_map);
                
                subplot(4,2,4);
                plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod));
                tmp.corr1=corrcoef(part_tc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr2=corrcoef(part_wc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr3=corrcoef(part_ec, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                title([tmp.varn, ' pcs, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%, ', ...
                    'r(tc):', num2str(round(tmp.corr1(1,2),2)), ...
                    ', r(wc):', num2str(round(tmp.corr2(1,2),2)), ...
                    ', r(ec):', num2str(round(tmp.corr3(1,2),2))]);
                grid on;
                
                amod=3;
                subplot(4,2,5)
                lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
                m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
                m_grid;  
                m_gshhs_i('color',[1 1 1]);
                m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
                title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
                mcax=max(abs(lv_l(:)));
                caxis([-mcax mcax]);
                colormap(fig_cfg.c_map);
                
                subplot(4,2,6);
                plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod));
                tmp.corr1=corrcoef(part_tc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr2=corrcoef(part_wc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr3=corrcoef(part_ec, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                title([tmp.varn, ' pcs, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%, ', ...
                    'r(tc):', num2str(round(tmp.corr1(1,2),2)), ...
                    ', r(wc):', num2str(round(tmp.corr2(1,2),2)), ...
                    ', r(ec):', num2str(round(tmp.corr3(1,2),2))]);
                grid on;
                
                amod=4;
                subplot(4,2,7)
                lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
                m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
                m_grid;  
                m_gshhs_i('color',[1 1 1]);
                m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
                title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
                mcax=max(abs(lv_l(:)));
                caxis([-mcax mcax]);
                colormap(fig_cfg.c_map);
                
                subplot(4,2,8);
                plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod));
                tmp.corr1=corrcoef(part_tc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr2=corrcoef(part_wc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                tmp.corr3=corrcoef(part_ec, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
                title([tmp.varn, ' pcs, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%, ', ...
                    'r(tc):', num2str(round(tmp.corr1(1,2),2)), ...
                    ', r(wc):', num2str(round(tmp.corr2(1,2),2)), ...
                    ', r(ec):', num2str(round(tmp.corr3(1,2),2))]);
                grid on;
                
                
                dir.figtgdir=[dir.figdir, '/', 'EOF_mov', '/', 'par', num2str(pmon)];
                mkdir(dir.figtgdir);
                
                cfg.figname=[dir.figtgdir, '/', 'EOF_mov_',tmp.varn,'_season_',num2str(si), '_', str_d, '_p_', num2str(pmon) '.tif'];
                print(fig_h, cfg.figname, '-dpng');
                RemoveWhiteSpace([], 'file', cfg.figname);
            end
        end
    end
end


%% timeseries
for vi=1:length(cfg.varnames)
    close all;
    fig_cfg.fig_size=[0 0 10 6];
    fig_h = figure('name', 'EOF','PaperUnits','inches', ...
            'PaperPosition',(fig_cfg.fig_size),'position', fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','on');

    tmp.varn=cfg.varnames{vi};
    loninfo=129;
    latinfo=30;
    [indw, inde, inds, indn]=Func_0012_findind_Y(1/10,[loninfo latinfo],lon,lat, 1);
    plot(cfg.years,squeeze(comb_data_an.(tmp.varn)(indw,inds,:)))
    tmp.corr=corrcoef(part, squeeze(comb_data_an.(tmp.varn)(indw,inds,:)));
    title([tmp.varn, ' ts, lon:', num2str(loninfo), ' lat: ', num2str(latinfo), ' ', num2str(round(tmp.corr(1,2),2))]);
    grid on;
    
    dir.figtgdir=[dir.figdir, '/', 'timeseries'];
        mkdir(dir.figtgdir);
    cfg.figname=[dir.figtgdir, '/', 'ts_',tmp.varn,'_', ...
        str_d, '_', num2str(loninfo),'E_', num2str(latinfo), 'N_', ...
        num2str(min(monrange)),'_',num2str(max(monrange)) '.tif'];
        print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);
end




















%% fig_zosto
fig_cfg.fig_size=[0 0 10 10];
fig_h = figure('name', 'EOF','PaperUnits','inches', ...
        'PaperPosition',(fig_cfg.fig_size),'position', fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','on');

for loni=1:size(comb_data_an.zeta,1)
    for lati=1:size(comb_data_an.zeta,2)
        tmp.corr=corrcoef(comb_data_an.zeta(loni,lati,:), comb_data_an.zosto(loni,lati,:));
        corr_zeta(loni,lati)=tmp.corr(1,2);
    end
end
m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
m_pcolor(lon',lat',corr_zeta'); shading flat; colorbar; caxis([-1 1]); colormap(fig_cfg.c_map);
m_grid;  
m_gshhs_i('color',[1 1 1])  
m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
dir.figtgdir=[dir.figdir, '/', 'check_steric'];
        mkdir(dir.figtgdir);
cfg.figname=[dir.figtgdir, '/', 'corr_zosto',  '.tif'];
    print(fig_h, cfg.figname, '-dpng');

%% fig_zosto_thermo
fig_cfg.fig_size=[0 0 10 10];
fig_h = figure('name', 'EOF','PaperUnits','inches', ...
        'PaperPosition',(fig_cfg.fig_size),'position', fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','on');

for loni=1:size(comb_data_an.zeta,1)
    for lati=1:size(comb_data_an.zeta,2)
        tmp.corr=corrcoef(comb_data_an.zeta(loni,lati,:), comb_data_an.zosto_thermo(loni,lati,:));
        corr_zeta(loni,lati)=tmp.corr(1,2);
    end
end
m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
m_pcolor(lon',lat',corr_zeta'); shading flat; colorbar; caxis([-1 1]); colormap(fig_cfg.c_map);
m_grid;  
m_gshhs_i('color',[1 1 1])  
m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
dir.figtgdir=[dir.figdir, '/', 'check_steric'];
        mkdir(dir.figtgdir);
cfg.figname=[dir.figtgdir, '/', 'corr_zosto_thermo',  '.tif'];
    print(fig_h, cfg.figname, '-dpng');
close all;

%% fig_zosto_halo
fig_cfg.fig_size=[0 0 10 10];
fig_h = figure('name', 'EOF','PaperUnits','inches', ...
        'PaperPosition',(fig_cfg.fig_size),'position', fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','on');

for loni=1:size(comb_data_an.zeta,1)
    for lati=1:size(comb_data_an.zeta,2)
        tmp.corr=corrcoef(comb_data_an.zeta(loni,lati,:), comb_data_an.zosto_halo(loni,lati,:));
        corr_zeta(loni,lati)=tmp.corr(1,2);
    end
end
m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
m_pcolor(lon',lat',corr_zeta'); shading flat; colorbar; caxis([-1 1]); colormap(fig_cfg.c_map);
m_grid;  
m_gshhs_i('color',[1 1 1])  
m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
dir.figtgdir=[dir.figdir, '/', 'check_steric'];
        mkdir(dir.figtgdir);
cfg.figname=[dir.figtgdir, '/', 'corr_zosto_halo',  '.tif'];
    print(fig_h, cfg.figname, '-dpng');
close all;





% % [lv, pc, var_exp] = Func_0024_EOF_3d(data,X);

var_l='zosto_thermo';
% var_l='Vwind';
var_r='zeta';
amod=4;
[mca_result.lv_left, ...
        mca_result.pcs_left, ...
        mca_result.lv_right, ...
        mca_result.pcs_right, ...
        mca_result.lambda, ...
        mca_result.scf] = ...
        mca( double(comb_data_an.(var_l).*mask_model), double(comb_data_an.(var_r).*mask_model), tmp.svd_modes);


close all;
subplot(2,2,1)
mca_result.scf(amod)
lv_l=mca_result.lv_left(:,:,amod);
m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
% m_grid('fontsize', m_grid_fontsize-5, 'tickdir', m_grid_tickdir_type, 'xticklabels', [], 'xtick',(120:20:160), 'parent', ax1_1);  
m_grid;  
m_gshhs_i('color',[1 1 1])  
m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
title([var_l, ' lv, ', num2str(round(mca_result.scf(amod).*100,2)), '%']);
mcax=max(abs(lv_l(:)));
caxis([-mcax mcax])
colormap(fig_cfg.c_map)

subplot(2,2,2);
plot(cfg.years,mca_result.pcs_left(amod,:))
title([var_l, ' pcs, ', num2str(round(mca_result.scf(amod).*100,2)), '%']);
grid on


subplot(2,2,3);
lv_r=mca_result.lv_right(:,:,amod);
m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
m_pcolor(lon',lat',lv_r'); shading flat; colorbar;
% m_grid('fontsize', m_grid_fontsize-5, 'tickdir', m_grid_tickdir_type, 'xticklabels', [], 'xtick',(120:20:160), 'parent', ax1_1);  
m_grid;  
m_gshhs_i('color',[1 1 1]);
m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
title([var_r, ' lv, ', num2str(round(mca_result.scf(amod).*100,2)), '%']);
mcax=max(abs(lv_r(:)));
caxis([-mcax mcax])
colormap(fig_cfg.c_map)

subplot(2,2,4);
plot(cfg.years,mca_result.pcs_right(amod,:))
title([var_r, ' pcs, ', num2str(round(mca_result.scf(amod).*100,2)), '%']);
grid on



% pcolor(lon',lat',mca_result.lv_right(:,:,1)'); shading flat; colorbar;
% pcolor(mask_model'); shading flat; colorbar;