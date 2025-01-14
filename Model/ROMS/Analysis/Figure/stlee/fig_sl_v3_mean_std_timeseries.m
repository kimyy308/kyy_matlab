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
cfg.varnames={'zeta','zosto','zosto_thermo','zosto_halo', 'ubar','vbar','Uwind','Vwind', 'sustr', 'svstr', 'wcurl', 'wscurl'};


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
[fig_cfg.c_map, tmp.err_stat] = Func_0009_get_colormaps('bwg_20', tmp.dropboxpath);

[fig_cfg.c_map2, tmp.err_stat] = Func_0009_get_colormaps('bwg_20', tmp.dropboxpath);



%% depth masked
% tmp.varn='Vwind'

for lm=1:12
    for ly=1:37
        for depthi=1:length(depths)
            dd=depths(depthi);
            str_d=num2str(dd);
            masks.(['d', str_d]);
            for vi=1:length(cfg.varnames)
                tmp.varn=cfg.varnames{vi};  
                comb_data_an.(tmp.varn)(lm,:,:,ly)= squeeze(comb_data_an.(tmp.varn)(lm,:,:,ly)).*mask_model.*masks.(['d', str_d]);
                comb_data_an_mov.(tmp.varn)(lm,:,:,ly)= squeeze(comb_data_an_mov.(tmp.varn)(lm,:,:,ly)).*mask_model.*masks.(['d', str_d]);
            end
        end
    end
end



% %% figure raw
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
si=[3,6,9,12];
% si=[2,5,8,11];
season_particle=7;


    tmp.intv=5;
    tmp.amp_size=5;
    tmp.intv_w=10;
    tmp.amp_size_w=20;

    
%% basic figure
    close all;
    fig_cfg.fig_size=[0 0 14 16];
    fig_h = figure('name', 'mean&std','PaperUnits','inches', ...
            'PaperPosition',(fig_cfg.fig_size),'position', fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','on');
    
    tmp.season=si(1);
    ax_s1=subplot(4,3,1); %Row, Column, order
    lv_u=squeeze(mean(comb_data_an_mov.ubar(tmp.season,:,:,:),4));
    lv_u(lv_u==lv_u(13,110))=NaN; %% make land grids NaN;
    lv_v=squeeze(mean(comb_data_an_mov.vbar(tmp.season,:,:,:),4));
    lv_v(lv_v==lv_v(13,110))=NaN; %% make land grids NaN;

    m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
    m_quiver(lon(1:tmp.intv:end, 1:tmp.intv:end)', ...
                    lat(1:tmp.intv:end, 1:tmp.intv:end)', ...
                    lv_u(1:tmp.intv:end, 1:tmp.intv:end)' * tmp.amp_size, ...
                    lv_v(1:tmp.intv:end, 1:tmp.intv:end)' * tmp.amp_size, ...
                    'AutoScale','off','LineWidth', 1);     m_grid;  
    m_gshhs_i('color',[1 1 1]);
    m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
    title(['vec', ' mean, ', num2str(tmp.season), '-', num2str(tmp.season+2), 'M']);
    
    ax_s2=subplot(4,3,2);
    lv_u=squeeze(comb_data_an_mov.ubar(tmp.season,:,:,:));
    lv_u(lv_u==lv_u(13,110,1))=NaN; %% make land grids NaN;
    [lv_u, error_status] = Func_0011_get_area_weighted_mean(lv_u, lon, lat);
    yyaxis left;
    plot(cfg.years,lv_u);
    yyaxis right;
    plot(cfg.years,part3_mov.tccount(:,season_particle));
    title([num2str(tmp.season), '-', num2str(tmp.season+2), 'M', ' btv-u, ', 'vs JAS par']);
    tmp.corr=corrcoef(lv_u,part3_mov.tccount(:,season_particle));
    text(cfg.years(1), ax_s2.YLim(1)+(ax_s2.YLim(2)-ax_s2.YLim(1))*0.04, ['r: ', num2str(round(tmp.corr(1,2),2))]);

    ax_s2=subplot(4,3,3);
    lv_v=squeeze(comb_data_an_mov.vbar(tmp.season,:,:,:));
    lv_v(lv_v==lv_v(13,110,1))=NaN; %% make land grids NaN;
    [lv_v, error_status] = Func_0011_get_area_weighted_mean(lv_v, lon, lat);
    yyaxis left;
    plot(cfg.years,lv_v);
    yyaxis right;
    plot(cfg.years,part3_mov.tccount(:,season_particle));
    title([num2str(tmp.season), '-', num2str(tmp.season+2), 'M', ' btv-v, ', 'vs JAS par']);
    tmp.corr=corrcoef(lv_v,part3_mov.tccount(:,season_particle));
    text(cfg.years(1), ax_s2.YLim(1)+(ax_s2.YLim(2)-ax_s2.YLim(1))*0.04, ['r: ', num2str(round(tmp.corr(1,2),2))]);



    tmp.season=si(2);
    ax_s1=subplot(4,3,4); %Row, Column, order
    lv_u=squeeze(mean(comb_data_an_mov.ubar(tmp.season,:,:,:),4));
    lv_u(lv_u==lv_u(13,110))=NaN; %% make land grids NaN;
    lv_v=squeeze(mean(comb_data_an_mov.vbar(tmp.season,:,:,:),4));
    lv_v(lv_v==lv_v(13,110))=NaN; %% make land grids NaN;

    m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
    m_quiver(lon(1:tmp.intv:end, 1:tmp.intv:end)', ...
                    lat(1:tmp.intv:end, 1:tmp.intv:end)', ...
                    lv_u(1:tmp.intv:end, 1:tmp.intv:end)' * tmp.amp_size, ...
                    lv_v(1:tmp.intv:end, 1:tmp.intv:end)' * tmp.amp_size, ...
                    'AutoScale','off','LineWidth', 1);     m_grid;  
    m_gshhs_i('color',[1 1 1]);
    m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
    title(['vec', ' mean, ', num2str(tmp.season), '-', num2str(tmp.season+2), 'M']);
    
    ax_s2=subplot(4,3,5);
    lv_u=squeeze(comb_data_an_mov.ubar(tmp.season,:,:,:));
    lv_u(lv_u==lv_u(13,110,1))=NaN; %% make land grids NaN;
    [lv_u, error_status] = Func_0011_get_area_weighted_mean(lv_u, lon, lat);
    yyaxis left;
    plot(cfg.years,lv_u);
    yyaxis right;
    plot(cfg.years,part3_mov.tccount(:,season_particle));
    title([num2str(tmp.season), '-', num2str(tmp.season+2), 'M', ' btv-u, ', 'vs JAS par']);
    tmp.corr=corrcoef(lv_u,part3_mov.tccount(:,season_particle));
    text(cfg.years(1), ax_s2.YLim(1)+(ax_s2.YLim(2)-ax_s2.YLim(1))*0.04, ['r: ', num2str(round(tmp.corr(1,2),2))]);

    ax_s2=subplot(4,3,6);
    lv_v=squeeze(comb_data_an_mov.vbar(tmp.season,:,:,:));
    lv_v(lv_v==lv_v(13,110,1))=NaN; %% make land grids NaN;
    [lv_v, error_status] = Func_0011_get_area_weighted_mean(lv_v, lon, lat);
    yyaxis left;
    plot(cfg.years,lv_v);
    yyaxis right;
    plot(cfg.years,part3_mov.tccount(:,season_particle));
    title([num2str(tmp.season), '-', num2str(tmp.season+2), 'M', ' btv-v, ', 'vs JAS par']);
    tmp.corr=corrcoef(lv_v,part3_mov.tccount(:,season_particle));
    text(cfg.years(1), ax_s2.YLim(1)+(ax_s2.YLim(2)-ax_s2.YLim(1))*0.04, ['r: ', num2str(round(tmp.corr(1,2),2))]);



    tmp.season=si(3);
    ax_s1=subplot(4,3,7); %Row, Column, order
    lv_u=squeeze(mean(comb_data_an_mov.ubar(tmp.season,:,:,:),4));
    lv_u(lv_u==lv_u(13,110))=NaN; %% make land grids NaN;
    lv_v=squeeze(mean(comb_data_an_mov.vbar(tmp.season,:,:,:),4));
    lv_v(lv_v==lv_v(13,110))=NaN; %% make land grids NaN;

    m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
    m_quiver(lon(1:tmp.intv:end, 1:tmp.intv:end)', ...
                    lat(1:tmp.intv:end, 1:tmp.intv:end)', ...
                    lv_u(1:tmp.intv:end, 1:tmp.intv:end)' * tmp.amp_size, ...
                    lv_v(1:tmp.intv:end, 1:tmp.intv:end)' * tmp.amp_size, ...
                    'AutoScale','off','LineWidth', 1);     m_grid;  
    m_gshhs_i('color',[1 1 1]);
    m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
    title(['vec', ' mean, ', num2str(tmp.season), '-', num2str(tmp.season+2), 'M']);
    
    ax_s2=subplot(4,3,8);
    lv_u=squeeze(comb_data_an_mov.ubar(tmp.season,:,:,:));
    lv_u(lv_u==lv_u(13,110,1))=NaN; %% make land grids NaN;
    [lv_u, error_status] = Func_0011_get_area_weighted_mean(lv_u, lon, lat);
    yyaxis left;
    plot(cfg.years,lv_u);
    yyaxis right;
    plot(cfg.years,part3_mov.tccount(:,season_particle));
    title([num2str(tmp.season), '-', num2str(tmp.season+2), 'M', ' btv-u, ', 'vs JAS par']);
    tmp.corr=corrcoef(lv_u(1:end-1),part3_mov.tccount(2:end,season_particle));
    text(cfg.years(1), ax_s2.YLim(1)+(ax_s2.YLim(2)-ax_s2.YLim(1))*0.04, ['r: ', num2str(round(tmp.corr(1,2),2))]);

    ax_s2=subplot(4,3,9);
    lv_v=squeeze(comb_data_an_mov.vbar(tmp.season,:,:,:));
    lv_v(lv_v==lv_v(13,110,1))=NaN; %% make land grids NaN;
    [lv_v, error_status] = Func_0011_get_area_weighted_mean(lv_v, lon, lat);
    yyaxis left;
    plot(cfg.years,lv_v);
    yyaxis right;
    plot(cfg.years,part3_mov.tccount(:,season_particle));
    title([num2str(tmp.season), '-', num2str(tmp.season+2), 'M', ' btv-v, ', 'vs JAS par']);
    tmp.corr=corrcoef(lv_v(1:end-1),part3_mov.tccount(2:end,season_particle));
    text(cfg.years(1), ax_s2.YLim(1)+(ax_s2.YLim(2)-ax_s2.YLim(1))*0.04, ['r: ', num2str(round(tmp.corr(1,2),2))]);

        

    tmp.season=si(4);
    ax_s1=subplot(4,3,10); %Row, Column, order
    if si(4)>10
        lv_u=squeeze(mean(comb_data_an_mov.ubar(tmp.season,:,:,1:end-1),4));
        lv_u(lv_u==lv_u(13,110))=NaN; %% make land grids NaN;
        lv_v=squeeze(mean(comb_data_an_mov.vbar(tmp.season,:,:,1:end-1),4));
        lv_v(lv_v==lv_v(13,110))=NaN; %% make land grids NaN;
    else
        lv_u=squeeze(mean(comb_data_an_mov.ubar(tmp.season,:,:,:),4));
        lv_u(lv_u==lv_u(13,110))=NaN; %% make land grids NaN;
        lv_v=squeeze(mean(comb_data_an_mov.vbar(tmp.season,:,:,:),4));
        lv_v(lv_v==lv_v(13,110))=NaN; %% make land grids NaN;
    end

    m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
    m_quiver(lon(1:tmp.intv:end, 1:tmp.intv:end)', ...
                    lat(1:tmp.intv:end, 1:tmp.intv:end)', ...
                    lv_u(1:tmp.intv:end, 1:tmp.intv:end)' * tmp.amp_size, ...
                    lv_v(1:tmp.intv:end, 1:tmp.intv:end)' * tmp.amp_size, ...
                    'AutoScale','off','LineWidth', 1);     m_grid;  
    m_gshhs_i('color',[1 1 1]);
    m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
    title(['vec', ' mean, ', num2str(tmp.season), '-', num2str(tmp.season+2), 'M']);
    
    ax_s2=subplot(4,3,11);
    if si(4)>10
        lv_u=squeeze(comb_data_an_mov.ubar(tmp.season,:,:,1:end-1));
    else
        lv_u=squeeze(comb_data_an_mov.ubar(tmp.season,:,:,:));
    end
    lv_u(lv_u==lv_u(13,110,1))=NaN; %% make land grids NaN;
    [lv_u, error_status] = Func_0011_get_area_weighted_mean(lv_u, lon, lat);
    yyaxis left;
    plot(cfg.years(1:end-1),lv_u);
    yyaxis right;
    plot(cfg.years,part3_mov.tccount(:,season_particle));
    title([num2str(tmp.season), '-', num2str(tmp.season+2), 'M', ' btv-u, ', 'vs JAS par']);
    tmp.corr=corrcoef(lv_u,part3_mov.tccount(2:end,season_particle));
    text(cfg.years(1), ax_s2.YLim(1)+(ax_s2.YLim(2)-ax_s2.YLim(1))*0.04, ['r: ', num2str(round(tmp.corr(1,2),2))]);

    ax_s2=subplot(4,3,12);
    if si(4)>10
        lv_v=squeeze(comb_data_an_mov.vbar(tmp.season,:,:,1:end-1));
    else
        lv_v=squeeze(comb_data_an_mov.vbar(tmp.season,:,:,:));
    end
    lv_v(lv_v==lv_v(13,110,1))=NaN; %% make land grids NaN;
    [lv_v, error_status] = Func_0011_get_area_weighted_mean(lv_v, lon, lat);
    yyaxis left;
    plot(cfg.years(1:end-1),lv_v);
    yyaxis right;
    plot(cfg.years,part3_mov.tccount(:,season_particle));
    title([num2str(tmp.season), '-', num2str(tmp.season+2), 'M', ' btv-v, ', 'vs JAS par']);
    tmp.corr=corrcoef(lv_v,part3_mov.tccount(2:end,season_particle));
    text(cfg.years(1), ax_s2.YLim(1)+(ax_s2.YLim(2)-ax_s2.YLim(1))*0.04, ['r: ', num2str(round(tmp.corr(1,2),2))]);




    dir.figtgdir=[dir.figdir, '/', 'mean_std_ts'];
    mkdir(dir.figtgdir);
    
    cfg.figname=[dir.figtgdir, '/', 'm_','vec','_season_',num2str(si), '_', 'JAS_', str_d '.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);




%% basic figure (wstr)
    close all;
    fig_cfg.fig_size=[0 0 14 16];
    fig_h = figure('name', 'mean&std','PaperUnits','inches', ...
            'PaperPosition',(fig_cfg.fig_size),'position', fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','on');
    
    tmp.season=si(1);
    ax_s1=subplot(4,3,1); %Row, Column, order
    lv_u=squeeze(mean(comb_data_an_mov.sustr(tmp.season,:,:,:),4));
    lv_u(lv_u==lv_u(13,110))=NaN; %% make land grids NaN;
    lv_v=squeeze(mean(comb_data_an_mov.svstr(tmp.season,:,:,:),4));
    lv_v(lv_v==lv_v(13,110))=NaN; %% make land grids NaN;

    m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
    m_quiver(lon(1:tmp.intv_w:end, 1:tmp.intv_w:end)', ...
                    lat(1:tmp.intv_w:end, 1:tmp.intv_w:end)', ...
                    lv_u(1:tmp.intv_w:end, 1:tmp.intv_w:end)' * tmp.amp_size_w, ...
                    lv_v(1:tmp.intv_w:end, 1:tmp.intv_w:end)' * tmp.amp_size_w, ...
                    'AutoScale','off','LineWidth', 1);     m_grid;  
    m_gshhs_i('color',[1 1 1]);
    m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
    title(['vec', ' mean, ', num2str(tmp.season), '-', num2str(tmp.season+2), 'M']);
    
    ax_s2=subplot(4,3,2);
    lv_u=squeeze(comb_data_an_mov.sustr(tmp.season,:,:,:));
    lv_u(lv_u==lv_u(13,110,1))=NaN; %% make land grids NaN;
    [lv_u, error_status] = Func_0011_get_area_weighted_mean(lv_u, lon, lat);
    yyaxis left;
    plot(cfg.years,lv_u);
    yyaxis right;
    plot(cfg.years,part3_mov.tccount(:,season_particle));
    title([num2str(tmp.season), '-', num2str(tmp.season+2), 'M', ' btv-u, ', 'vs JAS par']);
    tmp.corr=corrcoef(lv_u,part3_mov.tccount(:,season_particle));
    text(cfg.years(1), ax_s2.YLim(1)+(ax_s2.YLim(2)-ax_s2.YLim(1))*0.04, ['r: ', num2str(round(tmp.corr(1,2),2))]);

    ax_s2=subplot(4,3,3);
    lv_v=squeeze(comb_data_an_mov.svstr(tmp.season,:,:,:));
    lv_v(lv_v==lv_v(13,110,1))=NaN; %% make land grids NaN;
    [lv_v, error_status] = Func_0011_get_area_weighted_mean(lv_v, lon, lat);
    yyaxis left;
    plot(cfg.years,lv_v);
    yyaxis right;
    plot(cfg.years,part3_mov.tccount(:,season_particle));
    title([num2str(tmp.season), '-', num2str(tmp.season+2), 'M', ' btv-v, ', 'vs JAS par']);
    tmp.corr=corrcoef(lv_v,part3_mov.tccount(:,season_particle));
    text(cfg.years(1), ax_s2.YLim(1)+(ax_s2.YLim(2)-ax_s2.YLim(1))*0.04, ['r: ', num2str(round(tmp.corr(1,2),2))]);



    tmp.season=si(2);
    ax_s1=subplot(4,3,4); %Row, Column, order
    lv_u=squeeze(mean(comb_data_an_mov.sustr(tmp.season,:,:,:),4));
    lv_u(lv_u==lv_u(13,110))=NaN; %% make land grids NaN;
    lv_v=squeeze(mean(comb_data_an_mov.svstr(tmp.season,:,:,:),4));
    lv_v(lv_v==lv_v(13,110))=NaN; %% make land grids NaN;

    m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
    m_quiver(lon(1:tmp.intv_w:end, 1:tmp.intv_w:end)', ...
                    lat(1:tmp.intv_w:end, 1:tmp.intv_w:end)', ...
                    lv_u(1:tmp.intv_w:end, 1:tmp.intv_w:end)' * tmp.amp_size_w, ...
                    lv_v(1:tmp.intv_w:end, 1:tmp.intv_w:end)' * tmp.amp_size_w, ...
                    'AutoScale','off','LineWidth', 1);     m_grid;  
    m_gshhs_i('color',[1 1 1]);
    m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
    title(['vec', ' mean, ', num2str(tmp.season), '-', num2str(tmp.season+2), 'M']);
    
    ax_s2=subplot(4,3,5);
    lv_u=squeeze(comb_data_an_mov.sustr(tmp.season,:,:,:));
    lv_u(lv_u==lv_u(13,110,1))=NaN; %% make land grids NaN;
    [lv_u, error_status] = Func_0011_get_area_weighted_mean(lv_u, lon, lat);
    yyaxis left;
    plot(cfg.years,lv_u);
    yyaxis right;
    plot(cfg.years,part3_mov.tccount(:,season_particle));
    title([num2str(tmp.season), '-', num2str(tmp.season+2), 'M', ' btv-u, ', 'vs JAS par']);
    tmp.corr=corrcoef(lv_u,part3_mov.tccount(:,season_particle));
    text(cfg.years(1), ax_s2.YLim(1)+(ax_s2.YLim(2)-ax_s2.YLim(1))*0.04, ['r: ', num2str(round(tmp.corr(1,2),2))]);

    ax_s2=subplot(4,3,6);
    lv_v=squeeze(comb_data_an_mov.svstr(tmp.season,:,:,:));
    lv_v(lv_v==lv_v(13,110,1))=NaN; %% make land grids NaN;
    [lv_v, error_status] = Func_0011_get_area_weighted_mean(lv_v, lon, lat);
    yyaxis left;
    plot(cfg.years,lv_v);
    yyaxis right;
    plot(cfg.years,part3_mov.tccount(:,season_particle));
    title([num2str(tmp.season), '-', num2str(tmp.season+2), 'M', ' btv-v, ', 'vs JAS par']);
    tmp.corr=corrcoef(lv_v,part3_mov.tccount(:,season_particle));
    text(cfg.years(1), ax_s2.YLim(1)+(ax_s2.YLim(2)-ax_s2.YLim(1))*0.04, ['r: ', num2str(round(tmp.corr(1,2),2))]);



    tmp.season=si(3);
    ax_s1=subplot(4,3,7); %Row, Column, order
    lv_u=squeeze(mean(comb_data_an_mov.sustr(tmp.season,:,:,:),4));
    lv_u(lv_u==lv_u(13,110))=NaN; %% make land grids NaN;
    lv_v=squeeze(mean(comb_data_an_mov.svstr(tmp.season,:,:,:),4));
    lv_v(lv_v==lv_v(13,110))=NaN; %% make land grids NaN;

    m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
    m_quiver(lon(1:tmp.intv_w:end, 1:tmp.intv_w:end)', ...
                    lat(1:tmp.intv_w:end, 1:tmp.intv_w:end)', ...
                    lv_u(1:tmp.intv_w:end, 1:tmp.intv_w:end)' * tmp.amp_size_w, ...
                    lv_v(1:tmp.intv_w:end, 1:tmp.intv_w:end)' * tmp.amp_size_w, ...
                    'AutoScale','off','LineWidth', 1);     m_grid;  
    m_gshhs_i('color',[1 1 1]);
    m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
    title(['vec', ' mean, ', num2str(tmp.season), '-', num2str(tmp.season+2), 'M']);
    
    ax_s2=subplot(4,3,8);
    lv_u=squeeze(comb_data_an_mov.sustr(tmp.season,:,:,:));
    lv_u(lv_u==lv_u(13,110,1))=NaN; %% make land grids NaN;
    [lv_u, error_status] = Func_0011_get_area_weighted_mean(lv_u, lon, lat);
    yyaxis left;
    plot(cfg.years,lv_u);
    yyaxis right;
    plot(cfg.years,part3_mov.tccount(:,season_particle));
    title([num2str(tmp.season), '-', num2str(tmp.season+2), 'M', ' btv-u, ', 'vs JAS par']);
    tmp.corr=corrcoef(lv_u(1:end-1),part3_mov.tccount(2:end,season_particle));
    text(cfg.years(1), ax_s2.YLim(1)+(ax_s2.YLim(2)-ax_s2.YLim(1))*0.04, ['r: ', num2str(round(tmp.corr(1,2),2))]);

    ax_s2=subplot(4,3,9);
    lv_v=squeeze(comb_data_an_mov.svstr(tmp.season,:,:,:));
    lv_v(lv_v==lv_v(13,110,1))=NaN; %% make land grids NaN;
    [lv_v, error_status] = Func_0011_get_area_weighted_mean(lv_v, lon, lat);
    yyaxis left;
    plot(cfg.years,lv_v);
    yyaxis right;
    plot(cfg.years,part3_mov.tccount(:,season_particle));
    title([num2str(tmp.season), '-', num2str(tmp.season+2), 'M', ' btv-v, ', 'vs JAS par']);
    tmp.corr=corrcoef(lv_v(1:end-1),part3_mov.tccount(2:end,season_particle));
    text(cfg.years(1), ax_s2.YLim(1)+(ax_s2.YLim(2)-ax_s2.YLim(1))*0.04, ['r: ', num2str(round(tmp.corr(1,2),2))]);

        

    tmp.season=si(4);
    ax_s1=subplot(4,3,10); %Row, Column, order
    if si(4)>10
        lv_u=squeeze(mean(comb_data_an_mov.sustr(tmp.season,:,:,1:end-1),4));
        lv_u(lv_u==lv_u(13,110))=NaN; %% make land grids NaN;
        lv_v=squeeze(mean(comb_data_an_mov.svstr(tmp.season,:,:,1:end-1),4));
        lv_v(lv_v==lv_v(13,110))=NaN; %% make land grids NaN;
    else
        lv_u=squeeze(mean(comb_data_an_mov.sustr(tmp.season,:,:,:),4));
        lv_u(lv_u==lv_u(13,110))=NaN; %% make land grids NaN;
        lv_v=squeeze(mean(comb_data_an_mov.svstr(tmp.season,:,:,:),4));
        lv_v(lv_v==lv_v(13,110))=NaN; %% make land grids NaN;
    end

    m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
    m_quiver(lon(1:tmp.intv_w:end, 1:tmp.intv_w:end)', ...
                    lat(1:tmp.intv_w:end, 1:tmp.intv_w:end)', ...
                    lv_u(1:tmp.intv_w:end, 1:tmp.intv_w:end)' * tmp.amp_size_w, ...
                    lv_v(1:tmp.intv_w:end, 1:tmp.intv_w:end)' * tmp.amp_size_w, ...
                    'AutoScale','off','LineWidth', 1);     m_grid;  
    m_gshhs_i('color',[1 1 1]);
    m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
    title(['vec', ' mean, ', num2str(tmp.season), '-', num2str(tmp.season+2), 'M']);
    
    ax_s2=subplot(4,3,11);
    if si(4)>10
        lv_u=squeeze(comb_data_an_mov.sustr(tmp.season,:,:,1:end-1));
    else
        lv_u=squeeze(comb_data_an_mov.sustr(tmp.season,:,:,:));
    end
    lv_u(lv_u==lv_u(13,110,1))=NaN; %% make land grids NaN;
    [lv_u, error_status] = Func_0011_get_area_weighted_mean(lv_u, lon, lat);
    yyaxis left;
    plot(cfg.years(1:end-1),lv_u);
    yyaxis right;
    plot(cfg.years,part3_mov.tccount(:,season_particle));
    title([num2str(tmp.season), '-', num2str(tmp.season+2), 'M', ' btv-u, ', 'vs JAS par']);
    tmp.corr=corrcoef(lv_u,part3_mov.tccount(2:end,season_particle));
    text(cfg.years(1), ax_s2.YLim(1)+(ax_s2.YLim(2)-ax_s2.YLim(1))*0.04, ['r: ', num2str(round(tmp.corr(1,2),2))]);

    ax_s2=subplot(4,3,12);
    if si(4)>10
        lv_v=squeeze(comb_data_an_mov.svstr(tmp.season,:,:,1:end-1));
    else
        lv_v=squeeze(comb_data_an_mov.svstr(tmp.season,:,:,:));
    end
    lv_v(lv_v==lv_v(13,110,1))=NaN; %% make land grids NaN;
    [lv_v, error_status] = Func_0011_get_area_weighted_mean(lv_v, lon, lat);
    yyaxis left;
    plot(cfg.years(1:end-1),lv_v);
    yyaxis right;
    plot(cfg.years,part3_mov.tccount(:,season_particle));
    title([num2str(tmp.season), '-', num2str(tmp.season+2), 'M', ' btv-v, ', 'vs JAS par']);
    tmp.corr=corrcoef(lv_v,part3_mov.tccount(2:end,season_particle));
    text(cfg.years(1), ax_s2.YLim(1)+(ax_s2.YLim(2)-ax_s2.YLim(1))*0.04, ['r: ', num2str(round(tmp.corr(1,2),2))]);




    dir.figtgdir=[dir.figdir, '/', 'mean_std_ts'];
    mkdir(dir.figtgdir);
    
    cfg.figname=[dir.figtgdir, '/', 'm_','wstr','_season_',num2str(si), '_', 'JAS_', str_d '.tif'];
    print(fig_h, cfg.figname, '-dpng');
    RemoveWhiteSpace([], 'file', cfg.figname);

    close all;