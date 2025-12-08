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
% cfg.varnames={'vec', 'wvec', 'zeta','zosto','zosto_thermo','zosto_halo', 'wsvec', 'ubar','vbar','Uwind','Vwind', 'sustr', 'svstr', 'wcurl', 'wscurl'};

% cfg.varnames={'svstr', 'ubar','vbar', 'vec', 'zeta','zosto','zosto_thermo','zosto_halo', 'sustr',  'wsvec','wscurl'};
cfg.varnames={'svstr', 'zeta'};

% dir.figdir='/Volumes/kyy_raid/kimyy/Research/Postdoc/03_IBS/2024_TS_KS_connectivity/Figure_fin';
dir.figdir='/Volumes/kyy_raid/kimyy/Research/Postdoc/03_IBS/2024_TS_KS_connectivity/Figure_fin_250708';

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


% seasonal mean (5~7, summer)
for vi=1:length(cfg.varnames)
    tmp.varn=cfg.varnames{vi};
    for lm=1:6
        monrange=lm:lm+2;
        tmp.var=reshape(comb_data.(tmp.varn), [xn yn 12 tn/12]);
        comb_data_an_mov.(tmp.varn)(lm,:,:,:)=squeeze(mean(tmp.var(:,:,monrange,:),3));
    end
end



part=load('/Volumes/kyy_raid/kimyy/Research/Postdoc/03_IBS/2024_TS_KS_connectivity/layer_particle_month_2.mat');

% part.years=1982:2020;
% part.months=1:12;

part_names= {'tcage', 'wcage', 'ecage', 'tcdepth', 'wcdepth', 'ecdepth', 'tccount', 'wccount', 'eccount'};

for ni=1:length(part_names)
    tmp.data=squeeze(sum(part.(part_names{ni}), 1, 'omitnan'));
    part2.(part_names{ni})=tmp.data(2:end-1,:); % subsample 1983~2019
    for lm=1:6
        monrange=lm:lm;
        part3.(part_names{ni})(:,lm)=squeeze(mean(part2.(part_names{ni})(:,monrange),2));
    end
    for lm=1:6
        monrange=lm:lm+2;
        part3_mov.(part_names{ni})(:,lm)=squeeze(mean(part2.(part_names{ni})(:,monrange),2));
    end
end

% part_stlee
% stleemat='/Volumes/kyy_raid/kimyy/Research/Postdoc/03_IBS/2024_TS_KS_connectivity/num_particle_ks_stlee.mat';
% stleemat='/Volumes/kyy_raid/kimyy/Research/Postdoc/03_IBS/2024_TS_KS_connectivity/num_particle_ks_arrive_stlee2.mat';
stleemat='/Volumes/kyy_raid/kimyy/Research/Postdoc/03_IBS/2024_TS_KS_connectivity/manuscript/Revision_round1/num_particle_KS_hor_RW_250708.mat';

part_stlee=load(stleemat); % tt_t, num_A_t; 1982-01 ~ 2019-12;
fileInfo = matfile(stleemat);
whos(fileInfo);


% part_stlee.num_A_t=part_stlee.num_A_o_m(13:end,:);
% part_stlee.tt_t=part_stlee.tt_o_m(13:end,:);
% tn=size(part_stlee.tt_t, 1);
% part_stlee2.num_A_t=reshape(part_stlee.num_A_t, [12, tn/12, 12]);
% part_stlee2.tt_t=reshape(part_stlee.tt_t, [12, tn/12, 12]);
% part_stlee2.num_A_t=squeeze(sum(part_stlee2.num_A_t, 1));
% part_stlee2.tt_t=squeeze(sum(part_stlee2.tt_t, 1));
% part_stlee2.num_A_t = permute(part_stlee2.num_A_t, [2, 1]);
% part_stlee2.tt_t = permute(part_stlee2.tt_t, [2, 1]);

clear part3_mov_stlee
% for lm=1:9 % 1:12, month
%     monrange=lm:lm+3;
%     part3_mov_stlee.num_A_t(lm,:)=squeeze(mean(part_stlee2.num_A_t(monrange,:),1));
% end


part3_mov_stlee.num_A_t(6,:)=part_stlee.num_ptc_mean;


% wrong
% clear part3_mov_stlee
% part_stlee2.num_A_t=part_stlee.num_A_t(13:end);
% part_stlee2.tt_t=part_stlee.tt_t(13:end);
% tn=length(part_stlee2.num_A_t);
% part_stlee2.num_A_t=reshape(part_stlee2.num_A_t, [tn/12 12]);
% part_stlee2.tt_t=reshape(part_stlee2.tt_t, [tn/12 12]);
% for lm=1:9 % 1:12, month
%     monrange=lm:lm+3;
%     part3_mov_stlee.num_A_t(:,lm)=squeeze(mean(part_stlee2.num_A_t(:,monrange),2));
% end


plot(1983:2019,part3_mov_stlee.num_A_t(6,:))


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
    130, 38;
    120, 38;
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
    131, 38;
    122, 38;
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

for lm=6:6
    for depthi=1:length(depths)
        dd=depths(depthi);
        str_d=num2str(dd);
        masks.(['d', str_d]);
    
        tmp.svd_modes=4;
        for vi=1:length(cfg.varnames)
            tmp.varn=cfg.varnames{vi};            

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
    
        %% sign check
        for vi=1:length(cfg.varnames)
            tmp.varn=cfg.varnames{vi};
            for amodi=1:tmp.svd_modes
                if mean(mean(squeeze(EOF_result_mov.(tmp.varn).lv(lm,:,:,amodi)),'omitnan'),'omitnan')<0
                    EOF_result_mov.(tmp.varn).lv(lm,:,:,amodi)=-EOF_result_mov.(tmp.varn).lv(lm,:,:,amodi);
                    EOF_result_mov.(tmp.varn).pcs(lm,:,amodi)=-EOF_result_mov.(tmp.varn).pcs(lm,:,amodi);
                end
                
                if amodi==1
                    tmp.data_mov=squeeze(comb_data_an_mov.(tmp.varn)(lm,:,:,:));
                    tmp.data_mov(abs(tmp.data_mov)>10e6)=NaN;
                    tmp.data_mov_m=Func_0011_get_area_weighted_mean( ...
                        tmp.data_mov, lon, lat);                    
                end
            end
        end
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

% Spring: 3~5;
% Summer: 6~8;
% Autumn: 9~11;
% Winter: 12~2;

%% EOF LV & PCT (1-4 modes) - seasonal figures
for pmon= [6]
%     for si=[2, 3, 5, 6, 11, 12]
%     for si=[5, 2, 11]
    for si=[6]
    % for si=[9, 12]
        for vi=1:length(cfg.varnames)
            tmp.varn=cfg.varnames{vi};
            if si<=pmon
                part_tc = part3_mov.tccount(:,pmon); 
                part_wc = part3_mov.wccount(:,pmon); 
                part_ec = part3_mov.eccount(:,pmon); 
            else
                part_tc = part3_mov.tccount(2:end,pmon); 
                part_tc = [part_tc; NaN];
                part_wc = part3_mov.wccount(2:end,pmon); 
                part_wc = [part_wc; NaN];
                part_ec = part3_mov.eccount(2:end,pmon); 
                part_ec = [part_ec; NaN];
            end
        
            
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
            
            subplot(4,2,2); % Figure 4b
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



% %% key figure template
%             close all;
%             fig_cfg.fig_size=[0 0 14 16];
%             fig_h = figure('name', 'EOF','PaperUnits','inches', ...
%                     'PaperPosition',(fig_cfg.fig_size),'position', fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
%             
%             tmp.varn='svstr';
%             amod=1;
%             subplot(2,2,1); 
%             lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
%             m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
%             m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
%             % m_grid('fontsize', m_grid_fontsize-5, 'tickdir', m_grid_tickdir_type, 'xticklabels', [], 'xtick',(120:20:160), 'parent', ax1_1);  
%             m_grid;  
%             m_gshhs_i('color',[1 1 1]);
%             m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
%             title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
%             mcax=max(abs(lv_l(:)));
%             caxis([-mcax mcax]);
%             colormap(fig_cfg.c_map);
%             
%             subplot(2,2,2); 
%             plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod));
%             tmp.corr1=corrcoef(part_tc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
%             tmp.corr2=corrcoef(part_wc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
%             tmp.corr3=corrcoef(part_ec, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
%             title([tmp.varn, ' pcs, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%, ', ...
%                 'r(tc):', num2str(round(tmp.corr1(1,2),2)), ...
%                 ', r(wc):', num2str(round(tmp.corr2(1,2),2)), ...
%                 ', r(ec):', num2str(round(tmp.corr3(1,2),2))]);
%             grid on;
%             
%             
%             tmp.varn='zeta';
%             amod=2;
%             subplot(2,2,3); 
%             lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
%             m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
%             m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
%             m_grid;  
%             m_gshhs_i('color',[1 1 1]);
%             m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
%             title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
%             mcax=max(abs(lv_l(:)));
%             caxis([-mcax mcax]);
%             colormap(fig_cfg.c_map);
%             
%             subplot(2,2,4); 
%             plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod));
%             tmp.corr1=corrcoef(part_tc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
%             tmp.corr2=corrcoef(part_wc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
%             tmp.corr3=corrcoef(part_ec, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
%             title([tmp.varn, ' pcs, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%, ', ...
%                 'r(tc):', num2str(round(tmp.corr1(1,2),2)), ...
%                 ', r(wc):', num2str(round(tmp.corr2(1,2),2)), ...
%                 ', r(ec):', num2str(round(tmp.corr3(1,2),2))]);
%             grid on;
% 
% 
%             dir.figtgdir=[dir.figdir, '/', 'EOF_mov', '/', 'par', num2str(pmon)];
%             mkdir(dir.figtgdir);
%             
%             cfg.figname=[dir.figtgdir, '/', 'Fig_figure_EOF_mov_', '_season_',num2str(si), '_', str_d, '_p_', num2str(pmon) '.tif'];
%             print(fig_h, cfg.figname, '-dpng');
%             RemoveWhiteSpace([], 'file', cfg.figname);


% %% key figure template (par_stlee)
%             close all;
%             fig_cfg.fig_size=[0 0 14 16];
%             fig_h = figure('name', 'EOF','PaperUnits','inches', ...
%                     'PaperPosition',(fig_cfg.fig_size),'position', fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
%             
%             tmp.varn='svstr';
%             amod=1;
%             subplot(2,2,1);
%             lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
%             m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
%             m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
%             % m_grid('fontsize', m_grid_fontsize-5, 'tickdir', m_grid_tickdir_type, 'xticklabels', [], 'xtick',(120:20:160), 'parent', ax1_1);  
%             m_grid;  
%             m_gshhs_i('color',[1 1 1]);
%             m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
%             title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
%             mcax=max(abs(lv_l(:)));
%             caxis([-mcax mcax]);
%             colormap(fig_cfg.c_map);
%             
%             subplot(2,2,2);
%             abc=part3_mov_stlee.num_A_t(6,:);
% %             abc=part3_mov_stlee.num_A_t(:,6);
%             plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod));
%             tmp.corr1=corrcoef(abc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
%             tmp.corr2=corrcoef(abc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
%             tmp.corr3=corrcoef(abc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
%             title([tmp.varn, ' pcs, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%, ', ...
%                 'r(tc):', num2str(round(tmp.corr1(1,2),2)), ...
%                 ', r(wc):', num2str(round(tmp.corr2(1,2),2)), ...
%                 ', r(ec):', num2str(round(tmp.corr3(1,2),2))]);
%             grid on;
%             
%             
%             tmp.varn='zeta';
%             amod=2;
%             subplot(2,2,3);
%             lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
%             m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
%             m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
%             m_grid;  
%             m_gshhs_i('color',[1 1 1]);
%             m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
%             title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
%             mcax=max(abs(lv_l(:)));
%             caxis([-mcax mcax]);
%             colormap(fig_cfg.c_map);
%             
%             subplot(2,2,4);
%             plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod));
%             tmp.corr1=corrcoef(abc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
%             tmp.corr2=corrcoef(abc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
%             tmp.corr3=corrcoef(abc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
%             title([tmp.varn, ' pcs, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%, ', ...
%                 'r(tc):', num2str(round(tmp.corr1(1,2),2)), ...
%                 ', r(wc):', num2str(round(tmp.corr2(1,2),2)), ...
%                 ', r(ec):', num2str(round(tmp.corr3(1,2),2))]);
%             grid on;
% 
% 
%             dir.figtgdir=[dir.figdir, '/', 'EOF_mov', '/', 'par', num2str(pmon)];
%             mkdir(dir.figtgdir);
%             
%             cfg.figname=[dir.figtgdir, '/', 'Fig_figure_stleepar_EOF_mov_', '_season_',num2str(si), '_', str_d, '_p_', num2str(pmon) '.tif'];
%             print(fig_h, cfg.figname, '-dpng');
%             RemoveWhiteSpace([], 'file', cfg.figname);



% %% key figure draft
%             fntsize=15;
%             close all;
%             fig_cfg.fig_size=[0 0 14 16];
%             fig_h = figure('name', 'EOF','PaperUnits','inches', ...
%                     'PaperPosition',(fig_cfg.fig_size),'position', fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
%             
%             tmp.varn='svstr';
%             amod=1;
%             subplot(2,2,1); 
%             lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
%             m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
%             m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
%             % m_grid('fontsize', m_grid_fontsize-5, 'tickdir', m_grid_tickdir_type, 'xticklabels', [], 'xtick',(120:20:160), 'parent', ax1_1);  
%             m_grid;  
%             m_gshhs_i('color',[1 1 1]);
%             m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
%             title('(a) 1st mode LV of meridional wind stress');
%             mcax=max(abs(lv_l(:)));
%             caxis([-mcax mcax]);
%             colormap(fig_cfg.c_map);
%             set(gca, 'fontsize', fntsize);
%             
%             subplot(2,2,2); 
%             plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'k', 'linewidth', 2);
%             ylabel('wind stress (N/m^2)')
%             yyaxis right
%             plot(cfg.years,part_tc, 'k--', 'linewidth', 1);
%             ylabel('Particle #')
%             title('(b) 1st mode PCT of meridional wind stress');
%             grid on;
%             grid minor;
%             set(gca, 'fontsize', fntsize);
%             xlabel('Year')
%             legend('V wind str', 'Particle #', 'Location', 'NorthWest')
%             
%             
%             tmp.varn='zeta';
%             amod=2;
%             subplot(2,2,3); 
%             lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
%             m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
%             m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
%             m_grid;  
%             m_gshhs_i('color',[1 1 1]);
%             m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
% %             title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
%             title('(c) 2nd mode LV of sea level');
%             mcax=max(abs(lv_l(:)));
%             caxis([-mcax mcax]);
%             colormap(fig_cfg.c_map);
%             set(gca, 'fontsize', fntsize);
%             
%             subplot(2,2,4); 
%             plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'k', 'linewidth', 2);
%             ylabel('zeta (m)')
%             yyaxis right
%             plot(cfg.years,part_tc, 'k--', 'linewidth', 1);
%             ylabel('Particle #')
%             title('(d) 2nd mode PCT of sea level');
%             grid on;
%             grid minor;
%             set(gca, 'fontsize', fntsize);
%             xlabel('Year')
%             legend('V wind str', 'Particle #', 'Location', 'NorthWest')
%             
% 
%             dir.figtgdir=[dir.figdir, '/', 'EOF_mov', '/', 'par', num2str(pmon)];
%             mkdir(dir.figtgdir);
%             
%             cfg.figname=[dir.figtgdir, '/', 'Fin_figure_EOF_mov_', '_season_',num2str(si), '_', str_d, '_p_', num2str(pmon) '.tif'];
%             print(fig_h, cfg.figname, '-dpng');
%             RemoveWhiteSpace([], 'file', cfg.figname);
% 
% 
% 
%             corrcoef(EOF_result_mov.svstr.pcs(si,:,1), EOF_result_mov.zeta.pcs(si,:,2))
% 
% 
% 
% %% key figure draft (different aspect)
% 
%             fntsize=15;
%             close all;
%             fig_cfg.fig_size=[0 0 16 16];
%             fig_h = figure('name', 'EOF','PaperUnits','inches', ...
%                     'PaperPosition',(fig_cfg.fig_size),'position', fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
%             
%             tmp.varn='svstr';
%             amod=1;
%             sbp1=subplot(2,2,1);
%             lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
%             m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
%             m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
%             % m_grid('fontsize', m_grid_fontsize-5, 'tickdir', m_grid_tickdir_type, 'xticklabels', [], 'xtick',(120:20:160), 'parent', ax1_1);  
%             m_grid;  
%             m_gshhs_i('color',[1 1 1]);
%             m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
%             title('(a) 1st mode LV of meridional wind stress');
%             mcax=max(abs(lv_l(:)));
%             caxis([-mcax mcax]);
%             colormap(fig_cfg.c_map);
%             set(gca, 'fontsize', fntsize);
%             pos = get(sbp1, 'Position');  % Get current position
%             pos(1)=pos(1)+pos(3)/20;
%             set(sbp1, 'Position', pos);
%             
%             sbp3=subplot(2,2,3);
%             plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'k', 'linewidth', 2);
%             ylabel('wind stress (N/m^2)')
%             yyaxis right
%             plot(cfg.years,part_tc, 'k--', 'linewidth', 1);
%             ylabel('Particle #')
%             title('(c) 1st mode PCT of meridional wind stress');
%             grid on;
%             grid minor;
%             set(gca, 'fontsize', fntsize);
%             xlabel('Year')
%             legend('V wind str', 'Particle #', 'Location', 'NorthWest')
%             pos = get(sbp3, 'Position');  % Get current position
%             pos(4)=pos(4)/2;
%             pos(2)=pos(2)+pos(4)*4/3;
%             set(sbp3, 'Position', pos);
% 
%             tmp.varn='zeta';
%             amod=2;
%             sbp2=subplot(2,2,2);
%             lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
%             m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
%             m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
%             m_grid;  
%             m_gshhs_i('color',[1 1 1]);
%             m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
% %             title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
%             title('(b) 2nd mode LV of sea level');
%             mcax=max(abs(lv_l(:)));
%             caxis([-mcax mcax]);
%             colormap(fig_cfg.c_map);
%             set(gca, 'fontsize', fntsize);
%             pos = get(sbp2, 'Position');  % Get current position
%             pos(1)=pos(1)+pos(3)/20;
%             set(sbp2, 'Position', pos);
%             
%             sbp4=subplot(2,2,4);
%             plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'k', 'linewidth', 2);
%             ylabel('zeta (m)')
%             yyaxis right
%             plot(cfg.years,part_tc, 'k--', 'linewidth', 1);
%             ylabel('Particle #')
%             title('(d) 2nd mode PCT of sea level');
%             grid on;
%             grid minor;
%             set(gca, 'fontsize', fntsize);
%             xlabel('Year')
%             legend('V wind str', 'Particle #', 'Location', 'NorthWest')
%             pos = get(sbp4, 'Position');  % Get current position
%             pos(4)=pos(4)/2;
%             pos(2)=pos(2)+pos(4)*4/3;
%             set(sbp4, 'Position', pos);
%             
% 
%             dir.figtgdir=[dir.figdir, '/', 'EOF_mov', '/', 'par', num2str(pmon)];
%             mkdir(dir.figtgdir);
%             
%             cfg.figname=[dir.figtgdir, '/', 'Fin_figure2_EOF_mov_', '_season_',num2str(si), '_', str_d, '_p_', num2str(pmon) '.tif'];
%             print(fig_h, cfg.figname, '-dpng');
%             RemoveWhiteSpace([], 'file', cfg.figname);
% 
% 
% 
%             corrcoef(EOF_result_mov.svstr.pcs(si,:,1), EOF_result_mov.zeta.pcs(si,:,2))



%% key figure fin (different aspect, stlee particle. Figure 4)

            fntsize=15;
            close all;
            fig_cfg.fig_size=[0 0 16 16];
            fig_h = figure('name', 'EOF','PaperUnits','inches', ...
                    'PaperPosition',(fig_cfg.fig_size),'position', fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
            
            tmp.varn='svstr';
            amod=1;
            sbp1=subplot(2,2,1); % Fig 4a
            lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
            m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
            % m_grid('fontsize', m_grid_fontsize-5, 'tickdir', m_grid_tickdir_type, 'xticklabels', [], 'xtick',(120:20:160), 'parent', ax1_1);  
            m_grid;  
            m_gshhs_i('color',[1 1 1]);
            m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
            title('(a) 1st mode LV of meridional wind stress');
            mcax=max(abs(lv_l(:)));
            caxis([-mcax mcax]);
            colormap(fig_cfg.c_map);
            set(gca, 'fontsize', fntsize);
            pos = get(sbp1, 'Position');  % Get current position
            pos(1)=pos(1)+pos(3)/20;
            set(sbp1, 'Position', pos);
            
            sbp3=subplot(2,2,3); % Fig 4c
            plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'k', 'linewidth', 2);
            ylabel('wind stress (N/m^2)')
            yyaxis right
%             plot(cfg.years,part_tc, 'k--', 'linewidth', 1);
            plot(cfg.years,part3_mov_stlee.num_A_t(6,:), 'k--', 'linewidth', 1);
%             plot(cfg.years,part3_mov_stlee.num_A_t(:,6), 'k--', 'linewidth', 1);
            abcd=corrcoef(part3_mov_stlee.num_A_t(6,:), ...
                EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
            disp(['4c corr:', num2str(abcd(1,2))])
            ylabel('Particle #')
            title('(c) 1st mode PCT of meridional wind stress');
            grid on;
            grid minor;
            set(gca, 'fontsize', fntsize);
            xlabel('Year')
            legend('V wind str', 'Particle #', 'Location', 'NorthWest')
            pos = get(sbp3, 'Position');  % Get current position
            pos(4)=pos(4)/2;
            pos(2)=pos(2)+pos(4)*4/3;
            set(sbp3, 'Position', pos);

            tmp.varn='zeta';
            amod=2;
            sbp2=subplot(2,2,2); % Fig 4b
            lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
            m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
            m_grid;  
            m_gshhs_i('color',[1 1 1]);
            m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
%             title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
            title('(b) 2nd mode LV of sea level');
            mcax=max(abs(lv_l(:)));
            caxis([-mcax mcax]);
            colormap(fig_cfg.c_map);
            set(gca, 'fontsize', fntsize);
            pos = get(sbp2, 'Position');  % Get current position
            pos(1)=pos(1)+pos(3)/20;
            set(sbp2, 'Position', pos);
            
            sbp4=subplot(2,2,4); % Fig 4d
            plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'k', 'linewidth', 2);
            ylabel('SSH (m)')
            yyaxis right
%             plot(cfg.years,part_tc, 'k--', 'linewidth', 1);
            plot(cfg.years,part3_mov_stlee.num_A_t(6,:), 'k--', 'linewidth', 1);
%             plot(cfg.years,part3_mov_stlee.num_A_t(:,6), 'k--', 'linewidth', 1);
            
            abcd=corrcoef(part3_mov_stlee.num_A_t(6,:), ...
                EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
            disp(['4d corr:', num2str(abcd(1,2))])

            ylabel('Particle #')
            title('(d) 2nd mode PCT of sea level');
            grid on;
            grid minor;
            set(gca, 'fontsize', fntsize);
            xlabel('Year')
            legend('SSH', 'Particle #', 'Location', 'NorthWest')
            pos = get(sbp4, 'Position');  % Get current position
            pos(4)=pos(4)/2;
            pos(2)=pos(2)+pos(4)*4/3;
            set(sbp4, 'Position', pos);
            

            dir.figtgdir=[dir.figdir, '/', 'EOF_mov', '/', 'par', num2str(pmon)];
            mkdir(dir.figtgdir);
            
            cfg.figname=[dir.figtgdir, '/', 'Fin_figure2_stleepar_EOF_mov_', '_season_',num2str(si), '_', str_d, '_p_', num2str(pmon) '.tif'];
            print(fig_h, cfg.figname, '-dpng');
            RemoveWhiteSpace([], 'file', cfg.figname);








% %% key figure (supplementary) template
%             close all;
%             fig_cfg.fig_size=[0 0 14 16];
%             fig_h = figure('name', 'EOF','PaperUnits','inches', ...
%                     'PaperPosition',(fig_cfg.fig_size),'position', fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
%             
%             tmp.varn='svstr';
%             amod=2;
%             subplot(2,2,1);
%             lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
%             m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
%             m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
%             % m_grid('fontsize', m_grid_fontsize-5, 'tickdir', m_grid_tickdir_type, 'xticklabels', [], 'xtick',(120:20:160), 'parent', ax1_1);  
%             m_grid;  
%             m_gshhs_i('color',[1 1 1]);
%             m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
%             title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
%             mcax=max(abs(lv_l(:)));
%             caxis([-mcax mcax]);
%             colormap(fig_cfg.c_map);
%             
%             subplot(2,2,2);
%             plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod));
%             tmp.corr1=corrcoef(part_tc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
%             tmp.corr2=corrcoef(part_wc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
%             tmp.corr3=corrcoef(part_ec, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
%             title([tmp.varn, ' pcs, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%, ', ...
%                 'r(tc):', num2str(round(tmp.corr1(1,2),2)), ...
%                 ', r(wc):', num2str(round(tmp.corr2(1,2),2)), ...
%                 ', r(ec):', num2str(round(tmp.corr3(1,2),2))]);
%             grid on;
%             
%             
%             tmp.varn='zeta';
%             amod=1;
%             subplot(2,2,3);
%             lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
%             m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
%             m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
%             m_grid;  
%             m_gshhs_i('color',[1 1 1]);
%             m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
%             title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
%             mcax=max(abs(lv_l(:)));
%             caxis([-mcax mcax]);
%             colormap(fig_cfg.c_map);
%             
%             subplot(2,2,4);
%             plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod));
%             tmp.corr1=corrcoef(part_tc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
%             tmp.corr2=corrcoef(part_wc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
%             tmp.corr3=corrcoef(part_ec, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
%             title([tmp.varn, ' pcs, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%, ', ...
%                 'r(tc):', num2str(round(tmp.corr1(1,2),2)), ...
%                 ', r(wc):', num2str(round(tmp.corr2(1,2),2)), ...
%                 ', r(ec):', num2str(round(tmp.corr3(1,2),2))]);
%             grid on;
% 
% 
%             dir.figtgdir=[dir.figdir, '/', 'EOF_mov', '/', 'par', num2str(pmon)];
%             mkdir(dir.figtgdir);
%             
%             cfg.figname=[dir.figtgdir, '/', 'Fig_sup_figure_EOF_mov_', '_season_',num2str(si), '_', str_d, '_p_', num2str(pmon) '.tif'];
%             print(fig_h, cfg.figname, '-dpng');
%             RemoveWhiteSpace([], 'file', cfg.figname);
% 
% 
% %% key figure (supplementary) template
%             close all;
%             fig_cfg.fig_size=[0 0 14 16];
%             fig_h = figure('name', 'EOF','PaperUnits','inches', ...
%                     'PaperPosition',(fig_cfg.fig_size),'position', fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
%             
%             tmp.varn='svstr';
%             amod=2;
%             subplot(2,2,1);
%             lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
%             m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
%             m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
%             % m_grid('fontsize', m_grid_fontsize-5, 'tickdir', m_grid_tickdir_type, 'xticklabels', [], 'xtick',(120:20:160), 'parent', ax1_1);  
%             m_grid;  
%             m_gshhs_i('color',[1 1 1]);
%             m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
%             title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
%             mcax=max(abs(lv_l(:)));
%             caxis([-mcax mcax]);
%             colormap(fig_cfg.c_map);
%             
%             subplot(2,2,2);
%             plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod));
%             tmp.corr1=corrcoef(abc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
%             tmp.corr2=corrcoef(abc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
%             tmp.corr3=corrcoef(abc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
%             title([tmp.varn, ' pcs, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%, ', ...
%                 'r(tc):', num2str(round(tmp.corr1(1,2),2)), ...
%                 ', r(wc):', num2str(round(tmp.corr2(1,2),2)), ...
%                 ', r(ec):', num2str(round(tmp.corr3(1,2),2))]);
%             grid on;
%             
%             
%             tmp.varn='zeta';
%             amod=1;
%             subplot(2,2,3);
%             lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
%             m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
%             m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
%             m_grid;  
%             m_gshhs_i('color',[1 1 1]);
%             m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
%             title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
%             mcax=max(abs(lv_l(:)));
%             caxis([-mcax mcax]);
%             colormap(fig_cfg.c_map);
%             
%             subplot(2,2,4);
%             plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod));
%             tmp.corr1=corrcoef(abc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
%             tmp.corr2=corrcoef(abc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
%             tmp.corr3=corrcoef(abc, EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
%             title([tmp.varn, ' pcs, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%, ', ...
%                 'r(tc):', num2str(round(tmp.corr1(1,2),2)), ...
%                 ', r(wc):', num2str(round(tmp.corr2(1,2),2)), ...
%                 ', r(ec):', num2str(round(tmp.corr3(1,2),2))]);
%             grid on;
% 
% 
%             dir.figtgdir=[dir.figdir, '/', 'EOF_mov', '/', 'par', num2str(pmon)];
%             mkdir(dir.figtgdir);
%             
%             cfg.figname=[dir.figtgdir, '/', 'Fig_sup_figure_stleepar_EOF_mov_', '_season_',num2str(si), '_', str_d, '_p_', num2str(pmon) '.tif'];
%             print(fig_h, cfg.figname, '-dpng');
%             RemoveWhiteSpace([], 'file', cfg.figname);



%% key figure fin (supplementary; wind 2mode, ssh 1mode. Supplementary Figure 3)
            fntsize=15;
            close all;
            fig_cfg.fig_size=[0 0 6 8];
            fig_h = figure('name', 'EOF','PaperUnits','inches', ...
                    'PaperPosition',(fig_cfg.fig_size),'position', fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
                               
            tmp.varn='zeta';
            amod=1;
            subplot(2,1,1); % FigS 3a
            lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
            m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
            m_grid;  
            m_gshhs_i('color',[1 1 1]);
            m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
%             title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
            title('(a) 1st mode LV of sea level');
            mcax=max(abs(lv_l(:)));
            caxis([-mcax mcax]);
            colormap(fig_cfg.c_map);
            set(gca, 'fontsize', fntsize-5);
            
            sbp2=subplot(2,1,2); % FigS 3b
            plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'k', 'linewidth', 2);
            ylim([-0.4 0.6])
            ylabel('SSH (m)')
            yyaxis right
            abc=part3_mov_stlee.num_A_t(6,:);
            plot(cfg.years,abc, 'k--', 'linewidth', 1);

            abcd=corrcoef(part3_mov_stlee.num_A_t(6,:), ...
                EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
            disp(['supp corr:', num2str(abcd(1,2))]);

            ylabel('Particle #')
            title('(b) 1st mode PCT of sea level');
            grid on;
            grid minor;
            set(gca, 'fontsize', fntsize-5);
            xlabel('Year')
            legend('SSH', 'Particle #', 'Location', 'NorthWest')
            pos = get(sbp2, 'Position');  % Get current position
            pos(4)=pos(4)/2;
            pos(2)=pos(2)+pos(4)*4/3;
            set(sbp2, 'Position', pos);

            dir.figtgdir=[dir.figdir, '/', 'EOF_mov', '/', 'par', num2str(pmon)];
            mkdir(dir.figtgdir);
            
            cfg.figname=[dir.figtgdir, '/', 'Fin_sup_figure_stleepar_EOF_mov_', '_season_',num2str(si), '_', str_d, '_p_', num2str(pmon) '.tif'];
            print(fig_h, cfg.figname, '-dpng');
%             RemoveWhiteSpace([], 'file', cfg.figname);


%% key figure fin (supplementary; wind 2mode,)
            fntsize=15;
            close all;
            fig_cfg.fig_size=[0 0 6 8];
            fig_h = figure('name', 'EOF','PaperUnits','inches', ...
                    'PaperPosition',(fig_cfg.fig_size),'position', fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
                               
            tmp.varn='svstr';
            amod=2;
            subplot(2,1,1); % FigS 3a
            lv_l=squeeze(EOF_result_mov.(tmp.varn).lv(si,:,:,amod));
            m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            m_pcolor(lon',lat',lv_l'); shading flat; colorbar;
            m_grid;  
            m_gshhs_i('color',[1 1 1]);
            m_gshhs_i('patch',[0.7 0.7 0.7]);   % gray colored land
%             title([tmp.varn, ' lv, ', num2str(round(EOF_result_mov.(tmp.varn).var_exp(si,amod),2)), '%']);
            title('(a) 2nd mode LV of meridional wind stress');
            mcax=max(abs(lv_l(:)));
            caxis([-mcax mcax]);
            colormap(fig_cfg.c_map);
            set(gca, 'fontsize', fntsize-5);
            
            sbp2=subplot(2,1,2); % FigS 3b
            plot(cfg.years,EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'k', 'linewidth', 2);
%             ylim([-0.4 0.6])
            ylabel('wind stress (N/m^2)')
            yyaxis right
            abc=part3_mov_stlee.num_A_t(6,:);
            plot(cfg.years,abc, 'k--', 'linewidth', 1);

            abcd=corrcoef(part3_mov_stlee.num_A_t(6,:), ...
                EOF_result_mov.(tmp.varn).pcs(si,:,amod), 'Rows', 'complete');
            disp(['supp corr:', num2str(abcd(1,2))]);

            ylabel('Particle #')
            title('(b) 2nd mode PCT of meridional wind stress');
            grid on;
            grid minor;
            set(gca, 'fontsize', fntsize-5);
            xlabel('Year')
            legend('svstr', 'Particle #', 'Location', 'NorthWest')
            pos = get(sbp2, 'Position');  % Get current position
            pos(4)=pos(4)/2;
            pos(2)=pos(2)+pos(4)*4/3;
            set(sbp2, 'Position', pos);

            dir.figtgdir=[dir.figdir, '/', 'EOF_mov', '/', 'par', num2str(pmon)];
            mkdir(dir.figtgdir);
            
            cfg.figname=[dir.figtgdir, '/', 'Fin_sup_figure_svstr2_stleepar_EOF_mov_', '_season_',num2str(si), '_', str_d, '_p_', num2str(pmon) '.tif'];
            print(fig_h, cfg.figname, '-dpng');
%             RemoveWhiteSpace([], 'file', cfg.figname);

save('/Volumes/kyy_raid/kimyy/Research/Postdoc/03_IBS/2024_TS_KS_connectivity/manuscript/TS_KS_Fig4_Figs3.mat');