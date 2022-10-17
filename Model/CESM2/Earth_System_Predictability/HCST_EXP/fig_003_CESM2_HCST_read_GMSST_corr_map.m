% %  Created 05-Oct-2022 by Yong-Yub Kim
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
% dirs.root='/mnt/lustre/proj/earth.system.predictability/HCST_EXP';
% dirs.yoshi_root='/proj/yoshi/DATA/CESM2_ODA';
% dirs.archive=[dirs.root, filesep, 'archive'];
% dirs.saveroot='/mnt/lustre/proj/kimyy/Model/CESM2/ESP/HCST_EXP';
dirs.hcstroot='/Volumes/kyy_raid/kimyy/Model/CESM2/ESP/HCST_EXP';
dirs.figroot='/Volumes/kyy_raid/kimyy/Figure/CESM2/ESP/HCST_EXP';

config.iyears=1968:2021;
config.months=1:12;
config.scenname='HIST';
config.gridname='f09_g17';
config.assm_factor='10';
config.ens_member='1';
config.proj_year=10;

config.obsnames={'en4.2_ba'};
config.ensnames={'ba-10p1'};

config.component='ocn';
% config.varnames={'TEMP','SALT'};
config.varnames={'temp', 'salt'};
config.len_t_y = length(config.iyears);
config.len_t_m = length(config.months);
config.len_t = config.len_t_y *config.len_t_m;
config.len_obs= length(config.obsnames);
config.len_ens= length(config.ensnames);




%% grid set(mask from model)
tmp.obsname=config.obsnames{1};
iyear=min(config.iyears);
config.casename_m=[config.gridname, '.hcst.', tmp.obsname, '-', config.assm_factor, 'p', config.ens_member];
config.casename=[config.casename_m, '_i', num2str(iyear)];
dirs.datadir= [dirs.hcstroot, filesep, config.casename_m, filesep, 'GMSV'];

% f09_g17.hcst.en4.2_ba-10p1_i2021.pop.h.once.nc

% files = dir('/Volumes/kyy_raid/kimyy/bongan/example/data/*.nc');
% for i_model = 1:length(files)
% file_parts = split(files(i_model).name,{'_','.'});
% Table = file_parts{2};
% Model = file_parts{3};
% Grd   = file_parts{6};
% Time  = split(file_parts{end-1},'-');
% Time  = Time{1};
% file_in = fullfile(files(i_model).folder,files(i_model).name);
% end

% finfo   = ncinfo( fname );
% kk = find(ismember({finfo.Variables.Name},varnames),1,'last');

[tmp.error_status, tmp.value]=system(['ls ', dirs.datadir, '/*once*']);  % b.e21.BHISTsmbb.f09_g17.assm.oras4_ba-10p1.pop.h.once.nc
tmp.gridname = [tmp.value(1:end-1)];
grid.region_mask=ncread(tmp.gridname, 'REGION_MASK'); 
grid.ocean_mask=NaN(size(grid.region_mask));
grid.ocean_mask(grid.region_mask>0)=1;
grid.tarea = ncread(tmp.gridname, 'TAREA');

grid.tlong=ncread(tmp.gridname, 'TLONG');
grid.tlat=ncread(tmp.gridname, 'TLAT');
grid.nlon=size(grid.tlong,1);
grid.nlat=size(grid.tlong,2);
grid.ntime=config.proj_year.*12;


% % model filename example
% /mnt/lustre/proj/earth.system.predictability/HCST_EXP/archive/
% f09_g17.hcst.en4.2_ba-10p1/f09_g17.hcst.en4.2_ba-10p1_i1993/ocn/hist/
% f09_g17.hcst.en4.2_ba-10p1_i1993.pop.h.1993-01.nc

S = shaperead('landareas.shp');


%% read & plot data
for obsind=1:length(config.obsnames)
    tmp.obsname=config.obsnames{obsind};
    tmp.obsname_simple= f_obs_simple(tmp.obsname);
    for varind=1:length(config.varnames)
        tmp.varname=config.varnames{varind};
        tmp.varname_C= f_varname_C(tmp.varname);
        tmp.var_unit=f_varname_unit(tmp.varname);
        for lyear=0:config.proj_year-1
            tmp.lyear_str=num2str(lyear, '%02i');
            config.casename_m=[config.gridname, '.hcst.', tmp.obsname, '-', config.assm_factor, 'p', config.ens_member];
            dirs.datadir= [dirs.hcstroot, filesep, config.casename_m, filesep, 'GMSV'];
            fprintf('%02d_%s_%s  ',lyear,config.casename_m, tmp.varname); lap_time = tic;

            
%             grid.(['time_i', tmp.iyear_str])=NaN(1,grid.ntime);
            %% variables initialization
            data.([tmp.varname, '_model_', tmp.obsname_simple, '_l', tmp.lyear_str])=NaN(grid.nlon, grid.nlat, config.len_t_y+config.proj_year-1);
            data.([tmp.varname, '_bias_', tmp.obsname_simple, '_l', tmp.lyear_str])=NaN(grid.nlon, grid.nlat, config.len_t_y+config.proj_year-1);
            data.([tmp.varname, '_obs_', tmp.obsname_simple])=NaN(grid.nlon, grid.nlat, config.len_t_y);
%             data.([tmp.varname, '_sq_err_', tmp.obsname_simple, '_l', tmp.iyear_str])=NaN(grid.nlon, grid.nlat, config.len_t);
%             data.([tmp.varname, '_rmse_', tmp.obsname_simple, '_l', tmp.iyear_str])=NaN(grid.nlon, grid.nlat);

%             data.(['gm_', tmp.varname, '_model_', tmp.obsname_simple, '_i', tmp.iyear_str])=ncread(config.datafilename, ['gm_', tmp.varname]);
%             data.(['gm_', tmp.varname, '_bias_', tmp.obsname_simple, '_i', tmp.iyear_str])=ncread(config.datafilename, ['gm_bias_', tmp.varname]);
%             data.(['gm_', tmp.varname, '_obs_', tmp.obsname_simple, '_i', tmp.iyear_str])=ncread(config.datafilename, ['gm_obs_', tmp.varname]);
%             data.(['gm_', tmp.varname, '_sq_err_', tmp.obsname_simple, '_i', tmp.iyear_str])=ncread(config.datafilename, ['gm_sq_err_', tmp.varname]);

%             tmp.corrcoef= ...
%                 corrcoef(data.(['gm_', tmp.varname, '_model_', tmp.obsname_simple, '_i', tmp.iyear_str]), ...
%                 data.(['gm_', tmp.varname, '_obs_', tmp.obsname_simple, '_i', tmp.iyear_str]));
%             data.(['corr_gm_', tmp.varname])(iyear-min(config.iyears)+1)= tmp.corrcoef(1,2);
%             data.(['gm_', tmp.varname, '_rmse_', tmp.obsname_simple, '_i', tmp.iyear_str])= ncread(config.datafilename, ['gm_rmse_', tmp.varname]);

            for iyear=min(config.iyears):max(config.iyears)
                tmp.iyear_str=num2str(iyear, '%04i');
                config.casename=[config.casename_m, '_i', tmp.iyear_str];
                config.datafilename=[dirs.datadir, filesep, 'GMSV_', config.casename, '.nc'];
                data.time=ncread(config.datafilename, 'time');
                tmp.ymean= mean(ncread(config.datafilename, [tmp.varname_C], [1, 1, (lyear)*12+1], [inf, inf, 12]),3);
                data.([tmp.varname, '_model_', tmp.obsname_simple, '_l', tmp.lyear_str])(:,:,iyear-min(config.iyears)+1+lyear)= tmp.ymean;
%                 data.([tmp.varname, '_model_', tmp.obsname_simple, '_l', tmp.lyear_str])(data.([tmp.varname, '_model_', tmp.obsname_simple, '_l', tmp.lyear_str])==0)=NaN;
                tmp.ymean= mean(ncread(config.datafilename, ['obs_', tmp.varname], [1, 1, 1], [inf, inf, 12]), 3);
                data.([tmp.varname, '_obs_', tmp.obsname_simple])(:,:,iyear-min(config.iyears)+1)= tmp.ymean;
            end
            
            for loni=1:grid.nlon
                for lati=1:grid.nlat
                     if (isnan(data.([tmp.varname, '_obs_', tmp.obsname_simple])(loni,lati,1))~=1)
                         [tmp.corr, tmp.corr_p]=corrcoef(squeeze(data.([tmp.varname, '_model_', tmp.obsname_simple, '_l', tmp.lyear_str])(loni,lati,1+lyear:end-config.proj_year+1)), ...
                             squeeze(data.([tmp.varname, '_obs_', tmp.obsname_simple])(loni,lati,1+lyear:end)));
                         data.([tmp.varname, '_corr_', tmp.obsname_simple, '_l', tmp.lyear_str])(loni,lati)=tmp.corr(1,2);
                         data.([tmp.varname, '_corr_p_', tmp.obsname_simple, '_l', tmp.lyear_str])(loni,lati)=tmp.corr_p(1,2);
                     else
                         data.([tmp.varname, '_corr_', tmp.obsname_simple, '_l', tmp.lyear_str])(loni,lati)=NaN;
                         data.([tmp.varname, '_corr_p_', tmp.obsname_simple, '_l', tmp.lyear_str])(loni,lati)=NaN;                         
                     end
                end
            end
            disp('abc')


%%example of bongan
% name_rgn = 'Glob';  map_proj = 'robinson';  x_lim = [0 360];    y_lim = [-80 89];   fig_size = [0,0,6,2.8]; ax_size = [-0.15,-0.2,5.4,3];   cb_size = [5.15,0.20,0.15,2.3]; title_pos = [0.5,1.0];
% name_rgn = 'Glob';  map_proj = 'eqdcylin';  x_lim = [0 360];    y_lim = [-80 89];   fig_size = [0,0,6,2.8]; ax_size = [-0.15,0.15,5.4,2.7];	cb_size = [5.15,0.20,0.15,2.3]; title_pos = [0.5,0.93];

            fig_config.name_rgn = 'Glob';
            fig_config.map_proj = 'eqdcylin';  % robinson, eqdcylin
            fig_config.x_lim = [-180 180];
            fig_config.y_lim = [-80 89];
            fig_config.fig_size = [0,0,6,3.5];
            fig_config.ax_size = [0.3,0.7,5.4,2.7];
            fig_config.cb_size = [5.15,0.8,0.15,2.3];
            fig_config.title_pos = [0.5,0.93];
            fig_config.p_lim =0.1;
            fig_config.c_lim = [-1 1];
            [fig_config.c_map, tmp.err_stat] = Func_0009_get_colormaps('byr', tmp.dropboxpath);
            
            
            tmp.X=grid.tlong([end, 1:end],:);
%             tmp.X(tmp.X>180)= tmp.X(tmp.X>180)-360;
            tmp.Y=grid.tlat([end, 1:end],:);
            tmp.C=data.([tmp.varname, '_corr_', tmp.obsname_simple, '_l', tmp.lyear_str]);
            tmp.C=tmp.C([end, 1:end],:);
            tmp.C_H=data.([tmp.varname, '_corr_p_', tmp.obsname_simple, '_l', tmp.lyear_str]);
            tmp.C_H=tmp.C_H([end, 1:end],:);
            tmp.C_2=tmp.C;
            tmp.C_2(tmp.C_H<fig_config.p_lim)=NaN;
        
            fig_config.fig_name=['l',tmp.lyear_str, ', corr, ', tmp. varname];
            fig_h = figure('name',fig_config.fig_name,'PaperUnits','inches', ...
                'PaperPosition',fig_config.fig_size,'position',fig_config.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','off');
            %% map setting
            ax_m = axesm('MapProjection',fig_config.map_proj,'grid','on','fontsize',14, ...
                'fontname','freeserif'); 

            axis off; 
            hold on;
%             setm(ax_m, 'MapLatLimit',fig_config.y_lim);
%             setm(ax_m, 'MapLatLimit',fig_config.y_lim,'MapLonLimit',fig_config.x_lim);
%             setm(ax_m,'origin',[0,mean(fig_config.x_lim)],'MapLatLimit',fig_config.y_lim,'MapLonLimit',fig_config.x_lim);
            setm(ax_m,'origin',[0,180],'MapLatLimit',fig_config.y_lim);

%             setm(ax_m, 'MeridianLabel','on','ParallelLabel','on','MLineLocation',5,'PLineLocation',5,'MLabelLocation',5,'PLabelLocation',5)
%             setm(ax_m,'ParallelLabel','on','MeridianLabel','on','PLabelMeridian','west','MLabelParallel','south');
%             setm(ax_m,'ParallelLabel','on','MeridianLabel','on','PLabelMeridian','west','MLabelParallel','south', ...
%                 'PLineLocation',10, 'PLabelLocation',20, 'MLineLocation',20, 'MLabelLocation',60);
            

    
            set(ax_m,'Units','inches','Position',fig_config.ax_size);
            text(ax_m,fig_config.title_pos(1),fig_config.title_pos(2),fig_config.fig_name, ...
                'units','normalized','horizontalalignment','center','verticalalignment','middle','fontsize',14,'fontname','freeserif','interpreter','none')
            %% caxis & colorbar
            caxis(ax_m, fig_config.c_lim); 
            colormap(fig_config.c_map);
            cb = colorbar(ax_m,'units','inches','position',fig_config.cb_size);
            set(cb,'fontsize',12,'fontname','freeserif','TickDir','both');
%             ylabel(cb,[replace(tmp.varname,'_',' '), '(',tmp.var_unit,')']);
            title(cb,'R','fontsize',12);
            %% draw on ax_m
            
            
            h_pc = pcolorm(tmp.Y,tmp.X,tmp.C,'parent',ax_m); 
            shading flat;
            geoshow(ax_m,[S.Y],[S.X],'color','k','linewidth',0.5);

            pp2 = pcolorm(tmp.Y,tmp.X,tmp.C_2, 'parent', ax_m);
            set(pp2,'linestyle','none','Tag','HatchingRegion');
            hp = findobj(pp2,'Tag','HatchingRegion');
            hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','k','HatchLineWidth',0.5);

            %% frame and label setting
            setm(ax_m,'frame','on','FLineWidth',1);

            label_y=plabel('PlabelMeridian', 'west', 'PLineLocation',10, 'PLabelLocation',20, 'labelrotation','on');
            label_x=mlabel('MLabelParallel','south', 'MLineLocation',20, 'MLabelLocation',60, 'labelrotation','on');
            mlabel; plabel;
            label_y=plabel; label_x=mlabel;
            for lxi=1:length(label_x)
                tmp.tmppos=label_x(lxi,1).Position;
                tmp.tmppos(2)=-fig_config.ax_size(4)+1.55;
                label_x(lxi,1).Position=tmp.tmppos;
                label_x(lxi,1).String{2}=replace(label_x(lxi,1).String{2}, ' ','');
            end
            for lyi=1:length(label_y)
                label_y(lyi,1).String=replace(label_y(lyi,1).String, ' ','');
            end

            %% save
            dirs.figdir= [dirs.figroot, filesep, config.casename_m, filesep, 'GMSV', filesep, 'Corr_map'];
            if ~exist(dirs.figdir,'dir'), mkdir(dirs.figdir); end
            config.figname=[dirs.figdir, filesep, 'corr_map_', tmp.varname, '_l', tmp.lyear_str, '.tif'];
            print(fig_h, config.figname, '-dpng');
            RemoveWhiteSpace([], 'file', config.figname);
            fprintf('%7.1f sec\n', toc(lap_time) );

%             pcolor(data.([tmp.varname, '_corr_', tmp.obsname_simple, '_l', tmp.lyear_str])'); shading flat; colorbar;
            
%             set(gcf,'Renderer','opengl')
%             set(gcf,'Renderer','zbuffer')
%             set(gcf,'Renderer','painters')


%% sort by longitude - fail
%             var_src=data.([tmp.varname, '_corr_', tmp.obsname_simple, '_l', tmp.lyear_str]);
%             X = grid.tlong;
%             Y = NaN(size(grid.tlong)); C = Y;
%             X( grid.tlong < 0 ) = X( grid.tlong < 0 ) + 360;
%             [X, I] = sort( X );
%             for mm = 1:size(grid.tlat,2)
%                 Y(:,mm) = grid.tlat( I(:,mm), mm );
%                 C(:,mm) = squeeze( var_src( I(:,mm), mm ) );
%             end 
%             pcolor(X',Y',C')


            %% hatch
%             pp2 = pcolor(ref_lon,ref_lat,hatch_tr1');% shading flat
%             set(pp2,'linestyle','none','Tag','HatchingRegion');
%             hp = findobj(pp2,'Tag','HatchingRegion');
%             hh = hatchfill2(hp,'hatchstyle','single','HatchAngle',45,'HatchDensity',150,'HatchColor','k','HatchLineWidth',0.5);
            
%             m_pcolor(double(cut_lon_rho)',cut_lat_rho',squeeze(cmems_trend_yearly(:,:)'));
%             shading(gca, 'flat');
%             m_gshhs_c('color',[0.8 0.8 0.8]);
%             m_gshhs_c('patch',[0.8 0.8 0.8]);   % gray colored land
%             m_grid('fontsize', 20, 'box', 'fancy', 'tickdir', m_grid_tickdir_type);
%             hold off

        end
%% save ncfile

%     [f_exp, gof_exp] = fit(trendtime_yearly',squeeze(cmems_sla_yearly(i,j,:)),'poly1');
%                             cmems_sla_yearly_poly1_fit(i,j,1:length(trendtime_yearly))=f_exp(trendtime_yearly);

    end
    
    data.time_y=config.iyears;
    data.time_y_extended=[config.iyears, max(config.iyears)+1:max(config.iyears)+config.proj_year-1];

end

function obsname_simple = f_obs_simple(obsname)
    switch obsname
        case 'en4.2_ba'
            obsname_simple='en4';
        case 'projdv7.3'
            obsname_simple='projd';
    end
end

function gmsst = f_gm_var(var_2d, area)
    mask_var=NaN(size(var_2d));
    mask_var(isfinite(var_2d))=1;
    area=area.*mask_var;
    var_2d=squeeze(var_2d);
    var_2d_sum=var_2d.*area;
    var_2d_sum=sum(var_2d_sum(:), 'omitnan');
    gmsst=var_2d_sum ./ sum(area(:), 'omitnan');
end

function varname_C = f_varname_C(varname)
    switch varname
        case 'temp'
            varname_C='TEMP';
        case 'salt'
            varname_C='SALT';
    end
end

function varname_unit = f_varname_unit(varname)
    switch varname
        case 'temp'
            varname_unit='\circC';
        case 'salt'
            varname_unit='g/kg';
    end
end


% %                 % % %         cmems poly1 fitting yearly
% %                 tic;
% %                 for i=1:cmems_lonsize_cut
% %                     for j=1:cmems_latsize_cut
% %                         if isfinite(cmems_sla_yearly(i,j,1))
% %                             [f_exp, gof_exp] = fit(trendtime_yearly',squeeze(cmems_sla_yearly(i,j,:)),'poly1');
% %                             cmems_sla_yearly_poly1_fit(i,j,1:length(trendtime_yearly))=f_exp(trendtime_yearly);
% %                             cmems_sla_yearly_poly1_fit_rsquare(i,j)=gof_exp.rsquare;
% %                             cmems_sla_yearly_poly1_fit_rmse(i,j)=gof_exp.rmse;
% %                         else
% %                             cmems_sla_yearly_poly1_fit(i,j,1:length(trendtime_yearly))=NaN;
% %                             cmems_sla_yearly_poly1_fit_rsquare(i,j)=NaN;
% %                             cmems_sla_yearly_poly1_fit_rmse(i,j)=NaN;
% %                         end
% % %                         p=polyfit(trendtime_yearly,squeeze(cmems_sla_yearly(i,j,:))',1);
% % %                         cmems_trend_yearly(i,j)=p(1) * 1000.0 ;
% %                     end
% %                 end
% %                 disp('cmems poly1 fitting yearly complete') 
% %                 toc;
% %                 
% %                 for i=1:cmems_lonsize_cut
% %                     for j=1:cmems_latsize_cut
% %                         for k=1:size(cmems_sla_yearly,3)
% %                             cmems_sla_yearly_poly1_fit_detrended(i,j,k)=cmems_sla_yearly(i,j,k)-cmems_sla_yearly_poly1_fit(i,j,k);
% %                         end
% %                     end
% %                 end