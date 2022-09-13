disp('subroutine v_2p_2nd_001_sub_003_cal_EKWC') 
%% calculate EKWC


%% get data
RCM_info.savedir_v = [RCM_info.saveroot, RCM_info.testname, tmp.fs, ...
            'v', tmp.fs];
RCM_info.matname_v = [RCM_info.savedir_v,RCM_info.testname,'_',RCM_info.regionname, '_RCM_','v', '_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), '_', ...
    num2str(RCM_info.months(1),'%04i'),'_',num2str(RCM_info.months(end),'%04i'), '.mat'];
load(RCM_info.matname_v, 'RCM_data')
RCM_data_v=RCM_data;

RCM_info.savedir_u = [RCM_info.saveroot, RCM_info.testname, tmp.fs, ...
            'u', tmp.fs];
RCM_info.matname_u = [RCM_info.savedir_u,RCM_info.testname,'_',RCM_info.regionname, '_RCM_', 'u', '_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), '_', ...
    num2str(RCM_info.months(1),'%04i'),'_',num2str(RCM_info.months(end),'%04i'), '.mat'];
load(RCM_info.matname_u, 'RCM_data', 'RCM_grid')
RCM_data_u=RCM_data;

RCM_info_wind=RCM_info;
RCM_grid_wind=RCM_grid;
RCM_info_wind.matname = [RCM_info_wind.windsavedir,RCM_info_wind.testname,'_',RCM_info_wind.regionname, '_RCM_data_wind_', ...
    num2str(min(RCM_info_wind.years),'%04i'),'_',num2str(max(RCM_info_wind.years),'%04i'), '_', ...
    num2str(RCM_info_wind.months(1),'%04i'),'_',num2str(RCM_info_wind.months(end),'%04i'), '.mat'];
load(RCM_info_wind.matname, 'RCM_data_wind')

RCM_info.gos_matname = [RCM_info.atmfiledir, tmp.fs, 'zeta', tmp.fs, RCM_info.testname,'_',RCM_info.regionname, '_RCM_gos_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'_', ...
    num2str(RCM_info.months(1),'%04i'),'_',num2str(RCM_info.months(end),'%04i'), '.mat'];
load(RCM_info.gos_matname, 'RCM_data_gos')

RCM_info.size_t=length(RCM_info.years)*length(RCM_info.months);
RCM_info.size_years=length(RCM_info.years);
RCM_info.size_months=length(RCM_info.months);

%% get transport in the Korea Strait
RCM_EKWC.transport_w = RCM_data_v.all_transport_w;
RCM_EKWC.transport_e = RCM_data_v.all_transport_e;
RCM_EKWC.transport = RCM_data_v.all_transport_w + RCM_data_v.all_transport_e;
RCM_EKWC.yearly_transport_w = RCM_data_v.yearly_mean_transport_w;
RCM_EKWC.yearly_transport_e = RCM_data_v.yearly_mean_transport_e;
RCM_EKWC.yearly_transport = RCM_data_v.yearly_mean_transport_w + RCM_data_v.yearly_mean_transport_e;

%% get maximum velocity of the EKWC
RCM_EKWC.lat_limit_ind=149;  % avoid new WBC north of the EKB
[RCM_EKWC.maxval, RCM_EKWC.maxind] = max(RCM_data_v.maxval(1:RCM_EKWC.lat_limit_ind,:), [], 1, 'omitnan');
[RCM_EKWC.yearly_maxval, RCM_EKWC.yearly_maxind] = max(RCM_data_v.yearly_maxval(1:RCM_EKWC.lat_limit_ind,:), [], 1, 'omitnan');

%% get separation latitude (1/20 of maxv_EKWC)
for ti=1:RCM_info.size_t
    tmp.val = min(find(squeeze(RCM_data_v.maxval(:,ti))<=RCM_EKWC.maxval(ti)/20.0));
    if RCM_data_v.maxval(tmp.val,ti)<0
        tmp.val=tmp.val-1;
    end
    if isempty(tmp.val)
        RCM_EKWC.sep_lat_ind(ti)=RCM_grid.size_lat_rho-2;
    else
        RCM_EKWC.sep_lat_ind(ti)=tmp.val;
    end
    RCM_EKWC.sep_lat(ti)=RCM_grid.cut_lat_rho(1,RCM_EKWC.sep_lat_ind(ti));
end
for ti=1:RCM_info.size_years
    tmp.val = min(find(squeeze(RCM_data_v.yearly_maxval(:,ti))<=RCM_EKWC.yearly_maxval(ti)/10.0));
    if RCM_data_v.yearly_maxval(tmp.val,ti)<0
        tmp.val=tmp.val-1;
    end
    if isempty(tmp.val)
        RCM_EKWC.yearly_sep_lat_ind(ti)= RCM_grid.size_lat_rho-2;
    else
        RCM_EKWC.yearly_sep_lat_ind(ti)=tmp.val;
    end
    RCM_EKWC.yearly_sep_lat(ti)=RCM_grid.cut_lat_rho(1,RCM_EKWC.yearly_sep_lat_ind(ti));
end

%% get spatial mean value of Uwind and Vwind
%% EKWC coastal area only (<130.4E)
RCM_EKWC.lon_limit_ind= 70;
[RCM_EKWC.yearly_v_sp_mean, tmp.error_status] = Func_0011_get_area_weighted_mean(squeeze(RCM_data_v.yearly_mean(1:RCM_EKWC.lon_limit_ind,1:RCM_EKWC.lat_limit_ind,1,:)), ...
    RCM_grid.cut_lon_rho(1:RCM_EKWC.lon_limit_ind,1:RCM_EKWC.lat_limit_ind), RCM_grid.cut_lat_rho(1:RCM_EKWC.lon_limit_ind,1:RCM_EKWC.lat_limit_ind));
[RCM_EKWC.yearly_Vwind_sp_mean, tmp.error_status] = Func_0011_get_area_weighted_mean(RCM_data_wind.yearly_mean_Vwind(1:RCM_EKWC.lon_limit_ind,1:RCM_EKWC.lat_limit_ind,:), ...
    RCM_grid.cut_lon_rho(1:RCM_EKWC.lon_limit_ind,1:RCM_EKWC.lat_limit_ind), RCM_grid.cut_lat_rho(1:RCM_EKWC.lon_limit_ind,1:RCM_EKWC.lat_limit_ind));
[RCM_EKWC.yearly_Uwind_sp_mean, tmp.error_status] = Func_0011_get_area_weighted_mean(RCM_data_wind.yearly_mean_Uwind(1:RCM_EKWC.lon_limit_ind,1:RCM_EKWC.lat_limit_ind,:), ...
    RCM_grid.cut_lon_rho(1:RCM_EKWC.lon_limit_ind,1:RCM_EKWC.lat_limit_ind), RCM_grid.cut_lat_rho(1:RCM_EKWC.lon_limit_ind,1:RCM_EKWC.lat_limit_ind));
RCM_data_wind.all_NWwind = RCM_data_wind.all_Uwind * cosd(45) ...
    - RCM_data_wind.all_Vwind * cosd(45); 
RCM_data_wind.yearly_mean_NWwind = RCM_data_wind.yearly_mean_Uwind * cosd(45) ...
    - RCM_data_wind.yearly_mean_Vwind * cosd(45); 
[RCM_EKWC.yearly_NWwind_sp_mean, tmp.error_status] = Func_0011_get_area_weighted_mean(RCM_data_wind.yearly_mean_NWwind(1:RCM_EKWC.lon_limit_ind,1:RCM_EKWC.lat_limit_ind,:), ...
    RCM_grid.cut_lon_rho(1:RCM_EKWC.lon_limit_ind,1:RCM_EKWC.lat_limit_ind), RCM_grid.cut_lat_rho(1:RCM_EKWC.lon_limit_ind,1:RCM_EKWC.lat_limit_ind));

%% Northwesterly, v correlation
for loni=1:RCM_grid.size_lon_rho
    for lati=1:RCM_grid.size_lat_rho
        if RCM_grid.mask_ocean(loni,lati)==1
            [RCM_EKWC.corr_NWwind_v(loni,lati), RCM_EKWC.p_NWwind_v(loni,lati)]= ...
                corr(squeeze(RCM_data_wind.all_NWwind(loni,lati,:)), squeeze(RCM_data_v.all(loni,lati,1,:)));
%                 corr(squeeze(RCM_EKWC.yearly_NWwind_sp_mean), squeeze(RCM_data_v.yearly_mean(loni,lati,1,:)));
%                 corr(squeeze(RCM_data_wind.yearly_mean_NWwind(loni,lati,:)), squeeze(RCM_data_v.yearly_mean(loni,lati,1,:)));

        else
            RCM_EKWC.corr_NWwind_v(loni,lati)=NaN;
            RCM_EKWC.p_NWwind_v(loni,lati)=NaN;
        end
    end
end
% pcolor(RCM_EKWC.corr_NWwind_v'); shading flat; colorbar; colormap(cmap.bwr_map); caxis([-1 1])
% pcolor(RCM_EKWC.p_NWwind_v'); shading flat; colorbar;colormap(parula);
RCM_EKWC.corr_sig_NWwind_v=RCM_EKWC.corr_NWwind_v;
RCM_EKWC.corr_sig_NWwind_v(RCM_EKWC.p_NWwind_v>0.1)=NaN;

RCM_EKWC.mean_u_rho=mean(RCM_data_u.yearly_mean(:,:,1,:),4);
RCM_EKWC.mean_v_rho=mean(RCM_data_v.yearly_mean(:,:,1,:),4);

run(tmp.param_script)
m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);
m_pcolor(RCM_grid.cut_lon_rho', RCM_grid.cut_lat_rho', RCM_EKWC.corr_sig_NWwind_v'); shading flat; colorbar; colormap(cmap.bwr_map); caxis([-1 1])
m_gshhs_i('color',param.m_gshhs_line_color); m_gshhs_i('patch',param.m_gshhs_land_color);
m_grid('fontsize', param.m_grid_fontsize, 'box', param.m_grid_box_type, 'tickdir', param.m_grid_tickdir_type);
hold on
uvplot=m_quiver(RCM_grid.cut_lon_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)', ...
                RCM_grid.cut_lat_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)', ...
                RCM_EKWC.mean_u_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)' * param.m_quiver_vector_size, ...
                RCM_EKWC.mean_v_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)' * param.m_quiver_vector_size, ...
                'color', param.m_quiver_vector_color, 'AutoScale','off','LineWidth', param.m_quiver_LineWidth);

tmp.figdir = [RCM_info.figroot, tmp.fs, 'corr'];
if (exist(tmp.figdir , 'dir') ~= 7)
    mkdir(tmp.figdir);
end
tmp.tifname = [tmp.figdir, tmp.fs, 'corr_NWwind_v_',RCM_info.testname, '.tif'];
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [param.hor_paper_size_x, param.hor_paper_size_y]);
set(gcf,'PaperPosition', [param.paper_position_hor param.paper_position_ver param.paper_position_width param.paper_position_height]) 
hold off
saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);

%% Northwesterly, geostrophic v correlation
for loni=1:RCM_grid.size_lon_rho
    for lati=1:RCM_grid.size_lat_rho
        if RCM_grid.mask_ocean(loni,lati)==1
            [RCM_EKWC.corr_NWwind_vgos(loni,lati), RCM_EKWC.p_NWwind_vgos(loni,lati)]= ...
                corr(squeeze(RCM_data_wind.yearly_mean_NWwind(loni,lati,:)), squeeze(RCM_data_gos.yearly_vgos(loni,lati,:)));
%                 corr(squeeze(RCM_EKWC.yearly_NWwind_sp_mean), squeeze(RCM_data_v.yearly_mean(loni,lati,1,:)));
%                 corr(squeeze(RCM_data_wind.yearly_mean_NWwind(loni,lati,:)), squeeze(RCM_data_v.yearly_mean(loni,lati,1,:)));

        else
            RCM_EKWC.corr_NWwind_vgos(loni,lati)=NaN;
            RCM_EKWC.p_NWwind_vgos(loni,lati)=NaN;
        end
    end
end
RCM_EKWC.corr_sig_NWwind_vgos=RCM_EKWC.corr_NWwind_vgos;
RCM_EKWC.corr_sig_NWwind_vgos(RCM_EKWC.p_NWwind_vgos>0.1)=NaN;

RCM_EKWC.mean_u_rho=mean(RCM_data_gos.yearly_ugos(:,:,:),3);
RCM_EKWC.mean_v_rho=mean(RCM_data_gos.yearly_vgos(:,:,:),3);

run(tmp.param_script)
m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);
m_pcolor(RCM_grid.cut_lon_rho', RCM_grid.cut_lat_rho', RCM_EKWC.corr_sig_NWwind_vgos'); shading flat; colorbar; colormap(cmap.bwr_map); caxis([-1 1])
m_gshhs_i('color',param.m_gshhs_line_color); m_gshhs_i('patch',param.m_gshhs_land_color);
m_grid('fontsize', param.m_grid_fontsize, 'box', param.m_grid_box_type, 'tickdir', param.m_grid_tickdir_type);
hold on
uvplot=m_quiver(RCM_grid.cut_lon_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)', ...
                RCM_grid.cut_lat_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)', ...
                RCM_EKWC.mean_u_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)' * param.m_quiver_vector_size, ...
                RCM_EKWC.mean_v_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)' * param.m_quiver_vector_size, ...
                'color', param.m_quiver_vector_color, 'AutoScale','off','LineWidth', param.m_quiver_LineWidth);

tmp.figdir = [RCM_info.figroot, tmp.fs, 'corr'];
if (exist(tmp.figdir , 'dir') ~= 7)
    mkdir(tmp.figdir);
end
tmp.tifname = [tmp.figdir, tmp.fs, 'corr_NWwind_vgos_',RCM_info.testname, '.tif'];
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [param.hor_paper_size_x, param.hor_paper_size_y]);
set(gcf,'PaperPosition', [param.paper_position_hor param.paper_position_ver param.paper_position_width param.paper_position_height]) 
hold off
saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);


%% Northwesterly, ageostrophic v correlation
for loni=1:RCM_grid.size_lon_rho
    for lati=1:RCM_grid.size_lat_rho
        if RCM_grid.mask_ocean(loni,lati)==1
            [RCM_EKWC.corr_NWwind_vagos(loni,lati), RCM_EKWC.p_NWwind_vagos(loni,lati)]= ...
                corr(squeeze(RCM_data_wind.yearly_mean_NWwind(loni,lati,:)), squeeze(RCM_data_gos.yearly_vagos(loni,lati,:)));
%                 corr(squeeze(RCM_EKWC.yearly_NWwind_sp_mean), squeeze(RCM_data_v.yearly_mean(loni,lati,1,:)));
%                 corr(squeeze(RCM_data_wind.yearly_mean_NWwind(loni,lati,:)), squeeze(RCM_data_v.yearly_mean(loni,lati,1,:)));
        else
            RCM_EKWC.corr_NWwind_vagos(loni,lati)=NaN;
            RCM_EKWC.p_NWwind_vagos(loni,lati)=NaN;
        end
    end
end
mean(RCM_EKWC.corr_NWwind_vagos(:), 'omitnan')
RCM_EKWC.corr_sig_NWwind_vagos=RCM_EKWC.corr_NWwind_vagos;
RCM_EKWC.corr_sig_NWwind_vagos(RCM_EKWC.p_NWwind_vagos>0.1)=NaN;

RCM_EKWC.mean_u_rho=mean(RCM_data_gos.yearly_uagos(:,:,:),3);
RCM_EKWC.mean_v_rho=mean(RCM_data_gos.yearly_vagos(:,:,:),3);

run(tmp.param_script)
m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);
m_pcolor(RCM_grid.cut_lon_rho', RCM_grid.cut_lat_rho', RCM_EKWC.corr_sig_NWwind_vagos'); shading flat; colorbar; colormap(cmap.bwr_map); caxis([-1 1])
m_gshhs_i('color',param.m_gshhs_line_color); m_gshhs_i('patch',param.m_gshhs_land_color);
m_grid('fontsize', param.m_grid_fontsize, 'box', param.m_grid_box_type, 'tickdir', param.m_grid_tickdir_type);
hold on
uvplot=m_quiver(RCM_grid.cut_lon_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)', ...
                RCM_grid.cut_lat_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)', ...
                RCM_EKWC.mean_u_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)' * param.m_quiver_vector_size, ...
                RCM_EKWC.mean_v_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)' * param.m_quiver_vector_size, ...
                'color', param.m_quiver_vector_color, 'AutoScale','off','LineWidth', param.m_quiver_LineWidth);

tmp.figdir = [RCM_info.figroot, tmp.fs, 'corr'];
if (exist(tmp.figdir , 'dir') ~= 7)
    mkdir(tmp.figdir);
end
tmp.tifname = [tmp.figdir, tmp.fs, 'corr_NWwind_vagos_',RCM_info.testname, '.tif'];
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [param.hor_paper_size_x, param.hor_paper_size_y]);
set(gcf,'PaperPosition', [param.paper_position_hor param.paper_position_ver param.paper_position_width param.paper_position_height]) 
hold off
saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);


%% uwind v correlation
for loni=1:RCM_grid.size_lon_rho
    for lati=1:RCM_grid.size_lat_rho
        if RCM_grid.mask_ocean(loni,lati)==1
            [RCM_EKWC.corr_Uwind_v(loni,lati), RCM_EKWC.p_Uwind_v(loni,lati)]= ...
                corr(squeeze(RCM_data_wind.all_Uwind(loni,lati,:)), squeeze(RCM_data_v.all(loni,lati,1,:)));
%                 corr(squeeze(RCM_EKWC.yearly_Uwind_sp_mean), squeeze(RCM_data_v.yearly_mean(loni,lati,1,:)));
%                 corr(squeeze(RCM_data_wind.yearly_mean_Uwind(loni,lati,:)), squeeze(RCM_data_v.yearly_mean(loni,lati,1,:)));

        else
            RCM_EKWC.corr_Uwind_v(loni,lati)=NaN;
            RCM_EKWC.p_Uwind_v(loni,lati)=NaN;
        end
    end
end
% pcolor(RCM_EKWC.corr_NWwind_v'); shading flat; colorbar; colormap(cmap.bwr_map); caxis([-1 1])
% pcolor(RCM_EKWC.p_NWwind_v'); shading flat; colorbar;colormap(parula);
RCM_EKWC.corr_sig_Uwind_v=RCM_EKWC.corr_Uwind_v;
RCM_EKWC.corr_sig_Uwind_v(RCM_EKWC.p_Uwind_v>0.1)=NaN;

RCM_EKWC.mean_u_rho=mean(RCM_data_u.yearly_mean(:,:,1,:),4);
RCM_EKWC.mean_v_rho=mean(RCM_data_v.yearly_mean(:,:,1,:),4);

run(tmp.param_script)
m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);
m_pcolor(RCM_grid.cut_lon_rho', RCM_grid.cut_lat_rho', RCM_EKWC.corr_sig_Uwind_v'); shading flat; colorbar; colormap(cmap.bwr_map); caxis([-1 1])
m_gshhs_i('color',param.m_gshhs_line_color); m_gshhs_i('patch',param.m_gshhs_land_color);
m_grid('fontsize', param.m_grid_fontsize, 'box', param.m_grid_box_type, 'tickdir', param.m_grid_tickdir_type);
hold on
uvplot=m_quiver(RCM_grid.cut_lon_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)', ...
                RCM_grid.cut_lat_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)', ...
                RCM_EKWC.mean_u_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)' * param.m_quiver_vector_size, ...
                RCM_EKWC.mean_v_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)' * param.m_quiver_vector_size, ...
                'color', param.m_quiver_vector_color, 'AutoScale','off','LineWidth', param.m_quiver_LineWidth);

tmp.figdir = [RCM_info.figroot, tmp.fs, 'corr'];
if (exist(tmp.figdir , 'dir') ~= 7)
    mkdir(tmp.figdir);
end
tmp.tifname = [tmp.figdir, tmp.fs, 'corr_Uwind_v_',RCM_info.testname, '.tif'];
set(gcf, 'PaperUnits', 'points');
set(gcf, 'PaperSize', [param.hor_paper_size_x, param.hor_paper_size_y]);
set(gcf,'PaperPosition', [param.paper_position_hor param.paper_position_ver param.paper_position_width param.paper_position_height]) 
hold off
saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);







% corrcoef(RCM_EKWC.yearly_Uwind_sp_mean, RCM_EKWC.yearly_v_sp_mean)
% corrcoef(RCM_EKWC.yearly_Vwind_sp_mean, RCM_EKWC.yearly_v_sp_mean)
% corrcoef(RCM_EKWC.yearly_NWwind_sp_mean, RCM_EKWC.yearly_v_sp_mean)
% corrcoef(RCM_EKWC.yearly_NWwind_sp_mean, RCM_EKWC.yearly_maxval)
% corrcoef(RCM_EKWC.yearly_NWwind_sp_mean, RCM_EKWC.yearly_sep_lat)
% corrcoef(RCM_data_v.yearly_mean_transport_w, RCM_EKWC.yearly_v_sp_mean)
% corrcoef(RCM_data_v.yearly_mean_transport_w, RCM_EKWC.yearly_sep_lat)

% %% get spatial mean value of Uwind and Vwind
% %% EKWC2 all area
% [RCM_EKWC.yearly_v_sp_mean, tmp.error_status] = Func_0011_get_area_weighted_mean(squeeze(RCM_data_v.yearly_mean(:,:,1,:)), ...
%     RCM_grid.cut_lon_rho(:,:), RCM_grid.cut_lat_rho(:,:));
% [RCM_EKWC.yearly_Vwind_sp_mean, tmp.error_status] = Func_0011_get_area_weighted_mean(RCM_data_wind.yearly_mean_Vwind(:,:,:), ...
%     RCM_grid.cut_lon_rho(:,:), RCM_grid.cut_lat_rho(:,:));
% [RCM_EKWC.yearly_Uwind_sp_mean, tmp.error_status] = Func_0011_get_area_weighted_mean(RCM_data_wind.yearly_mean_Uwind(:,:,:), ...
%     RCM_grid.cut_lon_rho(:,:), RCM_grid.cut_lat_rho(:,:));
% RCM_data_wind.yearly_mean_NWwind = RCM_data_wind.yearly_mean_Uwind * cosd(45) ...
%     - RCM_data_wind.yearly_mean_Vwind * cosd(45); 
% [RCM_EKWC.yearly_NWwind_sp_mean, tmp.error_status] = Func_0011_get_area_weighted_mean(RCM_data_wind.yearly_mean_NWwind(:,:,:), ...
%     RCM_grid.cut_lon_rho(:,:), RCM_grid.cut_lat_rho(:,:));

%% get thickness
for ti=1:RCM_info.size_t
    for lati=1:RCM_grid.size_lat_rho
        tmp.val=RCM_data_v.all(RCM_data_v.maxind(lati,ti), lati, :, ti);
        if isnan(tmp.val(:,:,1))  % if all value is NaN;
            RCM_EKWC.thickness_ind(lati,ti)=NaN;
            RCM_EKWC.thickness(lati,ti)=NaN;
        else
            tmp.val2=max(find(tmp.val<=RCM_data_v.maxval(lati,ti)/5.0));
            if isempty(tmp.val2)  % if thickness of the EKWC = water column thickness
                RCM_EKWC.thickness_ind(lati,ti)=sum(isfinite(tmp.val));
            else % find thickness of the EKWC in the deep sea
                RCM_EKWC.thickness_ind(lati,ti)=min(find(tmp.val<=RCM_data_v.maxval(lati,ti)/2.0));
            end
            RCM_EKWC.thickness(lati,ti)=RCM_grid.depth(RCM_EKWC.thickness_ind(lati,ti)) + ...
                (RCM_grid.depth(RCM_EKWC.thickness_ind(lati,ti)+1)-RCM_grid.depth(RCM_EKWC.thickness_ind(lati,ti)))/2;
        end
%         disp(num2str(tmp.val2))
    end
    RCM_EKWC.thickness(RCM_EKWC.sep_lat_ind(ti)+1:end,ti)=NaN; % remove thickness of deep mixed water in the northern sea -> deeper than 2000, weird value
    
    RCM_EKWC.thickness_maxv(ti)=RCM_EKWC.thickness(RCM_EKWC.maxind(ti),ti);
end

%% get yearly thickness
for ti=1:RCM_info.size_years
    for lati=1:RCM_grid.size_lat_rho
        tmp.val=RCM_data_v.yearly_mean(RCM_data_v.yearly_maxind(lati,ti), lati, :, ti);
        if isnan(tmp.val(:,:,1))  % if all value is NaN;
            RCM_EKWC.yearly_thickness_ind(lati,ti)=NaN;
            RCM_EKWC.yearly_thickness(lati,ti)=NaN;
        else
            tmp.val2=max(find(tmp.val<=RCM_data_v.yearly_maxval(lati,ti)/5.0));
            if isempty(tmp.val2)  % if thickness of the EKWC = water column thickness
                RCM_EKWC.yearly_thickness_ind(lati,ti)=sum(isfinite(tmp.val));
            else % find thickness of the EKWC in the deep sea
                RCM_EKWC.yearly_thickness_ind(lati,ti)=min(find(tmp.val<=RCM_data_v.yearly_maxval(lati,ti)/2.0));
            end
            RCM_EKWC.yearly_thickness(lati,ti)=RCM_grid.depth(RCM_EKWC.yearly_thickness_ind(lati,ti)) + ...
                (RCM_grid.depth(RCM_EKWC.yearly_thickness_ind(lati,ti)+1)-RCM_grid.depth(RCM_EKWC.yearly_thickness_ind(lati,ti)))/2;
        end
%         disp(num2str(tmp.val2))
    end
    RCM_EKWC.yearly_thickness(RCM_EKWC.yearly_sep_lat_ind(ti)+1:end,ti)=NaN; % remove thickness of deep mixed water in the northern sea -> deeper than 2000, weird value
    
    RCM_EKWC.yearly_thickness_maxv(ti)=RCM_EKWC.yearly_thickness(RCM_EKWC.yearly_maxind(ti),ti);
end

% corrcoef(RCM_EKWC.yearly_thickness_maxv, RCM_EKWC.yearly_sep_lat)


% % %% get width
% % %% WSC weighted mean latitude based on latitudes where WSC was greater than 50 percent of its maximum.
% % RCM_EKWC.w_bndy_ind=NaN(RCM_grid.size_lat_rho, RCM_info.size_t);
% % RCM_EKWC.e_bndy_ind=NaN(RCM_grid.size_lat_rho, RCM_info.size_t);
% % RCM_EKWC.width_ind=NaN(RCM_grid.size_lat_rho, RCM_info.size_t);
% % RCM_EKWC.width=NaN(RCM_grid.size_lat_rho, RCM_info.size_t);
% % for ti=1:RCM_info.size_t
% % %     for lati=1:RCM_EKWC.sep_lat_ind(ti)
% %     for lati=1:RCM_EKWC.lat_limit_ind
% %         tmp.val=RCM_data_v.all(:, lati, 1, ti);
% %         if sum(isnan(tmp.val))==RCM_grid.size_lon_rho
% %             RCM_EKWC.w_bndy_ind(lati,ti)=NaN;
% %             RCM_EKWC.e_bndy_ind(lati,ti)=NaN;
% % %         elseif max(max(tmp.val))<0
% % %             RCM_EKWC.w_bndy_ind(lati,ti)=NaN;
% % %             RCM_EKWC.e_bndy_ind(lati,ti)=NaN;
% %         else
% % %             tmp.val(tmp.val<0)=NaN;
% %             tmp.val(tmp.val<RCM_data_v.maxval(lati,ti)/20.0)=NaN;
% %             RCM_EKWC.w_bndy_ind(lati,ti)=min(find(diff(isfinite(tmp.val))==1));
% %             RCM_EKWC.e_bndy_ind(lati,ti)=min(find(diff(isfinite(tmp.val))==-1));
% %             RCM_EKWC.width_ind(lati,ti)=RCM_EKWC.e_bndy_ind(lati,ti) - RCM_EKWC.w_bndy_ind(lati,ti) +1;
% %         end
% %     end
% %     RCM_EKWC.width_ind_maxv(ti)=RCM_EKWC.width_ind(RCM_EKWC.maxind(ti),ti);
% %     RCM_EKWC.width_maxv(ti)=m_lldist([RCM_grid.cut_lon_rho(1,lati) RCM_grid.cut_lon_rho(1+RCM_EKWC.width_ind_maxv(ti),lati)], ...
% %         [RCM_grid.cut_lat_rho(1,lati) RCM_grid.cut_lat_rho(1+RCM_EKWC.width_ind_maxv(ti),lati)]);
% % end

% % %% get yearly width
% % RCM_EKWC.yearly_w_bndy_ind=NaN(RCM_grid.size_lat_rho, RCM_info.size_years);
% % RCM_EKWC.yearly_e_bndy_ind=NaN(RCM_grid.size_lat_rho, RCM_info.size_years);
% % RCM_EKWC.yearly_width_ind=NaN(RCM_grid.size_lat_rho, RCM_info.size_years);
% % RCM_EKWC.yearly_width=NaN(RCM_grid.size_lat_rho, RCM_info.size_years);
% % for ti=1:RCM_info.size_years
% %     for lati=1:RCM_EKWC.yearly_sep_lat_ind(ti)
% %         tmp.val=RCM_data_v.yearly_mean(:, lati, 1, ti);
% %         if sum(isnan(tmp.val))==RCM_grid.size_lon_rho
% %             RCM_EKWC.yearly_w_bndy_ind(lati,ti)=NaN;
% %             RCM_EKWC.yearly_e_bndy_ind(lati,ti)=NaN;
% %         else
% %             tmp.val(tmp.val<RCM_data_v.yearly_maxval(lati,ti)/5.0)=NaN;
% %             RCM_EKWC.yearly_w_bndy_ind(lati,ti)=min(find(diff(isfinite(tmp.val))==1));
% %             RCM_EKWC.yearly_e_bndy_ind(lati,ti)=min(find(diff(isfinite(tmp.val))==-1));
% %             RCM_EKWC.yearly_width_ind(lati,ti)=RCM_EKWC.yearly_e_bndy_ind(lati,ti) - RCM_EKWC.yearly_w_bndy_ind(lati,ti) +1;
% %         end
% %     end
% %     RCM_EKWC.yearly_width_ind_maxv(ti)=RCM_EKWC.yearly_width_ind(RCM_EKWC.yearly_maxind(ti),ti);
% %     RCM_EKWC.yearly_width_maxv(ti)=m_lldist([RCM_grid.cut_lon_rho(1,lati) RCM_grid.cut_lon_rho(1+RCM_EKWC.yearly_width_ind_maxv(ti),lati)], ...
% %         [RCM_grid.cut_lat_rho(1,lati) RCM_grid.cut_lat_rho(1+RCM_EKWC.yearly_width_ind_maxv(ti),lati)]);
% % end

% RCM_EKWC.yearly_comb_val(:,1) = RCM_EKWC.yearly_maxval;
% RCM_EKWC.yearly_comb_val(:,2) = RCM_EKWC.yearly_NWwind_sp_mean;
% RCM_EKWC.yearly_comb_val(:,3) = RCM_EKWC.yearly_Uwind_sp_mean;
% RCM_EKWC.yearly_comb_val(:,4) = RCM_EKWC.yearly_sep_lat;
% RCM_EKWC.yearly_comb_val(:,5) = RCM_EKWC.yearly_transport_w;
% RCM_EKWC.yearly_comb_val(:,6) = RCM_EKWC.yearly_transport;
% RCM_EKWC.yearly_comb_val(:,7) = RCM_EKWC.yearly_width_maxv;
% RCM_EKWC.yearly_comb_val(:,8) = RCM_EKWC.yearly_thickness_maxv;
% RCM_EKWC.yearly_comb_val(:,9) = RCM_EKWC.yearly_maxind;
% 
% RCM_EKWC.analysis=array2table(RCM_EKWC.yearly_comb_val);
% RCM_EKWC.analysis.Properties.VariableNames = ...
% {'max-v', 'NWwind', 'Uwind', 'sep-lat', 'tr-w', 'tr', 'wid', 'thick', 'max-lat'};

% % % % % RCM_EKWC.yearly_comb_val(:,1) = RCM_EKWC.yearly_maxval;
% % % % % RCM_EKWC.yearly_comb_val(:,2) = RCM_EKWC.yearly_NWwind_sp_mean;
% % % % % RCM_EKWC.yearly_comb_val(:,3) = RCM_EKWC.yearly_sep_lat;
% % % % % RCM_EKWC.yearly_comb_val(:,4) = RCM_EKWC.yearly_transport;
% % % % % RCM_EKWC.yearly_comb_val(:,5) = RCM_EKWC.yearly_width_maxv;
% % % % % RCM_EKWC.yearly_comb_val(:,6) = RCM_EKWC.yearly_thickness_maxv;
% % % % % RCM_EKWC.yearly_comb_val(:,7) = RCM_EKWC.yearly_maxind;
% % % % % 
% % % % % RCM_EKWC.analysis=array2table(RCM_EKWC.yearly_comb_val);
% % % % % RCM_EKWC.analysis.Properties.VariableNames = ...
% % % % % {'max-v', 'NWwind', 'sep-lat', 'tr', 'wid', 'thick', 'max-lat'};
% % % % % 
% % % % % corrplot(RCM_EKWC.analysis, 'testR', 'on')
% % % % % 
% % % % % tmp.figdir = [RCM_info.figroot, tmp.fs, 'corr'];
% % % % % if (exist(tmp.figdir , 'dir') ~= 7)
% % % % %     mkdir(tmp.figdir);
% % % % % end
% % % % % tmp.tifname = [tmp.figdir, tmp.fs, 'corr_mat_',RCM_info.testname, '.tif'];
% % % % % 
% % % % % set(gcf, 'PaperUnits', 'points');
% % % % % set(gcf, 'PaperSize', [param.hor_paper_size_x, param.hor_paper_size_y]);
% % % % % set(gcf,'PaperPosition', [param.paper_position_hor param.paper_position_ver param.paper_position_width param.paper_position_height]) 
% % % % % saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
% % % % %          
% % % % % close all;
% % % % % % disp('abc')
% % % % % 
% % % % % pcolor(RCM_grid.cut_lon_rho', RCM_grid.cut_lat_rho',RCM_data_v.all(:,:,1,24)'); 
% % % % % shading flat; colorbar; colormap(cmap.bwr_map); caxis([-0.6, 0.6]);
% % % % % hold on
% % % % % 
% % % % % tmp.val=RCM_data_v.all(:,:,1,29);
% % % % % tmp.maxval=max(RCM_data_v.all(:,24));
% % % % % tmp.val(tmp.val<tmp.maxval*(1/10))=NaN;
% % % % % contour(RCM_grid.cut_lon_rho', RCM_grid.cut_lat_rho',tmp.val', 'y'); 
% % % % % hold off;

% mean(RCM_EKWC.yearly_thickness_maxv)
% mean(RCM_data_v.yearly_mean_transport_w)
% mean(RCM_EKWC.yearly_width_maxv)
% plot(RCM_EKWC.Vwind_sp_mean/mean(RCM_EKWC.Vwind_sp_mean))
% hold on
% plot(RCM_EKWC.v_sp_mean/mean(RCM_EKWC.v_sp_mean))
% hold off

% plot(RCM_EKWC.NWwind_sp_mean/mean(RCM_EKWC.NWwind_sp_mean))
% hold on
% plot(RCM_EKWC.yearly_sep_lat/mean(RCM_EKWC.yearly_sep_lat))
% hold off