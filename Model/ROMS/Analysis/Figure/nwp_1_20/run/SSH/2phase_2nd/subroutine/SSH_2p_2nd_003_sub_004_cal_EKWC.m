disp('subroutine v_2p_2nd_001_sub_003_cal_EKWC') 
%% calculate EKWC

%% get data
RCM_info.savedir_v = [RCM_info.saveroot, RCM_info.testname, tmp.fs, ...
            'v', tmp.fs];
RCM_info.matname_v = [RCM_info.savedir_v,RCM_info.testname,'_',RCM_info.regionname, '_RCM_','v', '_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), '_', ...
   RCM_info.season, '.mat'];
load(RCM_info.matname_v, 'RCM_data')
RCM_data_v=RCM_data;

RCM_info.savedir_u = [RCM_info.saveroot, RCM_info.testname, tmp.fs, ...
            'u', tmp.fs];
RCM_info.matname_u = [RCM_info.savedir_u,RCM_info.testname,'_',RCM_info.regionname, '_RCM_', 'u', '_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), '_', ...
   RCM_info.season, '.mat'];
load(RCM_info.matname_u, 'RCM_data', 'RCM_grid')
RCM_data_u=RCM_data;

RCM_info_wind=RCM_info;
RCM_grid_wind=RCM_grid;
RCM_info_wind.matname = [RCM_info_wind.windsavedir,RCM_info_wind.testname,'_',RCM_info_wind.regionname, '_RCM_data_wind_', ...
    num2str(min(RCM_info_wind.years),'%04i'),'_',num2str(max(RCM_info_wind.years),'%04i'), '_', ...
   RCM_info.season, '.mat'];
load(RCM_info_wind.matname, 'RCM_data_wind')

RCM_info.savedir_salt = [RCM_info.saveroot, RCM_info.testname, tmp.fs, ...
            'salt', tmp.fs];
RCM_info.matname_salt = [RCM_info.savedir_salt,RCM_info.testname,'_',RCM_info.regionname, '_RCM_','salt', '_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), '_', ...
   RCM_info.season, '.mat'];
load(RCM_info.matname_salt, 'RCM_data')
RCM_data_salt=RCM_data;

RCM_info.savedir_temp = [RCM_info.saveroot, RCM_info.testname, tmp.fs, ...
            'temp', tmp.fs];
RCM_info.matname_temp = [RCM_info.savedir_temp,RCM_info.testname,'_',RCM_info.regionname, '_RCM_','temp', '_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), '_', ...
   RCM_info.season, '.mat'];
load(RCM_info.matname_temp, 'RCM_data')
RCM_data_temp=RCM_data;

% RCM_info.gos_matname = [RCM_info.atmfiledir, tmp.fs, 'zeta', tmp.fs, RCM_info.testname,'_',RCM_info.regionname, '_RCM_gos_', ...
%     num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'),'_', ...
%    RCM_info.season, '.mat'];
% load(RCM_info.gos_matname, 'RCM_data_gos')

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
% RCM_EKWC.lat_limit_ind= find(RCM_grid.lat_rho(1,:)==[39.24]);
[tmp.indw, tmp.inde, tmp.inds, tmp.indn]=Func_0012_findind_Y(1/20,[130 130 41 41],RCM_grid.cut_lon_rho,RCM_grid.cut_lat_rho, 1);
RCM_EKWC.lat_limit_ind= round((tmp.inds+tmp.indn)/2);

[tmp.indw, tmp.inde, tmp.inds, tmp.indn]=Func_0012_findind_Y(1/20,[130 130 36.5 36.5],RCM_grid.cut_lon_rho,RCM_grid.cut_lat_rho, 1);
RCM_EKWC.lat_lowlimit_ind= round((tmp.inds+tmp.indn)/2);

% RCM_EKWC.lat_limit_ind=149;  % avoid separated current north of the EKB
% [RCM_EKWC.maxval_v, RCM_EKWC.maxind] = max(RCM_data_v.maxval(1:RCM_EKWC.lat_limit_ind,:), [], 1, 'omitnan');
% [RCM_EKWC.yearly_maxval_v, RCM_EKWC.yearly_maxind] = max(RCM_data_v.yearly_maxval(1:RCM_EKWC.lat_limit_ind,:), [], 1, 'omitnan');

RCM_EKWC.lon_limit=130;  % for defining separation latitude
RCM_EKWC.lon_limit_ind=find(RCM_grid.cut_lon_rho(:,1)==RCM_EKWC.lon_limit);

%% get separation latitude (1/20 of maxv_EKWC)
for ti=1:RCM_info.size_t
    tmp.checkudata=RCM_data_u.all;
    tmp.checkudata(:,1:RCM_EKWC.lat_lowlimit_ind-1,1,ti)=NaN;
    tmp.checkudata(:,RCM_EKWC.lat_limit_ind+1:end,1,ti)=NaN;

    [RCM_EKWC.maxval_u(ti), RCM_EKWC.maxlatind_u(ti)] = max(tmp.checkudata(RCM_EKWC.lon_limit_ind,:,1,ti), [], 2, 'omitnan');
    RCM_EKWC.maxlat_u(ti)=RCM_grid.cut_lat_rho(RCM_EKWC.lon_limit_ind, RCM_EKWC.maxlatind_u(ti));

%     tmp.val = min(find(squeeze(RCM_data_v.maxval(:,ti))<=RCM_EKWC.maxval(ti)/20.0));
%     if RCM_data_v.maxval(tmp.val,ti)<0
%         tmp.val=tmp.val-1;
%     end
%     if isempty(tmp.val)
%         RCM_EKWC.sep_lat_ind(ti)=RCM_grid.size_lat_rho-2;
%     else
%         RCM_EKWC.sep_lat_ind(ti)=tmp.val;
%     end
%     RCM_EKWC.sep_lat(ti)=RCM_grid.cut_lat_rho(1,RCM_EKWC.sep_lat_ind(ti));
end
for ti=1:RCM_info.size_years
    tmp.checkudata=RCM_data_u.yearly_mean;
    tmp.checkudata(:,1:RCM_EKWC.lat_lowlimit_ind-1,1,ti)=NaN;
    tmp.checkudata(:,RCM_EKWC.lat_limit_ind+1:end,1,ti)=NaN;
    
    [RCM_EKWC.yearly_maxval_u(ti), RCM_EKWC.yearly_maxlatind_u(ti)] = max(tmp.checkudata(RCM_EKWC.lon_limit_ind,:,1,ti), [], 2, 'omitnan');
    RCM_EKWC.yearly_maxlat_u(ti)=RCM_grid.cut_lat_rho(RCM_EKWC.lon_limit_ind, RCM_EKWC.yearly_maxlatind_u(ti));
    
%     tmp.val = min(find(squeeze(RCM_data_v.yearly_maxval(:,ti))<=RCM_EKWC.yearly_maxval(ti)/10.0));
%     if RCM_data_v.yearly_maxval(tmp.val,ti)<0
%         tmp.val=tmp.val-1;
%     end
%     if isempty(tmp.val)
%         RCM_EKWC.yearly_sep_lat_ind(ti)= RCM_grid.size_lat_rho-2;
%     else
%         RCM_EKWC.yearly_sep_lat_ind(ti)=tmp.val;
%     end
%     RCM_EKWC.yearly_sep_lat(ti)=RCM_grid.cut_lat_rho(1,RCM_EKWC.yearly_sep_lat_ind(ti));
end

%% get spatial mean value of Uwind and Vwind
%% EKWC coastal area only (<130.4E)
% RCM_EKWC.lon_limit_ind= 70;
[tmp.indw, tmp.inde, tmp.inds, tmp.indn]=Func_0012_findind_Y(1/20,[130.0 130.0 33 40],RCM_grid.lon_rho,RCM_grid.lat_rho, 1);
RCM_EKWC.lon_limit_ind= round((tmp.indw+tmp.inde)/2);

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

run(tmp.param_script)


RCM_EKWC.yearly_mean_surf_temp=mean(RCM_data_temp.yearly_mean(:,:,1),4);
RCM_EKWC.yearly_mean_surf_temp_100m=mean(RCM_data_temp.yearly_mean(:,:,11),4);
RCM_EKWC.yearly_mean_surf_salt=mean(RCM_data_salt.yearly_mean(:,:,1),4);
RCM_EKWC.yearly_mean_surf_salt_100m=mean(RCM_data_salt.yearly_mean(:,:,11),4);
RCM_EKWC.yearly_mean_surf_u=mean(RCM_data_u.yearly_mean(:,:,1),4);
RCM_EKWC.yearly_mean_surf_u_100m=mean(RCM_data_u.yearly_mean(:,:,11),4);
RCM_EKWC.yearly_mean_surf_v=mean(RCM_data_v.yearly_mean(:,:,1),4);
RCM_EKWC.yearly_mean_surf_v_100m=mean(RCM_data_v.yearly_mean(:,:,11),4);
% RCM_EKWC.yearly_mean_surf_zeta=mean(RCM_data_zeta.yearly_mean(:,:,1),4);

RCM_EKWC.yearly_mean_surf_CT=gsw_CT_from_pt(RCM_EKWC.yearly_mean_surf_salt, RCM_EKWC.yearly_mean_surf_temp);
RCM_EKWC.yearly_mean_surf_CT_100m=gsw_CT_from_pt(RCM_EKWC.yearly_mean_surf_salt_100m, RCM_EKWC.yearly_mean_surf_temp_100m);

RCM_EKWC.ref_pressure=10.1325;
RCM_EKWC.yearly_mean_surf_prho= gsw_rho(RCM_EKWC.yearly_mean_surf_salt,RCM_EKWC.yearly_mean_surf_CT,RCM_EKWC.ref_pressure);
RCM_EKWC.yearly_mean_surf_prho_100m= gsw_rho(RCM_EKWC.yearly_mean_surf_salt_100m,RCM_EKWC.yearly_mean_surf_CT_100m,RCM_EKWC.ref_pressure);

RCM_info.savedir_ekwc = [RCM_info.saveroot, RCM_info.testname, tmp.fs, ...
            'EKWC', tmp.fs];
if (exist(RCM_info.savedir_ekwc , 'dir') ~= 7)
    mkdir(RCM_info.savedir_ekwc);
end
RCM_info.matname_ekwc = [RCM_info.savedir_ekwc,RCM_info.testname,'_',RCM_info.regionname, '_RCM_','EKWC', '_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), '_', ...
   RCM_info.season, '.mat'];
save(RCM_info.matname_ekwc, 'RCM_EKWC', 'RCM_grid', 'RCM_info')


% for ti=1:RCM_info.size_years
%     if (isfield(tmp, 'ref_vec_x_range') ~= 1)
%         tmp.ref_vec_x_ind = find(abs(RCM_grid.cut_lon_rho(:,1)-param.m_quiver_ref_text_x_location) ...
%             == min(abs(RCM_grid.cut_lon_rho(:,1)-param.m_quiver_ref_text_x_location)));
%         tmp.ref_vec_y_ind = find(abs(RCM_grid.cut_lat_rho(1,:)-param.m_quiver_ref_text_y_location) ...
%             == min(abs(RCM_grid.cut_lat_rho(1,:)-param.m_quiver_ref_text_y_location)))+param.m_quiver_y_interval*2;
%         tmp.ref_vec_x_range = round(tmp.ref_vec_x_ind-(param.m_quiver_x_interval/2)) : ...
%             round(tmp.ref_vec_x_ind-(param.m_quiver_x_interval/2))+param.m_quiver_x_interval-4;
%         tmp.ref_vec_y_range = round(tmp.ref_vec_y_ind-(param.m_quiver_y_interval/2)) : ...
%             round(tmp.ref_vec_y_ind-(param.m_quiver_y_interval/2))+param.m_quiver_y_interval-1;
%     end
%     RCM_data_u.yearly_mean(tmp.ref_vec_x_range,tmp.ref_vec_y_range)=param.m_quiver_ref_u_value;
%     RCM_data_v.yearly_mean(tmp.ref_vec_x_range,tmp.ref_vec_y_range)=param.m_quiver_ref_v_value;     
%     
%     m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);
%     m_gshhs_i('color',param.m_gshhs_line_color); m_gshhs_i('patch',param.m_gshhs_land_color);
%     m_grid('fontsize', param.m_grid_fontsize, 'box', param.m_grid_box_type, 'tickdir', param.m_grid_tickdir_type);
%     hold on
%     uvplot=m_quiver(RCM_grid.cut_lon_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)', ...
%                     RCM_grid.cut_lat_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)', ...
%                     RCM_data_u.yearly_mean(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end,1,ti)' * param.m_quiver_vector_size, ...
%                     RCM_data_v.yearly_mean(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end,1,ti)' * param.m_quiver_vector_size, ...
%                     'color', param.m_quiver_vector_color, 'AutoScale','off','LineWidth', param.m_quiver_LineWidth);
%     m_text(param.m_quiver_ref_text_x_location, param.m_quiver_ref_text_y_location, param.m_quiver_ref_text, 'FontSize', param.m_quiver_ref_text_fontsize); 
% end






% %% get thickness
% for ti=1:RCM_info.size_t
%     for lati=1:RCM_grid.size_lat_rho
%         tmp.val=RCM_data_v.all(RCM_data_v.maxind(lati,ti), lati, :, ti);
%         if isnan(tmp.val(:,:,1))  % if all value is NaN;
%             RCM_EKWC.thickness_ind(lati,ti)=NaN;
%             RCM_EKWC.thickness(lati,ti)=NaN;
%         else
%             tmp.val2=max(find(tmp.val<=RCM_data_v.maxval(lati,ti)/5.0));
%             if isempty(tmp.val2)  % if thickness of the EKWC = water column thickness
%                 RCM_EKWC.thickness_ind(lati,ti)=sum(isfinite(tmp.val));
%             else % find thickness of the EKWC in the deep sea
%                 RCM_EKWC.thickness_ind(lati,ti)=min(find(tmp.val<=RCM_data_v.maxval(lati,ti)/2.0));
%             end
%             RCM_EKWC.thickness(lati,ti)=RCM_grid.depth(RCM_EKWC.thickness_ind(lati,ti)) + ...
%                 (RCM_grid.depth(RCM_EKWC.thickness_ind(lati,ti)+1)-RCM_grid.depth(RCM_EKWC.thickness_ind(lati,ti)))/2;
%         end
% %         disp(num2str(tmp.val2))
%     end
%     RCM_EKWC.thickness(RCM_EKWC.sep_lat_ind(ti)+1:end,ti)=NaN; % remove thickness of deep mixed water in the northern sea -> deeper than 2000, weird value
%     
%     RCM_EKWC.thickness_maxv(ti)=RCM_EKWC.thickness(RCM_EKWC.maxind(ti),ti);
% end
% 
% %% get yearly thickness
% for ti=1:RCM_info.size_years
%     for lati=1:RCM_grid.size_lat_rho
%         tmp.val=RCM_data_v.yearly_mean(RCM_data_v.yearly_maxind(lati,ti), lati, :, ti);
%         if isnan(tmp.val(:,:,1))  % if all value is NaN;
%             RCM_EKWC.yearly_thickness_ind(lati,ti)=NaN;
%             RCM_EKWC.yearly_thickness(lati,ti)=NaN;
%         else
%             tmp.val2=max(find(tmp.val<=RCM_data_v.yearly_maxval(lati,ti)/5.0));
%             if isempty(tmp.val2)  % if thickness of the EKWC = water column thickness
%                 RCM_EKWC.yearly_thickness_ind(lati,ti)=sum(isfinite(tmp.val));
%             else % find thickness of the EKWC in the deep sea
%                 RCM_EKWC.yearly_thickness_ind(lati,ti)=min(find(tmp.val<=RCM_data_v.yearly_maxval(lati,ti)/2.0));
%             end
%             RCM_EKWC.yearly_thickness(lati,ti)=RCM_grid.depth(RCM_EKWC.yearly_thickness_ind(lati,ti)) + ...
%                 (RCM_grid.depth(RCM_EKWC.yearly_thickness_ind(lati,ti)+1)-RCM_grid.depth(RCM_EKWC.yearly_thickness_ind(lati,ti)))/2;
%         end
% %         disp(num2str(tmp.val2))
%     end
%     RCM_EKWC.yearly_thickness(RCM_EKWC.yearly_sep_lat_ind(ti)+1:end,ti)=NaN; % remove thickness of deep mixed water in the northern sea -> deeper than 2000, weird value
%     
%     RCM_EKWC.yearly_thickness_maxv(ti)=RCM_EKWC.yearly_thickness(RCM_EKWC.yearly_maxind(ti),ti);
% end