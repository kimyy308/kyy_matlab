% %  Updated 19-May-2022 by Yong-Yub Kim, 

% start-------------------- earlier decadal SST, SSS plot
%% set variable name & figure directory
tmp.variable='UV';
dirs.figdir=[dirs.figrawdir,'surface', tmp.fs, tmp.regionname, tmp.fs, 'vec', tmp.fs, ...
    num2str(min(RCM_info.years)), '_', num2str(max(RCM_info.years)), tmp.fs];
if (exist(strcat(dirs.figdir) , 'dir') ~= 7)
    mkdir(strcat(dirs.figdir));
end 

%% set figure file name
tmp.tifname=strcat(dirs.figdir, tmp.testname, '_surf_', 'vec', '_', ...
    num2str(min(RCM_info.years),'%04i'), '_',num2str(max(RCM_info.years),'%04i'), ...
    '_', RCM_info.season,'.tif'); %% ~_year_month.jpg
if (exist(tmp.tifname , 'file') ~= 2 || fig_flag==2)      
    run(tmp.param_script);
%% set data file name (mat) for reading
    
    switch tmp.testname
            case 'ENS4_hist'
               tmp.tset= {'test2117', 'test2118','test2119','test2120'};
               for titi=1:4
                    tmp.ttname=tmp.tset{titi};
                    dirs.ttmatdir = strcat('/Volumes/kyy_raid/Data/Model/ROMS/', RCM_info.model, '/', tmp.ttname, '/run/mean/');
                    tmp.matname_u = [dirs.ttmatdir, tmp.ttname, '_', tmp.regionname, '_', 'u',...
                        '_mean_', num2str(min(RCM_info.years),'%04i'), '-', num2str(max(RCM_info.years),'%04i'),...
                    '_', RCM_info.season, '.mat'];
                    if (exist(tmp.matname_u , 'file') ~= 2)
                        disp('please get data from subroutine_004_001 first')
                    else
                        load(tmp.matname_u);
                    end 
                    RCM_data_u=RCM_data;
                    tmp.tmean_u(:,:,titi)=RCM_data.mean;

                    tmp.matname_v = [dirs.ttmatdir, tmp.ttname, '_', tmp.regionname, '_', 'v',...
                        '_mean_', num2str(min(RCM_info.years),'%04i'), '-', num2str(max(RCM_info.years),'%04i'),...
                    '_', RCM_info.season, '.mat'];
                    if (exist(tmp.matname_v , 'file') ~= 2)
                        disp('please get data from subroutine_004_001 first')
                    else
                        load(tmp.matname_v);
                    end 
                    RCM_data_v=RCM_data;
                    tmp.tmean_v(:,:,titi)=RCM_data.mean;
               end
               RCM_data_u.mean=mean(tmp.tmean_u,3);
               RCM_data_v.mean=mean(tmp.tmean_v,3);
           case 'ENS4_fut'
               tmp.tset= {'test2127', 'test2128','test2129','test2130'};
               for titi=1:4
                    tmp.ttname=tmp.tset{titi};
                    dirs.ttmatdir = strcat('/Volumes/kyy_raid/Data/Model/ROMS/', RCM_info.model, '/', tmp.ttname, '/run/mean/');
                    tmp.matname_u = [dirs.ttmatdir, tmp.ttname, '_', tmp.regionname, '_', 'u',...
                        '_mean_', num2str(min(RCM_info.years),'%04i'), '-', num2str(max(RCM_info.years),'%04i'),...
                    '_', RCM_info.season, '.mat'];
                    if (exist(tmp.matname_u , 'file') ~= 2)
                        disp('please get data from subroutine_004_001 first')
                    else
                        load(tmp.matname_u);
                    end 
                    RCM_data_u=RCM_data;
                    tmp.tmean_u(:,:,titi)=RCM_data.mean;

                    tmp.matname_v = [dirs.ttmatdir, tmp.ttname, '_', tmp.regionname, '_', 'v',...
                        '_mean_', num2str(min(RCM_info.years),'%04i'), '-', num2str(max(RCM_info.years),'%04i'),...
                    '_', RCM_info.season, '.mat'];
                    if (exist(tmp.matname_v , 'file') ~= 2)
                        disp('please get data from subroutine_004_001 first')
                    else
                        load(tmp.matname_v);
                    end 
                    RCM_data_v=RCM_data;
                    tmp.tmean_v(:,:,titi)=RCM_data.mean;
               end
               RCM_data_u.mean=mean(tmp.tmean_u,3);
               RCM_data_v.mean=mean(tmp.tmean_v,3);
            otherwise
                tmp.matname_u = [dirs.matdir, tmp.testname, '_', tmp.regionname, '_', 'u',...
                    '_mean_', num2str(min(RCM_info.years),'%04i'), '-', num2str(max(RCM_info.years),'%04i'),...
                '_', RCM_info.season, '.mat'];
                if (exist(tmp.matname_u , 'file') ~= 2)
                    disp('please get data from subroutine_004_001 first')
                else
                    load(tmp.matname_u);
                    RCM_data_u=RCM_data;
                end 
                
                tmp.matname_v = [dirs.matdir, tmp.testname, '_', tmp.regionname, '_', 'v',...
                    '_mean_', num2str(min(RCM_info.years),'%04i'), '-', num2str(max(RCM_info.years),'%04i'),...
                '_', RCM_info.season, '.mat'];
                if (exist(tmp.matname_v , 'file') ~= 2)
                    disp('please get data from subroutine_004_001 first')
                else
                    load(tmp.matname_v);
                    RCM_data_v=RCM_data;
                end 
        end

    

    RCM_grid.mask_model = double(inpolygon(RCM_grid.cut_lon_rho,RCM_grid.cut_lat_rho,RCM_grid.refpolygon(:,1),RCM_grid.refpolygon(:,2)));
    RCM_grid.mask_model(RCM_grid.mask_model==0)=NaN;
    RCM_data_u.mean=RCM_data_u.mean.*RCM_grid.mask_model;
    RCM_data_v.mean=RCM_data_v.mean.*RCM_grid.mask_model;
    
    if (isfield(tmp, 'ref_vec_x_range') ~= 1)
        tmp.ref_vec_x_ind = find(abs(RCM_grid.cut_lon_rho(:,1)-param.m_quiver_ref_text_x_location) ...
            == min(abs(RCM_grid.cut_lon_rho(:,1)-param.m_quiver_ref_text_x_location)));
        tmp.ref_vec_y_ind = find(abs(RCM_grid.cut_lat_rho(1,:)-param.m_quiver_ref_text_y_location) ...
            == min(abs(RCM_grid.cut_lat_rho(1,:)-param.m_quiver_ref_text_y_location)))+param.m_quiver_y_interval*2;
        tmp.ref_vec_x_range = round(tmp.ref_vec_x_ind-(param.m_quiver_x_interval/2)) : ...
            round(tmp.ref_vec_x_ind-(param.m_quiver_x_interval/2))+param.m_quiver_x_interval*1;
        switch RCM_info.model
            case 'nwp_1_20'
                tmp.ref_vec_y_range = round(tmp.ref_vec_y_ind-(param.m_quiver_y_interval/2)) : ...
                    round(tmp.ref_vec_y_ind-(param.m_quiver_y_interval/2))+param.m_quiver_y_interval*1;
            case 'nwp_1_10'
                 tmp.ref_vec_y_range = round(tmp.ref_vec_y_ind-(param.m_quiver_y_interval/2)) : ...
                    round(tmp.ref_vec_y_ind-(param.m_quiver_y_interval/2))+param.m_quiver_y_interval/2;
        end
    end
    RCM_data_u.mean(tmp.ref_vec_x_range,tmp.ref_vec_y_range)=param.m_quiver_ref_u_value;
    RCM_data_v.mean(tmp.ref_vec_x_range,tmp.ref_vec_y_range)=param.m_quiver_ref_v_value;     
    

    %% READ OBS, RMSE, CORR
    OBS_info.filename=['/Volumes/kyy_raid/Data/Observation/CMEMS/', 'AKP4cmems_gos_1995_2014_JanFeb.nc'];
    OBS_grid.lon=ncread(OBS_info.filename, 'lon');
    OBS_grid.lat=ncread(OBS_info.filename, 'lat');
    OBS_data.raw_u=ncread(OBS_info.filename, 'cmems_u');
    OBS_data.raw_v=ncread(OBS_info.filename, 'cmems_v');
     [OBS_grid.lat_rho, OBS_grid.lon_rho] = meshgrid(OBS_grid.lat, OBS_grid.lon);

    [OBS_grid.refpolygon, OBS_grid.domain, tmp.error_status] = Func_0007_get_polygon_data_from_regionname(tmp.regionname);

    [OBS_grid.lon_min, OBS_grid.lon_max, OBS_grid.lat_min, OBS_grid.lat_max] = ...
                                findind_Y(1/20, OBS_grid.domain(1:4), OBS_grid.lon_rho, OBS_grid.lat_rho);
    OBS_grid.cut_lon_rho = ...
        OBS_grid.lon_rho(OBS_grid.lon_min(1):OBS_grid.lon_max(1), OBS_grid.lat_min(1):OBS_grid.lat_max(1));
    OBS_grid.cut_lat_rho = ...
        OBS_grid.lat_rho(OBS_grid.lon_min(1):OBS_grid.lon_max(1), OBS_grid.lat_min(1):OBS_grid.lat_max(1));
    OBS_data.cut_raw_u=OBS_data.raw_u(OBS_grid.lon_min(1):OBS_grid.lon_max(1), OBS_grid.lat_min(1):OBS_grid.lat_max(1),:);
        OBS_data.cut_raw_v=OBS_data.raw_v(OBS_grid.lon_min(1):OBS_grid.lon_max(1), OBS_grid.lat_min(1):OBS_grid.lat_max(1),:);
    
    OBS_grid.mask_model = double(inpolygon(OBS_grid.cut_lon_rho,OBS_grid.cut_lat_rho,OBS_grid.refpolygon(:,1),OBS_grid.refpolygon(:,2)));
    OBS_grid.mask_model(OBS_grid.mask_model==0)=NaN;
    OBS_data.mean_u=mean(OBS_data.cut_raw_u, 3).*OBS_grid.mask_model;
    OBS_data.mean_v=mean(OBS_data.cut_raw_v, 3).*OBS_grid.mask_model;

    for i=1:size(OBS_data.cut_raw_u,3)
        RCM_data_u.interp_all(:,:,i)=griddata(RCM_grid.cut_lon_rho, RCM_grid.cut_lat_rho, RCM_data_u.all(:,:,ceil(i/2),2-mod(i,2)), ...
            OBS_grid.cut_lon_rho, OBS_grid.cut_lat_rho);
        RCM_data_v.interp_all(:,:,i)=griddata(RCM_grid.cut_lon_rho, RCM_grid.cut_lat_rho, RCM_data_v.all(:,:,ceil(i/2),2-mod(i,2)), ...
            OBS_grid.cut_lon_rho, OBS_grid.cut_lat_rho);
        RCM_data_u.sqe(:,:,i)=(RCM_data_u.interp_all(:,:,i)-OBS_data.cut_raw_u(:,:,i)).^2;
        RCM_data_v.sqe(:,:,i)=(RCM_data_v.interp_all(:,:,i)-OBS_data.cut_raw_v(:,:,i)).^2;
    end
    RCM_data_u.rmse=sqrt(mean(RCM_data_u.sqe,3));
    RCM_data_v.rmse=sqrt(mean(RCM_data_v.sqe,3));
    RCM_data_u.interp_mean=mean(RCM_data_u.interp_all,3);
    RCM_data_v.interp_mean=mean(RCM_data_v.interp_all,3);
    
    %% mean of RMSE
    [tmp.m_value_u, tmp.error_status] = Func_0011_get_area_weighted_mean(RCM_data_u.rmse, OBS_grid.cut_lon_rho, OBS_grid.cut_lat_rho);
    [tmp.m_value_v, tmp.error_status] = Func_0011_get_area_weighted_mean(RCM_data_v.rmse, OBS_grid.cut_lon_rho, OBS_grid.cut_lat_rho);

    [tmp.indw, tmp.inde, tmp.inds, tmp.indn]=Func_0012_findind_Y(1/10,[127, 130, 37, 41],OBS_grid.cut_lon_rho,OBS_grid.cut_lat_rho); % southern EKB
    EKB_lon=RCM_grid.lon_rho(tmp.indw:tmp.inde,tmp.inds:tmp.indn);
    EKB_lat=RCM_grid.lat_rho(tmp.indw:tmp.inde,tmp.inds:tmp.indn);
    EKB_data_u=RCM_data_u.rmse(tmp.indw:tmp.inde,tmp.inds:tmp.indn);
    [EKB_mean_u, error_status] = Func_0011_get_area_weighted_mean(EKB_data_u, EKB_lon, EKB_lat);
    EKB_data_v=RCM_data_v.rmse(tmp.indw:tmp.inde,tmp.inds:tmp.indn);
    [EKB_mean_v, error_status] = Func_0011_get_area_weighted_mean(EKB_data_v, EKB_lon, EKB_lat);

    %% pattern corr
    tmp.mask=RCM_data_u.interp_mean.*OBS_data.mean_u;
    tmp.a=RCM_data_u.interp_mean(isfinite(tmp.mask));
    tmp.b=OBS_data.mean_u(isfinite(tmp.mask));
    RCM_data.pattern_corr_u=corrcoef(tmp.a(:), tmp.b(:));

    tmp.mask=RCM_data_v.interp_mean.*OBS_data.mean_v;
    tmp.a=RCM_data_v.interp_mean(isfinite(tmp.mask));
    tmp.b=OBS_data.mean_v(isfinite(tmp.mask));
    RCM_data.pattern_corr_v=corrcoef(tmp.a(:), tmp.b(:));


    
    m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);
    hold on;


    m_quiver(RCM_grid.cut_lon_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)', ...
                    RCM_grid.cut_lat_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)', ...
                    RCM_data_u.mean(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)' * param.m_quiver_vector_size, ...
                    RCM_data_v.mean(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)' * param.m_quiver_vector_size, ...
                    'color', param.m_quiver_vector_color, 'AutoScale','off','LineWidth', param.m_quiver_LineWidth);


    m_gshhs_i('color',param.m_gshhs_line_color)  
    m_gshhs_i('patch',param.m_gshhs_land_color);   % gray colored land
    
    m_text(param.m_quiver_ref_text_x_location, param.m_quiver_ref_text_y_location, param.m_quiver_ref_text, 'FontSize', param.m_quiver_ref_text_fontsize); 

    m_grid('fontsize', param.m_grid_fontsize, 'box', param.m_grid_box_type, 'tickdir', param.m_grid_tickdir_type);
    if min(RCM_info.years) == max(RCM_info.years)
        tmp.titlename = strcat('vec', ', ', RCM_info.season(1:3), ', ', tmp.abb,',(',num2str(max(RCM_info.years),'%04i'),') ');                        
    else
        tmp.titlename = strcat('vec', ', ', RCM_info.season(1:3), ', ', tmp.abb,',(',num2str(min(RCM_info.years),'%04i'),'-',num2str(max(RCM_info.years),'%04i'),') ');
    end
    title(tmp.titlename,'fontsize',param.m_pcolor_title_fontsize);  %%title

%     disp(['M = ', num2str(tmp.m_value)]);
%         m_text(param.m_pcolor_ref_text_x_location, param.m_pcolor_ref_text_y_location, ['M = ', num2str(tmp.m_value)], 'FontSize', param.m_quiver_ref_text_fontsize); 

    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperSize', [param.hor_paper_size_x, param.hor_paper_size_y]);
    set(gcf,'PaperPosition', [param.paper_position_hor param.paper_position_ver param.paper_position_width param.paper_position_height]) 
    saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
    close all;
    clear RCM_data.mean
    RCM_grid=rmfield(RCM_grid, 'lon_rho');
end
