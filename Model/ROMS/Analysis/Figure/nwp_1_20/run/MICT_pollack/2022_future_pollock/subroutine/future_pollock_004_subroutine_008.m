% %  Updated 21-May-2021 by Yong-Yub Kim, 


RCM_info.years_ssp =RCM_info.years;
tmp.testname_ssp = tmp.testname;
dirs.filedir_ssp = strcat('D:\Data\Model\ROMS\', RCM_info.model, '\backup_surf\', tmp.testname_ssp, '\run\'); % % where data files are          
dirs.matdir_ssp = strcat('D:\Data\Model\ROMS\', RCM_info.model, '\', tmp.testname_ssp, '\run\mean\');

[tmp.testname_his, tmp.error_status] = Func_0023_RCM_CMIP6_testname_his(tmp.testname_ssp);
dirs.filedir_his = strcat('D:\Data\Model\ROMS\', RCM_info.model, '\backup_surf\', tmp.testname_his, '\run\'); % % where data files are          
dirs.matdir_his = strcat('D:\Data\Model\ROMS\', RCM_info.model, '\', tmp.testname_his, '\run\mean\');


%% set variable name & figure directory
tmp.variable='wind';

dirs.figdir=[dirs.figrawdir,'surface', tmp.fs, tmp.regionname, tmp.fs, tmp.variable, tmp.fs, ...
    num2str(min(RCM_info.years)), '_', num2str(max(RCM_info.years)), tmp.fs];
if (exist(strcat(dirs.figdir) , 'dir') ~= 7)
    mkdir(strcat(dirs.figdir));
end 

%% set figure file name
tmp.tifname=strcat(dirs.figdir, 'diff_', tmp.testname_ssp, '_', tmp.testname_his, ...
    '_surf_', tmp.variable, '_', num2str(min(RCM_info.years),'%04i'), '_', ...
    num2str(max(RCM_info.years),'%04i'), '_', RCM_info.season,'.tif'); %% ~_year_month.jpg
if (exist(tmp.tifname , 'file') ~= 2 || fig_flag==2)      
    run(tmp.param_script);

   %% set data file name (mat) and reading (ssp)
    tmp.matname_Uwind_ssp = [dirs.matdir_ssp, tmp.testname_ssp, '_', tmp.regionname, '_', 'Uwind',...
        '_mean_', num2str(min(RCM_info.years_ssp),'%04i'), '-', num2str(max(RCM_info.years_ssp),'%04i'),...
    '_', RCM_info.season, '.mat'];
    if (exist(tmp.matname_Uwind_ssp , 'file') ~= 2)
        disp('please get data from subroutine_004_001 first')
    else
        load(tmp.matname_Uwind_ssp, 'RCM_data', 'RCM_grid');
    end 
    RCM_data_Uwind_ssp = RCM_data;
    
    tmp.matname_Vwind_ssp = [dirs.matdir_ssp, tmp.testname_ssp, '_', tmp.regionname, '_', 'Vwind',...
        '_mean_', num2str(min(RCM_info.years_ssp),'%04i'), '-', num2str(max(RCM_info.years_ssp),'%04i'),...
    '_', RCM_info.season, '.mat'];
    if (exist(tmp.matname_Vwind_ssp , 'file') ~= 2)
        disp('please get data from subroutine_004_001 first')
    else
        load(tmp.matname_Vwind_ssp, 'RCM_data', 'RCM_grid');
    end 
    RCM_data_Vwind_ssp = RCM_data;

   %% set data file name (mat) and reading (his)
    tmp.matname_Uwind_his = [dirs.matdir_his, tmp.testname_his, '_', tmp.regionname, '_', 'Uwind',...
        '_mean_', num2str(min(RCM_info.years_his),'%04i'), '-', num2str(max(RCM_info.years_his),'%04i'),...
    '_', RCM_info.season, '.mat'];
    if (exist(tmp.matname_Uwind_his , 'file') ~= 2)
        disp('please get data from subroutine_004_001 first')
    else
        load(tmp.matname_Uwind_his, 'RCM_data', 'RCM_grid');
    end 
    RCM_data_Uwind_his = RCM_data;
    
    tmp.matname_Vwind_his = [dirs.matdir_his, tmp.testname_his, '_', tmp.regionname, '_', 'Vwind',...
        '_mean_', num2str(min(RCM_info.years_his),'%04i'), '-', num2str(max(RCM_info.years_his),'%04i'),...
    '_', RCM_info.season, '.mat'];
    if (exist(tmp.matname_Vwind_his , 'file') ~= 2)
        disp('please get data from subroutine_004_001 first')
    else
        load(tmp.matname_Vwind_his, 'RCM_data', 'RCM_grid');
    end 
    RCM_data_Vwind_his = RCM_data;
    
    RCM_data_Uwind.mean=RCM_data_Uwind_ssp.mean - RCM_data_Uwind_his.mean;
    RCM_data_Vwind.mean=RCM_data_Vwind_ssp.mean - RCM_data_Vwind_his.mean;
    
    RCM_grid.mask_model = double(inpolygon(RCM_grid.cut_lon_rho,RCM_grid.cut_lat_rho,RCM_grid.refpolygon(:,1),RCM_grid.refpolygon(:,2)));
    RCM_grid.mask_model(RCM_grid.mask_model==0)=NaN;
    RCM_data_Uwind.mean=RCM_data_Uwind.mean.*RCM_grid.mask_model;
    RCM_data_Vwind.mean=RCM_data_Vwind.mean.*RCM_grid.mask_model;
    
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
    RCM_data_Uwind.mean(tmp.ref_vec_x_range,tmp.ref_vec_y_range)=param.m_quiver_ref_uwind_value;
    RCM_data_Vwind.mean(tmp.ref_vec_x_range,tmp.ref_vec_y_range)=param.m_quiver_ref_vwind_value;     
    
    
    m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);
    hold on;


    m_quiver(RCM_grid.cut_lon_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)', ...
                    RCM_grid.cut_lat_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)', ...
                    RCM_data_Uwind.mean(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)' * param.m_quiver_wind_vector_size, ...
                    RCM_data_Vwind.mean(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)' * param.m_quiver_wind_vector_size, ...
                    'color', param.m_quiver_vector_color, 'AutoScale','off','LineWidth', param.m_quiver_LineWidth);


    m_gshhs_i('color',param.m_gshhs_line_color)  
    m_gshhs_i('patch',param.m_gshhs_land_color);   % gray colored land
    
    m_text(param.m_quiver_ref_text_x_location, param.m_quiver_ref_text_y_location, param.m_quiver_wind_ref_text, 'FontSize', param.m_quiver_ref_text_fontsize); 

    m_grid('fontsize', param.m_grid_fontsize, 'box', param.m_grid_box_type, 'tickdir', param.m_grid_tickdir_type);
    tmp.titlename = strcat('vec', ', ', RCM_info.season(1:3), ', ', tmp.abb,',(','diff',') ');                        
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
    
    disp([tmp.testname, ', future_pollock_004_subroutine_007']);
end


        
        
