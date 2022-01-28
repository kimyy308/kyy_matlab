% start-------------------- decadal current plot
   tmp.variable='UV';
    dirs.figdir=[dirs.figrawdir,'surface', tmp.fs, tmp.regionname, tmp.fs, 'vec', tmp.fs, ...
        num2str(min(RCM_info.years)), '_', num2str(max(RCM_info.years)), tmp.fs];
    if (exist(strcat(dirs.figdir) , 'dir') ~= 7)
        mkdir(strcat(dirs.figdir));
    end 
    tmp.tifname=strcat(dirs.figdir, 'diff_RCM_', RCM_info.ensname, '_clim_uv_',num2str(min(RCM_info.years),'%04i'), ...
    '_',num2str(max(RCM_info.years),'%04i'), ...
                    '_', num2str(RCM_info.months(1),'%04i'), '-', num2str(RCM_info.months(end),'%04i'),'.tif'); %% ~_year_month.jpg
    
    tmp.early_enssavefilename = [dirs.enssavedir, RCM_info.ensname, '_', tmp.regionname, '_', tmp.variable, ...
        '_mean_', num2str(RCM_info.inputyear1(1),'%04i'), '-', num2str(RCM_info.inputyear1(end),'%04i'),...
                    '_', num2str(RCM_info.months(1),'%04i'), '-', num2str(RCM_info.months(end),'%04i'), '.mat'];
    load(tmp.early_enssavefilename)
    early_u_rho=u_rho;
    early_v_rho=v_rho;
    
    tmp.late_enssavefilename = [dirs.enssavedir, RCM_info.ensname, '_', tmp.regionname, '_', tmp.variable, ...
        '_mean_', num2str(RCM_info.inputyear2(1),'%04i'), '-', num2str(RCM_info.inputyear2(end),'%04i'), ...
                    '_', num2str(RCM_info.months(1),'%04i'), '-', num2str(RCM_info.months(end),'%04i'), '.mat'];
    load(tmp.late_enssavefilename)
    late_u_rho=u_rho;
    late_v_rho=v_rho;
    
    u_rho = late_u_rho - early_u_rho;
    v_rho = late_v_rho - early_v_rho;
    
    if (exist('ref_vec_x_range' , 'var') ~= 1)
        ref_vec_x_ind = find(abs(cut_lon_rho(:,1)-param.m_quiver_ref_text_x_location) == min(abs(cut_lon_rho(:,1)-param.m_quiver_ref_text_x_location)));
        ref_vec_y_ind = find(abs(cut_lat_rho(1,:)-param.m_quiver_ref_text_y_location) == min(abs(cut_lat_rho(1,:)-param.m_quiver_ref_text_y_location)))+param.m_quiver_y_interval*2;
        ref_vec_x_range = round(ref_vec_x_ind-(param.m_quiver_x_interval/2)) : round(ref_vec_x_ind-(param.m_quiver_x_interval/2))+param.m_quiver_x_interval;
        ref_vec_y_range = round(ref_vec_y_ind-(param.m_quiver_y_interval/2)) : round(ref_vec_y_ind-(param.m_quiver_y_interval/2))+param.m_quiver_y_interval;
    end
    u_rho(ref_vec_x_range,ref_vec_y_range)=0.5;
    v_rho(ref_vec_x_range,ref_vec_y_range)=0.00001;     

    m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);
    hold on;
    m_gshhs_i('color',param.m_gshhs_line_color)  
    m_gshhs_i('patch',param.m_gshhs_land_color);   % gray colored land

    uvplot=m_quiver(cut_lon_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)', ...
                cut_lat_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)', ...
                u_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)' * param.m_quiver_vector_size, ...
                v_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)' * param.m_quiver_vector_size, ...
                'color', param.m_quiver_vector_color, 'AutoScale','off','LineWidth', param.m_quiver_LineWidth);


    m_text(param.m_quiver_ref_text_x_location, param.m_quiver_ref_text_y_location, '0.5 m/s', 'FontSize', param.m_quiver_ref_text_fontsize); 
    m_grid('fontsize', param.m_grid_fontsize, 'box', param.m_grid_box_type, 'tickdir', param.m_grid_tickdir_type);
    tmp.titlename = strcat('UV diff, ', season(1:3), ', ', RCM_info.ensname,',(',num2str(RCM_info.inputyear2(1),'%04i'),'-',num2str(RCM_info.inputyear1(end),'%04i'),') ');  %% + glacier contribution
    title(tmp.titlename,'fontsize',param.m_pcolor_title_fontsize);  %%title

    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperSize', [param.hor_paper_size_x, param.hor_paper_size_y]);
    set(gcf,'PaperPosition', [param.paper_position_hor param.paper_position_ver param.paper_position_width param.paper_position_height]) 
    saveas(gcf, tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
    close all;
    clear lon_rho ens_u_rho ref_vec_x_range