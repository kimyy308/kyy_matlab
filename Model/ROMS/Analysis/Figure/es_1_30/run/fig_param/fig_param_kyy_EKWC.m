% % this script is for EKWC region figure
% %  Updated 27-Apr-2018 by Yong-Yub Kim
% %  Updated 22-Jun-2018 by Yong-Yub Kim

if (exist('refdepth' , 'var') ~= 1)
% % set part of the m_quiver parameter 
    m_quiver_vector_size = 1;
    m_quiver_ref_u_value = 0.5;
elseif (refdepth >= -100)
    m_quiver_vector_size = 1;
    m_quiver_ref_u_value = 0.5;
elseif (refdepth >= -300  && refdepth < -100)
    m_quiver_vector_size = 2;
    m_quiver_ref_u_value = 0.2;
elseif (refdepth < -300)
    m_quiver_vector_size = 6;
    m_quiver_ref_u_value = 0.1;
end

% % set colorbar parameter
load C:\Users\KYY\Dropbox\source\matlab\Common\Figure\jet_mod  % % set colormap (jet_modified)
run('C:\Users\kyy\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
% for i = 1:3
%     bwrmap2(1:9901,i)=interp1((1:100)',bwrmap(1:100,i),(1:0.01:100)');
% end
colorbar_fontsize = 15;
colorbar_title_fontsize = 15;
colormap_style = jet_mod;  % % default

% % set varname
switch(var)
    case('SST')
        varname = 'temp'
        colorbar_title = '(^oC)';
    case('YSBCW')
        varname = 'temp'
        colorbar_title = '(^oC)';
    case('SSS')
        varname = 'salt'
        colorbar_title = ' ';
    case('SSH')
        varname = 'zeta'
        colorbar_title = '(m)';
        colormap_style = bwrmap;
    case('vert_temp')
        varname = 'temp'
        colorbar_title = '(^oC)';
%         colormap_style = bwrmap2;
        colormap_style = jet;
    case('vert_salt')
        varname = 'salt'
        colorbar_title = ' ';
        colormap_style = jet;
    case('vert_u')
        varname = 'u'
        colorbar_title = '(m/s)';
        colormap_style = bwrmap;
    case('vert_v')
        varname = 'v'
        colorbar_title = '(m/s)';
        colormap_style = bwrmap;
    case('vert_rho')
        varname = 'zeta'
        colorbar_title = '(kg/m^3)';
    case('vorticity')
        varname = 'relative vorticity'
        colorbar_title = '(s^-^1)';    
        colormap_style = bwrmap;
    otherwise
        varname = var;
        colorbar_title = ' ';
%         colormap_style =bwrmap2;  
end

% % set calendar name
calendarname=cell(1,12); calendarname{1} = 'January'; calendarname{2} = 'February'; calendarname{3} = 'March'; calendarname{4} = 'April'; calendarname{5} = 'May'; calendarname{6} = 'June';
calendarname{7} = 'July'; calendarname{8} = 'August'; calendarname{9} = 'September'; calendarname{10} = 'October'; calendarname{11} = 'November'; calendarname{12} = 'December';

% % set m_proj parameter
m_proj_name='mercator';

% % set m_grid parameter
m_grid_fontsize = 15;
m_grid_box_type = 'fancy';
m_grid_tickdir_type = 'in';

% % set m_pcolor parameter
m_pcolor_title_fontsize = 15;
m_pcolor_shading_method = 'interp';

% % set m_gshhs parameter
m_gshhs_line_color = 'k';
m_gshhs_land_color = [0.8 0.8 0.8]; % % gray
% %  Updated 27-Apr-2018 by Yong-Yub Kim


% % set m_contour parameter
m_contour_color = 'k';
m_contour_title_fontsize = 15;
m_contour_label_fontsize = 10;
m_contour_label_color = 'k';
m_contour_linewidth = 1.5;
m_contour_rotation = 0;
m_contour_labelspacing = 100000;
m_contour_fontweight = 'bold';

% % set m_quiver parameter
m_quiver_title_fontsize = 15;
m_quiver_x_interval = 3;
m_quiver_y_interval = 3;
m_quiver_vector_size = 0.5;
m_quiver_vector_color = 'k';
m_quiver_LineWidth = 0.5;
m_quiver_AutoScale = 'off';
m_quiver_ref_u_value = 0.5;
m_quiver_ref_v_value = 0.0001;
m_quiver_hor_ref_vec_x_range = 20:23;
m_quiver_hor_ref_vec_y_range = 72:73;
m_quiver_ref_vec_x_range = 281:285;
m_quiver_ref_vec_y_range = 468:470;
m_quiver_ref_text = '0.5 m/s';
m_quiver_ref_text_fontsize = 9;
m_quiver_ref_text_x_location = 129.1;
m_quiver_ref_text_y_location = 35.8;

% % set vert_grid parameter
vert_grid_fontsize = 15;
vert_grid_tickdir_type = 'in';
vert_grid_box = 'on';
vert_grid_box_linewidth = 1.5;

% % set vert_pcolor parameter
vert_pcolor_title_fontsize = 15;
vert_pcolor_shading_method = 'interp';

% % set vert_contour parameter
vert_contour_color = 'k';
vert_contour_title_fontsize = 15;
vert_contour_label_fontsize = 15;
vert_contour_label_color = 'k';
vert_contour_linewidth = 2;
vert_contour_rotation = 0;
vert_contour_labelspacing = 100000;
vert_contour_fontweight = 'bold';

% % set gcf parameter
if (exist('lonlat' , 'var') == 1)
    hor_paper_size_x = lonlat(2)-lonlat(1);
    hor_paper_size_y = lonlat(4)-lonlat(3);
    halt = 1;
    while(halt)
        if (hor_paper_size_x > 500 || hor_paper_size_y > 500)
            halt = 0;
        else
            hor_paper_size_x = hor_paper_size_x * 1.2; 
            hor_paper_size_y = hor_paper_size_y * 1.2;
        end
    end

    paper_position_width = hor_paper_size_x ;
    paper_position_height = hor_paper_size_y ;
end

paper_position_hor = 0; % % distance from left
paper_position_ver = 0; % % distance from bottom
vert_paper_size_y = 400;
