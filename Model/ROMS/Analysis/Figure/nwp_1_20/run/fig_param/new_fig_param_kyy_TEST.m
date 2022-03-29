% % this script is for north western pacific region figure
% %  Updated 27-Apr-2018 by Yong-Yub Kim
% %  Updated 22-Jun-2018 by Yong-Yub Kim
% %  Updated 05-Jul-2021 by Yong-Yub Kim  (structure)

if (exist('refdepth' , 'var') ~= 1)
% % set part of the m_quiver parameter 
    param.m_quiver_vector_size = 3;
    param.m_quiver_ref_u_value = 1;
elseif (refdepth >= -100)
    param.m_quiver_vector_size = 1.5;
    param.m_quiver_ref_u_value = 2;
elseif (refdepth >= -300  && refdepth < -100)
    param.m_quiver_vector_size = 3;
    param.m_quiver_ref_u_value = 1;
elseif (refdepth < -300)
    param.m_quiver_vector_size = 6;
    param.m_quiver_ref_u_value = 0.1;
end

% % set colorbar parameter
% load C:\Users\User\Dropbox\source\matlab\Common\Figure\jet_mod  % % set colormap (jet_modified)
% run('C:\Users\User\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
param.colorbar_fontsize = 20;
param.colorbar_title_fontsize = 15;
param.colormap_style = jet;  % % default

% % set varname
if (isfield(tmp, 'variable') == 1)
    switch(tmp.variable)
        case('SST')
            param.varname = 'temp';
            param.colorbar_title = '(^oC)';
            param.colorbar_lev = [-2 33];
            param.colorbar_lev_meanplot = [12 22];
            param.colorbar_diff_lev = [0 6];
        case('BT')
            param.varname = 'temp';
            param.colorbar_title = '(^oC)';
            param.colorbar_lev = [-2 33];
        case('SSS')
            param.varname = 'salt';
            param.colorbar_title = ' ';
            param.colorbar_lev = [25 35];
        case('SSH')
            param.varname = 'zeta';
            param.colorbar_title = '(m)';
            param.colormap_style = bwrmap;
%             colorbar_lev = [-0.2 1.0];
            param.colorbar_lev = [0 1.3];
            param.colorbar_diff_lev = [0 0.8];
        case('vert_temp')
            param.varname = 'temp';
            param.colorbar_title = '(^oC)';
            param.colormap_style = jet;
        case('vert_salt')
            param.varname = 'salt';
            param.colorbar_title = ' ';
            param.colormap_style = jet;
        case('vert_u')
            param.varname = 'u';
            param.colorbar_title = '(m/s)';
            param.colormap_style = bwrmap;
        case('vert_v')
            param.varname = 'v';
            param.colorbar_title = '(m/s)';
            param.colormap_style = bwrmap;
        case('vert_rho')
            param.varname = 'zeta';
            param.colorbar_title = '(kg/m^3)';
        case('speed')
            param.varname = 'speed';
            param.colorbar_title = '(m/s)';
            param.colorbar_lev_meanplot = [0.03 0.18];
        case('H')
            param.varname = 'h';
            param.colorbar_title = '(m)';
            param.colorbar_lev = [-3000 0];
        otherwise
            param.varname = tmp.variable;
            param.colorbar_title = ' ';
    end
end
% % set calendar name
param.calendarname=cell(1,12); 
param.calendarname{1} = 'January'; 
param.calendarname{2} = 'February'; 
param.calendarname{3} = 'March'; 
param.calendarname{4} = 'April'; 
param.calendarname{5} = 'May'; 
param.calendarname{6} = 'June';
param.calendarname{7} = 'July'; 
param.calendarname{8} = 'August'; 
param.calendarname{9} = 'September'; 
param.calendarname{10} = 'October'; 
param.calendarname{11} = 'November'; 
param.calendarname{12} = 'December';

% % set m_proj parameter
param.m_proj_name='mercator';

% % set m_grid parameter
param.m_grid_fontsize =20;
param.m_grid_box_type = 'fancy';
param.m_grid_tickdir_type = 'in';

% % set m_pcolor parameter
param.m_pcolor_title_fontsize = 20;
param.m_pcolor_shading_method = 'flat';

% % set m_gshhs parameter
param.m_gshhs_line_color = 'k';
param.m_gshhs_land_color = [0.8 0.8 0.8]; % % gray

% % set m_contour parameter
param.m_contour_color = 'k';
param.m_contour_title_fontsize = 20;
param.m_contour_label_fontsize = 13;
param.m_contour_label_color = 'k';
param.m_contour_linewidth = 1.5;
param.m_contour_rotation = 0;
param.m_contour_labelspacing = 100000;
param.m_contour_fontweight = 'bold';

% % set m_quiver parameter
param.m_quiver_title_fontsize = 20;
param.m_quiver_x_interval = 10;
param.m_quiver_y_interval = 10;

param.m_quiver_vector_color = 'k';
param.m_quiver_LineWidth = 0.5;
param.m_quiver_AutoScale = 'off';
param.m_quiver_ref_v_value = param.m_quiver_ref_u_value/10000.0;
param.m_quiver_ref_text = [num2str(param.m_quiver_ref_u_value),' m/s'];
param.m_quiver_ref_text_fontsize = 15;
param.m_quiver_ref_text_x_location = 124.5;
param.m_quiver_ref_text_y_location = 45.6;


% % set gcf parameter
param.hor_paper_size_x = RCM_grid.domain(2)-RCM_grid.domain(1);
param.hor_paper_size_y = RCM_grid.domain(4)-RCM_grid.domain(3);
tmp.halt = 1;
while(tmp.halt)
    if (param.hor_paper_size_x > 500 || param.hor_paper_size_y > 500)
        tmp.halt = 0;
    else
        param.hor_paper_size_x = param.hor_paper_size_x * 1.2; 
        param.hor_paper_size_y = param.hor_paper_size_y * 1.2;
    end
end

param.paper_position_hor = 0; % % distance from left
param.paper_position_ver = 0; % % distance from bottom
param.paper_position_width = param.hor_paper_size_x ;
param.paper_position_height = param.hor_paper_size_y ;
param.vert_paper_size_y = 400;