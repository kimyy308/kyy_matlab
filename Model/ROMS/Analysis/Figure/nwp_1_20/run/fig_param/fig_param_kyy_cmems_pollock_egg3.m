% % this script is for pollock_egg3 region figure
% %  Updated 15-Jun-2023 by Yong-Yub Kim

param.m_quiver_vector_size = 3;
param.m_quiver_ref_u_value = 0.2;


% % set colorbar parameter
% % set colorbar parameter
[cmap.bwr_map, error_status] = Func_0009_get_colormaps('bwr', tmp.dropboxpath);
[cmap.jet_mod, error_status] = Func_0009_get_colormaps('jet_mod', tmp.dropboxpath);
[cmap.byr_map, error_status] = Func_0009_get_colormaps('byr', tmp.dropboxpath);

colorbar_fontsize = 20;
colorbar_title_fontsize = 15;
colormap_style = cmap.jet_mod;  % % default

% % set varname
if (isfield(tmp, 'variable') == 1)
    switch(tmp.variable)
        case('SST')
            param.varname = 'temp';
            param.colorbar_title = '(^oC)';
            param.colorbar_lev = [-2 33];
            param.colorbar_lev_meanplot = [12 22];
            param.colorbar_diff_lev = [0 6];
            if isfield(RCM_info, 'season')
                switch RCM_info.season
                    case 'JF-'
                        param.colorbar_lev = [-2 20];
                    case 'February'
                        param.colorbar_lev = [-2 20];
                    case 'January'
                        param.colorbar_lev = [-2 20];
                end
            end
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
            param.colorbar_lev = [-0.5 0.4];
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
            param.colormap_style = cmap.bwr_map;
        case('vert_v')
            param.varname = 'v';
            param.colorbar_title = '(m/s)';
            param.colormap_style = cmap.bwr_map;
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
        case('Uwind')
            param.varname = 'Uwind';
            param.colorbar_title = '(m/s)';
%             param.colorbar_lev = [-4 4]; % annual
            param.colorbar_lev = [0 10]; % annual
            param.colormap_style = cmap.bwr_map;
        case('Vwind')
            param.varname = 'Vwind';
            param.colorbar_title = '(m/s)';
%             param.colorbar_lev = [-4 4]; % annual
            param.colorbar_lev = [-10 0]; % annual
            param.colormap_style = cmap.bwr_map;
        case('wstrcurl')
            param.varname = 'Vwind';
            param.colorbar_title = '/s';
            param.colorbar_lev = [-2e-7 2e-7]; % annual
            param.colormap_style = cmap.bwr_map;
        case('swrad')
            param.varname = 'swrad';
            param.colorbar_title = 'W/m^2';
            param.colorbar_lev = [0 200]; % winter
            param.colormap_style = colormap(flip(autumn));
        case('shflux')
            param.varname = 'shflux';
            param.colorbar_title = 'W/m^2';
            param.colorbar_lev = [-300 0]; % winter
            param.colormap_style = colormap(cool);
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
param.m_grid_fontsize = 20;
param.m_grid_box_type = 'fancy';
param.m_grid_tickdir_type = 'in';

% % set m_pcolor parameter
param.m_pcolor_title_fontsize = 20;
param.m_pcolor_shading_method = 'flat';
param.m_pcolor_ref_text_x_location = 127.5;
param.m_pcolor_ref_text_y_location = 37.5;

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
param.m_quiver_x_interval = 1;
param.m_quiver_y_interval = 1;

param.m_quiver_vector_color = 'k';
param.m_quiver_LineWidth = 0.5;
param.m_quiver_AutoScale = 'off';
param.m_quiver_ref_uwind_value = 1;
param.m_quiver_ref_vwind_value = param.m_quiver_ref_uwind_value/10000.0;
param.m_quiver_ref_v_value = param.m_quiver_ref_u_value/10000.0;
param.m_quiver_ref_text = [num2str(param.m_quiver_ref_u_value),' m/s'];
param.m_quiver_wind_ref_text = [num2str(param.m_quiver_ref_uwind_value),' m/s'];
param.m_quiver_ekman_ref_text = [num2str(param.m_quiver_ref_u_value),' m^3/s'];
param.m_quiver_ref_text_fontsize = 15;
param.m_quiver_ref_text_x_location = 127.5;
param.m_quiver_ref_text_y_location = 37.5;

% % set vert_grid parameter
param.vert_grid_fontsize = 15;
param.vert_grid_tickdir_type = 'in';
param.vert_grid_box = 'on';
param.vert_grid_box_linewidth = 1.5;

% % set vert_pcolor parameter
param.vert_pcolor_title_fontsize = 15;
param.vert_pcolor_shading_method = 'interp';

% % set vert_contour parameter
param.vert_contour_color = 'k';
param.vert_contour_title_fontsize = 15;
param.vert_contour_label_fontsize = 10;
param.vert_contour_label_color = 'k';
param.vert_contour_linewidth = 1;
param.vert_contour_rotation = 0;
param.vert_contour_labelspacing = 100000;
param.vert_contour_fontweight = 'bold';


% % set gcf parameter
param.hor_paper_size_x = RCM_grid.domain(2)-RCM_grid.domain(1);
param.hor_paper_size_y = RCM_grid.domain(4)-RCM_grid.domain(3);
param.halt = 1;
while(param.halt)
    if (param.hor_paper_size_x > 500 || param.hor_paper_size_y > 500)
        param.halt = 0;
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