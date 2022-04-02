disp('subroutine SSH_2p_2nd_010_sub_004_plot_bndy_RCM_GLORYS') 

tmp.matsavefile=[dirs.matdir, filesep, tmp.testname, '_', tmp.variable, '_bndy_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), '_', ...
    RCM_info.season, '.mat'];
load(tmp.matsavefile, 'RCM_all', 'RCM_all_bndy', 'RCM_bndy_mean', 'RCM_grid', 'RCM_eval')

tmp.gcmmatsavefile=[dirs.glorysdir, filesep, 'GLORYS', '_', tmp.variable_GCM, '_bndy_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), '_', ...
    RCM_info.season, '.mat'];
load(tmp.gcmmatsavefile, 'GCM_all', 'GCM_all_bndy', 'GCM_bndy_mean', 'GCM_grid')

abc=1;

for bndyij=1:length(RCM_info.bndy_directions)
    %% figdir configuration
    tmp.direction=RCM_info.bndy_directions{bndyij};

    dirs.figdir=[dirs.figrawdir, tmp.fs, 'bndy', tmp.fs, tmp.variable, tmp.fs, ...
        num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), tmp.fs, tmp.direction];
    if (exist(strcat(dirs.figdir) , 'dir') ~= 7)
        mkdir(strcat(dirs.figdir));
    end 
    
    %% nograd(basic configuration)
    tmp.tifname=[dirs.figdir, tmp.fs, tmp.testname, '_', tmp.variable, '_', ...
        num2str(min(RCM_info.years),'%04i'), '_',num2str(max(RCM_info.years),'%04i'), ...
        '_', tmp.direction, '_alldep_nograd_', RCM_info.season, '.tif'];
    switch tmp.direction
        case {'north', 'south'}
            tmp.lonlat='lon';
        case {'west', 'east'}
            tmp.lonlat='lat';
    end
    
    run(tmp.param_script);

    hold on
    pcolor(RCM_grid.(['lonlat_', tmp.direction]), RCM_grid.(['stddepth_', tmp.direction]), ...
        RCM_bndy_mean.([tmp.variable, '_stddepth_', tmp.direction])); 
    shading(gca,param.m_pcolor_shading_method);   
    switch tmp.variable
        case 'salt'
            [C,h2]=contour(RCM_grid.(['lonlat_', tmp.direction]), RCM_grid.(['stddepth_', tmp.direction]), ...
                    RCM_bndy_mean.([tmp.variable, '_stddepth_', tmp.direction]), [30:0.5:35], 'color','k', ...
                'linewidth', 1.5, 'linestyle', '-');
            clabel(C,h2,'FontSize',13,'Color','k', ...
                'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
        case 'temp'
            [C,h2]=contour(RCM_grid.(['lonlat_', tmp.direction]), RCM_grid.(['stddepth_', tmp.direction]), ...
                    RCM_bndy_mean.([tmp.variable, '_stddepth_', tmp.direction]), [0:5:35], 'color','k', ...
                'linewidth', 1.5, 'linestyle', '-');
            clabel(C,h2,'FontSize',13,'Color','k', ...
                'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
        otherwise
    end
    set(gca, 'box', 'on', 'fontsize', param.m_grid_fontsize);
    xlabel(tmp.lonlat)
    ylabel('depth')
    if min(RCM_info.years) == max(RCM_info.years)
        tmp.titlename = strcat(tmp.variable, ', ', RCM_info.season(1:3), ', ', tmp.abb,',(',num2str(max(RCM_info.years),'%04i'),') ');                        
    else
        tmp.titlename = strcat(tmp.variable, ', ', RCM_info.season(1:3), ', ', tmp.abb,',(',num2str(min(RCM_info.years),'%04i'),'-',num2str(max(RCM_info.years),'%04i'),') ');
    end
    title(tmp.titlename,'fontsize',param.m_pcolor_title_fontsize);  %%title
    h = colorbar;
    colormap(jet);
    set(h,'fontsize',param.colorbar_fontsize);
    title(h,param.colorbar_title,'fontsize',param.colorbar_title_fontsize);
    
%     disp(['M = ', num2str(tmp.m_value)]);
% param.text_pos_x = (RCM_info.vert_section(2)-RCM_info.vert_section(1))/20+RCM_info.vert_section(1);
% param.text_pos_y = (RCM_info.vert_section(6)-RCM_info.vert_section(5))/20+RCM_info.vert_section(5);
% 
% text(param.text_pos_x, param.text_pos_y, ['M = ', num2str(round(tmp.m_value,2))], 'FontSize', param.m_quiver_ref_text_fontsize); 
    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperSize', [param.hor_paper_size_x, param.hor_paper_size_y]);
    set(gcf,'PaperPosition', [param.paper_position_hor param.paper_position_ver param.paper_position_width param.paper_position_height]) 
    saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
%     close all;
    
    switch tmp.variable
        case {'temp', 'salt', 'rho'}
            caxis(param.colorbar_lev);
        case {'vert_u', 'vert_v'}
            caxis([-max(abs(mean_data(:))), max(abs(mean_data(:)))]);
            colormap(cmaps.byrmap);
        case{'vert_temp_nogrd', 'vert_salt_nogrd', 'vert_rho_nogrd'}
            colormap(jet)
    end
    %% fixed caxis configuration
    tmp.tifname=[dirs.figdir, tmp.fs, tmp.testname, '_', tmp.variable, '_', ...
        num2str(min(RCM_info.years),'%04i'), '_',num2str(max(RCM_info.years),'%04i'), ...
        '_', tmp.direction, '_alldep_', RCM_info.season, '.tif'];
    saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
    hold off
    close all;
    
    %% 0~refdep configuration
    tmp.refdep =200;
    tmp.refdep_str=num2str(tmp.refdep);
    tmp.tifname=[dirs.figdir, tmp.fs, tmp.testname, '_', tmp.variable, '_', ...
        num2str(min(RCM_info.years),'%04i'), '_',num2str(max(RCM_info.years),'%04i'), ...
        '_', tmp.direction, '_0-', tmp.refdep_str, 'm_nograd_', RCM_info.season, '.tif'];
    hold on
    tmp.ind_dep=find(RCM_grid.stddepth == -tmp.refdep); %% 1:21 -> 0:200m
    
    pcolor(RCM_grid.(['lonlat_', tmp.direction])(:,1:tmp.ind_dep), RCM_grid.(['stddepth_', tmp.direction])(:,1:tmp.ind_dep), ...
        RCM_bndy_mean.([tmp.variable, '_stddepth_', tmp.direction])(:,1:tmp.ind_dep)); 
    shading(gca,param.m_pcolor_shading_method);   
    switch tmp.variable
        case 'salt'
            [C,h2]=contour(RCM_grid.(['lonlat_', tmp.direction])(:,1:tmp.ind_dep), RCM_grid.(['stddepth_', tmp.direction])(:,1:tmp.ind_dep), ...
                    RCM_bndy_mean.([tmp.variable, '_stddepth_', tmp.direction])(:,1:tmp.ind_dep), [30:0.5:35], 'color','k', ...
                'linewidth', 1.5, 'linestyle', '-');
            clabel(C,h2,'FontSize',13,'Color','k', ...
                'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
        case 'temp'
            [C,h2]=contour(RCM_grid.(['lonlat_', tmp.direction])(:,1:tmp.ind_dep), RCM_grid.(['stddepth_', tmp.direction])(:,1:tmp.ind_dep), ...
                    RCM_bndy_mean.([tmp.variable, '_stddepth_', tmp.direction])(:,1:tmp.ind_dep), [0:5:35], 'color','k', ...
                'linewidth', 1.5, 'linestyle', '-');
            clabel(C,h2,'FontSize',13,'Color','k', ...
                'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
        otherwise
    end
    set(gca, 'box', 'on', 'fontsize', param.m_grid_fontsize);
    xlabel(tmp.lonlat)
    ylabel('depth')
    if min(RCM_info.years) == max(RCM_info.years)
        tmp.titlename = strcat(tmp.variable, ', ', RCM_info.season(1:3), ', ', tmp.abb,',(',num2str(max(RCM_info.years),'%04i'),') ');                        
    else
        tmp.titlename = strcat(tmp.variable, ', ', RCM_info.season(1:3), ', ', tmp.abb,',(',num2str(min(RCM_info.years),'%04i'),'-',num2str(max(RCM_info.years),'%04i'),') ');
    end
    title(tmp.titlename,'fontsize',param.m_pcolor_title_fontsize);  %%title
    h = colorbar;
    colormap(jet);
    set(h,'fontsize',param.colorbar_fontsize);
    title(h,param.colorbar_title,'fontsize',param.colorbar_title_fontsize);

    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperSize', [param.hor_paper_size_x, param.hor_paper_size_y]);
    set(gcf,'PaperPosition', [param.paper_position_hor param.paper_position_ver param.paper_position_width param.paper_position_height]) 
    saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
    
    switch tmp.variable
        case {'temp', 'salt', 'rho'}
            caxis(param.colorbar_lev);
        case {'vert_u', 'vert_v'}
            caxis([-max(abs(mean_data(:))), max(abs(mean_data(:)))]);
            colormap(cmaps.byrmap);
        case{'vert_temp_nogrd', 'vert_salt_nogrd', 'vert_rho_nogrd'}
            colormap(jet)
    end
    %% fixed caxis configuration
    tmp.tifname=[dirs.figdir, tmp.fs, tmp.testname, '_', tmp.variable, '_', ...
        num2str(min(RCM_info.years),'%04i'), '_',num2str(max(RCM_info.years),'%04i'), ...
        '_', tmp.direction, '_0-', tmp.refdep_str, 'm_', RCM_info.season, '.tif'];
    saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
    hold off
    close all;
    
end
% [tmp.m_value, tmp.error_status] = ...
%     Func_0021_get_vert_area_weighted_mean(mean_data, RCM_vert_grid.cut_lon_rho2, RCM_vert_grid.cut_lat_rho2, ...
%     RCM_vert_grid.cut_thick2, RCM_vert_grid.lonlatflag);

% hold on;
% pcolor(RCM_vert_grid.lonlat,RCM_vert_grid.cut_depth2, squeeze(mean_data));
% 
% shading(gca,param.m_pcolor_shading_method);   
% if(strcmp(tmp.variable, 'vert_temp'))
%     [C,h2]=contour(RCM_vert_grid.lonlat, RCM_vert_grid.cut_depth2, squeeze(mean_data), [5, 10, 15, 20, 25, 30], 'color','k', ...
%             'linewidth', 1.5, 'linestyle', '-');
%         clabel(C,h2,'FontSize',13,'Color','k', ...
%             'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
% end
% set(gca, 'box', 'on', 'fontsize', param.m_grid_fontsize);
%         m_grid('fontsize', param.m_grid_fontsize, 'box', param.m_grid_box_type, 'tickdir', param.m_grid_tickdir_type);
% if min(RCM_info.years) == max(RCM_info.years)
%     tmp.titlename = strcat(tmp.variable(6:end), ', ', RCM_info.season(1:4), ', ', tmp.abb,',(',num2str(max(RCM_info.years),'%04i'),') ');                        
% else
%     tmp.titlename = strcat(tmp.variable(6:end), ', ', RCM_info.season(1:4), ', ', tmp.abb,',(',num2str(min(RCM_info.years),'%04i'),'-',num2str(max(RCM_info.years),'%04i'),') ');
% end
% title(tmp.titlename,'fontsize',param.m_pcolor_title_fontsize);  %%title

% % set colorbar 
% h = colorbar;
% %                     if (strcmp(tmp.variable,'SSH')==1)
% %                         colormap(flip(cool));
% %                     else
%     colormap(jet);
% %                     end

% set(h,'fontsize',param.colorbar_fontsize);
% title(h,param.colorbar_title,'fontsize',param.colorbar_title_fontsize);
% switch tmp.variable
%     case {'vert_temp', 'vert_salt', 'vert_rho'}
%         caxis(param.colorbar_lev);
%     case {'vert_u', 'vert_v'}
%         caxis([-max(abs(mean_data(:))), max(abs(mean_data(:)))]);
%         colormap(cmaps.byrmap);
%     case{'vert_temp_nogrd', 'vert_salt_nogrd', 'vert_rho_nogrd'}
%         colormap(jet)
% %                 caxis([min(mean_data(:)) max(mean_data(:))])
% %                 max(abs(mean_data(:)))-mean(abs(mean_data(:)), 'omitnan')
% end

% disp(['M = ', num2str(tmp.m_value)]);
% param.text_pos_x = (RCM_info.vert_section(2)-RCM_info.vert_section(1))/20+RCM_info.vert_section(1);
% param.text_pos_y = (RCM_info.vert_section(6)-RCM_info.vert_section(5))/20+RCM_info.vert_section(5);
% 
% text(param.text_pos_x, param.text_pos_y, ['M = ', num2str(round(tmp.m_value,2))], 'FontSize', param.m_quiver_ref_text_fontsize); 
% 
% set(gcf, 'PaperUnits', 'points');
% set(gcf, 'PaperSize', [param.hor_paper_size_x, param.hor_paper_size_y]);
% set(gcf,'PaperPosition', [param.paper_position_hor param.paper_position_ver param.paper_position_width param.paper_position_height]) 
% saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
% close all;