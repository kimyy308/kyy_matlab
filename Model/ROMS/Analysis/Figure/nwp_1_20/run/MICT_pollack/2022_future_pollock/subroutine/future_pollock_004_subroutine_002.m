% %  Updated 26-Apr-2021 by Yong-Yub Kim, 

% start-------------------- earlier decadal SST, SSS plot
for varind2=1:length(RCM_info.vars)
%% set variable name & figure directory
    tmp.variable=RCM_info.vars{varind2};
    dirs.figdir=[dirs.figrawdir,'surface', tmp.fs, tmp.regionname, tmp.fs, tmp.variable, tmp.fs, ...
        num2str(min(RCM_info.years)), '_', num2str(max(RCM_info.years)), tmp.fs];
    if (exist(strcat(dirs.figdir) , 'dir') ~= 7)
        mkdir(strcat(dirs.figdir));
    end 

%% set figure file name
    tmp.tifname=strcat(dirs.figdir, tmp.testname, '_surf_', tmp.variable, '_', ...
        num2str(min(RCM_info.years),'%04i'), '_',num2str(max(RCM_info.years),'%04i'), ...
        '_', RCM_info.season,'.tif'); %% ~_year_month.jpg
    if (exist(tmp.tifname , 'file') ~= 2 || fig_flag==2)      
        run(tmp.param_script);
%% set data file name (mat) for writing
        tmp.matname = [dirs.matdir, tmp.testname, '_', tmp.regionname, '_', tmp.variable,...
            '_mean_', num2str(min(RCM_info.years),'%04i'), '-', num2str(max(RCM_info.years),'%04i'),...
        '_', RCM_info.season, '.mat'];
        if (exist(tmp.matname , 'file') ~= 2)
            disp('please get data from subroutine_004_001 first')
        else
            load(tmp.matname, 'RCM_data', 'RCM_grid');
        end 
        
        if(strcmp(tmp.variable, 'SST'))
            RCM_data.mean=RCM_data.mean;
        elseif (strcmp(tmp.variable,'SSH')==1)
            ssh_correction_for_fig = Func_0017_SSH_correction_for_CMIP6_RMSE(tmp.testname);
            RCM_data.mean=RCM_data.mean-ssh_correction_for_fig;
        end

        RCM_grid.mask_model = double(inpolygon(RCM_grid.cut_lon_rho,RCM_grid.cut_lat_rho,RCM_grid.refpolygon(:,1),RCM_grid.refpolygon(:,2)));
        RCM_grid.mask_model(RCM_grid.mask_model==0)=NaN;
        RCM_data.mean=RCM_data.mean.*RCM_grid.mask_model;

        m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);
        hold on;

        [tmp.m_value, tmp.error_status] = Func_0011_get_area_weighted_mean(RCM_data.mean, RCM_grid.cut_lon_rho, RCM_grid.cut_lat_rho);
        m_pcolor(RCM_grid.cut_lon_rho',RCM_grid.cut_lat_rho',RCM_data.mean');
        shading(gca,param.m_pcolor_shading_method);   
        if(strcmp(tmp.variable, 'SST'))
            [C,h2]=m_contour(RCM_grid.cut_lon_rho',RCM_grid.cut_lat_rho', RCM_data.mean', [2, 5], 'color','k', ...
                    'linewidth', 1.5, 'linestyle', '-');
                clabel(C,h2,'FontSize',13,'Color','k', ...
                    'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
        end
 
        m_gshhs_i('color',param.m_gshhs_line_color)  
        m_gshhs_i('patch',param.m_gshhs_land_color);   % gray colored land

        m_grid('fontsize', param.m_grid_fontsize, 'box', param.m_grid_box_type, 'tickdir', param.m_grid_tickdir_type);
        if min(RCM_info.years) == max(RCM_info.years)
            tmp.titlename = strcat(tmp.variable, ', ', RCM_info.season(1:3), ', ', tmp.abb,',(',num2str(max(RCM_info.years),'%04i'),') ');                        
        else
            tmp.titlename = strcat(tmp.variable, ', ', RCM_info.season(1:3), ', ', tmp.abb,',(',num2str(min(RCM_info.years),'%04i'),'-',num2str(max(RCM_info.years),'%04i'),') ');
        end
        title(tmp.titlename,'fontsize',param.m_pcolor_title_fontsize);  %%title

        % set colorbar 
        h = colorbar;
%                     if (strcmp(tmp.variable,'SSH')==1)
%                         colormap(flip(cool));
%                     else
            colormap(jet);
%                     end

        set(h,'fontsize',param.colorbar_fontsize);
        title(h,param.colorbar_title,'fontsize',param.colorbar_title_fontsize);
        switch tmp.variable
            case {'SST', 'SSH', 'SSS', 'Uwind', 'Vwind', 'u', 'v', 'swrad', 'shflux'}
                caxis(param.colorbar_lev);
            case {'wstrcurl', 'wcurl' }
%                 caxis([-max(abs(RCM_data.mean(:))), max(abs(RCM_data.mean(:)))]);
                caxis(param.colorbar_lev);
                colormap(cmaps.byrmap);
        end

        disp(['M = ', num2str(tmp.m_value)]);
%         m_text(param.m_pcolor_ref_text_x_location, param.m_pcolor_ref_text_y_location, ['M = ', num2str(tmp.m_value)], 'FontSize', param.m_quiver_ref_text_fontsize); 

        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [param.hor_paper_size_x, param.hor_paper_size_y]);
        set(gcf,'PaperPosition', [param.paper_position_hor param.paper_position_ver param.paper_position_width param.paper_position_height]) 
        saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
        close all;
        clear RCM_data.mean
        RCM_grid=rmfield(RCM_grid, 'lon_rho');
    end
end