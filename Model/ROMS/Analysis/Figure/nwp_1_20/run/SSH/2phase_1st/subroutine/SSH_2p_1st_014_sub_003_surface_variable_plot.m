% start-------------------- decadal SST, SSH, SSS plot
        for varind2=1:length(RCM_info.vars)
            tmp.variable=RCM_info.vars{varind2};
            dirs.figdir=[dirs.figrawdir,'surface', tmp.fs, tmp.regionname, tmp.fs, tmp.variable, tmp.fs, ...
                num2str(min(RCM_info.years)), '_', num2str(max(RCM_info.years)), tmp.fs];
            if (exist(strcat(dirs.figdir) , 'dir') ~= 7)
                mkdir(strcat(dirs.figdir));
            end 
%             tmp.tifname=strcat(dirs.figdir, 'diff_RCM_ENSg', '_clim_', variable, '_',num2str(min(RCM_info.years),'%04i'), ...
%                          '_',num2str(max(RCM_info.years),'%04i'), ...
%                         '_', num2str(min(RCM_info.months),'%04i'), '-', num2str(max(RCM_info.months),'%04i'),'.tif'); %% ~_year_month.jpg
                    
            tmp.tifname=strcat(dirs.figdir, 'diff_RCM_', RCM_info.ensname, '_clim_', tmp.variable, '_',num2str(min(RCM_info.years),'%04i'), ...
                            '_', num2str(max(RCM_info.years),'%04i'), ...
                            '_', num2str(RCM_info.months(1),'%04i'), '-', num2str(RCM_info.months(end),'%04i'),'.tif'); %% ~_year_month.jpg
                
            run(tmp.param_script);

            early_enssavefilename = [dirs.enssavedir, RCM_info.ensname, '_', tmp.regionname, '_', tmp.variable, ...
                '_mean_', num2str(RCM_info.inputyear1(1),'%04i'), '-', num2str(RCM_info.inputyear1(end),'%04i'),...
                    '_', num2str(RCM_info.months(1),'%04i'), '-', num2str(RCM_info.months(end),'%04i'), '.mat'];
            load(early_enssavefilename)
            early_data=mean_data;

            late_enssavefilename = [dirs.enssavedir, RCM_info.ensname, '_', tmp.regionname, '_', tmp.variable, ...
                '_mean_', num2str(RCM_info.inputyear2(1),'%04i'), '-', num2str(RCM_info.inputyear2(end),'%04i'), ...
                    '_', num2str(RCM_info.months(1),'%04i'), '-', num2str(RCM_info.months(end),'%04i'), '.mat'];
            load(late_enssavefilename)
            late_data=mean_data;

            mean_data = late_data - early_data;

    %         if (strcmp(variable,'SSH')==1)
    %             ssh_correction_for_fig = Func_0014_SSH_correction_for_pcolor('RCM_ENSg_historical');
    %             mean_data=mean_data-ssh_correction_for_fig;
    %         end

            RCM_grid.mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,RCM_grid.refpolygon(:,1),RCM_grid.refpolygon(:,2)));
            RCM_grid.mask_model(RCM_grid.mask_model==0)=NaN;
            mean_data=mean_data.*RCM_grid.mask_model;

            m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);
            hold on;

            [tmp.m_value, tmp.error_status] = Func_0011_get_area_weighted_mean(mean_data, cut_lon_rho, cut_lat_rho);
            m_pcolor(cut_lon_rho',cut_lat_rho',mean_data');
            shading(gca,param.m_pcolor_shading_method);      
            
            if(strcmp(tmp.variable, 'SST'))
                [C,h2]=m_contour(cut_lon_rho',cut_lat_rho', mean_data', [1, 1.5, 2], 'color','k', ...
                        'linewidth', 1.5, 'linestyle', '-');
                    clabel(C,h2,'FontSize',13,'Color','k', ...
                        'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
            end

            m_gshhs_i('color',param.m_gshhs_line_color)  
            m_gshhs_i('patch',param.m_gshhs_land_color);   % gray colored land

            m_grid('fontsize', param.m_grid_fontsize, 'box', param.m_grid_box_type, 'tickdir', param.m_grid_tickdir_type);
            tmp.titlename = strcat(tmp.variable, ' diff, ', season(1:3), ', ', RCM_info.ensname, ',(', num2str(RCM_info.inputyear2(1),'%04i'),'-',num2str(RCM_info.inputyear1(end),'%04i'),') ');  %% + glacier contribution
            title(tmp.titlename,'fontsize',param.m_pcolor_title_fontsize);  %%title

            % set colorbar 
            h = colorbar;
            colormap(jet);

            set(h,'fontsize',param.colorbar_fontsize);
            title(h,param.colorbar_title,'fontsize',param.colorbar_title_fontsize);
            if(strcmp(tmp.variable, 'SST'))
                caxis([-1 3]);
                colormap(jet);
            elseif(strcmp(tmp.variable, 'SSS'))
                caxis([-1 1]);
                [cmaps.byrmap, tmp.error_status] = Func_0009_get_colormaps('byr3', tmp.dropboxpath);
                colormap(cmaps.byrmap);
            elseif(strcmp(tmp.variable, 'SSH'))   
                caxis([0 0.4] + (RCM_info.inputyear2-2050)*0.006);
                colormap(jet);
            end

            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [param.hor_paper_size_x, param.hor_paper_size_y]);
            set(gcf,'PaperPosition', [param.paper_position_hor param.paper_position_ver param.paper_position_width param.paper_position_height]) 
            saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
            close all;

            clear lon_rho mean_data ens_data
        end