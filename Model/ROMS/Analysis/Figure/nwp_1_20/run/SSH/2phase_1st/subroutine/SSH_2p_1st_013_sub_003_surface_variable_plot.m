% start-------------------- decadal SST, SSH, SSS plot
            for varind2=1:length(RCM_info.vars)
                tmp.variable=RCM_info.vars{varind2};
                dirs.figdir=[dirs.figrawdir,'surface', tmp.fs, tmp.regionname, tmp.fs, tmp.variable, tmp.fs, ...
                    num2str(min(RCM_info.years)), '_', num2str(max(RCM_info.years)), tmp.fs];
                if (exist(strcat(dirs.figdir) , 'dir') ~= 7)
                    mkdir(strcat(dirs.figdir));
                end 
                tmp.tifname=strcat(dirs.figdir, 'RCM_ENSg', '_clim_', tmp.variable, '_',num2str(min(RCM_info.years),'%04i'), ...
                             '_',num2str(max(RCM_info.years),'%04i'),  ...
                            '_', num2str(RCM_info.months(1),'%04i'), '-', num2str(RCM_info.months(end),'%04i'),'.tif'); %% ~_year_month.jpg
                run(tmp.param_script);
                for testnameind2=1:length(RCM_info.name)
                    tmp.testname=RCM_info.name{testnameind2};
                    dirs.filedir = strcat('E:\Data\Model\ROMS\nwp_1_20\output\', tmp.testname, '\run\packed_monthly\'); % % where data files are          
                    dirs.matdir = strcat('D:\Data\Model\ROMS\nwp_1_20\', tmp.testname, '\run\mean\');
%                     griddir = strcat('E:\Data\Model\ROMS\nwp_1_20\output\', testname, '\run\packed_monthly\'); % % where data files are            

                    tmp.matname = [dirs.matdir, tmp.testname, '_', tmp.regionname, '_', tmp.variable,...
                            '_mean_', num2str(min(RCM_info.years),'%04i'), '-', num2str(max(RCM_info.years),'%04i'), ...
                            '_', num2str(RCM_info.months(1),'%04i'), '-', num2str(RCM_info.months(end),'%04i'), '.mat'];
                    load(tmp.matname);
                    if (exist('ens_data' , 'var') ~= 1)
                        ens_data=mean_data;
                    else
                        ens_data=ens_data+mean_data;
                    end
                end
                ens_data=ens_data/length(RCM_info.name);
                mean_data=ens_data;
        %     [m_value, error_status] = Func_0011_get_area_weighted_mean(mean_data, cut_lon_rho, cut_lat_rho)

                if (strcmp(tmp.variable,'SSH')==1)
        %             ssh_correction_for_fig = Func_0014_SSH_correction_for_pcolor('RCM_ENSg_historical');
                    ssh_correction_for_fig = Func_0017_SSH_correction_for_CMIP6_RMSE('RCM_ENSg_historical');
                    mean_data=mean_data-ssh_correction_for_fig;
                end

                RCM_grid.mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,RCM_grid.refpolygon(:,1),RCM_grid.refpolygon(:,2)));
                RCM_grid.mask_model(RCM_grid.mask_model==0)=NaN;
                mean_data=mean_data.*RCM_grid.mask_model;

                m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);
                hold on;

                [tmp.m_value, tmp.error_status] = Func_0011_get_area_weighted_mean(mean_data, cut_lon_rho, cut_lat_rho);
                m_pcolor(cut_lon_rho',cut_lat_rho',mean_data');
                shading(gca,param.m_pcolor_shading_method);   
                if(strcmp(tmp.variable, 'SST'))
                    [C,h2]=m_contour(cut_lon_rho',cut_lat_rho', mean_data', [5, 10, 15, 20, 25, 30], 'color','k', ...
                            'linewidth', 1.5, 'linestyle', '-');
                        clabel(C,h2,'FontSize',13,'Color','k', ...
                            'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
                end

                m_gshhs_i('color',param.m_gshhs_line_color)  
                m_gshhs_i('patch',param.m_gshhs_land_color);   % gray colored land

                m_grid('fontsize', param.m_grid_fontsize, 'box', param.m_grid_box_type, 'tickdir', param.m_grid_tickdir_type);
                if min(RCM_info.years) == max(RCM_info.years)
                    tmp.titlename = strcat(tmp.variable, ', ', season(1:3), ', ', tmp.abb,',(',num2str(max(RCM_info.years),'%04i'),') ');                        
                else
                    tmp.titlename = strcat(tmp.variable, ', ', season(1:3), ', ', tmp.abb,',(',num2str(min(RCM_info.years),'%04i'),'-',num2str(max(RCM_info.years),'%04i'),') ');
                end
                title(tmp.titlename,'fontsize',param.m_pcolor_title_fontsize);  %%title

                % set colorbar 
                h = colorbar;
                colormap(jet);

                set(h,'fontsize',param.colorbar_fontsize);
                title(h,param.colorbar_title,'fontsize',param.colorbar_title_fontsize);
                switch tmp.variable
                    case {'SST', 'SSH', 'SSS'}
                        caxis(param.colorbar_lev);
                end

                disp(['M = ', num2str(tmp.m_value)]);
                m_text(param.m_pcolor_ref_text_x_location, param.m_pcolor_ref_text_y_location, ['M = ', num2str(tmp.m_value)], 'FontSize', param.m_quiver_ref_text_fontsize); 

                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [param.hor_paper_size_x, param.hor_paper_size_y]);
                set(gcf,'PaperPosition', [param.paper_position_hor param.paper_position_ver param.paper_position_width param.paper_position_height]) 
                saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
                enssavefilename = [dirs.enssavedir, 'ENSg_', tmp.regionname, '_', tmp.variable, '_mean_', num2str(RCM_info.years(1),'%04i'), '-', num2str(RCM_info.years(end),'%04i'),...
                            '_', num2str(RCM_info.months(1),'%04i'), '-', num2str(RCM_info.months(end),'%04i'), '.mat'];
                save(enssavefilename, 'cut_lon_rho', 'cut_lat_rho', 'mean_data');
                close all;

                clear lon_rho mean_data ens_data
            end