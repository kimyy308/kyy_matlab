% start-------------------- decadal current plot
            tmp.variable='UV';
            dirs.figdir=[dirs.figrawdir,'surface', tmp.fs, tmp.regionname, tmp.fs, 'vec', tmp.fs, ...
                num2str(min(RCM_info.years)), '_', num2str(max(RCM_info.years)), tmp.fs];
            if (exist(strcat(dirs.figdir) , 'dir') ~= 7)
                mkdir(strcat(dirs.figdir));
            end 
            tmp.tifname=strcat(dirs.figdir, 'RCM_ENSg', '_clim_uv_',num2str(min(RCM_info.years),'%04i'), ...
            '_',num2str(max(RCM_info.years),'%04i'),  ...
                            '_', num2str(RCM_info.months(1),'%04i'), '-', num2str(RCM_info.months(end),'%04i'),'.tif'); %% ~_year_month.jpg
                        
            for testnameind2=1:length(RCM_info.name)
                tmp.testname=RCM_info.name{testnameind2};
                dirs.filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\output\', tmp.testname, '\run\packed_monthly\'); % % where data files are  
                dirs.matdir = strcat('D:\Data\Model\ROMS\nwp_1_20\', tmp.testname, '\run\mean\'); % variables of each test
%                 matdir = strcat('D:\Data\Model\ROMS\nwp_1_20\', tmp.testname, '\run\mean\');
%                 griddir = strcat('D:\Data\Model\ROMS\nwp_1_20\output\', tmp.testname, '\run\packed_monthly\'); % % where data files are            
                
                tmp.matname = [dirs.matdir, tmp.testname, '_', tmp.regionname, '_', tmp.variable, ...
                        '_mean_', num2str(RCM_info.years(1),'%04i'), '-', num2str(RCM_info.years(end),'%04i'), ...
                                    '_', num2str(RCM_info.months(1),'%04i'), '-', num2str(RCM_info.months(end),'%04i'),'.mat'];
%                 matname = [matdir, tmp.testname, '_', regionname, '_', variable, '_mean_', num2str(min(inputyear1),'%04i'), '-', num2str(max(inputyear1),'%04i'),  ...
%                             '_', num2str(min(RCM_info.months),'%04i'), '-', num2str(max(RCM_info.months),'%04i'),'.mat'];
                load(tmp.matname);
                if (exist('ens_u_rho' , 'var') ~= 1)
                    ens_u_rho=u_rho;
                    ens_v_rho=v_rho;
                else
                    ens_u_rho=ens_u_rho+u_rho;
                    ens_v_rho=ens_v_rho+v_rho;
                end
            end
            ens_u_rho=ens_u_rho/length(RCM_info.name);
            ens_v_rho=ens_v_rho/length(RCM_info.name);
            u_rho=ens_u_rho;
            v_rho=ens_v_rho;

            RCM_grid.mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,RCM_grid.refpolygon(:,1),RCM_grid.refpolygon(:,2)));
            RCM_grid.mask_model(RCM_grid.mask_model==0)=NaN;
            u_rho=u_rho(1:size(RCM_grid.mask_model,1),:).*RCM_grid.mask_model;
            v_rho=v_rho(:,1:size(RCM_grid.mask_model,2)).*RCM_grid.mask_model; 

            m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);
            hold on;
            m_gshhs_i('color',param.m_gshhs_line_color)  
            m_gshhs_i('patch',param.m_gshhs_land_color);   % gray colored land

            uvplot=m_quiver(cut_lon_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)', ...
                        cut_lat_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)', ...
                        u_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)' * param.m_quiver_vector_size, ...
                        v_rho(1:param.m_quiver_x_interval:end, 1:param.m_quiver_y_interval:end)' * param.m_quiver_vector_size, ...
                        'color', param.m_quiver_vector_color, 'AutoScale','off','LineWidth', param.m_quiver_LineWidth);


            m_text(param.m_quiver_ref_text_x_location, param.m_quiver_ref_text_y_location, param.m_quiver_ref_text, 'FontSize', param.m_quiver_ref_text_fontsize); 
            m_grid('fontsize', param.m_grid_fontsize, 'box', param.m_grid_box_type, 'tickdir', param.m_grid_tickdir_type);
            if min(RCM_info.years) == max(RCM_info.years)
                param.titlename = strcat('UV, ', season(1:3), ', ', tmp.abb,',(',num2str(min(RCM_info.years),'%04i'),') ');  %% + glacier contribution
            else
                param.titlename = strcat('UV, ', season(1:3), ', ', tmp.abb,',(',num2str(RCM_info.years(1),'%04i'),'-',num2str(RCM_info.years(end),'%04i'),') ');  %% + glacier contribution
            end
            title(param.titlename,'fontsize',param.m_pcolor_title_fontsize);  %%title

            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [param.hor_paper_size_x, param.hor_paper_size_y]);
            set(gcf,'PaperPosition', [param.paper_position_hor param.paper_position_ver param.paper_position_width param.paper_position_height]) 
            saveas(gcf, tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
        %                 system(['magick ', pngname, ' -trim ', pngname]);
            enssavefilename = [dirs.enssavedir, 'ENSg_', tmp.regionname, '_', tmp.variable, '_mean_', num2str(RCM_info.years(1),'%04i'), '-', num2str(RCM_info.years(end),'%04i'),...
                            '_', num2str(RCM_info.months(1),'%04i'), '-', num2str(RCM_info.months(end),'%04i'), '.mat'];
            save(enssavefilename, 'cut_lon_rho', 'cut_lat_rho', 'u_rho', 'v_rho');
            close all;
            clear lon_rho ens_u_rho ref_vec_x_range