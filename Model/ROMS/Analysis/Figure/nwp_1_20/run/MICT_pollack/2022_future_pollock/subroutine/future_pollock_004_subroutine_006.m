% %  Updated 21-May-2021 by Yong-Yub Kim, 


RCM_info.years_ssp =RCM_info.years;
tmp.testname_ssp = tmp.testname;
dirs.filedir_ssp = strcat('/Volumes/kyy_raid/Data/Model/ROMS/', RCM_info.model, '/backup_surf/', tmp.testname_ssp, '/run/'); % % where data files are          
dirs.matdir_ssp = strcat('/Volumes/kyy_raid/Data/Model/ROMS/', RCM_info.model, '/', tmp.testname_ssp, '/run/mean/');

[tmp.testname_his, tmp.error_status] = Func_0023_RCM_CMIP6_testname_his(tmp.testname_ssp);
dirs.filedir_his = strcat('/Volumes/kyy_raid/Data/Model/ROMS/', RCM_info.model, '/backup_surf/', tmp.testname_his, '/run/'); % % where data files are          
dirs.matdir_his = strcat('/Volumes/kyy_raid/Data/Model/ROMS/', RCM_info.model, '/', tmp.testname_his, '/run/mean/');


for varind2=1:length(RCM_info.vars)
    %% set variable name & figure directory
    tmp.variable=RCM_info.vars{varind2};
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
        switch tmp.testname
           case 'ENS4_fut'
               %% hist
               tmp.tset= {'test2117', 'test2118','test2119','test2120'};
               for titi=1:4
                    tmp.ttname=tmp.tset{titi};
                    dirs.ttmatdir = strcat('/Volumes/kyy_raid/Data/Model/ROMS/', RCM_info.model, '/', tmp.ttname, '/run/mean/');
                    tmp.matname = [dirs.ttmatdir, tmp.ttname, '_', tmp.regionname, '_', tmp.variable,...
                        '_mean_', num2str(min(RCM_info.years_his),'%04i'), '-', num2str(max(RCM_info.years_his),'%04i'),...
                    '_', RCM_info.season, '.mat'];
                    if (exist(tmp.matname , 'file') ~= 2)
                        disp('please get data from subroutine_004_001 first')
                    else
                        load(tmp.matname, 'RCM_data', 'RCM_grid');
                    end 
                    tmp.tmean_hist(:,:,titi)=RCM_data.mean;
               end

            %% fut
               tmp.tset= {'test2127', 'test2128','test2129','test2130'};
               for titi=1:4
                    tmp.ttname=tmp.tset{titi};
                    dirs.ttmatdir = strcat('/Volumes/kyy_raid/Data/Model/ROMS/', RCM_info.model, '/', tmp.ttname, '/run/mean/');
                    tmp.matname = [dirs.ttmatdir, tmp.ttname, '_', tmp.regionname, '_', tmp.variable,...
                        '_mean_', num2str(min(RCM_info.years_ssp),'%04i'), '-', num2str(max(RCM_info.years_ssp),'%04i'),...
                    '_', RCM_info.season, '.mat'];
                    if (exist(tmp.matname , 'file') ~= 2)
                        disp('please get data from subroutine_004_001 first')
                    else
                        load(tmp.matname, 'RCM_data', 'RCM_grid');
                    end 
                    tmp.tmean_fut(:,:,titi)=RCM_data.mean;
               end

           %% calculation of difference
                RCM_data.mean=mean(tmp.tmean_fut, 3) - mean(tmp.tmean_hist, 3);
            otherwise
                %% set data file name (mat) and reading (ssp)
                tmp.matname_ssp = [dirs.matdir_ssp, tmp.testname_ssp, '_', tmp.regionname, '_', tmp.variable,...
                    '_mean_', num2str(min(RCM_info.years_ssp),'%04i'), '-', num2str(max(RCM_info.years_ssp),'%04i'),...
                '_', RCM_info.season, '.mat'];
                if (exist(tmp.matname_ssp , 'file') ~= 2)
                    disp('please get data from subroutine_004_001 first')
                else
                    load(tmp.matname_ssp, 'RCM_data', 'RCM_grid');
                end 
                RCM_data_ssp = RCM_data;
               
               %% set data file name (mat) and reading (his)
                tmp.matname_his = [dirs.matdir_his, tmp.testname_his, '_', tmp.regionname, '_', tmp.variable,...
                    '_mean_', num2str(min(RCM_info.years_his),'%04i'), '-', num2str(max(RCM_info.years_his),'%04i'),...
                '_', RCM_info.season, '.mat'];
                if (exist(tmp.matname_his , 'file') ~= 2)
                    disp('please get data from subroutine_004_001 first')
                else
                    load(tmp.matname_his, 'RCM_data', 'RCM_grid');
                end 
                RCM_data_his = RCM_data;
                
                %% calculation of difference
                RCM_data.mean=RCM_data_ssp.mean - RCM_data_his.mean;
        end

        [tmp.m_value, tmp.error_status] = Func_0011_get_area_weighted_mean(RCM_data.mean, RCM_grid.cut_lon_rho, RCM_grid.cut_lat_rho);

        RCM_grid.mask_model = double(inpolygon(RCM_grid.cut_lon_rho,RCM_grid.cut_lat_rho,RCM_grid.refpolygon(:,1),RCM_grid.refpolygon(:,2)));
        RCM_grid.mask_model(RCM_grid.mask_model==0)=NaN;
        RCM_data.mean=RCM_data.mean.*RCM_grid.mask_model;
       
        %% figure
        m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);
        hold on;
        m_pcolor(RCM_grid.cut_lon_rho',RCM_grid.cut_lat_rho',RCM_data.mean');
        shading(gca,param.m_pcolor_shading_method);   
        if(strcmp(tmp.variable, 'SST'))
            [C,h2]=m_contour(RCM_grid.cut_lon_rho',RCM_grid.cut_lat_rho', RCM_data.mean', [0, 2, 4], 'color','k', ...
                    'linewidth', 1.5, 'linestyle', '-');
                clabel(C,h2,'FontSize',13,'Color','k', ...
                    'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
        end
 
        m_gshhs_i('color',param.m_gshhs_line_color)  
        m_gshhs_i('patch',param.m_gshhs_land_color);   % gray colored land

        m_grid('fontsize', param.m_grid_fontsize, 'box', param.m_grid_box_type, 'tickdir', param.m_grid_tickdir_type);
%         if min(RCM_info.years) == max(RCM_info.years)
%             tmp.titlename = strcat(tmp.variable, ', ', RCM_info.season(1:3), ', ', tmp.abb,',(',num2str(max(RCM_info.years),'%04i'),') ');                        
%         else
%             tmp.titlename = strcat(tmp.variable, ', ', RCM_info.season(1:3), ', ', tmp.abb,',(',num2str(min(RCM_info.years),'%04i'),'-',num2str(max(RCM_info.years),'%04i'),') ');
%         end
        tmp.titlename = strcat(tmp.variable, ', ', RCM_info.season(1:3), ', ', tmp.abb,',(','diff',') ');

        title(tmp.titlename,'fontsize',param.m_pcolor_title_fontsize);  %%title

        h = colorbar;
        colormap(jet);

        set(h,'fontsize',param.colorbar_fontsize);
        title(h,param.colorbar_title,'fontsize',param.colorbar_title_fontsize);
        switch tmp.variable
            case {'SST'}
%                 caxis(param.colorbar_lev);
%                 caxis([4 5.5])
                caxis([3 5])                
                colormap(flip(autumn));
            case {'SSH'}
                caxis([0.4 0.8])
                colormap(flip(autumn));
            case{'SSS'}
                caxis([-0.5 0.5])
                colormap(cmaps.byrmap);
            case{'Uwind', 'Vwind'}
                caxis([-1 1])
                colormap(cmaps.byrmap);
            case{'u','v'}
                caxis([-0.1 0.1])
                colormap(cmaps.byrmap);
            case {'wstrcurl', 'wcurl'}
%                 caxis([-max(abs(RCM_data.mean(:))), max(abs(RCM_data.mean(:)))]);
                caxis(param.colorbar_lev);
                colormap(cmaps.byrmap);
            case{'swrad', 'lwrad'}
                caxis([-20 20])
                colormap(cmaps.byrmap);
            case{'shflux', 'latent', 'sensible'}
                caxis([-100 100])
                colormap(cmaps.byrmap);
        end

        disp(['M = ', num2str(tmp.m_value)]);
%         m_text(param.m_pcolor_ref_text_x_location, param.m_pcolor_ref_text_y_location, ['M = ', num2str(tmp.m_value)], 'FontSize', param.m_quiver_ref_text_fontsize); 
        
        tmp.mean_data=RCM_data.mean;
        [tmp.indw, tmp.inde, tmp.inds, tmp.indn]=Func_0012_findind_Y(1/10,[127, 131, 36, 39],RCM_grid.cut_lon_rho,RCM_grid.cut_lat_rho); % southern EKB
        tmp.sEKB_lon=RCM_grid.cut_lon_rho(tmp.indw:tmp.inde,tmp.inds:tmp.indn);
        tmp.sEKB_lat=RCM_grid.cut_lat_rho(tmp.indw:tmp.inde,tmp.inds:tmp.indn);
        tmp.sEKB_data=tmp.mean_data(tmp.indw:tmp.inde,tmp.inds:tmp.indn);
        [tmp.sEKB_mean, tmp.error_status] = Func_0011_get_area_weighted_mean(tmp.sEKB_data, tmp.sEKB_lon, tmp.sEKB_lat);
        
        [tmp.indw, tmp.inde, tmp.inds, tmp.indn]=Func_0012_findind_Y(1/10,[127, 131, 39, 41],RCM_grid.cut_lon_rho,RCM_grid.cut_lat_rho); % southern EKB
        tmp.nEKB_lon=RCM_grid.cut_lon_rho(tmp.indw:tmp.inde,tmp.inds:tmp.indn);
        tmp.nEKB_lat=RCM_grid.cut_lat_rho(tmp.indw:tmp.inde,tmp.inds:tmp.indn);
        tmp.nEKB_data=tmp.mean_data(tmp.indw:tmp.inde,tmp.inds:tmp.indn);
        [tmp.nEKB_mean, tmp.error_status] = Func_0011_get_area_weighted_mean(tmp.nEKB_data, tmp.nEKB_lon, tmp.nEKB_lat);
        
        [tmp.indw, tmp.inde, tmp.inds, tmp.indn]=Func_0012_findind_Y(1/10,[127, 131, 36, 41],RCM_grid.cut_lon_rho,RCM_grid.cut_lat_rho); % southern EKB
        tmp.EKB_lon=RCM_grid.cut_lon_rho(tmp.indw:tmp.inde,tmp.inds:tmp.indn);
        tmp.EKB_lat=RCM_grid.cut_lat_rho(tmp.indw:tmp.inde,tmp.inds:tmp.indn);
        tmp.EKB_data=tmp.mean_data(tmp.indw:tmp.inde,tmp.inds:tmp.indn);
        [tmp.EKB_mean, tmp.error_status] = Func_0011_get_area_weighted_mean(tmp.EKB_data, tmp.EKB_lon, tmp.EKB_lat);
        
        disp([tmp.variable, ', ', 'EKB = ', num2str(tmp.EKB_mean), ', ', 'NEKB = ', num2str(tmp.nEKB_mean), ', ', 'SEKB = ', num2str(tmp.sEKB_mean)]);

        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [param.hor_paper_size_x, param.hor_paper_size_y]);
        set(gcf,'PaperPosition', [param.paper_position_hor param.paper_position_ver param.paper_position_width param.paper_position_height]) 
        saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
        close all;
        clear RCM_data.mean
        RCM_grid=rmfield(RCM_grid, 'lon_rho');
        
        disp('future_pollock_004_subroutine_006');

    end
end


        
        
