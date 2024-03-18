% %  Updated    19-May-2022 by Yong-Yub Kim   % make

if testnameind==1

    %% historical plot
    tmp.tifname=strcat(dirs.figdir_allmean, 'all_mean','_',tmp.regionname, '_spawning_ground_probability_', ...
                num2str(min(RCM_info.years_his),'%04i'),'_',num2str(max(RCM_info.years_his),'%04i'), 'y_', ...
                RCM_info.season, '.tif'); 

    tmp.testname_prefix=tmp.testname_his(1:4); 
    tmp.testlen=length(RCM_info.name);
    tmp.tlen_his=0;

    if (exist(tmp.tifname , 'file') ~= 2 || flags.fig_switch(flagi)==2)

        run(tmp.param_script);

        for sub_testnameind=1:length(RCM_info.name)
            tmp_sub.testname_ssp=RCM_info.name{sub_testnameind};
            [tmp_sub.testname_his, tmp.error_status] = Func_0023_RCM_CMIP6_testname_his(tmp_sub.testname_ssp);
            dirs_sub.filedir_his = strcat('/Volumes/kyy_raid/Data/Model/ROMS/nwp_1_20/', tmp_sub.testname_his, '/pollock/'); % % where data files are
            dirs_sub.savedir_his = strcat('/Volumes/kyy_raid/Data/Model/ROMS/nwp_1_20/', tmp_sub.testname_his, '/pollock/');

            %%    initialization
            if sub_testnameind==1
                for yearij = 1:length(RCM_info.years_his)
                    tmp.tempyear = RCM_info.years_his(yearij);
                    for monthij = 1:length(RCM_info.months)
                        tmp.tempmonth = RCM_info.months(monthij);
                        tmp.ncname = [dirs_sub.savedir_his,tmp_sub.testname_his,'_',tmp.regionname,'model_pollock_',num2str(tmp.tempyear,'%04i'),'_',num2str(tmp.tempmonth,'%02i'),'.nc'];
                        tmp.tlen_his=tmp.tlen_his + length(ncread(tmp.ncname, 'time'));
                    end
                end
                tmp.xlen=size(ncread(tmp.ncname,'lon_rho'),1);
                tmp.ylen=size(ncread(tmp.ncname,'lon_rho'),2);
                tmp.comb_egg_mask=NaN(tmp.xlen,tmp.ylen,tmp.tlen_his.*tmp.testlen);
            end

        %% read grids, eggs    
            for yearij = 1:length(RCM_info.years_his)
                tmp.tempyear = RCM_info.years_his(yearij);
                for monthij = 1:length(RCM_info.months)
                    tmp.tempmonth = RCM_info.months(monthij);
                    tmp.ncname = [dirs_sub.savedir_his,tmp_sub.testname_his,'_',tmp.regionname,'model_pollock_',num2str(tmp.tempyear,'%04i'),'_',num2str(tmp.tempmonth,'%02i'),'.nc'];
                    disp([num2str(yearij), 'y_',num2str(monthij),'m'])

                    tmp.egg_mask=ncread(tmp.ncname, 'egg_mask');
                    tmp.lastday_m=size(tmp.egg_mask,3);
                    if yearij==1 && monthij==1 && sub_testnameind==1
                        tmp.comb_egg_mask(:,:,1:tmp.lastday_m)=tmp.egg_mask;
                        tmp.endij=tmp.lastday_m;
                    else
                        tmp.comb_egg_mask(:,:,tmp.endij+1:tmp.endij+tmp.lastday_m) = tmp.egg_mask;
                        tmp.endij=tmp.endij+tmp.lastday_m;
                    end
                end
            end
            RCM_grid.lon_rho = ncread(tmp.ncname, 'lon_rho');
            RCM_grid.lat_rho = ncread(tmp.ncname, 'lat_rho');


        end
        tmp.mean_data = mean(tmp.comb_egg_mask,3);

        tmp.mean_data(tmp.mean_data==0)=NaN;

        RCM_grid.mask_EKBcoast = double(inpolygon(RCM_grid.lon_rho,RCM_grid.lat_rho,[127 127 132 132],tmp.refpolygon(:,2)));
        RCM_grid.mask_EKBcoast(RCM_grid.mask_EKBcoast==0)=NaN;
        tmp.mean_data=tmp.mean_data.* RCM_grid.mask_EKBcoast;
        tmp.mean_data_his=tmp.mean_data;

    %% land fig    
        RCM_grid.mask_model2 = double(inpolygon(RCM_grid.lon_rho,RCM_grid.lat_rho,tmp.refpolygon(:,1),tmp.refpolygon(:,2)));
        RCM_grid.mask_model2(RCM_grid.mask_model2==0)=NaN;
        tmp.mean_data = tmp.mean_data .* RCM_grid.mask_model2;

        [tmp.mean_egg_mask, tmp.error_status] = Func_0011_get_area_weighted_mean(tmp.mean_data, RCM_grid.lon_rho, RCM_grid.lat_rho);


        tmp.testnameind=1;
        m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);

        ax{tmp.testnameind,1}=axes;

        if strcmp(tmp.testname_prefix, 'prob')
            tmp.temp_surf=ncread(tmp.ncname, 'temp_surf', [1 1], [inf inf]);
        else
            tmp.temp_surf=ncread(tmp.ncname, 'temp_surf', [1 1 1], [inf inf 1]);
        end

        tmp.temp_surf(isnan(tmp.temp_surf))=50000;
        tmp.temp_surf(tmp.temp_surf<50000)=NaN;
        RCM_grid.model_land=tmp.temp_surf;
        RCM_grid.model_land(tmp.temp_surf==50000)=1;
        pc{tmp.testnameind,1}=m_pcolor(RCM_grid.lon_rho',RCM_grid.lat_rho', RCM_grid.model_land','parent',ax{tmp.testnameind,1});
        colormap(ax{tmp.testnameind,1},[0.8 0.8 0.8]);
        shading(gca,param.m_pcolor_shading_method); 

        m_grid('fontsize', param.m_grid_fontsize, 'tickdir', param.m_grid_tickdir_type, 'box', param.m_grid_box_type,  ...
               'parent', ax{tmp.testnameind,1});

        pos_ax{tmp.testnameind}=get(ax{tmp.testnameind,1}, 'pos');
        col_bar{tmp.testnameind,1}=colorbar;
        set(col_bar{tmp.testnameind,1}, 'TickLabels', []);
        hold on

    %% data fig
        m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);
        ax{tmp.testnameind,2}=axes;
        pc{tmp.testnameind,2}=m_pcolor(RCM_grid.lon_rho',RCM_grid.lat_rho', tmp.mean_data','parent',ax{tmp.testnameind,2});
    %     colormap(ax{tmp.testnameind,2},jet);
%         colormap(ax{tmp.testnameind,2},parula);
        cmp.wr_08=Func_0009_get_colormaps('wr_08', tmp.dropboxpath);
        colormap(ax{tmp.testnameind,2},cmp.wr_08);

%         caxis([0, 1.0]);
        caxis([0, 0.8]);
        
        shading(gca,param.m_pcolor_shading_method);   

        m_grid('fontsize', param.m_grid_fontsize, 'tickdir', param.m_grid_tickdir_type, 'box', param.m_grid_box_type,  ...
               'backcolor', 'none','parent', ax{tmp.testnameind,2});
        col_bar{tmp.testnameind,2}=colorbar;
        set(col_bar{tmp.testnameind,2}, 'fontsize',param.m_grid_fontsize+5);

        tmp.titlename = strcat('sp ground, ','allmean', ',(', ...
            num2str(min(RCM_info.years_his),'%04i'),'-', num2str(max(RCM_info.years_his),'%04i'), ',',  ...
            RCM_info.season, ')'); 

        title(tmp.titlename,'fontsize',param.m_pcolor_title_fontsize, 'parent', ax{tmp.testnameind,1});  %%title for land figure
        title(tmp.titlename,'fontsize',param.m_pcolor_title_fontsize, 'parent', ax{tmp.testnameind,2});  %%title for data figure

        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [param.hor_paper_size_x, param.hor_paper_size_y]);
        set(gcf,'PaperPosition', [param.paper_position_hor param.paper_position_ver param.paper_position_width param.paper_position_height]) 

        [tmp.indw, tmp.inde, tmp.inds, tmp.indn]=Func_0012_findind_Y(1/10,[128, 130, 36, 39],RCM_grid.lon_rho,RCM_grid.lat_rho); % southern EKB
        sEKB_lon=RCM_grid.lon_rho(tmp.indw:tmp.inde,tmp.inds:tmp.indn);
        sEKB_lat=RCM_grid.lat_rho(tmp.indw:tmp.inde,tmp.inds:tmp.indn);
        sEKB_data=tmp.mean_data(tmp.indw:tmp.inde,tmp.inds:tmp.indn);
        [sEKB_mean, error_status] = Func_0011_get_area_weighted_mean(sEKB_data, sEKB_lon, sEKB_lat)
        
        [indw, inde, inds, indn]=Func_0012_findind_Y(1/10,[128, 130, 39, 41],RCM_grid.lon_rho,RCM_grid.lat_rho); % southern EKB
        nEKB_lon=RCM_grid.lon_rho(indw:inde,inds:indn);
        nEKB_lat=RCM_grid.lat_rho(indw:inde,inds:indn);
        nEKB_data=tmp.mean_data(indw:inde,inds:indn);
        [nEKB_mean, error_status] = Func_0011_get_area_weighted_mean(nEKB_data, nEKB_lon, nEKB_lat)
        
        [EKB_mean, error_status] = Func_0011_get_area_weighted_mean(tmp.mean_data, RCM_grid.lon_rho, RCM_grid.lat_rho)
        
        saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);

        disp(' ')
        disp([num2str(flagi), ' plot is created.'])
        disp(' ')
        disp([' File path is : ',tmp.tifname])
        disp(' ')
        close all;
    end



    %% future(ssp) plot
    tmp.tifname=strcat(dirs.figdir_allmean, 'all_mean','_',tmp.regionname, '_spawning_ground_probability_', ...
                num2str(min(RCM_info.years_ssp),'%04i'),'_',num2str(max(RCM_info.years_ssp),'%04i'), 'y_', ...
                RCM_info.season, '.tif'); 

    tmp.testname_prefix=tmp.testname_ssp(1:4); 
    tmp.testlen=length(RCM_info.name);
    tmp.tlen_ssp=0;

    if (exist(tmp.tifname , 'file') ~= 2 || flags.fig_switch(flagi)==2)

        run(tmp.param_script);

        for sub_testnameind=1:length(RCM_info.name)
            tmp_sub.testname_ssp=RCM_info.name{sub_testnameind};
            dirs_sub.filedir_ssp = strcat('/Volumes/kyy_raid/Data/Model/ROMS/nwp_1_20/', tmp_sub.testname_ssp, '/pollock/'); % % where data files are
            dirs_sub.savedir_ssp = strcat('/Volumes/kyy_raid/Data/Model/ROMS/nwp_1_20/', tmp_sub.testname_ssp, '/pollock/');

            %%    initialization
            if sub_testnameind==1
                for yearij = 1:length(RCM_info.years_ssp)
                    tmp.tempyear = RCM_info.years_ssp(yearij);
                    for monthij = 1:length(RCM_info.months)
                        tmp.tempmonth = RCM_info.months(monthij);
                        tmp.ncname = [dirs_sub.savedir_ssp,tmp_sub.testname_ssp,'_',tmp.regionname,'model_pollock_',num2str(tmp.tempyear,'%04i'),'_',num2str(tmp.tempmonth,'%02i'),'.nc'];
                        tmp.tlen_ssp=tmp.tlen_ssp + length(ncread(tmp.ncname, 'time'));
                    end
                end
                tmp.xlen=size(ncread(tmp.ncname,'lon_rho'),1);
                tmp.ylen=size(ncread(tmp.ncname,'lon_rho'),2);
                tmp.comb_egg_mask=NaN(tmp.xlen,tmp.ylen,tmp.tlen_ssp.*tmp.testlen);
            end

        %% read grids, eggs    
            for yearij = 1:length(RCM_info.years_ssp)
                tmp.tempyear = RCM_info.years_ssp(yearij);
                for monthij = 1:length(RCM_info.months)
                    tmp.tempmonth = RCM_info.months(monthij);
                    tmp.ncname = [dirs_sub.savedir_ssp,tmp_sub.testname_ssp,'_',tmp.regionname,'model_pollock_',num2str(tmp.tempyear,'%04i'),'_',num2str(tmp.tempmonth,'%02i'),'.nc'];
                    disp([num2str(yearij), 'y_',num2str(monthij),'m'])

                    tmp.egg_mask=ncread(tmp.ncname, 'egg_mask');
                    tmp.lastday_m=size(tmp.egg_mask,3);
                    if yearij==1 && monthij==1 && sub_testnameind==1
                        tmp.comb_egg_mask(:,:,1:tmp.lastday_m)=tmp.egg_mask;
                        tmp.endij=tmp.lastday_m;
                    else
                        tmp.comb_egg_mask(:,:,tmp.endij+1:tmp.endij+tmp.lastday_m) = tmp.egg_mask;
                        tmp.endij=tmp.endij+tmp.lastday_m;
                    end
                end
            end
            RCM_grid.lon_rho = ncread(tmp.ncname, 'lon_rho');
            RCM_grid.lat_rho = ncread(tmp.ncname, 'lat_rho');


        end
        tmp.mean_data = mean(tmp.comb_egg_mask,3);

        tmp.mean_data(tmp.mean_data==0)=NaN;
        tmp.mean_data_ssp=tmp.mean_data;
        
        RCM_grid.mask_EKBcoast = double(inpolygon(RCM_grid.lon_rho,RCM_grid.lat_rho,[127 127 132 132],tmp.refpolygon(:,2)));
        RCM_grid.mask_EKBcoast(RCM_grid.mask_EKBcoast==0)=NaN;
        tmp.mean_data=tmp.mean_data.* RCM_grid.mask_EKBcoast;

    %% land fig    
        RCM_grid.mask_model2 = double(inpolygon(RCM_grid.lon_rho,RCM_grid.lat_rho,tmp.refpolygon(:,1),tmp.refpolygon(:,2)));
        RCM_grid.mask_model2(RCM_grid.mask_model2==0)=NaN;
        tmp.mean_data = tmp.mean_data .* RCM_grid.mask_model2;

        [tmp.mean_egg_mask, tmp.error_status] = Func_0011_get_area_weighted_mean(tmp.mean_data, RCM_grid.lon_rho, RCM_grid.lat_rho);


        tmp.testnameind=1;
        m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);

        ax{tmp.testnameind,1}=axes;

        if strcmp(tmp.testname_prefix, 'prob')
            tmp.temp_surf=ncread(tmp.ncname, 'temp_surf', [1 1], [inf inf]);
        else
            tmp.temp_surf=ncread(tmp.ncname, 'temp_surf', [1 1 1], [inf inf 1]);
        end

        tmp.temp_surf(isnan(tmp.temp_surf))=50000;
        tmp.temp_surf(tmp.temp_surf<50000)=NaN;
        RCM_grid.model_land=tmp.temp_surf;
        RCM_grid.model_land(tmp.temp_surf==50000)=1;
        pc{tmp.testnameind,1}=m_pcolor(RCM_grid.lon_rho',RCM_grid.lat_rho', RCM_grid.model_land','parent',ax{tmp.testnameind,1});
        colormap(ax{tmp.testnameind,1},[0.8 0.8 0.8]);
        shading(gca,param.m_pcolor_shading_method); 

        m_grid('fontsize', param.m_grid_fontsize, 'tickdir', param.m_grid_tickdir_type, 'box', param.m_grid_box_type,  ...
               'parent', ax{tmp.testnameind,1});

        pos_ax{tmp.testnameind}=get(ax{tmp.testnameind,1}, 'pos');
        col_bar{tmp.testnameind,1}=colorbar;
        set(col_bar{tmp.testnameind,1}, 'TickLabels', []);
        hold on

    %% data fig
        m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);
        ax{tmp.testnameind,2}=axes;
        pc{tmp.testnameind,2}=m_pcolor(RCM_grid.lon_rho',RCM_grid.lat_rho', tmp.mean_data','parent',ax{tmp.testnameind,2});
    %     colormap(ax{tmp.testnameind,2},jet);
%         colormap(ax{tmp.testnameind,2},parula);
        cmp.wr_08=Func_0009_get_colormaps('wr_08', tmp.dropboxpath);
        colormap(ax{tmp.testnameind,2},cmp.wr_08);
%         caxis([0, 1.0]);
        caxis([0, 0.8]);
        shading(gca,param.m_pcolor_shading_method);   

        m_grid('fontsize', param.m_grid_fontsize, 'tickdir', param.m_grid_tickdir_type, 'box', param.m_grid_box_type,  ...
               'backcolor', 'none','parent', ax{tmp.testnameind,2});
        col_bar{tmp.testnameind,2}=colorbar;
        set(col_bar{tmp.testnameind,2}, 'fontsize',param.m_grid_fontsize+5);

        tmp.titlename = strcat('sp ground, ','allmean', ',(', ...
            num2str(min(RCM_info.years_ssp),'%04i'),'-', num2str(max(RCM_info.years_ssp),'%04i'), ',',  ...
            RCM_info.season, ')'); 

        title(tmp.titlename,'fontsize',param.m_pcolor_title_fontsize, 'parent', ax{tmp.testnameind,1});  %%title for land figure
        title(tmp.titlename,'fontsize',param.m_pcolor_title_fontsize, 'parent', ax{tmp.testnameind,2});  %%title for data figure

        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [param.hor_paper_size_x, param.hor_paper_size_y]);
        set(gcf,'PaperPosition', [param.paper_position_hor param.paper_position_ver param.paper_position_width param.paper_position_height]) 

        [tmp.indw, tmp.inde, tmp.inds, tmp.indn]=Func_0012_findind_Y(1/10,[128, 130, 36, 39],RCM_grid.lon_rho,RCM_grid.lat_rho); % southern EKB
        sEKB_lon=RCM_grid.lon_rho(tmp.indw:tmp.inde,tmp.inds:tmp.indn);
        sEKB_lat=RCM_grid.lat_rho(tmp.indw:tmp.inde,tmp.inds:tmp.indn);
        sEKB_data=tmp.mean_data(tmp.indw:tmp.inde,tmp.inds:tmp.indn);
        [sEKB_mean, error_status] = Func_0011_get_area_weighted_mean(sEKB_data, sEKB_lon, sEKB_lat)
        
        [indw, inde, inds, indn]=Func_0012_findind_Y(1/10,[128, 130, 39, 41],RCM_grid.lon_rho,RCM_grid.lat_rho); % southern EKB
        nEKB_lon=RCM_grid.lon_rho(indw:inde,inds:indn);
        nEKB_lat=RCM_grid.lat_rho(indw:inde,inds:indn);
        nEKB_data=tmp.mean_data(indw:inde,inds:indn);
        [nEKB_mean, error_status] = Func_0011_get_area_weighted_mean(nEKB_data, nEKB_lon, nEKB_lat)
    
        [EKB_mean, error_status] = Func_0011_get_area_weighted_mean(tmp.mean_data, RCM_grid.lon_rho, RCM_grid.lat_rho)

        saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);

        disp(' ')
        disp([num2str(flagi), ' plot is created.'])
        disp(' ')
        disp([' File path is : ',tmp.tifname])
        disp(' ')
        close all;
    end
    
    %% diff plot
    tmp.tifname=strcat(dirs.figdir_allmean, 'all_mean','_',tmp.regionname, '_sp_ground_', ...
                'diff_', num2str(min(RCM_info.years_ssp),'%04i'),'_',num2str(max(RCM_info.years_ssp),'%04i'), 'y_', ...
                num2str(min(RCM_info.years_his),'%04i'),'_',num2str(max(RCM_info.years_his),'%04i'), 'y_', ...
                RCM_info.season, '.tif'); 

    if (exist(tmp.tifname , 'file') ~= 2 || flags.fig_switch(flagi)==2)

        run(tmp.param_script);
        
        tmp.mean_data= tmp.mean_data_ssp-tmp.mean_data_his;
        
        RCM_grid.mask_EKBcoast = double(inpolygon(RCM_grid.lon_rho,RCM_grid.lat_rho,[127 127 132 132],tmp.refpolygon(:,2)));
        RCM_grid.mask_EKBcoast(RCM_grid.mask_EKBcoast==0)=NaN;
        tmp.mean_data=tmp.mean_data.* RCM_grid.mask_EKBcoast;

    %% land fig    
        RCM_grid.mask_model2 = double(inpolygon(RCM_grid.lon_rho,RCM_grid.lat_rho,tmp.refpolygon(:,1),tmp.refpolygon(:,2)));
        RCM_grid.mask_model2(RCM_grid.mask_model2==0)=NaN;
        tmp.mean_data = tmp.mean_data .* RCM_grid.mask_model2;
        tmp.mean_data_ssp=tmp.mean_data_ssp .* RCM_grid.mask_model2;
        tmp.mean_data_his=tmp.mean_data_his  .* RCM_grid.mask_model2;

        [tmp.mean_egg_mask, tmp.error_status] = Func_0011_get_area_weighted_mean(tmp.mean_data, RCM_grid.lon_rho, RCM_grid.lat_rho);


        tmp.testnameind=1;
        m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);

        ax{tmp.testnameind,1}=axes;

        if strcmp(tmp.testname_prefix, 'prob')
            tmp.temp_surf=ncread(tmp.ncname, 'temp_surf', [1 1], [inf inf]);
        else
            tmp.temp_surf=ncread(tmp.ncname, 'temp_surf', [1 1 1], [inf inf 1]);
        end

        tmp.temp_surf(isnan(tmp.temp_surf))=50000;
        tmp.temp_surf(tmp.temp_surf<50000)=NaN;
        RCM_grid.model_land=tmp.temp_surf;
        RCM_grid.model_land(tmp.temp_surf==50000)=1;
        pc{tmp.testnameind,1}=m_pcolor(RCM_grid.lon_rho',RCM_grid.lat_rho', RCM_grid.model_land','parent',ax{tmp.testnameind,1});
        colormap(ax{tmp.testnameind,1},[0.8 0.8 0.8]);
        shading(gca,param.m_pcolor_shading_method); 

        m_grid('fontsize', param.m_grid_fontsize, 'tickdir', param.m_grid_tickdir_type, 'box', param.m_grid_box_type,  ...
               'parent', ax{tmp.testnameind,1});

        pos_ax{tmp.testnameind}=get(ax{tmp.testnameind,1}, 'pos');
        col_bar{tmp.testnameind,1}=colorbar;
        set(col_bar{tmp.testnameind,1}, 'TickLabels', []);
        hold on

    %% data fig
        m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);
        ax{tmp.testnameind,2}=axes;
        pc{tmp.testnameind,2}=m_pcolor(RCM_grid.lon_rho',RCM_grid.lat_rho', tmp.mean_data','parent',ax{tmp.testnameind,2});
    %     colormap(ax{tmp.testnameind,2},jet);
        colormap(ax{tmp.testnameind,2},cmaps.bwr_10);

        caxis([-1, 1]);
        shading(gca,param.m_pcolor_shading_method);   

        m_grid('fontsize', param.m_grid_fontsize, 'tickdir', param.m_grid_tickdir_type, 'box', param.m_grid_box_type,  ...
               'backcolor', 'none','parent', ax{tmp.testnameind,2});
        col_bar{tmp.testnameind,2}=colorbar;
        set(col_bar{tmp.testnameind,2}, 'fontsize',param.m_grid_fontsize+5);

        tmp.titlename = strcat('sp ground, ','allmean', ',(', ...
            'diff', ',',  ...
            RCM_info.season, ')'); 

        title(tmp.titlename,'fontsize',param.m_pcolor_title_fontsize, 'parent', ax{tmp.testnameind,1});  %%title for land figure
        title(tmp.titlename,'fontsize',param.m_pcolor_title_fontsize, 'parent', ax{tmp.testnameind,2});  %%title for data figure

        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [param.hor_paper_size_x, param.hor_paper_size_y]);
        set(gcf,'PaperPosition', [param.paper_position_hor param.paper_position_ver param.paper_position_width param.paper_position_height]) 

        [tmp.indw, tmp.inde, tmp.inds, tmp.indn]=Func_0012_findind_Y(1/10,[127, 130, 36, 39],RCM_grid.lon_rho,RCM_grid.lat_rho); % southern EKB
        sEKB_lon=RCM_grid.lon_rho(tmp.indw:tmp.inde,tmp.inds:tmp.indn);
        sEKB_lat=RCM_grid.lat_rho(tmp.indw:tmp.inde,tmp.inds:tmp.indn);
        sEKB_data=tmp.mean_data(tmp.indw:tmp.inde,tmp.inds:tmp.indn);
        [sEKB_mean, error_status] = Func_0011_get_area_weighted_mean(sEKB_data, sEKB_lon, sEKB_lat)
        
        [indw, inde, inds, indn]=Func_0012_findind_Y(1/10,[127, 130, 39, 41],RCM_grid.lon_rho,RCM_grid.lat_rho); % southern EKB
        nEKB_lon=RCM_grid.lon_rho(indw:inde,inds:indn);
        nEKB_lat=RCM_grid.lat_rho(indw:inde,inds:indn);
        nEKB_data=tmp.mean_data(indw:inde,inds:indn);
        [nEKB_mean, error_status] = Func_0011_get_area_weighted_mean(nEKB_data, nEKB_lon, nEKB_lat)
        
        [EKB_mean, error_status] = Func_0011_get_area_weighted_mean(tmp.mean_data, RCM_grid.lon_rho, RCM_grid.lat_rho)
        [EKB_mean_ssp, error_status] = Func_0011_get_area_weighted_mean(tmp.mean_data_ssp, RCM_grid.lon_rho, RCM_grid.lat_rho)
        [EKB_mean_his, error_status] = Func_0011_get_area_weighted_mean(tmp.mean_data_his, RCM_grid.lon_rho, RCM_grid.lat_rho)

%         tmp.mean_data_ssp-tmp.mean_data_his;
        
        saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);

        disp(' ')
        disp([num2str(flagi), ' plot is created.'])
        disp(' ')
        disp([' File path is : ',tmp.tifname])
        disp(' ')
        close all;
    end
    
    %% diff plot (%)
    tmp.tifname=strcat(dirs.figdir_allmean, 'all_mean','_',tmp.regionname, '_sp_ground_', ...
                'diff_percent_', num2str(min(RCM_info.years_ssp),'%04i'),'_',num2str(max(RCM_info.years_ssp),'%04i'), 'y_', ...
                num2str(min(RCM_info.years_his),'%04i'),'_',num2str(max(RCM_info.years_his),'%04i'), 'y_', ...
                RCM_info.season, '.tif'); 

    if (exist(tmp.tifname , 'file') ~= 2 || flags.fig_switch(flagi)==2)

        run(tmp.param_script);

        tmp.mean_data= tmp.mean_data_ssp-tmp.mean_data_his;
        tmp.mean_data=(tmp.mean_data./tmp.mean_data_his).*100;
        
        RCM_grid.mask_EKBcoast = double(inpolygon(RCM_grid.lon_rho,RCM_grid.lat_rho,[127 127 132 132],tmp.refpolygon(:,2)));
        RCM_grid.mask_EKBcoast(RCM_grid.mask_EKBcoast==0)=NaN;
        tmp.mean_data=tmp.mean_data.* RCM_grid.mask_EKBcoast;

    %% land fig    
        RCM_grid.mask_model2 = double(inpolygon(RCM_grid.lon_rho,RCM_grid.lat_rho,tmp.refpolygon(:,1),tmp.refpolygon(:,2)));
        RCM_grid.mask_model2(RCM_grid.mask_model2==0)=NaN;
        tmp.mean_data = tmp.mean_data .* RCM_grid.mask_model2;

        [tmp.mean_egg_mask, tmp.error_status] = Func_0011_get_area_weighted_mean(tmp.mean_data, RCM_grid.lon_rho, RCM_grid.lat_rho);


        tmp.testnameind=1;
        m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);

        ax{tmp.testnameind,1}=axes;

        if strcmp(tmp.testname_prefix, 'prob')
            tmp.temp_surf=ncread(tmp.ncname, 'temp_surf', [1 1], [inf inf]);
        else
            tmp.temp_surf=ncread(tmp.ncname, 'temp_surf', [1 1 1], [inf inf 1]);
        end

        tmp.temp_surf(isnan(tmp.temp_surf))=50000;
        tmp.temp_surf(tmp.temp_surf<50000)=NaN;
        RCM_grid.model_land=tmp.temp_surf;
        RCM_grid.model_land(tmp.temp_surf==50000)=1;
        pc{tmp.testnameind,1}=m_pcolor(RCM_grid.lon_rho',RCM_grid.lat_rho', RCM_grid.model_land','parent',ax{tmp.testnameind,1});
        colormap(ax{tmp.testnameind,1},[0.8 0.8 0.8]);
        shading(gca,param.m_pcolor_shading_method); 

        m_grid('fontsize', param.m_grid_fontsize, 'tickdir', param.m_grid_tickdir_type, 'box', param.m_grid_box_type,  ...
               'parent', ax{tmp.testnameind,1});

        pos_ax{tmp.testnameind}=get(ax{tmp.testnameind,1}, 'pos');
        col_bar{tmp.testnameind,1}=colorbar;
        set(col_bar{tmp.testnameind,1}, 'TickLabels', []);
        hold on

    %% data fig
        m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);
        ax{tmp.testnameind,2}=axes;
        pc{tmp.testnameind,2}=m_pcolor(RCM_grid.lon_rho',RCM_grid.lat_rho', tmp.mean_data','parent',ax{tmp.testnameind,2});
    %     colormap(ax{tmp.testnameind,2},jet);
%         colormap(ax{tmp.testnameind,2},cmaps.bw_10);
        colormap(ax{tmp.testnameind,2},cmaps.bwr_10);

%         caxis([-100, 0]);
        caxis([-100, 100]);

        shading(gca,param.m_pcolor_shading_method);   

        m_grid('fontsize', param.m_grid_fontsize, 'tickdir', param.m_grid_tickdir_type, 'box', param.m_grid_box_type,  ...
               'backcolor', 'none','parent', ax{tmp.testnameind,2});
        col_bar{tmp.testnameind,2}=colorbar;
        set(col_bar{tmp.testnameind,2}, 'fontsize',param.m_grid_fontsize+5);
        title(col_bar{tmp.testnameind,2},'%','fontsize',param.colorbar_title_fontsize);
        
        tmp.titlename = strcat('sp ground, ','allmean', ',(', ...
            'diff', ',',  ...
            RCM_info.season, ')'); 

        title(tmp.titlename,'fontsize',param.m_pcolor_title_fontsize, 'parent', ax{tmp.testnameind,1});  %%title for land figure
        title(tmp.titlename,'fontsize',param.m_pcolor_title_fontsize, 'parent', ax{tmp.testnameind,2});  %%title for data figure

        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [param.hor_paper_size_x, param.hor_paper_size_y]);
        set(gcf,'PaperPosition', [param.paper_position_hor param.paper_position_ver param.paper_position_width param.paper_position_height]) 

        [tmp.indw, tmp.inde, tmp.inds, tmp.indn]=Func_0012_findind_Y(1/10,[128, 130, 36, 39],RCM_grid.lon_rho,RCM_grid.lat_rho); % southern EKB
        sEKB_lon=RCM_grid.lon_rho(tmp.indw:tmp.inde,tmp.inds:tmp.indn);
        sEKB_lat=RCM_grid.lat_rho(tmp.indw:tmp.inde,tmp.inds:tmp.indn);
        sEKB_data=tmp.mean_data(tmp.indw:tmp.inde,tmp.inds:tmp.indn);
        [sEKB_mean, error_status] = Func_0011_get_area_weighted_mean(sEKB_data, sEKB_lon, sEKB_lat)
        
        [tmp.indw, tmp.inde, tmp.inds, tmp.indn]=Func_0012_findind_Y(1/10,[128, 130, 36, 38.62],RCM_grid.lon_rho,RCM_grid.lat_rho); % ssouthern EKB
        ssEKB_lon=RCM_grid.lon_rho(tmp.indw:tmp.inde,tmp.inds:tmp.indn);
        ssEKB_lat=RCM_grid.lat_rho(tmp.indw:tmp.inde,tmp.inds:tmp.indn);
        ssEKB_data=tmp.mean_data(tmp.indw:tmp.inde,tmp.inds:tmp.indn);
        [ssEKB_mean, error_status] = Func_0011_get_area_weighted_mean(ssEKB_data, ssEKB_lon, ssEKB_lat)
        
        [indw, inde, inds, indn]=Func_0012_findind_Y(1/10,[128, 130, 39, 41],RCM_grid.lon_rho,RCM_grid.lat_rho); % southern EKB
        nEKB_lon=RCM_grid.lon_rho(indw:inde,inds:indn);
        nEKB_lat=RCM_grid.lat_rho(indw:inde,inds:indn);
        nEKB_data=tmp.mean_data(indw:inde,inds:indn);
        [nEKB_mean, error_status] = Func_0011_get_area_weighted_mean(nEKB_data, nEKB_lon, nEKB_lat)
        
        [EKB_mean, error_status] = Func_0011_get_area_weighted_mean(tmp.mean_data, RCM_grid.lon_rho, RCM_grid.lat_rho)
        
        
        saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);

        disp(' ')
        disp([num2str(flagi), ' plot is created.'])
        disp(' ')
        disp([' File path is : ',tmp.tifname])
        disp(' ')
        close all;
    end
    
    
    



end