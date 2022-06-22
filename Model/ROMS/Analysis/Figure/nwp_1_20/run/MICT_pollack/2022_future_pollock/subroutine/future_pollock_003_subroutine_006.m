% %  Updated    19-May-2022 by Yong-Yub Kim   % make

if testnameind==1
    
    for checkti=1:length(RCM_info.checktime)
        tmp.checktime=RCM_info.checktime(checkti);

        %% historical plot
        tmp.tifname=strcat(dirs.figdir_allmean, 'all_mean','_',tmp.regionname, '_loc_num_', num2str(tmp.checktime, '%02i'),'days', ...
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
                dirs_sub.filedir_his = strcat('D:\Data\Model\ROMS\nwp_1_20\', tmp_sub.testname_his, '\pollock\'); % % where data files are
                dirs_sub.savedir_his = strcat('D:\Data\Model\ROMS\nwp_1_20\', tmp_sub.testname_his, '\pollock\');

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
                        
                        tmp.all_checktime=ncread(tmp.ncname, 'checktime');
                        tmp.ind_checktime=find(tmp.all_checktime==tmp.checktime);
                    
                        tmp.egg_mask=squeeze(ncread(tmp.ncname, 'mask_par', [1 1 1 tmp.ind_checktime], [inf inf inf 1]));
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
            tmp.comb_egg_mask_testsep=reshape(tmp.comb_egg_mask, ...
                [size(tmp.comb_egg_mask,1) size(tmp.comb_egg_mask,2), 5, size(tmp.comb_egg_mask,3)/5]);
            tmp.comb_egg_mask_testsep=squeeze(sum(tmp.comb_egg_mask_testsep,4));
            tmp.mean_data = sum(tmp.comb_egg_mask,3)./tmp.testlen;
            tmp.mean_data_his=tmp.mean_data;
            tmp.std_data = 2*std(tmp.comb_egg_mask_testsep,0,3);
            tmp.mean_data(tmp.mean_data==0)=NaN;
            tmp.agreement_data=NaN(size(tmp.mean_data));
            tmp.agreement_data(tmp.mean_data>tmp.std_data)=1;

%             pcolor(tmp.agreement_data'); shading flat; colorbar;
%             pcolor(tmp.std_data'); shading flat; colorbar;
%             pcolor(tmp.mean_data'); shading flat; colorbar;
%             pcolor(tmp.comb_egg_mask_testsep(:,:,2)'); shading flat; colorbar;
            
%             RCM_grid.mask_EKBcoast = double(inpolygon(RCM_grid.lon_rho,RCM_grid.lat_rho,[127 127 132 132],tmp.refpolygon(:,2)));
%             RCM_grid.mask_EKBcoast(RCM_grid.mask_EKBcoast==0)=NaN;
%             tmp.mean_data=tmp.mean_data.* RCM_grid.mask_EKBcoast;

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
            colormap(ax{tmp.testnameind,2},jet);
%             colormap(ax{tmp.testnameind,2},parula);

            caxis([0, 200]);
            shading(gca,param.m_pcolor_shading_method);   

            m_grid('fontsize', param.m_grid_fontsize, 'tickdir', param.m_grid_tickdir_type, 'box', param.m_grid_box_type,  ...
                   'backcolor', 'none','parent', ax{tmp.testnameind,2});
            col_bar{tmp.testnameind,2}=colorbar;
            set(col_bar{tmp.testnameind,2}, 'fontsize',param.m_grid_fontsize+5);
            
            tmp.titlename = strcat(num2str(tmp.checktime, '%02i'), 'd, egg #, ','allmean', ',(', ...
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

            saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);

            disp(' ')
            disp([num2str(flagi), ' plot is created.'])
            disp(' ')
            disp([' File path is : ',tmp.tifname])
            disp(' ')
            close all;
        end



        %% future(ssp) plot
        tmp.tifname=strcat(dirs.figdir_allmean, 'all_mean','_',tmp.regionname, '_loc_num_', num2str(tmp.checktime, '%02i'),'days', ...
                    num2str(min(RCM_info.years_ssp),'%04i'),'_',num2str(max(RCM_info.years_ssp),'%04i'), 'y_', ...
                    RCM_info.season, '.tif'); 

        tmp.testname_prefix=tmp.testname_ssp(1:4); 
        tmp.testlen=length(RCM_info.name);
        tmp.tlen_ssp=0;

        if (exist(tmp.tifname , 'file') ~= 2 || flags.fig_switch(flagi)==2)

            run(tmp.param_script);

            for sub_testnameind=1:length(RCM_info.name)
                tmp_sub.testname_ssp=RCM_info.name{sub_testnameind};
                dirs_sub.filedir_ssp = strcat('D:\Data\Model\ROMS\nwp_1_20\', tmp_sub.testname_ssp, '\pollock\'); % % where data files are
                dirs_sub.savedir_ssp = strcat('D:\Data\Model\ROMS\nwp_1_20\', tmp_sub.testname_ssp, '\pollock\');

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
                        
                        tmp.all_checktime=ncread(tmp.ncname, 'checktime');
                        tmp.ind_checktime=find(tmp.all_checktime==tmp.checktime);
                    
                        tmp.egg_mask=squeeze(ncread(tmp.ncname, 'mask_par', [1 1 1 tmp.ind_checktime], [inf inf inf 1]));
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
            tmp.mean_data = sum(tmp.comb_egg_mask,3)./tmp.testlen;
            tmp.mean_data_ssp=tmp.mean_data;

            tmp.mean_data(tmp.mean_data==0)=NaN;
            
            tmp.explen=length(RCM_info.name);
            tmp.comb_egg_mask_testsep=reshape(tmp.comb_egg_mask, ...
                [size(tmp.comb_egg_mask,1) size(tmp.comb_egg_mask,2), tmp.explen, size(tmp.comb_egg_mask,3)/tmp.explen]);
            tmp.comb_egg_mask_testsep=squeeze(sum(tmp.comb_egg_mask_testsep,4));
            tmp.std_data = 2*std(tmp.comb_egg_mask_testsep,0,3);
            tmp.agreement_data=NaN(size(tmp.mean_data));
            tmp.agreement_data(tmp.mean_data>tmp.std_data)=1;
            
%             pcolor(tmp.agreement_data'); shading flat; colorbar;
%             pcolor(tmp.std_data'); shading flat; colorbar;
%             pcolor(tmp.mean_data'); shading flat; colorbar;
%             pcolor(tmp.comb_egg_mask_testsep(:,:,2)'); shading flat; colorbar;
            
%             RCM_grid.mask_EKBcoast = double(inpolygon(RCM_grid.lon_rho,RCM_grid.lat_rho,[127 127 132 132],tmp.refpolygon(:,2)));
%             RCM_grid.mask_EKBcoast(RCM_grid.mask_EKBcoast==0)=NaN;
%             tmp.mean_data=tmp.mean_data.* RCM_grid.mask_EKBcoast;

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
            colormap(ax{tmp.testnameind,2},jet);
%             colormap(ax{tmp.testnameind,2},parula);

%             caxis([0, 1.0]);
            shading(gca,param.m_pcolor_shading_method);   

            m_grid('fontsize', param.m_grid_fontsize, 'tickdir', param.m_grid_tickdir_type, 'box', param.m_grid_box_type,  ...
                   'backcolor', 'none','parent', ax{tmp.testnameind,2});
            col_bar{tmp.testnameind,2}=colorbar;
            set(col_bar{tmp.testnameind,2}, 'fontsize',param.m_grid_fontsize+5);
            
            tmp.titlename = strcat(num2str(tmp.checktime, '%02i'), 'd, egg #, ','allmean', ',(', ...
                num2str(min(RCM_info.years_ssp),'%04i'),'-', num2str(max(RCM_info.years_ssp),'%04i'), ',',  ...
                RCM_info.season, ')'); 
            
            
            caxis([0, 200]);

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

            saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);

            disp(' ')
            disp([num2str(flagi), ' plot is created.'])
            disp(' ')
            disp([' File path is : ',tmp.tifname])
            disp(' ')
            close all;
        end
        
       %% diff plot relative to spawning decrease
        tmp.tifname=strcat(dirs.figdir_allmean, 'all_mean','_',tmp.regionname, '_loc_num_', num2str(tmp.checktime, '%02i'),'days_',...
                    'diff_rel_sp_', num2str(min(RCM_info.years_ssp),'%04i'),'_',num2str(max(RCM_info.years_ssp),'%04i'), 'y_', ...
                    num2str(min(RCM_info.years_his),'%04i'),'_',num2str(max(RCM_info.years_his),'%04i'), 'y_', ...
                    RCM_info.season, '.tif'); 

        if (exist(tmp.tifname , 'file') ~= 2 || flags.fig_switch(flagi)==2)

            run(tmp.param_script);

            tmp.mean_data= tmp.mean_data_ssp-tmp.mean_data_his.*(100-76.4)./100;
            tmp.mean_data(tmp.mean_data_ssp==0 & tmp.mean_data_his==0)=NaN;


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
            colormap(ax{tmp.testnameind,2},cmaps.byrmap3);

            caxis([-50, 50]);
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

            saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);

            disp(' ')
            disp([num2str(flagi), ' plot is created.'])
            disp(' ')
            disp([' File path is : ',tmp.tifname])
            disp(' ')
            close all;
        end
        
        %% diff plot
        tmp.tifname=strcat(dirs.figdir_allmean, 'all_mean','_',tmp.regionname, '_loc_num_', num2str(tmp.checktime, '%02i'),'days_',...
                    'diff_', num2str(min(RCM_info.years_ssp),'%04i'),'_',num2str(max(RCM_info.years_ssp),'%04i'), 'y_', ...
                    num2str(min(RCM_info.years_his),'%04i'),'_',num2str(max(RCM_info.years_his),'%04i'), 'y_', ...
                    RCM_info.season, '.tif'); 

        if (exist(tmp.tifname , 'file') ~= 2 || flags.fig_switch(flagi)==2)

            run(tmp.param_script);

            tmp.mean_data= tmp.mean_data_ssp-tmp.mean_data_his;
            tmp.mean_data(tmp.mean_data_ssp==0 & tmp.mean_data_his==0)=NaN;


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
            colormap(ax{tmp.testnameind,2},cmaps.byrmap3);

            caxis([-200, 200]);
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

            saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);

            disp(' ')
            disp([num2str(flagi), ' plot is created.'])
            disp(' ')
            disp([' File path is : ',tmp.tifname])
            disp(' ')
            close all;
        end
        
        
       %% diff plot (%)
        tmp.tifname=strcat(dirs.figdir_allmean, 'all_mean','_',tmp.regionname, '_loc_num_', num2str(tmp.checktime, '%02i'),'days_',...
                    'diff_percent_', num2str(min(RCM_info.years_ssp),'%04i'),'_',num2str(max(RCM_info.years_ssp),'%04i'), 'y_', ...
                    num2str(min(RCM_info.years_his),'%04i'),'_',num2str(max(RCM_info.years_his),'%04i'), 'y_', ...
                    RCM_info.season, '.tif'); 

        if (exist(tmp.tifname , 'file') ~= 2 || flags.fig_switch(flagi)==2)

            run(tmp.param_script);

            tmp.mean_data= tmp.mean_data_ssp-tmp.mean_data_his;
            tmp.mean_data(tmp.mean_data_ssp==0 & tmp.mean_data_his==0)=NaN;
            tmp.mean_data =(tmp.mean_data./tmp.mean_data_his).*100;


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
%             colormap(ax{tmp.testnameind,2},cmaps.byrmap);
            
            colormap(ax{tmp.testnameind,2},cmaps.bymap3);
            caxis([-100, 0]);
            if strcmp(RCM_info.name{1}, 't17s_t27c')
                colormap(ax{tmp.testnameind,2},cmaps.byrmap);
                caxis([-100, 100]);
            end
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
        
        
       %% diff plot (%), anomaly relative with (EKB, 76.4%)
        tmp.tifname=strcat(dirs.figdir_allmean, 'all_mean','_',tmp.regionname, '_loc_num_', num2str(tmp.checktime, '%02i'),'days_',...
                    'diff_percent_rel_sp_', num2str(min(RCM_info.years_ssp),'%04i'),'_',num2str(max(RCM_info.years_ssp),'%04i'), 'y_', ...
                    num2str(min(RCM_info.years_his),'%04i'),'_',num2str(max(RCM_info.years_his),'%04i'), 'y_', ...
                    RCM_info.season, '.tif'); 

        if (exist(tmp.tifname , 'file') ~= 2 || flags.fig_switch(flagi)==2)

            run(tmp.param_script);

%             tmp.mean_data= tmp.mean_data_ssp./(100-76.4).*100-tmp.mean_data_his;
            tmp.mean_data= tmp.mean_data_ssp-tmp.mean_data_his.*(100-76.4)./100;
            tmp.mean_data(tmp.mean_data_ssp==0 & tmp.mean_data_his==0)=NaN;
            tmp.mean_data =(tmp.mean_data./tmp.mean_data_his).*100;


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
            tmp.spawn_diff_avg=76.4340;
%             pc{tmp.testnameind,2}=m_pcolor(RCM_grid.lon_rho',RCM_grid.lat_rho', tmp.mean_data'+tmp.spawn_diff_avg,'parent',ax{tmp.testnameind,2});
            pc{tmp.testnameind,2}=m_pcolor(RCM_grid.lon_rho',RCM_grid.lat_rho', tmp.mean_data','parent',ax{tmp.testnameind,2});

        %     colormap(ax{tmp.testnameind,2},jet);
%             colormap(ax{tmp.testnameind,2},cmaps.bymap3);
            colormap(ax{tmp.testnameind,2},cmaps.byrmap);

            caxis([-25, 25]);
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
        
        
    %% diff plot (%), anomaly relative with (korean fishing area, <38.6N, 89.6%)
        tmp.tifname=strcat(dirs.figdir_allmean, 'all_mean','_',tmp.regionname, '_loc_num_', num2str(tmp.checktime, '%02i'),'days_',...
                    'diff_percent_rel_sp2_', num2str(min(RCM_info.years_ssp),'%04i'),'_',num2str(max(RCM_info.years_ssp),'%04i'), 'y_', ...
                    num2str(min(RCM_info.years_his),'%04i'),'_',num2str(max(RCM_info.years_his),'%04i'), 'y_', ...
                    RCM_info.season, '.tif'); 

        if (exist(tmp.tifname , 'file') ~= 2 || flags.fig_switch(flagi)==2)

            run(tmp.param_script);

%             tmp.mean_data= tmp.mean_data_ssp./(100-76.4).*100-tmp.mean_data_his;
            tmp.mean_data= tmp.mean_data_ssp-tmp.mean_data_his.*(100-89.6)./100;
            tmp.mean_data(tmp.mean_data_ssp==0 & tmp.mean_data_his==0)=NaN;
            tmp.mean_data =(tmp.mean_data./tmp.mean_data_his).*100;


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
            tmp.spawn_diff_avg=76.4340;
%             pc{tmp.testnameind,2}=m_pcolor(RCM_grid.lon_rho',RCM_grid.lat_rho', tmp.mean_data'+tmp.spawn_diff_avg,'parent',ax{tmp.testnameind,2});
            pc{tmp.testnameind,2}=m_pcolor(RCM_grid.lon_rho',RCM_grid.lat_rho', tmp.mean_data','parent',ax{tmp.testnameind,2});

        %     colormap(ax{tmp.testnameind,2},jet);
%             colormap(ax{tmp.testnameind,2},cmaps.bymap3);
            colormap(ax{tmp.testnameind,2},cmaps.byrmap);

            caxis([-25, 25]);
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

end