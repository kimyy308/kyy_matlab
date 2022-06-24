% %  Updated    24-May-2022 by Yong-Yub Kim   % make

if testnameind==1
    
%     for checkti=1:length(RCM_info.checktime)
%         tmp.checktime=RCM_info.checktime(checkti);

        %% historical plot
        tmp.tifname=strcat(dirs.figdir_allmean, 'all_mean','_',tmp.regionname, '_spawning_std_', ...
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
                        
%                         tmp.all_checktime=ncread(tmp.ncname, 'checktime');
%                         tmp.ind_checktime=find(tmp.all_checktime==tmp.checktime);
                    
                        tmp.egg_mask=squeeze(ncread(tmp.ncname, 'egg_mask'));
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
                [size(tmp.comb_egg_mask,1) size(tmp.comb_egg_mask,2), size(tmp.comb_egg_mask,3)/5, 5]);
            tmp.comb_egg_mask_testsep=squeeze(mean(tmp.comb_egg_mask_testsep,3));
%             tmp.comb_egg_mask_testsep=squeeze(sum(tmp.comb_egg_mask_testsep,3));
            tmp.mean_data = mean(tmp.comb_egg_mask,3);
            tmp.std_data = 1*std(tmp.comb_egg_mask_testsep,0,3);
            tmp.mean_data(tmp.mean_data==0)=NaN;
            tmp.agreement_data=NaN(size(tmp.mean_data));
            tmp.agreement_data(tmp.mean_data>tmp.std_data)=1;
            tmp.std_data(isnan(tmp.mean_data))=NaN;
            
            RCM_grid.mask_EKBcoast = double(inpolygon(RCM_grid.lon_rho,RCM_grid.lat_rho,[127 127 132 132],tmp.refpolygon(:,2)));
            RCM_grid.mask_EKBcoast(RCM_grid.mask_EKBcoast==0)=NaN;
            tmp.mean_data=tmp.mean_data.* RCM_grid.mask_EKBcoast;
            tmp.std_data=tmp.std_data.* RCM_grid.mask_EKBcoast;
            tmp.agreement_data=tmp.agreement_data.* RCM_grid.mask_EKBcoast;
        
%% sig_mask            
%             tmp.sp_mean=mean(tmp.mean_data(:), 'omitnan');
%             tmp.sp_std=std(tmp.mean_data(:), 'omitnan');
%             tmp.sp_std_mean=mean(tmp.std_data(:), 'omitnan');
%             tmp.sig_mask=NaN(size(tmp.mean_data));
%             tmp.sig_mask(tmp.mean_data > (tmp.sp_mean-tmp.sp_std)) = 1;
%             
%             tmp.agreement_data=tmp.agreement_data.*tmp.sig_mask;
%             tmp.std_data=tmp.std_data.*tmp.sig_mask;
%             tmp.mean_data=tmp.mean_data.*tmp.sig_mask;
            
%             pcolor(tmp.sig_mask'); shading flat; colorbar;   

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
            pc{tmp.testnameind,2}=m_pcolor(RCM_grid.lon_rho',RCM_grid.lat_rho', tmp.std_data',...
                'parent',ax{tmp.testnameind,2});
%             colormap(ax{tmp.testnameind,2},jet);
%             colormap(ax{tmp.testnameind,2},parula);
            colormap(ax{tmp.testnameind,2},cmaps.yrmap);

            caxis([0, 0.5]);
            shading(gca,param.m_pcolor_shading_method);   
%             hold off


%% dot
            tmp.dot_lon_agr=RCM_grid.lon_rho(tmp.agreement_data==1);
            tmp.dot_lat_agr=RCM_grid.lat_rho(tmp.agreement_data==1);
            for agri=1:1:length(tmp.dot_lon_agr)
                hold on
                tmp.dot=m_plot(tmp.dot_lon_agr(agri), ...
                tmp.dot_lat_agr(agri), ...
                 'marker','o','color','k', ...
                  'markerfacecolor', 'none', 'parent',ax{tmp.testnameind,2}, 'markersize', 1.5); 
            end

%% contour
%             hold on
%             tmp.agreement_data2=tmp.agreement_data;
%             tmp.agreement_data2(isnan(tmp.agreement_data2))=0;
%             [C,h2]=m_contour(RCM_grid.lon_rho',RCM_grid.lat_rho', tmp.agreement_data2', [0, 1], 'color','k', ...
%                     'linewidth', 2, 'linestyle', '-', 'parent',ax{tmp.testnameind,2});

%% hatch
%             tmp.dot_lon_agr=RCM_grid.lon_rho(tmp.agreement_data==1);
%             tmp.dot_lat_agr=RCM_grid.lat_rho(tmp.agreement_data==1);
%             for agri=1:2:length(tmp.dot_lon_agr)
%                 hold on
%                 tmp.dot=m_hatch(tmp.dot_lon_agr(agri), ...
%                 tmp.dot_lat_agr(agri), ...
%                  'single', ...
%                   'parent',ax{tmp.testnameind,2}); 
%             end
%             tmp.dot=m_hatch(tmp.dot_lon_agr, ...
%                 tmp.dot_lat_agr, ...
%                  'single', ...
%                   'parent',ax{tmp.testnameind,2}, 'color', 'k'); 

            m_grid('fontsize', param.m_grid_fontsize, 'tickdir', param.m_grid_tickdir_type, 'box', param.m_grid_box_type,  ...
                   'backcolor', 'none','parent', ax{tmp.testnameind,2});
            col_bar{tmp.testnameind,2}=colorbar;
            set(col_bar{tmp.testnameind,2}, 'fontsize',param.m_grid_fontsize+5);
            
            %             tmp.dot=m_plot(RCM_grid.lon_rho(tmp.agreement_data==1), ...
%                 RCM_grid.lat_rho(tmp.agreement_data==1), ...
%                 tmp.agreement_data(tmp.agreement_data==1), 'marker','o','color','y', ...
%                   'markerfacecolor',bwrmap(colind,:), 'parent',ax{tmp.testnameind,2}); 
            
            

            
            
              
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

            saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);

            disp(' ')
            disp([num2str(flagi), ' plot is created.'])
            disp(' ')
            disp([' File path is : ',tmp.tifname])
            disp(' ')
            close all;
        end



        %% future(ssp) plot
        tmp.tifname=strcat(dirs.figdir_allmean, 'all_mean','_',tmp.regionname, '_spawning_std_', ...
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
                                            
                        tmp.egg_mask=squeeze(ncread(tmp.ncname, 'egg_mask'));
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

            tmp.mean_data(tmp.mean_data==0)=NaN;
            
            tmp.comb_egg_mask_testsep=reshape(tmp.comb_egg_mask, ...
                [size(tmp.comb_egg_mask,1) size(tmp.comb_egg_mask,2), size(tmp.comb_egg_mask,3)/5, 5]);
            tmp.comb_egg_mask_testsep=squeeze(mean(tmp.comb_egg_mask_testsep, 3));
%             tmp.comb_egg_mask_testsep=squeeze(sum(tmp.comb_egg_mask_testsep,3));
            tmp.std_data = 1*std(tmp.comb_egg_mask_testsep,0,3);
            tmp.mean_data(tmp.mean_data==0)=NaN;
            tmp.agreement_data=NaN(size(tmp.mean_data));
            tmp.agreement_data(tmp.mean_data>tmp.std_data)=1;
            tmp.std_data(isnan(tmp.mean_data))=NaN;
            
            
            RCM_grid.mask_EKBcoast = double(inpolygon(RCM_grid.lon_rho,RCM_grid.lat_rho,[127 127 132 132],tmp.refpolygon(:,2)));
            RCM_grid.mask_EKBcoast(RCM_grid.mask_EKBcoast==0)=NaN;
            tmp.mean_data=tmp.mean_data.* RCM_grid.mask_EKBcoast;
            tmp.std_data=tmp.std_data.* RCM_grid.mask_EKBcoast;
            tmp.agreement_data=tmp.agreement_data.* RCM_grid.mask_EKBcoast;
            
            
%% sig_mask            
%             tmp.sp_mean=mean(tmp.mean_data(:), 'omitnan');
%             tmp.sp_std=std(tmp.mean_data(:), 'omitnan');
%             tmp.sp_std_mean=mean(tmp.std_data(:), 'omitnan');
%             tmp.sig_mask=NaN(size(tmp.mean_data));
%             tmp.sig_mask(tmp.mean_data > (tmp.sp_mean-tmp.sp_std)) = 1;
%             
%             tmp.agreement_data=tmp.agreement_data.*tmp.sig_mask;
%             tmp.std_data=tmp.std_data.*tmp.sig_mask;
%             tmp.mean_data=tmp.mean_data.*tmp.sig_mask;
            
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
            pc{tmp.testnameind,2}=m_pcolor(RCM_grid.lon_rho',RCM_grid.lat_rho', tmp.std_data','parent',ax{tmp.testnameind,2});
%             colormap(ax{tmp.testnameind,2},jet);
%             colormap(ax{tmp.testnameind,2},parula);
            colormap(ax{tmp.testnameind,2},cmaps.yrmap);

            

            caxis([0, 0.5]);
            shading(gca,param.m_pcolor_shading_method);   
            
%% dot
            tmp.dot_lon_agr=RCM_grid.lon_rho(tmp.agreement_data==1);
            tmp.dot_lat_agr=RCM_grid.lat_rho(tmp.agreement_data==1);
            for agri=1:1:length(tmp.dot_lon_agr)
                hold on
                tmp.dot=m_plot(tmp.dot_lon_agr(agri), ...
                tmp.dot_lat_agr(agri), ...
                 'marker','o','color','k', ...
                  'markerfacecolor', 'none', 'parent',ax{tmp.testnameind,2}, 'markersize', 1.5); 
            end
            
            m_grid('fontsize', param.m_grid_fontsize, 'tickdir', param.m_grid_tickdir_type, 'box', param.m_grid_box_type,  ...
                   'backcolor', 'none','parent', ax{tmp.testnameind,2});
            col_bar{tmp.testnameind,2}=colorbar;
            set(col_bar{tmp.testnameind,2}, 'fontsize',param.m_grid_fontsize+5);
            
            tmp.titlename = strcat('sp ground, ','allmean', ',(', ...
                num2str(min(RCM_info.years_ssp),'%04i'),'-', num2str(max(RCM_info.years_ssp),'%04i'), ',',  ...
                RCM_info.season, ')'); 
            
            
%             caxis([0, 200]);

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

            %% line width
            
            disp(' ')
            disp([num2str(flagi), ' plot is created.'])
            disp(' ')
            disp([' File path is : ',tmp.tifname])
            disp(' ')
            close all;
        end

%     end

end