% %  Updated    18-Apr-2022 by Yong-Yub Kim   %
% %  Updated    03-May-2022 by Yong-Yub Kim   % structure
% %  Updated    16-May-2022 by Yong-Yub Kim   % add reading prob_ens???? test

for checkti=1:length(RCM_info.checktime)
    tmp.checktime=RCM_info.checktime(checkti);

    %% historical plot
    tmp.tifname=strcat(dirs.figdir, tmp.testname_his,'_',tmp.regionname, '_loc_num_', num2str(tmp.checktime, '%02i'),'days', ...
                num2str(min(RCM_info.years_his),'%04i'),'_',num2str(max(RCM_info.years_his),'%04i'), 'y', ...
                RCM_info.season, '.tif'); %% ~_year_month.jpg
          
            
    tmp.testname_prefix=tmp.testname_his(1:4);        

    if (exist(tmp.tifname , 'file') ~= 2 || flags.fig_switch(flagi)==2)
        run(tmp.param_script);
    %%    initialization
        if(isfield(tmp, 'tlen_his')~=1)
            tmp.tlen_his=0;
            if strcmp(tmp.testname_prefix, 'prob')
                for yearij = 1:length(RCM_info.years_his)
                    tmp.tempyear = RCM_info.years_his(yearij);
                    tmp.ncname = [dirs.savedir_his,tmp.testname_his,'_',tmp.regionname,'model_pollock_',num2str(tmp.tempyear,'%04i'),'.nc'];
                    tmp.tlen_his=tmp.tlen_his + length(ncread(tmp.ncname, 'time'));
                end
            else
                for yearij = 1:length(RCM_info.years_his)
                    tmp.tempyear = RCM_info.years_his(yearij);
                    for monthij = 1:length(RCM_info.months)
                        tmp.tempmonth = RCM_info.months(monthij);
                        tmp.ncname = [dirs.savedir_his,tmp.testname_his,'_',tmp.regionname,'model_pollock_',num2str(tmp.tempyear,'%04i'),'_',num2str(tmp.tempmonth,'%02i'),'.nc'];
                        tmp.tlen_his=tmp.tlen_his + length(ncread(tmp.ncname, 'time'));
                    end
                end
            end
            
            tmp.xlen=size(ncread(tmp.ncname,'lon_rho'),1);
            tmp.ylen=size(ncread(tmp.ncname,'lon_rho'),2);
        end
        tmp.comb_egg_mask=NaN(tmp.xlen,tmp.ylen,tmp.tlen_his);

    %% read grids, eggs    
        if strcmp(tmp.testname_prefix, 'prob')
            for yearij = 1:length(RCM_info.years_his)
                tmp.tempyear = RCM_info.years_his(yearij);
                tmp.ncname = [dirs.savedir_his,tmp.testname_his,'_',tmp.regionname,'model_pollock_',num2str(tmp.tempyear,'%04i'),'.nc'];
                disp([num2str(yearij), 'y'])
                tmp.all_checktime=ncread(tmp.ncname, 'checktime');
                tmp.ind_checktime=find(tmp.all_checktime==tmp.checktime);
                tmp.egg_mask=ncread(tmp.ncname, 'mask_par', [1 1 1 tmp.ind_checktime], [inf inf inf 1]);
                tmp.lastday_m=size(tmp.egg_mask,3);
                tmp.comb_egg_mask(:,:,1:tmp.lastday_m)=tmp.egg_mask;
                tmp.endij=tmp.lastday_m;
            end
        else
            for yearij = 1:length(RCM_info.years_his)
                tmp.tempyear = RCM_info.years_his(yearij);
                for monthij = 1:length(RCM_info.months)
                    tmp.tempmonth = RCM_info.months(monthij);
                    tmp.ncname = [dirs.savedir_his,tmp.testname_his,'_',tmp.regionname,'model_pollock_',num2str(tmp.tempyear,'%04i'),'_',num2str(tmp.tempmonth,'%02i'),'.nc'];
                    disp([num2str(yearij), 'y_',num2str(monthij),'m'])

                    tmp.all_checktime=ncread(tmp.ncname, 'checktime');
                    tmp.ind_checktime=find(tmp.all_checktime==tmp.checktime);

                    tmp.egg_mask=squeeze(ncread(tmp.ncname, 'mask_par', [1 1 1 tmp.ind_checktime], [inf inf inf 1]));
                    tmp.lastday_m=size(tmp.egg_mask,3);
                    if yearij==1 && monthij==1
                        tmp.comb_egg_mask(:,:,1:tmp.lastday_m)=tmp.egg_mask;
                        tmp.endij=tmp.lastday_m;
                    else
                        tmp.comb_egg_mask(:,:,tmp.endij+1:tmp.endij+tmp.lastday_m)=tmp.egg_mask;
                        tmp.endij=tmp.endij+tmp.lastday_m;
                    end
                end
            end
        end
        

        RCM_grid.lon_rho = ncread(tmp.ncname, 'lon_rho');
        RCM_grid.lat_rho = ncread(tmp.ncname, 'lat_rho');

        tmp.mean_data = sum(tmp.comb_egg_mask,3);

        tmp.mean_data(tmp.mean_data==0)=NaN;


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
        colormap(ax{tmp.testnameind,2},jet(10));
        caxis([0, 200]);
        shading(gca,param.m_pcolor_shading_method);   

        m_grid('fontsize', param.m_grid_fontsize, 'tickdir', param.m_grid_tickdir_type, 'box', param.m_grid_box_type,  ...
               'backcolor', 'none','parent', ax{tmp.testnameind,2});
        col_bar{tmp.testnameind,2}=colorbar;
        set(col_bar{tmp.testnameind,2}, 'fontsize',param.m_grid_fontsize+5);

        tmp.titlename = strcat(num2str(tmp.checktime, '%02i'), 'd, egg #, ',tmp.testname_his, ',(', ...
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

    tmp.tifname=strcat(dirs.figdir, tmp.testname_ssp,'_',tmp.regionname, '_loc_num_', num2str(tmp.checktime, '%02i'),'days', ...
                num2str(min(RCM_info.years_ssp),'%04i'),'_',num2str(max(RCM_info.years_ssp),'%04i'), 'y', ...
                RCM_info.season, '.tif'); %% ~_year_month.jpg
            
    if (exist(tmp.tifname , 'file') ~= 2 || flags.fig_switch(flagi)==2)
        run(tmp.param_script);
    %%    initialization
        if(isfield(tmp, 'tlen_ssp')~=1)
            tmp.tlen_ssp=0;
            if strcmp(tmp.testname_prefix, 'prob')
                for yearij = 1:length(RCM_info.years_ssp)
                    tmp.tempyear = RCM_info.years_ssp(yearij);
                        tmp.ncname = [dirs.savedir_ssp,tmp.testname_ssp,'_',tmp.regionname,'model_pollock_',num2str(tmp.tempyear,'%04i'),'.nc'];
                        tmp.tlen_ssp=tmp.tlen_ssp + length(ncread(tmp.ncname, 'time'));
                end
            else
                for yearij = 1:length(RCM_info.years_ssp)
                    tmp.tempyear = RCM_info.years_ssp(yearij);
                    for monthij = 1:length(RCM_info.months)
                        tmp.tempmonth = RCM_info.months(monthij);
                        tmp.ncname = [dirs.savedir_ssp,tmp.testname_ssp,'_',tmp.regionname,'model_pollock_',num2str(tmp.tempyear,'%04i'),'_',num2str(tmp.tempmonth,'%02i'),'.nc'];
                        tmp.tlen_ssp=tmp.tlen_ssp + length(ncread(tmp.ncname, 'time'));
                    end
                end
            end            
            
            tmp.xlen=size(ncread(tmp.ncname,'lon_rho'),1);
            tmp.ylen=size(ncread(tmp.ncname,'lon_rho'),2);
        end
        tmp.comb_egg_mask=NaN(tmp.xlen,tmp.ylen,tmp.tlen_ssp);

    %% read grids, eggs    
        if strcmp(tmp.testname_prefix, 'prob')
            for yearij = 1:length(RCM_info.years_ssp)
                tmp.tempyear = RCM_info.years_ssp(yearij);
                tmp.ncname = [dirs.savedir_ssp,tmp.testname_ssp,'_',tmp.regionname,'model_pollock_',num2str(tmp.tempyear,'%04i'),'.nc'];
                disp([num2str(yearij), 'y'])
                tmp.all_checktime=ncread(tmp.ncname, 'checktime');
                tmp.ind_checktime=find(tmp.all_checktime==tmp.checktime);
                tmp.egg_mask=ncread(tmp.ncname, 'mask_par', [1 1 1 tmp.ind_checktime], [inf inf inf 1]);
                tmp.lastday_m=size(tmp.egg_mask,3);
                tmp.comb_egg_mask(:,:,1:tmp.lastday_m)=tmp.egg_mask;
                tmp.endij=tmp.lastday_m;
            end
        else
            for yearij = 1:length(RCM_info.years_ssp)
                tmp.tempyear = RCM_info.years_ssp(yearij);
                for monthij = 1:length(RCM_info.months)
                    tmp.tempmonth = RCM_info.months(monthij);
                    tmp.ncname = [dirs.savedir_ssp,tmp.testname_ssp,'_',tmp.regionname,'model_pollock_',num2str(tmp.tempyear,'%04i'),'_',num2str(tmp.tempmonth,'%02i'),'.nc'];
                    disp([num2str(yearij), 'y_',num2str(monthij),'m'])

                    tmp.all_checktime=ncread(tmp.ncname, 'checktime');
                    tmp.ind_checktime=find(tmp.all_checktime==tmp.checktime);

                    tmp.egg_mask=squeeze(ncread(tmp.ncname, 'mask_par', [1 1 1 tmp.ind_checktime], [inf inf inf 1]));
                    tmp.lastday_m=size(tmp.egg_mask,3);
                    if yearij==1 && monthij==1
                        tmp.comb_egg_mask(:,:,1:tmp.lastday_m)=tmp.egg_mask;
                        tmp.endij=tmp.lastday_m;
                    else
                        tmp.comb_egg_mask(:,:,tmp.endij+1:tmp.endij+tmp.lastday_m)=tmp.egg_mask;
                        tmp.endij=tmp.endij+tmp.lastday_m;
                    end
                end
            end
        end

        RCM_grid.lon_rho = ncread(tmp.ncname, 'lon_rho');
        RCM_grid.lat_rho = ncread(tmp.ncname, 'lat_rho');

        tmp.mean_data = sum(tmp.comb_egg_mask,3);

        tmp.mean_data(tmp.mean_data==0)=NaN;


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
        colormap(ax{tmp.testnameind,2},jet(10));
%         caxis([0.1, 1.0]);
        caxis([0, 200]);

        shading(gca,param.m_pcolor_shading_method);   

        m_grid('fontsize', param.m_grid_fontsize, 'tickdir', param.m_grid_tickdir_type, 'box', param.m_grid_box_type,  ...
               'backcolor', 'none','parent', ax{tmp.testnameind,2});
        col_bar{tmp.testnameind,2}=colorbar;
        set(col_bar{tmp.testnameind,2}, 'fontsize',param.m_grid_fontsize+5);

        tmp.titlename = strcat(num2str(tmp.checktime, '%02i'), 'd, egg #, ',tmp.testname_ssp, ',(', ...
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

        saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);

        disp(' ')
        disp([num2str(flagi), ' plot is created.'])
        disp(' ')
        disp([' File path is : ',tmp.tifname])
        disp(' ')
        close all;
    end
    
end



% % % if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
% % %     run(param_script);
% % %     ind=1;
% % %     clear comb_egg_mask
% % %     for yearij = 1:length(inputyear)
% % %         tempyear = inputyear(yearij);
% % %         for monthij = 1:length(inputmonth)
% % %             tempmonth = inputmonth(monthij);
% % %             ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
% % %             disp([num2str(yearij), 'y_',num2str(monthij),'m'])
% % %             all_checktime=ncread(ncname, 'checktime');
% % %             ind_checktime=find(all_checktime==temp_checktime);
% % %             egg_mask=squeeze(ncread(ncname, 'mask_par', [1 1 1 ind_checktime], [inf inf inf 1]));
% % % %                         egg_mask=ncread(ncname, 'mask_15day');
% % %             lastday_m=size(egg_mask,3);
% % %             if (exist('comb_egg_mask')==0)
% % %                 comb_egg_mask=egg_mask;
% % %             else
% % %                 comb_egg_mask(:,:,end+1:end+lastday_m)=egg_mask;
% % %             end
% % %         end
% % %     end
% % %     lon_rho = ncread(ncname, 'lon_rho');
% % %     lat_rho = ncread(ncname, 'lat_rho');
% % % %                 pcolor(sum(comb_egg_mask,3)'/size(comb_egg_mask,3)); shading flat; colorbar
% % %     mean_data = sum(comb_egg_mask,3);
% % %     mean_data(mean_data==0)=NaN;
% % %     mask_model2 = double(inpolygon(lon_rho,lat_rho,refpolygon(:,1),refpolygon(:,2)));
% % %     mask_model2(mask_model2==0)=NaN;
% % %     mean_data = mean_data .* mask_model2;
% % %     testnameind=1;
% % %     m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
% % % %                 sb{testnameind}=subplot(8,4,[1+(testnameind-1)*8 2+(testnameind-1)*8 5+(testnameind-1)*8 6+(testnameind-1)*8]);  % Auto-fitted to the figure.
% % % %                 pos_sb{testnameind}=get(sb{testnameind}, 'pos'); % Get the position.
% % % %                 pos_sb{testnameind} = pos_sb{testnameind} + correction_large_fig;
% % % %                 delete(sb{testnameind}); % Delete the subplot axes
% % %     ax{testnameind,1}=axes;
% % %     pos_ax{testnameind}=get(ax{testnameind,1}, 'pos');
% % % %                 set(ax{testnameind,1},'pos',pos_sb{testnameind});
% % %     temp_surf=ncread(ncname, 'temp_surf', [1 1 1], [inf inf 1]);
% % %     temp_surf(isnan(temp_surf))=50000;
% % %     temp_surf(temp_surf<50000)=NaN;
% % %     model_land=temp_surf;
% % %     model_land(temp_surf==50000)=1;
% % %     pc{testnameind,1}=m_pcolor(lon_rho',lat_rho', model_land','parent',ax{testnameind,1});
% % %     colormap(ax{testnameind,1},[0.8 0.8 0.8]);
% % %     shading(gca,m_pcolor_shading_method); 
% % %     m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
% % %            'parent', ax{testnameind,1});
% % %     col_bar{testnameind,1}=colorbar;
% % %     set(col_bar{testnameind,1}, 'TickLabels', []);
% % %     hold on
% % % 
% % %     m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
% % %     ax{testnameind,2}=axes;
% % %     set(ax{testnameind,2},'pos',pos_ax{testnameind});
% % %     pc{testnameind,2}=m_pcolor(lon_rho',lat_rho', mean_data','parent',ax{testnameind,2});
% % %     colormap(ax{testnameind,2},jet(10));
% % %     caxis([1, 100]);
% % %     shading(gca,m_pcolor_shading_method);   
% % %     m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
% % %            'backcolor', 'none','parent', ax{testnameind,2});
% % %     col_bar{testnameind,2}=colorbar;
% % %     set(col_bar{testnameind,2}, 'fontsize',m_grid_fontsize+5);
% % % 
% % % %                if testnameind==4
% % % %                     m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
% % % %                         'xticklabels', [120,  130, 140], 'xtick',[120,  130, 140], ...
% % % %                         'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
% % % %                 else
% % % %                     m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
% % % %                         'xticklabels', [], 'xtick',[120,  130, 140], ...
% % % %                         'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
% % % %                 end
% % % %                 txt{testnameind,2}=m_text(119, 50, test_text{testnameind}, 'FontSize', m_grid_fontsize+4); 
% % % %                 txt{testnameind,3}=m_text(119, 48, test_text2{testnameind}, 'FontSize', m_grid_fontsize+4); 
% % % 
% % %     titlename = strcat(num2str(temp_checktime, '%02i'), 'd, egg density, ',testname, ',(', ...
% % %         num2str(min(inputyear),'%04i'),'-', num2str(max(inputyear),'%04i'), ',',  ...
% % %         num2str(min(inputmonth),'%02i'),'-', num2str(max(inputmonth),'%02i'), ')'); 
% % % 
% % %     title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
% % % 
% % %     set(gcf, 'PaperUnits', 'points');
% % %     set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
% % %     set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
% % % 
% % %     saveas(gcf,jpgname,'tif'); RemoveWhiteSpace([], 'file', jpgname);
% % % 
% % %     disp(' ')
% % %     disp([fig_flags{1,1}, ' plot is created.'])
% % %     disp(' ')
% % %     disp([' File path is : ',jpgname])
% % %     disp(' ')
% % %     close all;
% % % end