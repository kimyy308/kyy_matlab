% %  Updated    18-Apr-2022 by Yong-Yub Kim   % 
% %  Updated    03-May-2022 by Yong-Yub Kim   % structure
% %  Updated    16-May-2022 by Yong-Yub Kim   % add reading prob_ens???? test


%% historical plot
tmp.tifname=strcat(dirs.figdir, tmp.testname_his,'_',tmp.regionname, '_spawning_ground_probability_', ...
            num2str(min(RCM_info.years_his),'%04i'),'_',num2str(max(RCM_info.years_his),'%04i'), 'y_', ...
            RCM_info.season, '.tif'); 
        
tmp.testname_prefix=tmp.testname_his(1:4);        
if (exist(tmp.tifname , 'file') ~= 2 || flags.fig_switch(flagi)==2)
    run(tmp.param_script);
%%    initialization
    if(isfield(tmp, 'tlen_his')~=1)
        tmp.tlen_his=0;
        if strcmp(tmp.testname_prefix, 'prob')
            for yearij = 1:length(RCM_info.years_his)
                tmp.tempyear = RCM_info.years_his(yearij);
%                 for dayij = 1:length(RCM_info.days)
%                     tmp.tempday = RCM_info.days(dayij);
                    tmp.ncname = [dirs.savedir_his,tmp.testname_his,'_',tmp.regionname,'model_pollock_',num2str(tmp.tempyear,'%04i'),'.nc'];
                    tmp.tlen_his=tmp.tlen_his + length(ncread(tmp.ncname, 'time'));
%                 end
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
%             for monthij = 1:length(RCM_info.months)
%                 tmp.tempmonth = RCM_info.months(monthij);
                tmp.ncname = [dirs.savedir_his,tmp.testname_his,'_',tmp.regionname,'model_pollock_',num2str(tmp.tempyear,'%04i'),'.nc'];
                disp([num2str(yearij), 'y'])

                tmp.egg_mask=ncread(tmp.ncname, 'egg_mask');
                tmp.lastday_m=size(tmp.egg_mask,3);
                tmp.comb_egg_mask(:,:,1:tmp.lastday_m)=tmp.egg_mask;
                tmp.endij=tmp.lastday_m;
%             end
        end
    else
        for yearij = 1:length(RCM_info.years_his)
            tmp.tempyear = RCM_info.years_his(yearij);
            for monthij = 1:length(RCM_info.months)
                tmp.tempmonth = RCM_info.months(monthij);
                tmp.ncname = [dirs.savedir_his,tmp.testname_his,'_',tmp.regionname,'model_pollock_',num2str(tmp.tempyear,'%04i'),'_',num2str(tmp.tempmonth,'%02i'),'.nc'];
                disp([num2str(yearij), 'y_',num2str(monthij),'m'])

                tmp.egg_mask=ncread(tmp.ncname, 'egg_mask');
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

    tmp.mean_data = mean(tmp.comb_egg_mask,3);

    tmp.mean_data(tmp.mean_data==0)=NaN;
%     tmp.mean_data(tmp.mean_data<0.01)=NaN;
    
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
    colormap(ax{tmp.testnameind,2},parula);
    
    caxis([0, 1.0]);
    shading(gca,param.m_pcolor_shading_method);   

    m_grid('fontsize', param.m_grid_fontsize, 'tickdir', param.m_grid_tickdir_type, 'box', param.m_grid_box_type,  ...
           'backcolor', 'none','parent', ax{tmp.testnameind,2});
    col_bar{tmp.testnameind,2}=colorbar;
    set(col_bar{tmp.testnameind,2}, 'fontsize',param.m_grid_fontsize+5);

    tmp.titlename = strcat('sp ground, ',tmp.testname_his, ',(', ...
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

tmp.tifname=strcat(dirs.figdir, tmp.testname_ssp,'_',tmp.regionname, '_spawning_ground_probability_', ...
            num2str(min(RCM_info.years_ssp),'%04i'),'_',num2str(max(RCM_info.years_ssp),'%04i'), 'y_', ...
            RCM_info.season, '.tif'); 
        
if (exist(tmp.tifname , 'file') ~= 2 || flags.fig_switch(flagi)==2)
    run(tmp.param_script);
%     ind=1;
%     clear egg_mask comb_egg_mask
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

            tmp.egg_mask=ncread(tmp.ncname, 'egg_mask');
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

                tmp.egg_mask=ncread(tmp.ncname, 'egg_mask');
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
%                 pcolor(sum(comb_egg_mask,3)'/size(comb_egg_mask,3)); shading flat; colorbar
%                 mean_data = sum(comb_egg_mask,3)/size(comb_egg_mask,3);
    tmp.mean_data = mean(tmp.comb_egg_mask,3);

    tmp.mean_data(tmp.mean_data==0)=NaN;
%     tmp.mean_data(tmp.mean_data<0.01)=NaN;

    
%% land fig    
    RCM_grid.mask_model2 = double(inpolygon(RCM_grid.lon_rho,RCM_grid.lat_rho,tmp.refpolygon(:,1),tmp.refpolygon(:,2)));
    RCM_grid.mask_model2(RCM_grid.mask_model2==0)=NaN;
    tmp.mean_data = tmp.mean_data .* RCM_grid.mask_model2;
    
    [tmp.mean_egg_mask, tmp.error_status] = Func_0011_get_area_weighted_mean(tmp.mean_data, RCM_grid.lon_rho, RCM_grid.lat_rho);
    

    tmp.testnameind=1;
    m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);
%                 sb{tmp.testnameind}=subplot(8,4,[1+(tmp.testnameind-1)*8 2+(tmp.testnameind-1)*8 5+(tmp.testnameind-1)*8 6+(tmp.testnameind-1)*8]);  % Auto-fitted to the figure.
%                 pos_sb{tmp.testnameind}=get(sb{tmp.testnameind}, 'pos'); % Get the position.
%                 pos_sb{tmp.testnameind} = pos_sb{tmp.testnameind} + correction_large_fig;
%                 delete(sb{tmp.testnameind}); % Delete the subplot axes
    ax{tmp.testnameind,1}=axes;

%                 set(ax{tmp.testnameind,1},'pos',pos_sb{tmp.testnameind});
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
%                colorbar;
%                caxis([0,1]);
    pos_ax{tmp.testnameind}=get(ax{tmp.testnameind,1}, 'pos');
    col_bar{tmp.testnameind,1}=colorbar;
    set(col_bar{tmp.testnameind,1}, 'TickLabels', []);
    hold on

%% data fig
    m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);
    ax{tmp.testnameind,2}=axes;
    pc{tmp.testnameind,2}=m_pcolor(RCM_grid.lon_rho',RCM_grid.lat_rho', tmp.mean_data','parent',ax{tmp.testnameind,2});
%     colormap(ax{tmp.testnameind,2},jet);
    colormap(ax{tmp.testnameind,2},parula);
    caxis([0, 1.0]);
    shading(gca,param.m_pcolor_shading_method);   

    m_grid('fontsize', param.m_grid_fontsize, 'tickdir', param.m_grid_tickdir_type, 'box', param.m_grid_box_type,  ...
           'backcolor', 'none','parent', ax{tmp.testnameind,2});
    col_bar{tmp.testnameind,2}=colorbar;
    set(col_bar{tmp.testnameind,2}, 'fontsize',param.m_grid_fontsize+5);

%                if tmp.testnameind==4
%                     m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
%                         'xticklabels', [120,  130, 140], 'xtick',[120,  130, 140], ...
%                         'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{tmp.testnameind,2});  
%                 else
%                     m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
%                         'xticklabels', [], 'xtick',[120,  130, 140], ...
%                         'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{tmp.testnameind,2});  
%                 end
%                 txt{tmp.testnameind,2}=m_text(119, 50, test_text{tmp.testnameind}, 'FontSize', m_grid_fontsize+4); 
%                 txt{tmp.testnameind,3}=m_text(119, 48, test_text2{tmp.testnameind}, 'FontSize', m_grid_fontsize+4); 

    tmp.titlename = strcat('sp ground, ',tmp.testname_ssp, ',(', ...
        num2str(min(RCM_info.years_ssp),'%04i'),'-', num2str(max(RCM_info.years_ssp),'%04i'), ',',  ...
        RCM_info.season, ')'); 

    title(tmp.titlename,'fontsize',param.m_pcolor_title_fontsize, 'parent', ax{tmp.testnameind,1});  %%title for land figure
    title(tmp.titlename,'fontsize',param.m_pcolor_title_fontsize, 'parent', ax{tmp.testnameind,2});  %%title for data figure

    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperSize', [param.hor_paper_size_x, param.hor_paper_size_y]);
    set(gcf,'PaperPosition', [param.paper_position_hor param.paper_position_ver param.paper_position_width param.paper_position_height]) 
    
    [tmp.indw, tmp.inde, tmp.inds, tmp.indn]=Func_0012_findind_Y(1/10,[128, 130, 36, 39],RCM_grid.lon_rho,RCM_grid.lat_rho); % southern EKB
%     [indw, inde, inds, indn]=Func_0012_findind_Y(1/10,[127, 131, 37, 42],lon_rho,lat_rho);
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