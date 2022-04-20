 if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
    run(param_script);
    ind=1;
    clear comb_egg_mask
    tempyear = inputyear(1);
    tempmonth = inputmonth(1);
    ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
    ocean_mask=ncread(ncname, 'zeta',[1, 1, 1], [inf, inf, 1]);
    ocean_mask(isfinite(ocean_mask))=1;
    lon_rho = ncread(ncname, 'lon_rho');
    lat_rho = ncread(ncname, 'lat_rho');
    SK_EEZ_mask = double(inpolygon(lon_rho,lat_rho,SK_EEZ_polygon(:,1),SK_EEZ_polygon(:,2)));
    NK_EEZ_mask = double(inpolygon(lon_rho,lat_rho,NK_EEZ_polygon(:,1),NK_EEZ_polygon(:,2)));
    RU_PR_EEZ_mask = double(inpolygon(lon_rho,lat_rho,RU_PR_EEZ_polygon(:,1),RU_PR_EEZ_polygon(:,2)));
    RU_SH_EEZ_mask = double(inpolygon(lon_rho,lat_rho,RU_SH_EEZ_polygon(:,1),RU_SH_EEZ_polygon(:,2)));
    JP_EEZ_mask = double(inpolygon(lon_rho,lat_rho,JP_EEZ_polygon(:,1),JP_EEZ_polygon(:,2)));

    ocean_SK_EEZ=ocean_mask.*SK_EEZ_mask;
    ocean_SK_EEZ(ocean_SK_EEZ==0)=NaN;
    ocean_NK_EEZ=ocean_mask.*NK_EEZ_mask;
    ocean_NK_EEZ(ocean_NK_EEZ==0)=NaN;
    ocean_RU_PR_EEZ=ocean_mask.*RU_PR_EEZ_mask;
    ocean_RU_PR_EEZ(ocean_RU_PR_EEZ==0)=NaN;
    ocean_RU_SH_EEZ=ocean_mask.*RU_SH_EEZ_mask;
    ocean_RU_SH_EEZ(ocean_RU_SH_EEZ==0)=NaN;
    ocean_JP_EEZ=ocean_mask.*JP_EEZ_mask;
    ocean_JP_EEZ(ocean_JP_EEZ==0)=NaN;

    lon_rho = ncread(ncname, 'lon_rho');
    lat_rho = ncread(ncname, 'lat_rho');
%                 pcolor(sum(comb_egg_mask,3)'/size(comb_egg_mask,3)); shading flat; colorbar
%                 mean_data = sum(comb_egg_mask,3)/size(comb_egg_mask,3);
%                 mean_data(mean_data==0)=NaN;
    testnameind=1;
    m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
%                 sb{testnameind}=subplot(8,4,[1+(testnameind-1)*8 2+(testnameind-1)*8 5+(testnameind-1)*8 6+(testnameind-1)*8]);  % Auto-fitted to the figure.
%                 pos_sb{testnameind}=get(sb{testnameind}, 'pos'); % Get the position.
%                 pos_sb{testnameind} = pos_sb{testnameind} + correction_large_fig;
%                 delete(sb{testnameind}); % Delete the subplot axes
    ax{testnameind,1}=axes;
    pos_ax{testnameind}=get(ax{testnameind,1}, 'pos');
%                 set(ax{testnameind,1},'pos',pos_sb{testnameind});
    temp_surf=ncread(ncname, 'temp_surf', [1 1 1], [inf inf 1]);
    temp_surf(isnan(temp_surf))=50000;
    temp_surf(temp_surf<50000)=NaN;
    model_land=temp_surf;
    model_land(temp_surf==50000)=1;
    pc{testnameind,1}=m_pcolor(lon_rho',lat_rho', model_land','parent',ax{testnameind,1});
    colormap(ax{testnameind,1},[0.8 0.8 0.8]);
    shading(gca,m_pcolor_shading_method); 
%                 m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
%                         'xticklabels', [], 'xtick',[120, 130, 140], ...
%                         'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'parent', ax{testnameind,1});  
    m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
           'parent', ax{testnameind,1});
%                 col_bar{testnameind,1}=colorbar;
%                 set(col_bar{testnameind,1}, 'TickLabels', []);
    hold on

    m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
    ax{testnameind,2}=axes;
    set(ax{testnameind,2},'pos',pos_ax{testnameind});
    pc{testnameind,2}=m_pcolor(lon_rho',lat_rho', ocean_SK_EEZ','parent',ax{testnameind,2});
    colormap(ax{testnameind,2},[0 0 0]); % black
    caxis([0, 1]);
    shading(gca,m_pcolor_shading_method);   
    m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
           'backcolor', 'none','parent', ax{testnameind,2});
%                 col_bar{testnameind,2}=colorbar;
%                 set(col_bar{testnameind,2}, 'fontsize',m_grid_fontsize+5);

    m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
    ax{testnameind,3}=axes;
    set(ax{testnameind,3},'pos',pos_ax{testnameind});
    pc{testnameind,3}=m_pcolor(lon_rho',lat_rho', ocean_NK_EEZ','parent',ax{testnameind,3});
    colormap(ax{testnameind,3},[1 0 0]); % red
    caxis([0, 1]);
    shading(gca,m_pcolor_shading_method);   
    m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
           'backcolor', 'none','parent', ax{testnameind,3});
%                 col_bar{testnameind,3}=colorbar;
%                 set(col_bar{testnameind,3}, 'fontsize',m_grid_fontsize+5);

    m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
    ax{testnameind,4}=axes;
    set(ax{testnameind,4},'pos',pos_ax{testnameind});
    pc{testnameind,4}=m_pcolor(lon_rho',lat_rho', ocean_RU_PR_EEZ','parent',ax{testnameind,4});
    colormap(ax{testnameind,4},[0 0 1]); % blue
    caxis([0, 1]);
    shading(gca,m_pcolor_shading_method);   
    m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
           'backcolor', 'none','parent', ax{testnameind,4});

    m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
    ax{testnameind,5}=axes;
    set(ax{testnameind,5},'pos',pos_ax{testnameind});
    pc{testnameind,5}=m_pcolor(lon_rho',lat_rho', ocean_RU_SH_EEZ','parent',ax{testnameind,5});
    colormap(ax{testnameind,5},[1 0 1]); % magenta
    caxis([0, 1]);
    shading(gca,m_pcolor_shading_method);   
    m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
           'backcolor', 'none','parent', ax{testnameind,5});

    m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
    ax{testnameind,6}=axes;
    set(ax{testnameind,6},'pos',pos_ax{testnameind});
    pc{testnameind,6}=m_pcolor(lon_rho',lat_rho', ocean_JP_EEZ','parent',ax{testnameind,6});
    colormap(ax{testnameind,6},[0 1 0]); % green
    caxis([0, 1]);
    shading(gca,m_pcolor_shading_method);   
    m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
           'backcolor', 'none','parent', ax{testnameind,6});



%                if testnameind==4
%                     m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
%                         'xticklabels', [120,  130, 140], 'xtick',[120,  130, 140], ...
%                         'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
%                 else
%                     m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
%                         'xticklabels', [], 'xtick',[120,  130, 140], ...
%                         'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
%                 end
%                 txt{testnameind,2}=m_text(119, 50, test_text{testnameind}, 'FontSize', m_grid_fontsize+4); 
%                 txt{testnameind,3}=m_text(119, 48, test_text2{testnameind}, 'FontSize', m_grid_fontsize+4); 

    titlename = strcat('15d, egg density, ',testname, ',(', ...
        num2str(min(inputyear),'%04i'),'-', num2str(max(inputyear),'%04i'), ',',  ...
        num2str(min(inputmonth),'%02i'),'-', num2str(max(inputmonth),'%02i'), ')'); 

    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

    saveas(gcf,jpgname,'tif');

    disp(' ')
    disp([fig_flags{1,1}, ' plot is created.'])
    disp(' ')
    disp([' File path is : ',jpgname])
    disp(' ')
    close all;
end