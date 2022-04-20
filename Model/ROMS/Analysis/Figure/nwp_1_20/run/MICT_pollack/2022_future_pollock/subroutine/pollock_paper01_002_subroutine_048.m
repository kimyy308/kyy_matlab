if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
    run(param_script);
    clear ref_data comb_ref_data
    for yearij = 1:length(refyear)
        tempyear = refyear(yearij);
        for monthij = 1:length(inputmonth)
            tempmonth = inputmonth(monthij);
            ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
            disp([num2str(yearij), 'y_',num2str(monthij),'m'])
            ref_data=ncread(ncname, 'mask_30day');
            lastday_m=size(ref_data,3);
            if (exist('comb_ref_data')==0)
                comb_ref_data=ref_data;
            else
                comb_ref_data(:,:,end+1:end+lastday_m)=ref_data;
            end
        end
    end
    mean_ref_data=mean(comb_ref_data,3);

%                 mean_ref_data(mean_ref_data==0)=NaN;

%                 pcolor(mean_ref_data'); shading flat; colorbar; caxis([0 0.7]); colormap(jet)
    clear data comb_data
    for yearij = 1:length(inputyear)
        tempyear = inputyear(yearij);
        for monthij = 1:length(inputmonth)
            tempmonth = inputmonth(monthij);
            ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
            disp([num2str(yearij), 'y_',num2str(monthij),'m'])
            data=ncread(ncname, 'mask_30day');
            lastday_m=size(data,3);
            if (exist('comb_data')==0)
                comb_data=data;
            else
                comb_data(:,:,end+1:end+lastday_m)=data;
            end
        end
    end
    mean_later_data=mean(comb_data,3);

%                 mean_later_data(mean_later_data==0)=NaN;

    mean_data=mean_later_data-mean_ref_data;
    mean_data(mean_data==0)=NaN;

    lon_rho = ncread(ncname, 'lon_rho');
    lat_rho = ncread(ncname, 'lat_rho');
    testnameind=1;
    m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
%                 sb{testnameind}=subplot(8,4,[1+(testnameind-1)*8 2+(testnameind-1)*8 5+(testnameind-1)*8 6+(testnameind-1)*8]);  % Auto-fitted to the figure.
%                 pos_sb{testnameind}=get(sb{testnameind}, 'pos'); % Get the position.
%                 pos_sb{testnameind} = pos_sb{testnameind} + correction_large_fig;
%                 delete(sb{testnameind}); % Delete the subplot axes
    ax{testnameind,1}=axes;

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
%                colorbar;
%                caxis([0,1]);
   pos_ax{testnameind}=get(ax{testnameind,1}, 'pos');
    col_bar{testnameind,1}=colorbar;
    set(col_bar{testnameind,1}, 'TickLabels', []);
    hold on

    m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
    ax{testnameind,2}=axes;
%                 set(ax{testnameind,2},'pos',pos_ax{testnameind});
    pc{testnameind,2}=m_pcolor(lon_rho',lat_rho', mean_data','parent',ax{testnameind,2});
    colormap(ax{testnameind,2},byrmap);
    caxis([-0.2,0.2]);
    shading(gca,m_pcolor_shading_method);   

    m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
           'backcolor', 'none','parent', ax{testnameind,2});
    col_bar{testnameind,2}=colorbar;
    set(col_bar{testnameind,2}, 'fontsize',m_grid_fontsize+5);

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

    titlename = strcat('egg-prob(30d), ',testname, ',(80s-', ...
        num2str(min(inputyear),'%04i'),'-', num2str(max(inputyear),'%04i'), ',',  ...
        num2str(min(inputmonth),'%02i'),'-', num2str(max(inputmonth),'%02i'), ')'); 

    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

    saveas(gcf,jpgname,'tif');

    disp(' ')
    disp([fig_name, ' plot is created.'])
    disp(' ')
    disp([' File path is : ',jpgname])
    disp(' ')
    close all;
end