    if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
        run(param_script);
        ind=1;
        clear comb_data
        for yearij = 1:length(inputyear)
            tempyear = inputyear(yearij);
            for monthij = 1:length(inputmonth)
                tempmonth = inputmonth(monthij);
                ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
                disp([num2str(yearij), 'y_',num2str(monthij),'m'])

                data=ncread(ncname, 'temp_surf');
                lastday_m=size(data,3);
                if (exist('comb_data')==0)
                    comb_data=data;
                else
                    comb_data(:,:,end+1:end+lastday_m)=data;
                end
            end
        end
        lon_rho = ncread(ncname, 'lon_rho');
        lat_rho = ncread(ncname, 'lat_rho');
    %                 pcolor(sum(comb_data,3)'/size(comb_data,3)); shading flat; colorbar
        mean_data = sum(comb_data,3)/size(comb_data,3);
        mean_data(mean_data==0)=NaN;
        testnameind=1;
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);

        ax{testnameind,1}=axes;
        pos_ax{testnameind}=get(ax{testnameind,1}, 'pos');
        temp_surf=ncread(ncname, 'temp_surf', [1 1 1], [inf inf 1]);
        temp_surf(isnan(temp_surf))=50000;
        temp_surf(temp_surf<50000)=NaN;
        model_land=temp_surf;
        model_land(temp_surf==50000)=1;
%         avhrrname= [filedir,testname,regionname,'_rms_clim_',num2str(1983,'%04i'),'_',num2str(2018,'%02i'),'.nc'];
%         model_land=ncread(avhrrname, 'model_land',[1 1], [inf inf]);
        pc{testnameind,1}=m_pcolor(lon_rho',lat_rho', model_land','parent',ax{testnameind,1});
        colormap(ax{testnameind,1},[0.8 0.8 0.8]);
        shading(gca,m_pcolor_shading_method); 

        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
               'xticklabels', [], 'yticklabels', [], 'parent', ax{testnameind,1});
        col_bar{testnameind,1}=colorbar;
    %                 set(col_bar{testnameind,1}, 'fontsize',m_grid_fontsize);
        set(col_bar{testnameind,1}, 'TickLabels', []);

        hold on
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        ax{testnameind,2}=axes;
        set(ax{testnameind,2},'pos',pos_ax{testnameind});
        
        mask_model = double(inpolygon(lon_rho,lat_rho,refpolygon(:,1),refpolygon(:,2)));
        mask_model(mask_model==0)=NaN;
        
        pc{testnameind,2}=m_pcolor(lon_rho',lat_rho', mean_data'.*mask_model','parent',ax{testnameind,2});

        colormap(ax{testnameind,2},jet_mod);
        caxis([0, 15]);
        shading(gca,m_pcolor_shading_method);   

        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
               'xticklabels', [], 'yticklabels', [], 'backcolor', 'none','parent', ax{testnameind,2});
        col_bar{testnameind,2}=colorbar;
    %                 set(col_bar{testnameind,2}, 'TickLabels', []);
        set(col_bar{testnameind,2}, 'fontsize',m_grid_fontsize+5);


        ax{testnameind,3}=axes;
        set(ax{testnameind,3},'pos',pos_ax{testnameind});
        [con_C,con_h]=m_contour(lon_rho',lat_rho', mean_data'.*mask_model', [2, 5, 10], m_contour_color, 'linewidth', m_contour_linewidth, 'parent',ax{testnameind,3});
        clabel(con_C,con_h,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
            'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);
        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
               'backcolor', 'none','parent', ax{testnameind,3});
        col_bar{testnameind,3}=colorbar;
        caxis([0, 15]);
        colormap(ax{testnameind,3},jet_mod);
        set(col_bar{testnameind,3}, 'fontsize',m_grid_fontsize+5);
    %                 set(gca, 'fontsize', m_grid_fontsize)

        titlename = strcat('SST, ',testname, ',(', ...
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