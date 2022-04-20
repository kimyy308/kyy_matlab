   if (exist(figname , 'file') ~= 2 || fig_flag==2)
    run(param_script);
    ind=1;
    clear comb_data
    for yearij = 1:length(inputyear)
        tempyear = inputyear(yearij);
        for monthij = 1:length(inputmonth)
            tempmonth = inputmonth(monthij);
            ncname = [filedir,testname,regionname,'_sq_er_clim_',num2str(1983,'%04i'),'_',num2str(2018,'%02i'),'.nc'];
% model_land=ncread(ncname, 'avhrr_land');
            disp([num2str(yearij), 'y_',num2str(monthij),'m'])
            sst_ind=(inputyear(yearij)-1983)*12+inputmonth(monthij);
            data=ncread(ncname, 'sq_er',[1 1 sst_ind], [inf inf 1]);
            sst=ncread(ncname, 'sst',[1 1 sst_ind], [inf inf 1]);
            avhrr_sst=ncread(ncname, 'avhrr_sst',[1 1 sst_ind], [inf inf 1]);
            avhrr_err=ncread(ncname, 'avhrr_err',[1 1 sst_ind], [inf inf 1]);
            lon_rho=ncread(ncname, 'lon');
            lat_rho=ncread(ncname, 'lat');
%                         lastday_m=size(data,3);
            if (exist('comb_data')==0)
                comb_data=data;
                comb_sst=sst;
                comb_avhrr_sst=avhrr_sst;
                comb_avhrr_err=avhrr_err;
            else
                comb_data(:,:,end+1)=data;
                comb_sst(:,:,end+1)=sst;
                comb_avhrr_sst(:,:,end+1)=avhrr_sst;
                comb_avhrr_err(:,:,end+1)=avhrr_err;
            end
        end
    end
[lat_2d, lon_2d]=meshgrid(lat_rho,lon_rho);    
    
mean_sst=mean(comb_sst,3);
mean_avhrr_sst=mean(comb_avhrr_sst,3);
mean_avhrr_err=mean(comb_avhrr_err,3);
coast_mean_avhrr_err=mean_avhrr_err;
coast_mean_avhrr_err(mean_sst>2)=NaN;
coast_mean_avhrr_err(lon_2d>128.5)=NaN;
[mean_coast_mean_avhrr_err, error_status] = Func_0011_get_area_weighted_mean(coast_mean_avhrr_err, lon_rho, lat_rho);
[mean_mean_avhrr_err, error_status] = Func_0011_get_area_weighted_mean(mean_avhrr_err, lon_rho, lat_rho);
sst_pt_corr= corr2(mean_sst(logical(~isnan(mean_sst).*~isnan(mean_avhrr_sst))), ...
    mean_avhrr_sst(logical(~isnan(mean_sst).*~isnan(mean_avhrr_sst))));

    

%                 lon_rho = ncread(ncname, 'lon_rho');
%                 lat_rho = ncread(ncname, 'lat_rho');
    rms_data = sqrt(mean(comb_data,3));
    rms_data(rms_data==0)=NaN;
    testnameind=1;
    m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);

    ax{testnameind,1}=axes;
    pos_ax{testnameind}=get(ax{testnameind,1}, 'pos');
%                 temp_surf=ncread(ncname, 'temp_surf', [1 1 1], [inf inf 1]);
%                 temp_surf(isnan(temp_surf))=50000;
%                 temp_surf(temp_surf<50000)=NaN;

    model_land=ncread(ncname, 'avhrr_land',[1 1], [inf inf]);
    %     pcolor(model_land'); shading flat; colorbar
%     model_land=data;
%     model_land(isnan(model_land))=50000;
%     model_land(model_land<50000)=NaN;
%     model_land(model_land==50000)=1;
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

    pc{testnameind,2}=m_pcolor(lon_rho',lat_rho', rms_data','parent',ax{testnameind,2});

    colormap(ax{testnameind,2},jet_mod);
    caxis([0, 3]);
    shading(gca,m_pcolor_shading_method);   

    m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
           'xticklabels', [], 'yticklabels', [], 'backcolor', 'none','parent', ax{testnameind,2});
    col_bar{testnameind,2}=colorbar;
%                 set(col_bar{testnameind,2}, 'TickLabels', []);
    set(col_bar{testnameind,2}, 'fontsize',m_grid_fontsize+5);


    ax{testnameind,3}=axes;
    set(ax{testnameind,3},'pos',pos_ax{testnameind});
%     [con_C,con_h]=m_contour(lon_rho',lat_rho', rms_data', [2, 5, 10], m_contour_color, 'linewidth', m_contour_linewidth, 'parent',ax{testnameind,3});
%     clabel(con_C,con_h,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
%         'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);
    m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
           'backcolor', 'none','parent', ax{testnameind,3});
    col_bar{testnameind,3}=colorbar;
    caxis([0, 3]);
    colormap(ax{testnameind,3},jet_mod);
    set(col_bar{testnameind,3}, 'fontsize',m_grid_fontsize+5);
%                 set(gca, 'fontsize', m_grid_fontsize)
    

    [mean_rms, error_status] = Func_0011_get_area_weighted_mean(rms_data, lon_rho, lat_rho);
    titlename = strcat('RMSE,(', ...
        num2str(min(inputyear),'%04i'),'-', num2str(max(inputyear),'%04i'), ',',  ...
        num2str(min(inputmonth),'%02i'),'-', num2str(max(inputmonth),'%02i'), ', M=', num2str(round(mean_rms,2)), ')'); 

    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

    saveas(gcf,figname,'tif');

    disp(' ')
    disp([fig_name, ' plot is created.'])
    disp(' ')
    disp([' File path is : ',figname])
    disp(' ')
    close all;
end