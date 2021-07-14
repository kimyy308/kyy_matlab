% clear all;close all;
%==========================================================================
% %  Updated 19-Jul-2019 by Yong-Yub Kim
close all; clear all;
warning off;
dropboxpath='C:\Users\KYY\Dropbox';
addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));


testname='test06';
param_script='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\EAST\4th\EAST_fig_param_kyy.m';
section = [121, 128, 29.5, 35.5, -5, -5];
uvfig=0;
for fig_title = {'20190718', '20190719', '20190720'}
    disp(fig_title{1})
end
for fig_title = {'20190718', '20190719', '20190720'}

    outfile='D:\OneDrive - 서울대학교\MEPL\project\EAST\4th_year\Figure\fig';

    var='temp';
    shadlev=[17 24.0];
    conlev=[0:0];

    lonlat = section(1:4);
    run(param_script);
    %==========================================================================
    varname='analysed_sst';

        filedir=['E:\Data\Observation\OISST\daily\',fig_title{1}(1:8),'\'];
        filename = strcat(filedir, ...
                    fig_title{1}(1:8), '120000-NCEI-L4_GHRSST-SSTblend-AVHRR_OI-GLOB-v02.0-fv02.0', '.nc');
        avhrrfilename=filename;
        avhrrinfo=ncinfo(avhrrfilename);
        avhrr_lon = ncread(avhrrfilename,'lon');
        avhrr_lat = ncread(avhrrfilename,'lat');

        avhrr_lon_west = abs(avhrr_lon - (lonlat(1)));
        min_avhrr_lon_west=min(avhrr_lon_west);
        avhrr_lon_east = abs(avhrr_lon - (lonlat(2)));
        min_avhrr_lon_east=min(avhrr_lon_east);
        avhrr_lat_south = abs(avhrr_lat - (lonlat(3)));
        min_avhrr_lat_south=min(avhrr_lat_south);
        avhrr_lat_north = abs(avhrr_lat - (lonlat(4)));
        min_avhrr_lat_north=min(avhrr_lat_north);

        avhrr_lon_min = find(avhrr_lon_west == min_avhrr_lon_west);
        avhrr_lon_max = find(avhrr_lon_east == min_avhrr_lon_east);
        avhrr_lat_min = find(avhrr_lat_south == min_avhrr_lat_south);
        avhrr_lat_max = find(avhrr_lat_north == min_avhrr_lat_north);

    %         ncinfo('E:\Data\Observation\OISST\monthly\avhrr_monthly1983_11.nc');

        avhrr_lon = ncread(avhrrfilename,'lon', [avhrr_lon_min(1)], [avhrr_lon_max(1)-avhrr_lon_min(1)+1]);
        avhrr_lat = ncread(avhrrfilename,'lat', [avhrr_lat_min(1)], [avhrr_lat_max(1)-avhrr_lat_min(1)+1]);

%         [lon_min, lon_max, lat_min, lat_max] = findind_Y(dl, section(1:4), lon_rho', lat_rho');


        data_info = ncinfo(filename, varname); 

        data = ncread(filename,varname,[avhrr_lon_min(1) avhrr_lat_min(1) 1], [avhrr_lon_max(1)-avhrr_lon_min(1)+1 avhrr_lat_max(1)-avhrr_lat_min(1)+1 data_info.Size(3)]);  %% cut horizontal area [x,y,z] (wider than target area)
        data=squeeze(data);

        cut_data = permute(data, [3 2 1]);  %% permute [x y z] -> [z y x]
        data=cut_data;
        data=squeeze(data)-273.15;
        
        
         % % set m_quiver parameter
        m_quiver_vector_size = 1;
        m_quiver_ref_u_value = 0.2;
        m_quiver_title_fontsize = 15;
        m_quiver_x_interval = 3;
        m_quiver_y_interval = 3;
        m_quiver_vector_color = 'w';
        m_quiver_LineWidth = 0.5;
        m_quiver_AutoScale = 'off';
        m_quiver_ref_text_fontsize = 15;
        m_quiver_ref_text_x_location = 126.5;
        m_quiver_ref_text_y_location = 34.8;
        m_quiver_ref_v_value = m_quiver_ref_u_value/10000.0;
        m_quiver_ref_text = [num2str(m_quiver_ref_u_value),' m/s'];

        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
        set(gcf,'PaperPosition', [paper_position_hor, paper_position_ver, 500, vert_paper_size_y]) 

        % % %     plot
        m_proj(m_proj_name,'lon',[section(1) section(2)],'lat',[section(3) section(4)]);
        hold on;
        m_pcolor(avhrr_lon,avhrr_lat,data);
        shading(gca,m_pcolor_shading_method);
        m_gshhs_i('color',m_gshhs_line_color)  
        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        
        titlename = strcat( fig_title{1}(1:8),', ', var);
        title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

    %         % contour
    %         [C,h2]=m_contour(cut_lon_rho,cut_lat_rho, data, conlev, m_contour_color, 'linewidth', m_contour_linewidth);
    %         clabel(C,h2,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
    %             'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);

            % set colorbar 
        h = colorbar;
        colormap(colormap_style);
        set(h,'fontsize',colorbar_fontsize);
        title(h,colorbar_title,'fontsize',colorbar_title_fontsize);
        caxis(shadlev);

        % set grid
        m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

        jpgname=strcat(outfile, '_', 'avhrr', '_', fig_title{1}(1:8), '.jpg'); %% ~_year_month.jpg
    %     pause(2); %% prevent too short time between previous command and save command


%         hold on
%         for i= 1:jetfinind
%             m_line(data_ctd(i,2), data_ctd(i,1), 'marker', 'o','markerfacecolor', jetctd(i,:), 'markeredgecolor', 'k');
%         end
%         hold off
        drawnow; %% prevent too short time between previous command and save command
        saveas(gcf,jpgname,'jpg');

        disp(' ')
%         disp([num2str(tempyear), '_', num2str(tempmonth), '_', var, ' plot is created.'])
        disp(' ')
        disp([' File path is : ',jpgname])
        disp(' ')

        hold off
        close all;
%     end
%     end
end
