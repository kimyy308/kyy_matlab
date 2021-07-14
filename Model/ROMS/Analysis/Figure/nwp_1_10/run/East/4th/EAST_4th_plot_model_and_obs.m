% clear all;close all;
%==========================================================================
% %  Updated 19-Jul-2019 by Yong-Yub Kim
close all; clear all;
warning off;
dropboxpath='C:\Users\KYY\Dropbox';
addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));


testname='test06';
param_script='EAST_fig_param_kyy.m';

uvfig=0;
for fig_title = {'1804_SST','1804_SSS','1807_SST','1807_SSS'}
    disp(fig_title{1})
end
for fig_title = {'1804_SST','1804_SSS','1807_SST','1807_SSS'}

    outfile='D:\OneDrive - 서울대학교\MEPL\project\EAST\4th_year\Figure\fig';
    switch fig_title{1}(1:4)
            case '1804'
                load('1804_CTD_info.mat'); data_ctd = data18;
                tempyear=2018;
                tempmonth=4;
                title_str = {' 2018.04.21',' ~2018.04.29'};
                lat_lim = [31.0, 37.0];
                lon_lim = [123.0, 131.00];
                filedir = strcat(['E:\Data\Model\ROMS\nwp_1_10\',testname,'\DA\',num2str(tempyear,'%04i'),'\']);
                filename = strcat(filedir, ...
                    testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(04210429,'%08i'), '.nc');
                title_str = {' 2018.04.21',' ~2018.04.29'};
            case '1807'
                load('1807_CTD_info.mat');
                tempyear=2018;
                tempmonth=7;
                title_str = {' 2018.07.17',' ~2018.07.25'};;
                lat_lim = [34.5, 38.5];
                lon_lim = [127.0, 133.00];
                filedir = strcat(['E:\Data\Model\ROMS\nwp_1_10\',testname,'\DA\',num2str(tempyear,'%04i'),'\']);
                filename = strcat(filedir, ...
                    testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(07170725,'%08i'), '.nc');
                title_str = {' 2018.07.17',' ~2018.07.25'};
    end
    switch fig_title{1}(6:8)
            case 'SST'
%                 v=data_ctd(:,9);
                var='temp';
                shadlev=[9.5 19.0];
                if( strcmp(fig_title{1}(1:4),'1807') ), shadlev = [17.7, 27.2]; end
                conlev=[0:0];
                clear jetctd
                jet10000=parula(10000);
                jetfinind=sum(isfinite(data_ctd(:,9)));
                jetind=round((data_ctd(1:jetfinind,9)-shadlev(1))/(shadlev(2)-shadlev(1))*10000);
                jetind(jetind>=10000)=10000;
                jetind(jetind<1)=1;
                jetctd=jet10000(jetind,:);        
            case 'SSS'
%                 v=data_ctd(:,10);
                var='salt';
                shadlev=[32.45 34.65];
                conlev=[0:0];
                clear jetctd
                jet10000=parula(10000);
                jetfinind=sum(isfinite(data_ctd(:,10)));
                jetind=round((data_ctd(1:jetfinind,10)-shadlev(1))/(shadlev(2)-shadlev(1))*10000);
                jetind(jetind>=10000)=10000;
                jetind(jetind<1)=1;
                jetctd=jet10000(jetind,:);
    end
    section=[lon_lim, lat_lim -5 -5];
    lonlat = section(1:4);
    refdepth=(section(5)+section(6))/2.0   %% section(5) and section (6) should be equal, but just preparation for wrong case
    run(param_script);
    %==========================================================================

    % tempmonth=inputmonth(1);
    % filename = strcat(filedir, '\run\', num2str(tempyear,'%04i'), '\', ...
    %                 testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');

%     tempmonth=inputmonth(1);

    Vstretching = ncread(filename,'Vstretching')';
    Vtransform = ncread(filename,'Vtransform')';
    theta_s = ncread(filename,'theta_s')';
    theta_b = ncread(filename,'theta_b')';
    s_rho = ncread(filename,'s_rho')';
    N=length(s_rho);
    hc = ncread(filename,'hc')';

%     for monthij=1:length(inputmonth)
%         tempmonth = inputmonth(monthij);
%         tempmonth=04;
        filename = strcat(filedir, ...
                    testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');

%         if (exist('lon_rho' , 'var') ~= 1)
            gd = EAST_ROMS_read_grid(filename,Vtransform,Vstretching,theta_s,theta_b,hc,N);
            lon_rho  = gd.lon_rho;
            lat_rho  = gd.lat_rho; 
            lon_u  = gd.lon_u;
            lat_u  = gd.lat_u; 
            lon_v  = gd.lon_v;
            lat_v  = gd.lat_v; 
            mask_rho = gd.mask_rho;
            h = gd.h;
            N = gd.N;
            depth=gd.z_r;
            [Nlat, Nlon, Nz]=size(lon_rho);
            dl=(max(max(lon_rho))-min(min(lon_rho)))/Nlon;
            depth_u=0.5*(depth(:,:,1:Nlon-1)+depth(:,:,2:Nlon));
            depth_v=0.5*(depth(:,1:Nlat-1,:)+depth(:,2:Nlat,:));

            [lon_min, lon_max, lat_min, lat_max] = findind_Y(dl, section(1:4), lon_rho', lat_rho');
            if (uvfig==1)
                [lon_u_min, lon_u_max, lat_u_min, lat_u_max] = findind_Y(dl, section(1:4), lon_u', lat_u');
                [lon_v_min, lon_v_max, lat_v_min, lat_v_max] = findind_Y(dl, section(1:4), lon_v', lat_v');
            end
%         end

        data_info = ncinfo(filename, varname); 
        if (length(data_info.Dimensions)==4)
            data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 data_info.Size(3) 1]);  %% cut horizontal area [x,y,z] (wider than target area)
            if (uvfig==1)
                u = ncread(filename,'u',[lon_u_min(1) lat_u_min(1) 1 1], [lon_u_max(1)-lon_u_min(1)+1 lat_u_max(1)-lat_u_min(1)+1 data_info.Size(3) 1]);  %% cut horizontal area [x,y,z] (wider than target area)
                v = ncread(filename,'v',[lon_v_min(1) lat_v_min(1) 1 1], [lon_v_max(1)-lon_v_min(1)+1 lat_v_max(1)-lat_v_min(1)+1 data_info.Size(3) 1]);  %% cut horizontal area [x,y,z] (wider than target area)
            end
        else
            data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 data_info.Size(3)]);  %% cut horizontal area [x,y,z] (wider than target area)
            if (uvfig==1)
                u = ncread(filename,'u',[lon_u_min(1) lat_u_min(1) 1], [lon_u_max(1)-lon_u_min(1)+1 lat_u_max(1)-lat_u_min(1)+1 data_info.Size(3)]);  %% cut horizontal area [x,y,z] (wider than target area)
                v = ncread(filename,'v',[lon_v_min(1) lat_v_min(1) 1], [lon_v_max(1)-lon_v_min(1)+1 lat_v_max(1)-lat_v_min(1)+1 data_info.Size(3)]);  %% cut horizontal area [x,y,z] (wider than target area)
            end
        end
        data=squeeze(data);
        if (uvfig==1)
            u=squeeze(u);
            v=squeeze(v);
        end

        cut_data = permute(data, [3 2 1]);  %% permute [x y z] -> [z y x]
        if (uvfig==1)
            cut_u = permute(u, [3 2 1]);  %% permute [x y z] -> [z y x]
            cut_v = permute(v, [3 2 1]);  %% permute [x y z] -> [z y x]
        end
        cut_depth = depth(:,lat_min(1):lat_max(1), lon_min(1):lon_max(1));
        if (uvfig==1)
            cut_depth_u = depth(:,lat_u_min(1):lat_u_max(1), lon_u_min(1):lon_u_max(1));
            cut_depth_v = depth(:,lat_v_min(1):lat_v_max(1), lon_v_min(1):lon_v_max(1));
        end

        cut_lon_rho = lon_rho(lat_min(1):lat_max(1), lon_min(1):lon_max(1));
        cut_lat_rho = lat_rho(lat_min(1):lat_max(1), lon_min(1):lon_max(1));
        if (uvfig==1)
            cut_lon_u = lon_u(lat_u_min(1):lat_u_max(1), lon_u_min(1):lon_u_max(1));
            cut_lat_u = lat_u(lat_u_min(1):lat_u_max(1), lon_u_min(1):lon_u_max(1));
            cut_lon_v = lon_v(lat_v_min(1):lat_v_max(1), lon_v_min(1):lon_v_max(1));
            cut_lat_v = lat_v(lat_v_min(1):lat_v_max(1), lon_v_min(1):lon_v_max(1));
        end

        cut_mask_rho = mask_rho(lat_min(1):lat_max(1), lon_min(1):lon_max(1));

        cut_data2 = EAST_roms_vinterp(cut_data,cut_depth,refdepth);
        if (uvfig==1)
            cut_u2 = EAST_roms_vinterp(cut_u,cut_depth_u,refdepth);
            cut_v2 = EAST_roms_vinterp(cut_v,cut_depth_v,refdepth);
            cut_u3=griddata(double(cut_lon_u), double(cut_lat_u), double(cut_u2(:,:)), double(cut_lon_rho), double(cut_lat_rho));
            cut_v3=griddata(double(cut_lon_v), double(cut_lat_v), double(cut_v2(:,:)), double(cut_lon_rho), double(cut_lat_rho));
            u_rho=cut_u3';
            v_rho=cut_v3';
        end
        data=cut_data2;
        
        if (uvfig==1)
    %         if (exist('ref_vec_x_range' , 'var') ~= 1)
                ref_vec_x_ind = find(abs(cut_lon_rho(1,:)-m_quiver_ref_text_x_location) == min(abs(cut_lon_rho(1,:)-m_quiver_ref_text_x_location)));
                ref_vec_y_ind = find(abs(cut_lat_rho(:,1)-m_quiver_ref_text_y_location) == min(abs(cut_lat_rho(:,1)-m_quiver_ref_text_y_location)))+m_quiver_y_interval*2;
                ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval/2;
                ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2))+m_quiver_y_interval/2;
    %         end
            u_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_u_value;
            v_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_v_value;       
        %     u_rho(m_quiver_hor_ref_vec_y_range,m_quiver_hor_ref_vec_x_range)=m_quiver_ref_u_value;
        %     v_rho(m_quiver_hor_ref_vec_y_range,m_quiver_hor_ref_vec_x_range)=m_quiver_ref_v_value;
        end
        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
        set(gcf,'PaperPosition', [paper_position_hor, paper_position_ver, 500, vert_paper_size_y]) 

    % % %     plot
        m_proj(m_proj_name,'lon',[section(1) section(2)],'lat',[section(3) section(4)]);
        hold on;
        m_pcolor(cut_lon_rho,cut_lat_rho,data);
        shading(gca,m_pcolor_shading_method);
        m_gshhs_i('color',m_gshhs_line_color)  
        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
    
        if (uvfig==1)
            uvplot=m_quiver(cut_lon_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                                cut_lat_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                                u_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end) * m_quiver_vector_size, ...
                                v_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end) * m_quiver_vector_size, ...
                                'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth);
        end
        if (refdepth>0  && refdepth<Nz)
            titlename = strcat('T. layer ', num2str(refdepth),' (',char(calendarname(tempmonth)), ', ', num2str(tempyear),')');
        elseif (refdepth ==Nz)
            titlename = strcat('T. Surface ', ' (',char(calendarname(tempmonth)), ', ', num2str(tempyear),')');
        else
%             titlename = strcat('T. ', num2str(refdepth),'m (',char(calendarname(tempmonth)), ', ', num2str(tempyear),')');
            titlename = strcat(var,', ', num2str(refdepth),'m ',title_str);
        end
        title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
        if (uvfig==1)
            m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_text, 'FontSize', m_quiver_ref_text_fontsize); 
        end

        % contour
        [C,h2]=m_contour(cut_lon_rho,cut_lat_rho, data, conlev, m_contour_color, 'linewidth', m_contour_linewidth);
        clabel(C,h2,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
            'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);

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
        
        if (uvfig==1)
            jpgname=strcat(outfile, '_', testname, '_', fig_title{1}, '_uv.jpg'); %% ~_year_month.jpg
        else
            jpgname=strcat(outfile, '_', testname, '_', fig_title{1}, '.jpg'); %% ~_year_month.jpg
        end
    %     pause(2); %% prevent too short time between previous command and save command
        
        
        hold on
        for i= 1:jetfinind
            m_line(data_ctd(i,2), data_ctd(i,1), 'marker', 'o','markerfacecolor', jetctd(i,:), 'markeredgecolor', 'k');
        end
        hold off
        drawnow; %% prevent too short time between previous command and save command
        saveas(gcf,jpgname,'jpg');

        disp(' ')
        disp([num2str(tempyear), '_', num2str(tempmonth), '_', var, ' plot is created.'])
        disp(' ')
        disp([' File path is : ',jpgname])
        disp(' ')

        hold off
        close all;
%     end
end
