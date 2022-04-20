if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
    run(param_script);
    clear ref_u_rho comb_ref_u_rho
    for yearij = 1:length(refyear)
        tempyear = refyear(yearij);
        for monthij = 1:length(inputmonth)
            tempmonth = inputmonth(monthij);
            ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
            disp([num2str(yearij), 'y_',num2str(monthij),'m'])
            ref_u_rho=ncread(ncname, 'uwind');
            lastday_m=size(ref_u_rho,3);
            if (exist('comb_ref_u_rho')==0)
                comb_ref_u_rho=ref_u_rho;
            else
                comb_ref_u_rho(:,:,end+1:end+lastday_m)=ref_u_rho;
            end
        end
    end
    mean_ref_u_rho=mean(comb_ref_u_rho,3);

    clear ref_v_rho comb_ref_v_rho
    for yearij = 1:length(refyear)
        tempyear = refyear(yearij);
        for monthij = 1:length(inputmonth)
            tempmonth = inputmonth(monthij);
            ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
            disp([num2str(yearij), 'y_',num2str(monthij),'m'])
            ref_v_rho=ncread(ncname, 'vwind');
            lastday_m=size(ref_v_rho,3);
            if (exist('comb_ref_v_rho')==0)
                comb_ref_v_rho=ref_v_rho;
            else
                comb_ref_v_rho(:,:,end+1:end+lastday_m)=ref_v_rho;
            end
        end
    end
    mean_ref_v_rho=mean(comb_ref_v_rho,3);

    clear u_rho comb_u_rho
    for yearij = 1:length(inputyear)
        tempyear = inputyear(yearij);
        for monthij = 1:length(inputmonth)
            tempmonth = inputmonth(monthij);
            ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
            disp([num2str(yearij), 'y_',num2str(monthij),'m'])
            u_rho=ncread(ncname, 'uwind');
            lastday_m=size(u_rho,3);
            if (exist('comb_u_rho')==0)
                comb_u_rho=u_rho;
            else
                comb_u_rho(:,:,end+1:end+lastday_m)=u_rho;
            end
        end
    end
    mean_later_u_rho=mean(comb_u_rho,3);
    mean_u_rho=mean_later_u_rho-mean_ref_u_rho;

    clear v_rho comb_v_rho
    for yearij = 1:length(inputyear)
        tempyear = inputyear(yearij);
        for monthij = 1:length(inputmonth)
            tempmonth = inputmonth(monthij);
            ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
            disp([num2str(yearij), 'y_',num2str(monthij),'m'])
            v_rho=ncread(ncname, 'vwind');
            lastday_m=size(v_rho,3);
            if (exist('comb_v_rho')==0)
                comb_v_rho=v_rho;
            else
                comb_v_rho(:,:,end+1:end+lastday_m)=v_rho;
            end
        end
    end
    mean_later_v_rho=mean(comb_v_rho,3);
    mean_v_rho=mean_later_v_rho-mean_ref_v_rho;  

    lon_rho = ncread(ncname, 'lon_rho');
    lat_rho = ncread(ncname, 'lat_rho');
%                 if (exist('ref_vec_x_range' , 'var') ~= 1)
        ref_vec_x_ind = find(abs(lon_rho(:,1)-m_quiver_ref_text_x_location) == min(abs(lon_rho(:,1)-m_quiver_ref_text_x_location)));
        ref_vec_y_ind = find(abs(lat_rho(1,:)-m_quiver_ref_text_y_location) == min(abs(lat_rho(1,:)-m_quiver_ref_text_y_location)))+m_quiver_y_interval*2;
        ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval;
        ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2))+m_quiver_y_interval;
%                 end
    mean_u_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_uwind_value;
    mean_v_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_vwind_value;

%                 mean_data(mean_data==0)=NaN;
    testnameind=1;
    m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
    ax{testnameind,1}=axes;

    temp_surf=ncread(ncname, 'temp_surf', [1 1 1], [inf inf 1]);
    temp_surf(isnan(temp_surf))=50000;
    temp_surf(temp_surf<50000)=NaN;
    model_land=temp_surf;
    model_land(temp_surf==50000)=1;
    pc{testnameind,1}=m_pcolor(lon_rho',lat_rho', model_land','parent',ax{testnameind,1});
    colormap(ax{testnameind,1},[0.8 0.8 0.8]);
    shading(gca,m_pcolor_shading_method); 
    m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
           'parent', ax{testnameind,1});
   pos_ax{testnameind}=get(ax{testnameind,1}, 'pos');

    hold on

    m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
    ax{testnameind,2}=axes;
    pc{testnameind,2}=m_quiver(lon_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                    lat_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                    mean_u_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_wind_vector_size*1.5, ...
                    mean_v_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_wind_vector_size*1.5, ...
                    'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth);

    m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_wind_text, 'FontSize', m_quiver_ref_text_fontsize); 

    m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
           'backcolor', 'none','parent', ax{testnameind,2});

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

    titlename = strcat('wind-diff, ',testname, ',(80s-', ...
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