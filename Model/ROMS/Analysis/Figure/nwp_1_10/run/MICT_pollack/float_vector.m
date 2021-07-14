close all; clear all; clc;
dropboxpath='C:\Users\KYY\Dropbox'; %% SNU_desktop, kyy_laptop
addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run']));
addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Roms_tools\Preprocessing_tools']));


% % set m_quiver parameter
m_quiver_title_fontsize = 15;
m_quiver_x_interval = 2;
m_quiver_y_interval = 2;
m_quiver_vector_color = 'c';
m_quiver_LineWidth = 0.5;
m_quiver_AutoScale = 'off';
% m_quiver_ref_vec_x_range = 390:395;
% m_quiver_ref_vec_y_range = 745:750;
m_quiver_ref_text_fontsize = 15;
m_quiver_ref_text_x_location = 127.2;
m_quiver_ref_text_y_location = 37.5;
m_quiver_ref_u_value = 0.2;
m_quiver_ref_v_value = m_quiver_ref_u_value/10000.0;
m_quiver_ref_text = [num2str(m_quiver_ref_u_value),' m/s'];
m_quiver_text_color = 'c';
m_quiver_vector_size = 6;


year= 1985;

filename = ['E:\Data\Reanalysis\nwp_1_10_seo\avg_ens_10km_ens01\out_',num2str(year,'%04i'),'\output.nc'];
lon = ncread(filename,'lon');
lat = ncread(filename,'lat');
gridname = ['E:\Data\Reanalysis\nwp_1_10_seo\avg_ens_10km_ens01\roms_grid02.nc'];
mask_rho = ncread(gridname,'mask_rho');
lon_rho = ncread(gridname,'lon_rho');
lat_rho = ncread(gridname,'lat_rho');





hold on;
% m_proj('mercator','lon',[127 130],'lat',[38 41]);
% m_grid('fontsize',20);
% m_gshhs_c('color','k')  
% m_gshhs_c('patch',[.8 .8 .8]);
% for i=1:90
% m_line(squeeze(lon(i,:)),squeeze(lat(i,:)),'marker', 'o', 'color','r', 'linewidth', 0.5, 'markersize', 2, 'markerfacecolor', 'r');
% end
% m_line(squeeze(lon(2,:)),squeeze(lat(2,:)),'marker', 'o', 'color','b', 'linewidth', 0.5, 'markersize', 2, 'markerfacecolor', 'r');
% m_line(squeeze(lon(3,:)),squeeze(lat(3,:)),'marker', 'o', 'color','g', 'linewidth', 0.5, 'markersize', 2, 'markerfacecolor', 'r');

if (exist(['D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\3rd year\figure\ens_10km_ens01\drifter_LTRANS\',num2str(year,'%04i')] , 'dir') ~= 7)
    mkdir(['D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\3rd year\figure\ens_10km_ens01\drifter_LTRANS\',num2str(year,'%04i')]);
end 
outfile = ['D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\3rd year\figure\ens_10km_ens01\drifter_LTRANS\',num2str(year,'%04i'),'\vector'];
% pcolor
% lonlat = [127 130 38 40];
lonlat = [127 140 35 43];
mask_rho_nan=mask_rho;
mask_rho_nan(mask_rho_nan==0)=NaN;
[indw, inde, inds, indn]=findind_Y(1/10,[lonlat(1) lonlat(2) lonlat(3) lonlat(4)],lon_rho',lat_rho');
% pcolor(lon_rho(indw:inde,inds:indn),lat_rho(indw:inde,inds:indn),mask_rho_nan(indw:inde,inds:indn));
% shading flat


hor_paper_size_x= lonlat(2)-lonlat(1);
hor_paper_size_y = lonlat(4)-lonlat(3);
halt = 1;
while(halt)
    if (hor_paper_size_x > 500 || hor_paper_size_y > 500)
        halt = 0;
    else
        hor_paper_size_x = hor_paper_size_x * 1.2; 
        hor_paper_size_y = hor_paper_size_y * 1.2;
    end
end
paper_position_hor = 0; % % distance from left
paper_position_ver = 0; % % distance from bottom
paper_position_width = hor_paper_size_x ;
paper_position_height = hor_paper_size_y ;
vert_paper_size_y = 400;

vid_fps=3;
vid_qty=100;
flag_record_mp4=1;
flag_record_gif=1;
plot_dir = ['D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\3rd year\figure\ens_10km_ens01\drifter_LTRANS\',num2str(year,'%04i'),'\'];
exp_name = 'pollack';
plot_type = 'vector';
k=90;
figure(100);
if(flag_record_mp4),    % prepare video object
    vidName = [plot_dir,'',exp_name,'_',plot_type,'_',num2str(year),'_', num2str(k), '.mp4'];
    vidObj = VideoWriter(vidName,'MPEG-4');
%         vidObj = VideoWriter(vidName,'wmv');
    vidObj.Quality = vid_qty;
    vidObj.FrameRate = vid_fps;
    open(vidObj);
end
if(flag_record_gif), vidName = [plot_dir,'',exp_name,'_',plot_type,'_',num2str(year),'_', num2str(k), '.gif']; end

jet96=jet(96);
for t=1:90
    figure(100);
    ylim([lonlat(3) lonlat(4)]);
    xlim([lonlat(1) lonlat(2)]);
    hold on;
    uvfilename=['E:\Data\Reanalysis\nwp_1_10_seo\avg_ens_10km_ens01\out_',num2str(year,'%04i'),'\avg_',num2str(year,'%04i'),'_',num2str(t,'%04i'),'.nc'];
    u = ncread(uvfilename,'u');
    v = ncread(uvfilename,'v');
    u_r = u2rho_2d(u(:,:,1)')'; %%[x y]
    v_r = v2rho_2d(v(:,:,1)')';
    u_rho=u_r(indw:inde,inds:indn);
    v_rho=v_r(indw:inde,inds:indn);
    cut_lon_rho=lon_rho(indw:inde,inds:indn);
    cut_lat_rho=lat_rho(indw:inde,inds:indn);
    if (exist('ref_vec_x_range' , 'var') ~= 1)
        ref_vec_x_ind = find(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location) == min(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location)));
        ref_vec_y_ind = find(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location) == min(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location)))+m_quiver_y_interval*2;
        ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2))+m_quiver_x_interval-1;
        ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2))+m_quiver_y_interval-1;
    end
    u_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_u_value;
    v_rho(ref_vec_x_range,ref_vec_y_range)=m_quiver_ref_v_value;
    
    uvplot=quiver(cut_lon_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
                        cut_lat_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end), ...
                        u_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end) * m_quiver_vector_size, ...
                        v_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end) * m_quiver_vector_size, ...
                        'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth);

    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
    set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
    jpgname=strcat(outfile, '_', num2str(t,'%04i'), '.jpg'); %% ~_year_month.jpg
    plot_google_map('MapType','hybrid', 'Alpha', 1, 'Scale', 2)
    titlename = strcat(num2str(t,'%04i'), 'day vectors');
    title(titlename,'fontsize',20);  %%title
    text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_text, 'FontSize', m_quiver_ref_text_fontsize, 'color', m_quiver_text_color); 
    if(flag_record_mp4), writeVideo(vidObj,getframe(gcf)); end
    if(flag_record_gif)
        f = getframe(gcf);
        if( ~exist('map','var') ),  [im,map] = rgb2ind(f.cdata,256,'nodither');
                else  im(:,:,1,t) = rgb2ind(f.cdata,map,'nodither');
        end
    end
    
    drawnow; %% prevent too short time between previous command and save command
    saveas(gcf,jpgname,'jpg');
    hold off;
    close figure 100;
    disp(t)
end
if(flag_record_mp4), close(vidObj); end
if(flag_record_gif), imwrite(im,map,vidName,'delaytime',0.33,'loopcount',inf);   end

% plot_google_map('MapType','satellite', 'Alpha', 1, 'Scale', 2)

hold off;
