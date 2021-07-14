close all; clear all; clc;

system_name=computer;
if (strcmp(system_name,'PCWIN64'))
    % % for windows
    dropboxpath='C:\Users\KYY\Dropbox'; %% SNU_desktop, kyy_laptop
    addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
    addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
    addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run']));

elseif (strcmp(system_name,'GLNXA64'))
    dropboxpath='/home/kimyy/Dropbox'; %% DAMO
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
    addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_20/run']));
end

warning off;

workdir='/data1/kimyy/Model/LTRANS/LTRANSv2b/LTRANSv2b/output/seo_ens_mean/';

totyear=1980:2010;
totmonth=1:3;
day=1;
totalday=90;
calendarname=cell(1,12); calendarname{1} = 'January'; calendarname{2} = 'February'; calendarname{3} = 'March'; calendarname{4} = 'April'; calendarname{5} = 'May'; calendarname{6} = 'June';
calendarname{7} = 'July'; calendarname{8} = 'August'; calendarname{9} = 'September'; calendarname{10} = 'October'; calendarname{11} = 'November'; calendarname{12} = 'December';


for nyear=1:length(totyear)
    for nmonth =1:length(totmonth)
        year=totyear(nyear);
        month=totmonth(nmonth);
        if month==12 
            year=totyear(nyear)-1;
        end
        if (strcmp(system_name,'PCWIN64'))
            filename = ['E:\Data\Reanalysis\nwp_1_10_seo\avg_ens_10km_ens01\out_',num2str(year,'%04i'),'\output.nc'];
        elseif (strcmp(system_name,'GLNXA64'))
            filename = [workdir, ...
                num2str(year,'%04i'),'/output_',num2str(year,'%04i'),'_',num2str(month,'%02i'),'_',num2str(day,'%02i'),'_',num2str(totalday,'%04i'),'.nc'];
        end
        lon = ncread(filename,'lon');
        lat = ncread(filename,'lat');
        if (strcmp(system_name,'PCWIN64'))
            gridname = ['E:\Data\Reanalysis\nwp_1_10_seo\avg_ens_10km_ens01\roms_grid02.nc'];
        elseif (strcmp(system_name,'GLNXA64'))
            gridname = ['/data2/kimyy/Reanalysis/nwp_1_10_seo/ens_mean_daily/' ...
                num2str(year,'%04i'),'/avg_',num2str(year,'%04i'),'_ens_mean_0001.nc'];
        end
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

        if (strcmp(system_name,'PCWIN64'))
            if (exist(['D:\OneDrive - ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Ð±ï¿?\MEPL\project\MICT_pollack\3rd year\figure\ens_10km_ens01\drifter_LTRANS\',num2str(year,'%04i')] , 'dir') ~= 7)
                mkdir(['D:\OneDrive - ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Ð±ï¿?\MEPL\project\MICT_pollack\3rd year\figure\ens_10km_ens01\drifter_LTRANS\',num2str(year,'%04i')]);
            end 
            outfile = ['D:\OneDrive - ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Ð±ï¿?\MEPL\project\MICT_pollack\3rd year\figure\ens_10km_ens01\drifter_LTRANS\',num2str(year,'%04i'),'\drifter'];
        elseif (strcmp(system_name,'GLNXA64'))
            if (exist([workdir,num2str(year,'%04i'),'/figures'] , 'dir') ~= 7)
                mkdir([workdir,num2str(year,'%04i'),'/figures']);
            end 
            outfile = [workdir,num2str(year,'%04i'),'/figures','/drifter_',num2str(year,'%04i'),'_',num2str(month,'%02i'),'_',num2str(day,'%02i'),'_',num2str(totalday,'%04i')];
        end
        % pcolor
        % lonlat = [127 130 38 40];
        lonlat = [128 131 37 41];
        % mask_rho_nan=mask_rho;
        % mask_rho_nan(mask_rho_nan==0)=NaN;
        % [indw, inde, inds, indn]=findind_Y(1/10,[127 130 38 40],lon_rho',lat_rho');
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
        flag_record_mp4=0;
        flag_record_gif=0;
        if (strcmp(system_name,'PCWIN64'))
            plot_dir = ['D:\OneDrive - ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Ð±ï¿?\MEPL\project\MICT_pollack\3rd year\figure\ens_10km_ens01\drifter_LTRANS\',num2str(year,'%04i'),'\'];
        elseif (strcmp(system_name,'GLNXA64'))
            plot_dir = [workdir,num2str(year,'%04i'),'/figures','/'];
        end
        exp_name = 'pollack';
        plot_type = 'passive'
        k=90;
        figure(100);
        if(flag_record_mp4),    % prepare video object
            vidName = [plot_dir,'',exp_name,'_',plot_type,'_',num2str(year,'%04i'),'_',num2str(month,'%02i'),'_',num2str(day,'%02i'),'_',num2str(totalday,'%04i'), '.mp4'];
            vidObj = VideoWriter(vidName,'MPEG-4');
        %         vidObj = VideoWriter(vidName,'wmv');
            vidObj.Quality = vid_qty;
            vidObj.FrameRate = vid_fps;
            open(vidObj);
        end
        if(flag_record_gif), vidName = [plot_dir,'',exp_name,'_',plot_type,'_',num2str(year,'%04i'),'_',num2str(month,'%02i'),'_',num2str(day,'%02i'),'_',num2str(totalday,'%04i'), '.gif']; end

        jet96=jet(96);
        for t=1:totalday
            figure(100);
            ylim([lonlat(3) lonlat(4)]);
            xlim([lonlat(1) lonlat(2)]);
            hold on;
            if t < 5
                for i=1:96-32
                    line(squeeze(lon(i,1:t*24)),squeeze(lat(i,1:t*24)),'marker', 'o', 'color',jet96(i,:), 'linewidth', 1, 'markersize', 1, 'markerfacecolor', 'k');
                end
            else
                for i=1:96-32
                    line(squeeze(lon(i,(t-4)*24:t*24)),squeeze(lat(i,(t-4)*24:t*24)),'marker', 'o', 'color',jet96(i,:), 'linewidth', 1, 'markersize', 1, 'markerfacecolor', 'k');
                end
            end
            xlabel('longitude');
            ylabel('latitude');
            set(gcf, 'PaperUnits', 'points');
            set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
            set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
            jpgname=strcat(outfile, '_', num2str(t,'%04i'), '.jpg'); %% ~_year_month.jpg
%             pause
            plot_google_map('MapType','hybrid', 'Alpha', 1, 'Scale', 2)
%             ylim([lonlat(3) lonlat(4)]);
%             xlim([lonlat(1) lonlat(2)]);
            set(gcf,'PaperPositionMode','auto');
            
%             titlename = strcat(num2str(year,'%04i'),'y, ',num2str(month,'%02i'),'m, ',num2str(t,'%04i'), 'day particles');
            titlename = strcat(num2str(year,'%04i'),'y, ',calendarname{month}(1:3),', ',num2str(t,'%04i'), 'day particles');
            title(titlename,'fontsize',20);  %%title
            
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
%             pause
            close figure 100;
            disp(t)
        end
        if(flag_record_mp4), close(vidObj); end
        if(flag_record_gif), imwrite(im,map,vidName,'delaytime',0.33,'loopcount',inf);   end

        % plot_google_map('MapType','satellite', 'Alpha', 1, 'Scale', 2)

        hold off;
    end
end

