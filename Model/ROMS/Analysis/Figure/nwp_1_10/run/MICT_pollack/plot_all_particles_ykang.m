close all; clear all; clc;

addpath(genpath('/home/auto/ext_hdi/ykang_pollack'));
testname = 'test06_DA';
totyear=1983:2019;
totmonth=1:3;
day=1;
totalday=30; 
col_lat = 40; % latitude 40N
calendarname=cell(1,12); calendarname{1} = 'January'; calendarname{2} = 'February'; calendarname{3} = 'March'; calendarname{4} = 'April'; calendarname{5} = 'May'; calendarname{6} = 'June';
calendarname{7} = 'July'; calendarname{8} = 'August'; calendarname{9} = 'September'; calendarname{10} = 'October'; calendarname{11} = 'November'; calendarname{12} = 'December';
folder_name = {'1_Jan','2_Feb','3_Mar'};

for nyear=1:length(totyear)
    for nmonth =1:length(totmonth)
        year=totyear(nyear);
        month=totmonth(nmonth);

        workdir=['/home/auto/ext_hdi/LTRANSv2b_auto/output/',testname,'/',num2str(year,'%04i'),'/',...
            num2str(totalday,'%04i'),'d_',num2str(year,'%04i'),num2str(month,'%02i'),num2str(day,'%02i'),'/'];
        filename = [workdir,'output_',num2str(totalday,'%04i'),'d_',num2str(year,'%04i'), ...
            num2str(month,'%02i'),num2str(day,'%02i'),'.nc'];
        
        lon = ncread(filename,'lon');
        lat = ncread(filename,'lat');
        gridname = ['/home/auto/ext_hdi/nwp_1_10/reanalysis/',  ...
            num2str(year,'%04i'),'/ocean_avg_0001.nc'];
        
        mask_rho = ncread(gridname,'mask_rho');
        lon_rho = ncread(gridname,'lon_rho');
        lat_rho = ncread(gridname,'lat_rho');
        
        figdir = ['/home/auto/ext_hdi/ykang_pollack/fig/all_plot/',num2str(year,'%04i'),...
            '/',folder_name{nmonth},'/'];
        if (exist([figdir] , 'dir') ~= 7)
            mkdir([figdir]);
        end
        outfile = [figdir,'drifter_',num2str(year,'%04i'),'_',...
            num2str(month,'%02i'),'_',num2str(day,'%02i'),'_',num2str(totalday,'%04i')];
        
        lonlat = [127 133 35.5 44.5];
        [size_p, size_t] = size(lat);

        
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
        
        plot_dir = [figdir,num2str(year,'%04i')];
        
        exp_name = 'pollack';
        plot_type = 'passive'
        k=90;
        %
        if(flag_record_mp4),    % prepare video object
            vidName = [plot_dir,'',exp_name,'_',plot_type,'_',num2str(year,'%04i'),'_',num2str(month,'%02i'),'_',num2str(day,'%02i'),'_',num2str(totalday,'%04i'), '.mp4'];
            vidObj = VideoWriter(vidName,'MPEG-4');
            %         vidObj = VideoWriter(vidName,'wmv');
            vidObj.Quality = vid_qty;
            vidObj.FrameRate = vid_fps;
            open(vidObj);
        end
        if(flag_record_gif), vidName = [plot_dir,'',exp_name,'_',plot_type,'_',num2str(year,'%04i'),'_',num2str(month,'%02i'),'_',num2str(day,'%02i'),'_',num2str(totalday,'%04i'), '.gif']; end
        
        
        for t = 1:totalday
            m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            m_grid('xtick',[127:2:133],'ytick',[36:1:44]);
            m_gshhs_c('color','k')
            m_gshhs_c('patch',[.8 .8 .8]);
            hold on;
            
            if t < 5
                for i=1:size_p
                    if lat(i,1) >= col_lat
                        col = 'r';
                    else col = 'b';
                    end
                    m_line(squeeze(lon(i,2*24+1:(2+t)*24)),squeeze(lat(i,2*24+1:(2+t)*24)),...
                        'marker', 'o', 'color',col, 'linewidth', 1, 'markersize', 1, 'markerfacecolor',col);
                end
            else
                for i=1:size_p
                    if lat(i,1) >= col_lat
                        col = 'r';
                    else col = 'b';
                    end
                    m_line(squeeze(lon(i,(t+2-4)*24+1:(t+2)*24)),squeeze(lat(i,(t+2-4)*24+1:(t+2)*24)),...
                        'marker', 'o', 'color',col, 'linewidth', 1, 'markersize', 1, 'markerfacecolor', 'k');
                end
            end
            
            m_text(127.2,44,['initial particles : ',num2str(size_p)],'fontweig','bold')
            %m_text(127.2,43.5,['surviving particles : ',num2str(survive_p)],'fontweig','bold')
            
                xlabel('longitude');
                ylabel('latitude');
                set(gca,'fontsize',15)
                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height])
                jpgname=strcat(outfile, '_', num2str(t,'%04i'), '.jpg'); %% ~_year_month.jpg
                %             pause
                %            plot_google_map('MapType','hybrid', 'Alpha', 1, 'Scale', 2)
                %             ylim([lonlat(3) lonlat(4)]);
                %             xlim([lonlat(1) lonlat(2)]);
                set(gcf,'PaperPositionMode','auto');
                
                %             titlename = strcat(num2str(year,'%04i'),'y, ',num2str(month,'%02i'),'m, ',num2str(t,'%04i'), 'day particles');
                titlename = strcat(num2str(year,'%04i'),'y, ',calendarname{month}(1:3),', ',num2str(t,'%04i'), 'day');
                title(titlename,'fontsize',17);  %%title
                
                if(flag_record_mp4), writeVideo(vidObj,getframe(gcf)); end
                if(flag_record_gif)
                    f = getframe(gcf);
                    if( ~exist('map','var') ),  [im,map] = rgb2ind(f.cdata,256,'nodither');
                    else  im(:,:,1,t) = rgb2ind(f.cdata,map,'nodither');
                    end
                end
                
                drawnow; %% prevent too short time between previous command and save command
                pause(2); %% cocolink server is really slow
                saveas(gcf,jpgname,'jpg');
                hold off;
                %             pause
                %             close figure 100;
                close all;
                disp(t)
            end
            if(flag_record_mp4), close(vidObj); end
            if(flag_record_gif), imwrite(im,map,vidName,'delaytime',0.33,'loopcount',inf);   end
            
            % plot_google_map('MapType','satellite', 'Alpha', 1, 'Scale', 2)
            
            hold off;
        end
    end
    
