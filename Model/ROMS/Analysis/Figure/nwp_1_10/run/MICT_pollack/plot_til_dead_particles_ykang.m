close all; clear all; clc;

addpath(genpath('/home/auto/ext_hdi/ykang_pollack/'));
min_temp = 2; max_temp = 5;
max_depth = 500;
dead_time = 3*24;

%     addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
%     addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
%     addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
%     addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
%     addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_20/run']));

testname = 'test06_DA';
totyear = 1983 : 2019;
totmonth = 1:2;
day=1; dd = num2str(day,'%02i');
totalday=30; len_day = num2str(totalday,'%04i');
col_lat = 40; % latitude 40N
calendarname=cell(1,12); calendarname{1} = 'January'; calendarname{2} = 'February'; calendarname{3} = 'March'; calendarname{4} = 'April'; calendarname{5} = 'May'; calendarname{6} = 'June';
calendarname{7} = 'July'; calendarname{8} = 'August'; calendarname{9} = 'September'; calendarname{10} = 'October'; calendarname{11} = 'November'; calendarname{12} = 'December';

for nyear=1:length(totyear)
    for nmonth =1:length(totmonth)
        year=totyear(nyear); yy = num2str(year);
        month=totmonth(nmonth); mm = num2str(month, '%02i');
        
        if month==12
            year=totyear(nyear)-1;
        end
        
        
        workdir=['/home/auto/ext_hdi/LTRANSv2b_auto/output/',testname,'/',yy,'/',...
            len_day,'d_',yy,mm,dd,'/'];
        filename = [workdir,'output_',len_day,'d_',yy,mm,dd,'.nc'];
        
        tempname = [workdir,'LTRANS_egg_depth_',yy,'_',mm,'_',len_day,'.mat'];
        load(tempname);
        
        lon = ncread(filename,'lon');
        lat = ncread(filename,'lat');
        gridname = ['/home/auto/ext_hdi/nwp_1_10/reanalysis/',yy,'/ocean_avg_0001.nc'];
        
        mask_rho = ncread(gridname,'mask_rho');
        lon_rho = ncread(gridname,'lon_rho');
        lat_rho = ncread(gridname,'lat_rho');
        
        figdir = '/home/auto/ext_hdi/ykang_pollack/fig/plot_til_dead/';
        if (exist([figdir,yy] , 'dir') ~= 7)
            mkdir([figdir,yy]);
        end
        outfile = [figdir,yy,'/drifter_til_dead_',yy,'_',mm,'_',dd,'_',len_day];
        
        lonlat = [127 133 37.5 44];
        [size_p, size_t] = size(lat);
        count_dead = 0;
        
        % find dead eggs
        for i = 1 : size_p
            for j = 1: size_t
                if temp(j,i)<=min_temp | temp(j,i)>=max_temp | egg_depth(j,i)>=max_depth
                    idx_dead(j,i) = 1;
                else idx_dead(j,i) = 0;
                end
            end
        end
        
        k = 0;
        
        for pp = 1: size_p
            k = k + 1;
            
            for tt = 24*2 + 2 : size_t
                if tt+(dead_time-1) >= 769
                    sum_dead = sum(idx_dead(tt:end,pp));
                else
                    sum_dead = sum(idx_dead(tt:tt+(dead_time-1),pp));
                end
                
                if sum_dead == dead_time
                    index_dead(1,k) = tt+1;
                    break
                elseif tt == size_t & sum_dead ~= dead_time
                    index_dead(1,k) = NaN;
                end
            end
            
        end
   
    
    for i = 1: size_p
        if ~isnan(index_dead(i)) == 1
            lon(i,index_dead(i):end) = NaN;
            lat(i,index_dead(i):end) = NaN;
        end
    end
    
    %%%%%
    
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
    
    plot_dir = [figdir,yy];
    
    exp_name = 'pollack';
    plot_type = 'passive'
    k=90;
    %
    if(flag_record_mp4),    % prepare video object
        vidName = [plot_dir,'',exp_name,'_',plot_type,'_',yy,'_',mm,'_',dd,'_',len_day, '.mp4'];
        vidObj = VideoWriter(vidName,'MPEG-4');
        %         vidObj = VideoWriter(vidName,'wmv');
        vidObj.Quality = vid_qty;
        vidObj.FrameRate = vid_fps;
        open(vidObj);
    end
    if(flag_record_gif), vidName = [plot_dir,'',exp_name,'_',plot_type,'_',...
            yy,'_',mm,'_',dd,'_',len_day, '.gif']; end
    
    
    for t = 15  %1:totalday
        figure('position',[200 400 350 500])
        m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        m_grid('xtick',[127:2:133],'ytick',[36:1:44]);
        m_gshhs_c('color','k')
        m_gshhs_c('patch',[.8 .8 .8]);
        hold on;
        
        find_survive = ~isnan(lon(:,24*(t+2) + 1));
        survive_p = sum(find_survive);
        survive_south = length(find(lat(:,24*(t+2)+1)<=40));
        survive_north = length(find(lat(:,24*(t+2)+1)>40));
        total_particle = size_p;
        surviving_particle = survive_p;
        ratio = surviving_particle/total_particle * 100;
        %             save([workdir,yy,mm,dd,'_',len_day,'_survival_ratio'],'total_particle','surviving_particle','ratio')
        
        for i=1:size_p
            if lat(i,1) >= col_lat
                col = 'r';
            else col = 'b';
            end
            m_line(squeeze(lon(i,(2+t)*24 +1)),squeeze(lat(i,(2+t)*24+1)),...
                'marker', 'o', 'color',col, 'markersize',3,'linewid',3,'markerfacecolor',col);
        end
        
        
        m_text(127.2,43.5,['initial particles : ',num2str(size_p)],'fontweig','bold','fontsize',13)
        m_text(127.2,43,['South : ',num2str(survive_south)],'fontweig','bold','fontsize',13)
        m_text(127.2,42.5,['North : ',num2str(survive_north)],'fontweig','bold','fontsize',13)
        m_text(127.2,42,['Total : ',num2str(survive_north+survive_south)],'fontweig','bold','fontsize',13)
        
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
        titlename = [yy,'y, ',calendarname{month}(1:3),', ',num2str(t,'%04i'), 'day'];
        title(titlename,'fontsize',20);  %%title
        
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
        close all;
        disp(t)
    end
    if(flag_record_mp4), close(vidObj); end
    if(flag_record_gif), imwrite(im,map,vidName,'delaytime',0.33,'loopcount',inf);   end
    
    % plot_google_map('MapType','satellite', 'Alpha', 1, 'Scale', 2)
    
end
end

