% % This code based on MATLAB R2016b.

clc;close all;clear all;
warning off;

linux=0; windows=1;

if (windows ==1)
    % % for windows
    dropboxpath='C:\Users\KYY\Dropbox'; %% SNU_desktop, kyy_laptop
    addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
    addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
    addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run']));
elseif (linux==1)
    % % for linux
    dropboxpath='/home01/kimyy/Dropbox'; %% ROMS server
%     dropboxpath='/home/kimyy/Dropbox'; %% DAMO server
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
    addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_20/run']));
end



for year=2013:2016;

    filename = ['C:\Users\kyy\Desktop\ocean_flt_', num2str(year), '.nc'];
    calendar{1} = ' January'; calendar{2} = ' February'; calendar{3} = ' March'; calendar{4} = ' April'; calendar{5} = ' May'; calendar{6} = ' June';
    calendar{7} = ' July'; calendar{8} = ' August'; calendar{9} = ' September'; calendar{10} = ' October'; calendar{11} = ' November'; calendar{12} = ' December';

    lon = ncread(filename,'lon');
    lat = ncread(filename,'lat');
    lon(find(lon==0))=NaN;
    lon(find(lat==0))=NaN;
    lonlat = [115 164 15 52];

    % set( gcf, 'position', [0,0,810,760] );
    daylimit=[4 9 14];
    selectind(1) =14;
    for i =1:12
        selectind(i+1) = 14 +sum(eomday(year,1:i));
    end
    for selectmonth = 1:3
        for k= 1:3
            for dmonth=1:3
                plot_dir = 'C:\Users\kyy\Desktop\snapshot\';
                exp_name = 'sanchi';
    %             plot_type = 'drifter'
    %             vid_fps=5;
    %             vid_qty=100;
    %             flag_record_mp4=1;
    %             flag_record_gif=1;
    %             if(flag_record_mp4),    % prepare video object
                jpgName = [plot_dir,'',exp_name,'_',num2str(year),'_',num2str(selectmonth),char(calendar(selectmonth)),'_',num2str(daylimit(k)+1),'_',num2str(selectmonth+dmonth),char(calendar(selectmonth+dmonth)),'.jpg'];
    %                 vidObj = VideoWriter(vidName,'MPEG-4');
    %         %         vidObj = VideoWriter(vidName,'wmv');
    %                 vidObj.Quality = vid_qty;
    %                 vidObj.FrameRate = vid_fps;
    %                 open(vidObj);
    %             end
    %             if(flag_record_gif), vidName = [plot_dir,'',exp_name,'_',plot_type,'_',num2str(year),'_',char(calendar(selectmonth)),'_', num2str(daylimit(k)), '.gif']; end

%                 ind=1;
%                 for month=selectmonth:12
                indstart=selectind(selectmonth); 
                indend = indstart + daylimit(k);
%                     if month==selectmonth
%                         dstart=15;
%                         dend=31;
%                     else
%                         dstart=1;
%                         dend=eomday(year,month);
%                     end
%                     for day=dstart:dend 
%                     figure(year+100); set(gcf,'renderer', 'zbuffer'); 
                m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
                m_grid('fontsize',20);
                m_gshhs_c('color','k')  
                m_gshhs_c('patch',[.8 .8 .8]);
                hold on
                for i=indstart:indend
            %             m_line(lon(i,1:24:24*ind),lat(i,1:24:24*ind));
                    m_line(lon(i,indstart*24+1:24:selectind(selectmonth+dmonth)*24),lat(i,indstart*24+1:24:selectind(selectmonth+dmonth)*24)); 
                    m_line(lon(i,selectind(selectmonth+dmonth)*24),lat(i,selectind(selectmonth+dmonth)*24), ...
                            'marker', 'o', 'color','r', 'linewi', 2, 'markersize', 2, 'markerfacecolor', 'r');%%
    %                     m_line(lon(i,indstart*24+1:24:indend*24),lat(i,indstart*24+1:24:indend*24));
                end
                hold off;
                titlename = strcat('(from ', char(calendar(selectmonth)),') ', num2str(year),', ', char(calendar(selectmonth+dmonth)),' ', num2str(15), '. (', num2str(daylimit(k)+1), ' particles)');
                title(titlename,'fontsize',20);  %%title
                set( gcf, 'position', [0,0,1300,920] );
%                 if(flag_record_mp4), writeVideo(vidObj,getframe(gcf)); end
%                 if(flag_record_gif)
%                     f = getframe(gcf);
%                     if( ~exist('map','var') ),  [im,map] = rgb2ind(f.cdata,256,'nodither');
%                             else  im(:,:,1,ind) = rgb2ind(f.cdata,map,'nodither');
%                     end
%                 end
                saveas(gcf,jpgName,'jpg');
                close all;
                
%                     ind=ind+1;
%                     end
%                 end
            end
%             if(flag_record_mp4), close(vidObj); end
%             if(flag_record_gif), imwrite(im,map,vidName,'delaytime',0.2,'loopcount',inf);   end
        end
    end
end
% 
% clear all; clc; close all;
% 
% rname='D:\MEPL\project\EAST\downwelling\downwelling_fig_2.nc';
% x=ncread(rname,'XT_OCEAN');
% y=ncread(rname,'YT_OCEAN');
% z=ncread(rname,'ST_OCEAN');
% xv=ncread(rname, 'XV');
% yu=ncread(rname, 'YU');
% xw=ncread(rname, 'XW');
% yw=ncread(rname, 'YW');
% xtemp=ncread(rname, 'XTEMP');
% ytemp=ncread(rname, 'YTEMP');
% time=1:1440;
% depth=-z
% ym=5:10:495;
% % % total x,y : 500m , dx=dy=dz=10m(~500m depth)
% dt=1; nz=200; nx=50; ny=50;
% t=1:dt:1440;
% 
% ninit=1:21;
% n=21;
% i=1;
% plot_dir='D:\MEPL\project\EAST\downwelling\';
% exp_name='downwelling_2';
% plot_type='shade';
% kk=1;
% exp_num(1)=1;
% vid_fps=5;
% vid_qty=100;
% flag_record_mp4=1;
% flag_record_gif=1;
% if(flag_record_mp4),    % prepare video object
%     vidName = [plot_dir,'',exp_name,'_',plot_type,'_',num2str(exp_num(kk))];
%     vidObj = VideoWriter(vidName,'MPEG-4');
%     vidObj.Quality = vid_qty;
%     vidObj.FrameRate = vid_fps;
%     open(vidObj);
% end
% if(flag_record_gif), vidName = [plot_dir,'',exp_name,'_',plot_type,'_',num2str(exp_num(kk)),'.gif']; end
% % figure(exp_num(kk)+100); set(gcf,'renderer', 'zbuffer'); set(gcf,'position',fig_size);
% figure(exp_num(kk)+100); set(gcf,'renderer', 'zbuffer'); 
% set( gcf, 'position', [680,190,580,780] );
% set( gca, 'fontsize', 15 );
% nfr = 0;
% % for ii=1size(time,1)
% for ii=1:120
%     nfr = nfr + 1;
% %     contourf(x,zu,pt_xz(,,ii)',clev,'edgecolor','none'); axis equal; axis tight;
%     pcolor(ym,depth,squeeze(ytemp(:,1,:,ii))'); shading interp;
%     ccc=colorbar;
%     caxis([-0.5,0]);
% %     ccc.Label.String= 'Temp (^oC)';
%     colormap jet;
%     xlabel('x(m)','fontsize',15); ylabel('depth(m)','fontsize',15); set(gca,'fontsize',15);
% %         bar = colorbar('fontsize',15);
% %     set(get(bar,'title'),'string','Temp. (กษ)','FontSize',15)
%     %     xlim(xlimit);
% %     hcb=colorbar; colormap(cmap); caxis(cb_lim); set(hcb,'fontsize',15);
%     set(get(ccc,'title'),'string','Temp (^oC)','fontsize',15);
% %     title( ['Temperature(theta) at ',num2str(time(ii),'%8.2f'),' sec'] );
%      title( [num2str(time(ii),'%8.2f'),' min.'] );
% %     drawnow;
%     if(flag_record_mp4), writeVideo(vidObj,getframe(gcf)); end
%     if(flag_record_gif)
%         f = getframe(gcf);
%         if( ~exist('map','var') ),  [im,map] = rgb2ind(f.cdata,256,'nodither');
% %               else im(,,1,nfr) = rgb2ind(f.cdata,map,'nodither');
%                 else  im(:,:,1,ii) = rgb2ind(f.cdata,map,'nodither');
%         end
%     end
% end
% if(flag_record_mp4), close(vidObj); end
% if(flag_record_gif), imwrite(im,map,vidName,'delaytime',0.2,'loopcount',inf);   end