clear all; clc; close all;

warning off;
testname='test23';
var='vert_temp';
linux=1; windows=0;
if (windows ==1)
    % % for windows
    dropboxpath='C:\Users\KYY\Dropbox'; %% SNU_desktop
    addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
    addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
    addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run']));
    rootdir='C:\Users\kyy\Desktop\conv\';
    filedir=[rootdir,testname,'\'];
    outfile=[filedir,'fig\'];
    run('C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\convection\fig_param_kyy_conv');
elseif (linux==1)
    % % for linux
    dropboxpath='/home01/kimyy/Dropbox'; %% ROMS server
%     dropboxpath='/home/kimyy/Dropbox'; %% DAMO server
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
    addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_20/run']));
    
    rootdir='/data/kimyy/todgist/conv/output/';
    filedir=[rootdir,testname,'/'];
    outfile=[filedir,'fig/'];
    if (exist(outfile,'dir')~=7)
        mkdir(outfile);
    end
    run([dropboxpath,'/source/matlab/Model/ROMS/Analysis/Figure/convection/fig_param_kyy_conv.m']);
end




inputminute=1:4;
% shadlev=[-1 -0.05];
shadlev=[-0.5 -0.025];
% shadlev=[-0.2 -0.01];
% shadlev=[-0.1 -0.005];
vid_fps=5;
vid_qty=100;
flag_record_mp4=0;
flag_record_gif=1;
if(flag_record_mp4),    % prepare video object
    vidName = [plot_dir,'',exp_name,'_',plot_type,'_',num2str(exp_num(kk))];
    vidObj = VideoWriter(vidName,'MPEG-4');
    vidObj.Quality = vid_qty;
    vidObj.FrameRate = vid_fps;
    open(vidObj);
end
if(flag_record_gif), vidName = strcat(outfile, testname, '_ani.gif'); end
conlev=10:10;

tempminute=inputminute(1);
filename = strcat(filedir,'ocean_his_', num2str(tempminute,'%04i'), '.nc');
Vstretching = ncread(filename,'Vstretching')';
Vtransform = ncread(filename,'Vtransform')';
theta_s = ncread(filename,'theta_s')';
theta_b = ncread(filename,'theta_b')';
s_rho = ncread(filename,'s_rho')';
N=length(s_rho);
hc = ncread(filename,'hc')';

jet_conv=jet;
close all;

nfr = 0;
for minuteij=1:length(inputminute)
    nfr = nfr + 1;
    tempminute = inputminute(minuteij);
    filename = strcat(filedir,'ocean_his_', num2str(tempminute,'%04i'), '.nc');
    
    grd_file=filename;
    if (exist('grd.lon_rho' , 'var') ~= 1)
        grd.lon_rho = ncread(grd_file,'x_rho')';
        grd.lat_rho = ncread(grd_file,'y_rho')';
        grd.mask_rho = ncread(grd_file,'mask_rho')';
        grd.lon_u = ncread(grd_file,'x_u')';
        grd.lat_v = ncread(grd_file,'y_v')';
        grd.h = ncread(grd_file,'h')';
        grd.mask_rho_nan = grd.mask_rho;
        land = find(grd.mask_rho_nan==0);
        grd.mask_rho_nan(land) = NaN;

        h = grd.h;
%         kgrid=0;
%         [sc_r,Cs_r]=stretching(Vstretching, theta_s, theta_b, hc, N, kgrid);
%         kgrid=1;
%         [sc_w,Cs_w]=stretching(Vstretching, theta_s, theta_b, hc, N, kgrid);

        % zeta  
          
        grd.N = N;
        lon_rho  = grd.lon_rho;
        lat_rho  = grd.lat_rho;
        lon_u  = grd.lon_u;
        lat_v  = grd.lat_v; 
        mask_rho = grd.mask_rho;
        h = grd.h;
        N = grd.N;
        pi=3.141592;
        R=6392;
        for j=2:511
            for i=2:511
                dA(j,i)= (lon_u(j,i)-lon_u(j,i-1))*(lat_v(j,i)-lat_v(j-1,i));
            end
        end
        dA(1,2:511)=dA(2,2:511);
        dA(512,2:511)=dA(511,2:511);
        dA(2:511,1)=dA(2:511,2);
        dA(2:511,512)=dA(2:511,511);
        dA(1,1)=dA(2,2); dA(1,512)=dA(2,511); dA(512,1)=dA(511,2); dA(512,512)=dA(511,511);
    end
    pt_eps = -1e-6;
    prho=1027.306;
    c_p=3.99e3;
    zeta = zeros(size(grd.h)); % default
      if (strcmp(grd_file,'E:\Data\Reanalysis\nwp_1_10_seo\avg_ens_10km_mean\roms_grid_final.nc')==1)
        grd.zeta = zeta;
      else
%         grd.zeta = ncread(grd_file,'zeta')'
        grd.zeta=zeros(size(h));
      end
%       grd.z_r=zlevs(Vtransform,Vstretching,h,grd.zeta,theta_s,theta_b,hc,N,'r');  
%     depth=grd.z_r;
    
    lon_min=1; lon_max=512; lat_min=256; lat_max=256;
    grd.z_r=zlevs(Vtransform,Vstretching,h(lat_min(1):lat_max(1), lon_min(1):lon_max(1)),grd.zeta(lat_min(1):lat_max(1), lon_min(1):lon_max(1)),theta_s,theta_b,hc,N,'r');  
%     depth=grd.z_r;
    data_info = ncinfo(filename, varname); 
    if (length(data_info.Dimensions)==4)
        data = ncread(filename,varname,[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 data_info.Size(3) 1]);  %% cut horizontal area [x,y,z] (wider than target area)
    else
        data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 data_info.Size(3)]);  %% cut horizontal area [x,y,z] (wider than target area)    end
    end
    data=squeeze(data);
    cut_data = permute(data, [3 2 1]);  %% permute [x y z] -> [z y x]
    cut_data=reshape(cut_data,[208,1,512]);
    
    cut_mask_rho = mask_rho(lat_min(1):lat_max(1), lon_min(1):lon_max(1));
%     cut_depth = depth(:,lat_min(1):lat_max(1), lon_min(1):lon_max(1));
    cut_depth = grd.z_r;


    x1=1; y1=1;
    x2=1; y2=512;
    lat1=lat_rho(x1(1),y1(1));  lon1=lon_rho(x1(1),y1(1));
    lat2=lat_rho(x2(1),y2(1));  lon2=lon_rho(x2(1),y2(1));

    lon_line = lon_rho(min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1)));  %% for 1/20^o horizontal resolution
    lat_line = (lon_line-lon1)/((lon2-lon1)/(lat2-lat1))+lat1; %% weight for latitude index (lat1 : 0.05/((lon2-lon1) * (lat2-lat1) : lat2 
    x=repmat(lon_line,grd.N,1);  %% copy lon_line (grd.N times) to make matrix 
    x_label='distance from western boundary (m)';
    titlename = [num2str(tempminute,'%04i'), 'min'];

    Temp(:,:) = squeeze(cut_data(:,min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1))));
    Yi(:,:)= squeeze(cut_depth(:,min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1))));

    data=Temp;
    
    close all;
    figure('position',[100,190,580,780],'PaperUnits','inches','PaperPosition',[0 0 5.8 7.8]); %%   %%figure window, figure file 
    set(gcf,'renderer', 'zbuffer');
%         figure('PaperUnits','inches','PaperPosition',[0 0 5.8 7.8]); %%   %%figure window, figure file 
    set(gca,'position',[0.2 0.25 0.65 0.5]);  %% figure
    axis tight;
    hold on
    pcolor(x,Yi,data)
    shading flat;
    caxis(shadlev)
    set(gca,'box','on','linewidth',1.5,'fontsize',17)
    xlabel(x_label,'color','k','FontSize',17,'fontweight','bold')
    ylabel('Depth(m)','color','k','FontSize',17,'fontweight','bold')
    title(titlename,'fontsize',17); % stlee
      hold on
%       [C,h2]=contour(x,Yi, data, conlev, m_contour_color, 'linewidth', m_contour_linewidth);
%     clabel(C,h2,'FontSize',m_contour_label_fontsize,'Color',m_contour_label_color, ...
%         'labelspacing',m_contour_labelspacing,'Rotation',m_contour_rotation,'fontweight',m_contour_fontweight);              
    ch = colorbar;
    colormap(jet);
    set(ch,'fontsize',colorbar_fontsize);
    title(ch,colorbar_title,'fontsize',colorbar_title_fontsize);    
    jpgname=strcat(outfile, testname, '_', num2str(tempminute,'%04i'), '.png'); %% ~_year_minute.jpg
%     set( gcf, 'position', [680,190,580,780] );
%     figure('position',[680,190,580,780],'PaperUnits','inches','PaperPosition',[0 0 5.8 7.8]); %%   %%figure window, figure file 
%     pause(2);  %% sleep 2 seconds for uniform figure printing   
    drawnow; %% prevent too short time between previous command and save command
    saveas(gcf,jpgname,'png');
    if(flag_record_mp4), writeVideo(vidObj,getframe(gcf)); end
    if(flag_record_gif)
        f = getframe(gcf);
        if( ~exist('map','var') ),  [im,map] = rgb2ind(f.cdata,256,'nodither');
%               else im(,,1,nfr) = rgb2ind(f.cdata,map,'nodither');
                else  im(:,:,1,minuteij) = rgb2ind(f.cdata,map,'nodither');
        end
    end
    if(flag_record_mp4), close(vidObj); end
    if(flag_record_gif), imwrite(im,map,vidName,'delaytime',0.2,'loopcount',inf);   end
    close all
end
% 
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
% %     set(get(bar,'title'),'string','Temp. (��)','FontSize',15)
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
% 
% 
% 
