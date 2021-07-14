%--------------------------------------------------------------------------
%
%
%                        Plotting ROMS OUT ( 365daily )
%
%                                             date : 2005. 11. 19
%                                             made by K.S. Moon
%
%                                             date : 2007. 10. 28
%                                             edited by Y.G.
%
%--------------------------------------------------------------------------
clear all;close all;
%==========================================================================
max_level= 30;
variation = 1;  %  1 = temperature , 2 = salinity ,3 = density

current   = 1;   skip_v = 3;     size_v = 2;    color_v = 'k';

plot_pcolorjw = 1;    temp_lim = [2 25];    salt_lim = [29.0 34.5]; den_lim = [25 27.5];

plot_contour  = 1;    color_c  ='-k' ;      temp_c  =[2:1:30];  salt_c  =[29:0.5:34.5];

% plot_geoshow  = 0;    color_g = [.7 .7 .7];'black';

switch_save   = 1;    out_type = 'png';

section = 1; % 0 => whole, 1=> yellowsea, 2=> eastsea

% grdfile       = 'd:\add2_ini_bry_grd\grid\roms_grid2_ADD_10_ep.nc';

yy = 2006;
start_mm=1;
end_mm=12;
time_step=1;

% file_dir=[num2str(yy),'\monthly_mean\'];
file_dir=['ocean_avg_mon_\'];

mm=start_mm;
%==========================================================================

gd = grd('NP');
lon_rho  = gd.lon_rho;
lat_rho  = gd.lat_rho; 
mask_rho = gd.mask_rho;
h=gd.h;
N = gd.N;
% depth=zlevs(h,gd.zeta,gd.theta_s,gd.theta_b,gd.hc,N,'r',1);
depth=zlevs(h,gd.zeta,gd.theta_s,gd.theta_b,gd.hc,N,'r');
% angle    = gd{'angle'}(:);
% mask_u = gd{'mask_u'}(:);
% mask_v = gd{'mask_v'}(:);
warning off
mask_rho = mask_rho./mask_rho;
% mask_u = mask_u./mask_u;
% mask_v = mask_v./mask_v;
warning on
% vname = {'temp','salt'};%,'zeta','ubar','vbar','u','v','omega'};
vname = {'TEMP_MONTHLY','SALT_MONTHLY'};%,'zeta','ubar','vbar','u','v','omega'};
     
for im=start_mm:time_step:end_mm
		mid=[num2char(im,2)];
%         file = [file_dir,'ocean_avg_mon_',mid,'.nc'];
        file = ['ocean_avg_mon_',mid,'.nc']; % original_stlee
%         file = ['ocean_avg_nwp_yyk.nc'];
        disp([file,' : ', num2char(im,2)])
            nc=netcdf(file);
            date=[num2str(yy),'. ',num2str(im)];
            switch variation
                case 1
                    value=nc{char(vname(variation))}(:);
                    value(find(value<-1000))=nan; % stlee
                    val_name='Temperature';
                    unit = '^oC';
                    out_name_1=['vertical_NS_Temp-',num2str(yy),'-'];
                    val_caxis=temp_lim;
                    level_c=temp_c;
                    data=value;
                   clear value;         
                case 2
                    value=nc{char(vname(variation))}(:);
                    val_name='Salinity';
                    unit = 'psu';
                    out_name_1=['vertical_NS_Salt-',num2str(yy),'-'];
                    val_caxis=salt_lim;
                    level_c=salt_c;
                    data=value;
                   clear value;     
            end
   

            if (plot_pcolorjw)
                for i=1:1:length(data(:,1))
                    for j=1:1:length(data(1,:))
                        if data(i,j) > 10000
                            data(i,j) = NaN;
                        end
                    end
                end
                
                switch section
                    case 1
%                 domaxis=[38 39.5 120 124 -100 0]; % 발해
%                 domaxis=[35.28 33.42 128.88 130.09 -100 0]; % 대한해협 along_stlee  
%                 domaxis=[35.28 33.42 128.88 130.09 -180 0]; % 대한해협 (NWP기준) along_stlee  
%                 domaxis=[41 42 140.8 141 -50 0]; % 쓰가루 stlee
%                 domaxis=[41 42 140.5 140.7 -190 0]; % 쓰가루 stlee
                  domaxis=[45.1872096987 46.4790236999 142.043821201 142.151473264  -100 0]; % 소야 stlee
%                 domaxis=[34.85 33.2265 128.6 129.5 -100 0]; % 대한해협 along_stlee2
%                 domaxis=[128.88 35.28 130.09 33.42] % 대한해협 along_stlee  
%                 domaxis=[30 36 122 125 -200 0]; % 발해
%                 domaxis=[23 30 118.3 124 -100 0]; % 대만해협 along
%                 domaxis=[25 24 118.5 121 -100 0]; % 대만해협 cross
                dist=sqrt((lon_rho-domaxis(3)).^2+(lat_rho-domaxis(1)).^2);
                min_dist=min(min(dist));
                dist2=sqrt((lon_rho-domaxis(4)).^2+(lat_rho-domaxis(2)).^2);
                min_dist2=min(min(dist2));                
                [x1,y1]=find(dist==min_dist);
                [x2,y2]=find(dist2==min_dist2);                
                lat1=lat_rho(x1(1),y1(1));lon1=lon_rho(x1(1),y1(1));
                lat2=lat_rho(x2(1),y2(1));lon2=lon_rho(x2(1),y2(1));
                if (lon2-lon1) >= (lat2-lat1)
                    lon_line=[lon1:0.05:lon2];
                    lat_line=(lon_line-lon1)/((lon2-lon1)/(lat2-lat1))+lat1;
                    x=repmat(lon_line,gd.N,1);
                    x_label='Longitude(^oE)';
                    domaxis=[domaxis(3) domaxis(4) domaxis(5) domaxis(6) domaxis(5) domaxis(6)];
                else
                    lat_line=[min(lat1,lat2):0.05:max(lat1,lat2)];
                    lon_line=(lat_line-lat1)*((lon2-lon1)/(lat2-lat1))+lon1;
                    x=repmat(lat_line,gd.N,1);
                    x_label='Latitude(^oN)';
                end
                Temp=zeros(gd.N,length(lat_line));
                for k=1:1:gd.N
                    lon_range=lon_rho(min(x1,x2):max(x1,x2),min(y1,y2):max(y1,y2));
                    lat_range=lat_rho(min(x1,x2):max(x1,x2),min(y1,y2):max(y1,y2));
                    data_range=squeeze(data(k,min(x1,x2):max(x1,x2),min(y1,y2):max(y1,y2)));
                    depth_range=squeeze(depth(k,min(x1,x2):max(x1,x2),min(y1,y2):max(y1,y2)));
                    Temp(k,:)=griddata(lon_range,lat_range,data_range,lon_line,lat_line);
                    Yi(k,:)=griddata(lon_range,lat_range,depth_range,lon_line,lat_line);
                end
                data=Temp;
                end
            end    
            
        figure('position',[400 100 550 550],'PaperUnits','inches','PaperPosition',[0 0 5.7 5.8]);
        set(gca,'Position',[0.2 0.15 0.73 0.75]);
        text_posi_x=(domaxis(2)-domaxis(1))/20+domaxis(1);
        text_posi_y1=(domaxis(6)-domaxis(5))/20+domaxis(5);
        text_posi_y2=2*(domaxis(6)-domaxis(5))/20+domaxis(5);
        text_posi_y3=3*(domaxis(6)-domaxis(5))/20+domaxis(5);
            switch section
                case 1
                hold on
                pcolor(x,Yi,data)
                axis([domaxis(1) domaxis(2) domaxis(5) domaxis(6)]);
                shading flat;caxis(val_caxis)
                set(gca,'box','on','linewidth',1.5,'fontsize',17)
                xlabel(x_label,'color',color_v,'FontSize',17,'fontweight','bold')
                ylabel('Depth(m)','color',color_v,'FontSize',17,'fontweight','bold')
                text(text_posi_x,text_posi_y2,val_name,'color',color_v,'FontSize',17,'fontweight','bold')
                text(text_posi_x,text_posi_y3,date,'color',color_v,'FontSize',17,'fontweight','bold')
                title('Vertical Temperature','fontsize',17); % stlee
                out_name_1=['YellowSea',out_name_1];
                if (plot_contour)
                  hold on
                  [C,h]=contour(x,Yi,data,level_c,color_c,'linewidth',1);
                  [C2,h2]=contour(x,Yi,data,level_c-0.5,color_c,'linewidth',1);
                  [C2,h2]=contour(x,Yi,data,[-1:2:1],'-w','linewidth',1);                  
                  clabel(C,h,'FontSize',15,'Color','k','labelspacing',100000,'Rotation',0,'fontweight','bold');
                end
            end
%             caxis(val_caxis);
            bar = colorbar('fontsize',17,'fontweight','bold');
            set(get(bar,'title'),'string',unit,'FontSize',17,'fontweight','bold')
            
            out_name=['Soya_vertical_temp_',num2char(im,2)];

            if (switch_save)
                saveas(gcf,out_name,out_type);
            end
            close all
end
%      close all
% % % end