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
clear all;close all;clc
%==========================================================================
addpath(genpath('D:\MEPL\project\SSH\중간보고\smooth13_vtvs\kyy_plot_subroutine'))

% 정선 307 - 36.925  // 308 - 36.33 // 309 - 35.855 // 310 - 35.335  // 311 - 34.716 // 312 - 33.975 ~ 34.0917
st_line=35;

section =1 ;  % 1 =yellow sea   //  2 =east sea

max_level= 40;

Vtransform=2;  Vstretching=4;

variation = 3;  %  1 = temperature , 2 = salinity , 3 = V velocity , 4 = E velocity , 5 = density , 6 = AKv, 7 = AKt

current   = 1;   skip_v = 3;     size_v = 2;    color_v = 'k';

plot_pcolorjw = 1; plot_contour  = 1;    color_c  ='-k' ;  

temp_lim = [0 20];  level_tem  = [ 0:2:35];
salt_lim = [31 35]; level_sal=[31:0.5:35];
vel_lim=[-5 5];   level_vel=[-10:1:10];
vele_lim=[-10 10];  level_vele=[-10:2:10];
den_lim= [23 28];  level_den=[23:0.2:28];
AKv_lim= [0 0.05];   level_AKv=[0:0.01:0.05];
AKt_lim= [0 0.05];   level_AKt=[0:0.01:0.05];
% plot_geoshow  = 0;    color_g = [.7 .7 .7]; 'black';

switch_save   = 1;    out_type = 'tif';

% grdfile       = 'd:\add2_ini_bry_grd\grid\roms_grid2_ADD_10_ep.nc';
yy = 1999;
% file_dir=[num2str(yy),'\monthly_mean4\'];
file_dir='D:\MEPL\project\SSH\중간보고\smooth13_vtvs\';
start_mm=1;
end_mm=2;
time_step=1;

mm=start_mm;
Vtransform = 2;  Vstretching = 4; theta_s = 7; theta_b = 2; hc = 250; N = max_level;
%==========================================================================

                
switch section
    case 1
    domaxis=[119 127 -100 0];
    case 2
    domaxis=[129.5 131.2 -500 0];    
end

gd = grd('NWP_1_20');
lon_rho  = gd.lon_rho;
lat_rho  = gd.lat_rho; 
mask_rho = gd.mask_rho;
h=gd.h;
N = gd.N;
dist=abs(st_line-lat_rho);
[num_stline,temp1,temp2]=find(dist==min(min(dist)));   %%수정요망 (거리가 최소가 되게끔)
num_stline=num_stline(1);

depth=zlevs(h,gd.zeta,gd.theta_s,gd.theta_b,gd.hc,N,'r');

% angle    = gd{'angle'}(:);
% mask_u = gd{'mask_u'}(:);
% mask_v = gd{'mask_v'}(:);
warning off
mask_rho = mask_rho./mask_rho;
% mask_u = mask_u./mask_u;
% mask_v = mask_v./mask_v;
warning on
vname = {'temp','salt','v','u','AKv','AKt'};%,'zeta','ubar','vbar','u','v','omega'};
str_month={'Jan.','Feb.','Mar.','Apr.','May','Jun.','Jul.','Aug.','Sep.','Oct.','Nov.','Dec.'};
for im=start_mm:time_step:end_mm
        mid=num2str(im+(yy-1992)*12,'%04i')
%         file = [file_dir,'ocean_avg_mon_',mid,'.nc'];
        file = [file_dir,'monthly_spinup_',num2str(mid,'%04i'),'.nc']; % kimyy
%         disp([file,' : ', num2char(im,2)])
        nc=netcdf(file);
%         date=[str_month{im},' ',num2str(yy)];
                    date=[num2str(yy),'. ',num2str(im)];

        zeta=nc{'zeta'}(:);
        switch variation
            case 1
                value=nc{char(vname(variation))}(:);
                val_name='Temperature';
                unit = '^oC';
                out_name_1=['vertical_Temp-',num2str(yy),'-'];
                val_caxis=temp_lim;
                level_c=level_tem;
                data=squeeze(value(:,num_stline,:)); 
            case 2
                value=nc{char(vname(variation))}(:);
                val_name='Salinity';
                unit = 'psu';
                out_name_1=['Sal',num2str(yy)];
                val_caxis=salt_lim;
                level_c=level_sal;
                data=squeeze(value(:,num_stline,:)); 
             case 3
                value=nc{char(vname(variation))}(:);
                val_name= ' '; % 'N-ward Velocity';
                unit = 'cm/s';
                out_name_1=['Nvel',num2str(yy)];
                val_caxis=vel_lim;
                level_c=level_vel;
                data=squeeze(value(:,num_stline,:))*100; 
             case 4
                value=nc{char(vname(variation))}(:);
                val_name='E-ward Velocity';
                unit = 'cm/s';
                out_name_1=['Evel',num2str(yy)];
                val_caxis=vele_lim;
                level_c=level_vele;
                data=squeeze(value(:,num_stline,:))*100;                 
             case 5
                temp=nc{'temp'}(:);
                salt=nc{'salt'}(:);
                val_name='Density';
                unit = 'σ';
                out_name_1=['Den',num2str(yy)];
                val_caxis=den_lim;
                level_c=level_den;
                data=sw_dens(squeeze(salt(:,num_stline,:)),squeeze(temp(:,num_stline,:)),0)-1000;
            case 6
                value=nc{char(vname(variation-1))}(:);
                val_name='Vertical eddy viscosity';
                unit = 'm/s';
                out_name_1=['AKv',num2str(yy)];
                val_caxis=AKv_lim;
                level_c=level_AKv;
                data=squeeze(value(2:end,num_stline,:)); 
            case 7
                value=nc{char(vname(variation-1))}(:);
                val_name='Vertical Tracer diffusion';
                unit = 'm/s';
                out_name_1=['AKt',num2str(yy)];
                val_caxis=AKt_lim;
                level_c=level_AKt;
                data=squeeze(value(2:end,num_stline,:)); 
        end
        Yi=squeeze(depth(:,num_stline,:));
        clear value;         

        if (plot_pcolorjw)
            for i=1:1:length(data(:,1))
                for j=1:1:length(data(1,:))
                    if data(i,j) > 100
                        data(i,j) = NaN;
                    end
                end
            end                       
            x_1=find(lon_rho(1,:)>=domaxis(1));
            x_2=find(lon_rho(1,:)>=domaxis(2));
            x=lon_rho(num_stline,x_1(1):x_2(1));x=repmat(x,40,1);
            data=data(:,x_1(1):x_2(1));
            Yi=Yi(:,x_1(1):x_2(1));
        end               
        
        figure('position',[87 200 500 500],'PaperUnits','inches','PaperPosition',[1.22 0.66 8 6.02]);
        set(gca,'Position',[0.13 0.11 0.825 0.81]);
        text_posi_x=(domaxis(2)-domaxis(1))/20+domaxis(1);
        text_posi_y1=(domaxis(4)-domaxis(3))/20+domaxis(3);
        text_posi_y2=2.5*(domaxis(4)-domaxis(3))/20+domaxis(3);
        text_posi_y3=4*(domaxis(4)-domaxis(3))/20+domaxis(3);
        text_posi_y21=21*(domaxis(4)-domaxis(3))/20+domaxis(3);
            switch section
                case 0
                hold on
                pcolor(Xi,Yi,data)
                axis([domaxis(1) domaxis(2) domaxis(3) domaxis(4)]);
                shading interp;caxis(val_caxis)
                set(gca,'box','on','linewidth',1.5,'fontsize',17)
                lab_lat=['Latitude',num2str(lat_rho(num_stline,1)),'(^oN)'];
                xlabel('Longitude(^oE)','color',color_v,'FontSize',17,'fontweight','bold')
                ylabel('Depth(m)','color',color_v,'FontSize',17,'fontweight','bold')
                text(text_posi_x,text_posi_y1,lab_lat,'color',color_v,'FontSize',17,'fontweight','bold')
                text(text_posi_x,text_posi_y2,val_name,'color',color_v,'FontSize',17,'fontweight','bold')
%                 text(text_posi_x,text_posi_y3,date,'color',color_v,'FontSize',17,'fontweight','bold')
                if (plot_contour)
                  hold on
                  [C,h]=contour(x,-Yi,Zi,level_c,color_c,'linewidth',1);
                  [C2,h2]=contour(x,-Yi,Zi,level_c-1,color_c,'linewidth',1);
                  clabel(C,h,'FontSize',10,'Color','k','labelspacing',100000,'Rotation',0,'fontweight','bold');
                end
                case 1
                hold on
                pcolor(x,Yi,data);
                shading interp;caxis(val_caxis)
                set(gca,'box','on','linewidth',1.5,'fontsize',17)
                axis([domaxis(1) domaxis(2) domaxis(3) domaxis(4)]);
                lab_lat=['Latitude ',num2str(st_line),'^oN'];
                xlabel('Longitude (^oE)','color',color_v,'FontSize',17,'fontweight','bold')
                ylabel('Depth (m)','color',color_v,'FontSize',17,'fontweight','bold')
                text(text_posi_x,text_posi_y2,lab_lat,'color',color_v,'FontSize',22,'fontweight','bold')
                text(text_posi_x,text_posi_y2,val_name,'color',color_v,'FontSize',20,'fontweight','bold')
%                 text(text_posi_x,text_posi_y3,date,'color',color_v,'FontSize',17,'fontweight','bold')
                text(text_posi_x,text_posi_y21,'(c)','color',color_v,'FontSize',22,'fontweight','bold')
                out_name_1=['YellowSea',out_name_1];
                if (plot_contour)
                  hold on
                  [C,h]=contour(x,Yi,data,level_c,color_c,'linewidth',1);
%                   [C2,h2]=contour(x,Yi,data,level_c+1,color_c,'linewidth',1);
                  clabel(C,h,'Rotation',0,'Color','k','fontsize',24);

%                   clabel(C,h,'FontSize',15,'Color','k','labelspacing',100000,'Rotation',0,'fontweight','bold');
                end
                case 2 
                hold on
                pcolor(x,Yi,data)
                shading interp;caxis(val_caxis)
                set(gca,'box','on','linewidth',1.5,'fontsize',17)
                axis([domaxis(1) domaxis(2) domaxis(3) domaxis(4)])
                lab_lat=['Latitude',num2str(st_line),'(^oN)'];
                xlabel('Longitude (^oE)','color',color_v,'FontSize',17,'fontweight','bold')
                ylabel('Depth(m)','color',color_v,'FontSize',17,'fontweight','bold')
                text(text_posi_x,text_posi_y1,lab_lat,'color',color_v,'FontSize',17,'fontweight','bold')
                text(text_posi_x,text_posi_y2,val_name,'color',color_v,'FontSize',17,'fontweight','bold')
                text(text_posi_x,text_posi_y3,date,'color',color_v,'FontSize',17,'fontweight','bold')
                out_name_1=['EastSea',out_name_1];
                if (plot_contour)
                  hold on
                  [C,h]=contour(x,Yi,data,level_c,color_c,'linewidth',1);
%                   [C2,h2]=contour(x,Yi,data,level_c-1,color_c,'linewidth',1);
                  clabel(C,h,'FontSize',15,'Color','k','labelspacing',100000,'Rotation',0,'fontweight','bold');
                end
            end
%             if variation ==3 ||variation ==4
%                 color=cbrewer('div','RdBu',50);color=flipud(color);
%                 colormap(color)
%             end
            bar = colorbar('fontsize',22,'fontweight','bold');
            set(get(bar,'title'),'string',unit,'FontSize',22,'fontweight','bold')
%             title(['CONTROL  ',date],'fontsize',18)            
            out_name=[file_dir,num2str(floor(st_line)),out_name_1,num2str(mid,'%02i')];
%                         out_name=['09_13_YellowSea_temperature_35',num2str(mid,'%02i')];
%             pts = ginput(1)
% point=[num2str(pts(1))]
% pts = ginput(2)
% point=[num2str(pts(2))]
            if (switch_save)
                saveas(gcf,out_name,out_type);
            end
            close all
end
%      close all
% % % end