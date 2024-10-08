% function status=plot_SST(outfile, workdir, lonlat, year, inputmonth, inputdir)


clc; clear all; close all;

addpath(genpath('/scratch/snu02/roms_nwp/output/test34'))
% % lonlat   : [lon_start lon_end lat_start lat_end]
% % year     : [year_start year_end]
% % month    : [first month of the year_start, last month of the year_end]
% % filename : optional input


figdir ='/scratch/snu02/roms_nwp/output/test34/';
workdir = '/scratch/snu02/roms_nwp/output/test34/';
inputdir = '/scratch/snu02/roms_nwp/output/test34/';
outfile =[figdir, 'model_nwp_KE_all_nwp_1_20']
lonlat = [115 164 15 52]; % [lon_start lon_end lat_start lat_end] (NWP)
year = [1992 2000] % [year_start year_end]
inputmonth = [1 12] % [first month of the year_start, last month of the year_end]



startmonth=(year(1)-1992)*12+inputmonth(1)
endmonth=(year(2)-1992)*12+inputmonth(2)
name = 'monthly_spinup_';
filename_suffix = '.nc';
% ex : ~workdir\monthly_spinup_0001.nc
filename = strcat(workdir,name,num2str(startmonth,'%04i'),filename_suffix);

lon = ncread(filename,'lon_rho');
lat = ncread(filename,'lat_rho');
depth = ncread(filename,'h');

theta_s=7.0; theta_b = 2.0; N = 40; hc = 250;
z_w = zlevs(depth,0,theta_s,theta_b,hc,N,'w');
for i=1:N
    layer_thickness(i,:,:)=z_w(i+1,:,:)-z_w(i,:,:);
end
[L M] = size(depth);
u_thick=0.5*(layer_thickness(:,1:L-1,:)+layer_thickness(:,2:L,:));
v_thick=0.5*(layer_thickness(:,:,1:M-1)+layer_thickness(:,:,2:M));

for i=1:M-1
    zon_dis(i)=m_lldist([lon(1,1) lon(2,1)], [lat(1,i) lat(1,i)])*1000.;
    mer_dis(i)=m_lldist([lon(1,i) lon(1,i)], [lat(1,i) lat(1,i+1)])*1000.;
end
zon_dis(M)=m_lldist([lon(1,1) lon(2,1)], [lat(1,M) lat(1,M)])*1000.;

for i=1:M-1
    temp_area(i)=(zon_dis(i)+zon_dis(i+1)) * sqrt(mer_dis(i)^2-((zon_dis(i)-zon_dis(i+1))/2.0)^2)/2.0;
end


% slope_index : avg(abs(center depth - adjacent depth))
% and normalize  (slope_index/center depth)

% % read depth and get slope_index
% nlon = length(lon(1,:));
% nlat = length(lat(:,1));
% sum_ind(1:nlat,1:nlon)=0;
% slp_ind(1:nlat,1:nlon)=0;
% for j= 2:nlon-1
%     for i=2:nlat-1
%         sum_ind(i,j)=(abs(depth(i,j)-depth(i+1,j)) + abs(depth(i,j)-depth(i-1,j)) ...
%                         + abs(depth(i,j)-depth(i,j+1)) + abs(depth(i,j)-depth(i,j-1)))/4.0;
% %         slp_ind(i,j)=sum_ind(i,j)/depth(i,j);
%         slp_ind(i,j)=sum_ind(i,j);
%     end
% end
% clear sum_ind;
% slp_ind=depth;

% load D:\MEPL\project\SSH\�߰�����\smooth13_vtvs\kyy_plot_subroutine\jet_mod


%     calendar=cell(1,12);
%     calendar{1} = ' January'; calendar{2} = ' February'; calendar{3} = ' March'; calendar{4} = ' April'; calendar{5} = ' May'; calendar{6} = ' June';
%     calendar{7} = ' July'; calendar{8} = ' August'; calendar{9} = ' September'; calendar{10} = ' October'; calendar{11} = ' November'; calendar{12} = ' December';
    for month = startmonth:endmonth
        tempyear = int32(fix(month/12) +1);
        tempmonth = mod(month,12);
        if (tempmonth==0) 
            tempmonth=12;
            tempyear=tempyear-1;
        end
        name = 'monthly_spinup_';
        filename_suffix = '.nc';
        % ex : ~workdir\monthly_spinup_0001.nc
        filename = strcat(workdir,name,num2str(month,'%04i'),filename_suffix);
        
        % read data
        if (exist('lon_u')==0)
            lon_u = ncread(filename,'lon_u');
            lat_u = ncread(filename,'lat_u');
            lon_v = ncread(filename,'lon_v');
            lat_v = ncread(filename,'lat_v');
            lon = ncread(filename,'lon_rho');
            lat = ncread(filename,'lat_rho');
        end
%         u = ncread(filename,'u',[1 1 40 1], [length(lon_u(:,1)) length(lat_u(1,:)) 1 1]);
%         v = ncread(filename,'v',[1 1 40 1], [length(lon_v(:,1)) length(lat_v(1,:)) 1 1]);
%         u = ncread(filename,'u',[1 1 1 1], [L-1 M 40 1]);
%         v = ncread(filename,'v',[1 1 1 1], [L M-1 40 1]);
        u = ncread(filename,'u');
        v = ncread(filename,'v');        
        
%         rho = ncread(filename,'rho',[1 1 40 1], [length(lon(:,1)) length(lat(1,:)) 1 1]) + 1000.;
        rho = ncread(filename,'rho',[1 1 1 1], [length(lon(:,1)) length(lat(1,:)) 40 1]) + 1000.;
%         u_rho=griddata(lon_u, lat_u, u, lon, lat);
%         v_rho=griddata(lon_v, lat_v, v, lon, lat);

% for i=1:40
%         u = ncread(filename,'u',[1 1 i 1], [length(lon_u(:,1)) length(lat_u(1,:)) 1 1]);
%         v = ncread(filename,'v',[1 1 i 1], [length(lon_v(:,1)) length(lat_v(1,:)) 1 1]);
%         u_rho(:,:,41-i)=griddata(lon_u, lat_u, u, lon, lat);
%         v_rho(:,:,41-i)=griddata(lon_v, lat_v, v, lon, lat);
% end
%         kinetic_E = 1/2 * rho .* u_rho.^2 .* v_rho.^2;
        kinetic_E_u=0.5*rho(1:L-1,:,:,1).*u.^2;
        kinetic_E_v=0.5*rho(:,1:M-1,:,1).*v.^2;
        for k=1:40
            for j=1:M-1
                for i=1:L-1
                    kinetic_E_u(i,j,k)=kinetic_E_u(i,j,k)*u_thick(k,i,j)*temp_area(j);
                    kinetic_E_v(i,j,k)=kinetic_E_v(i,j,k)*v_thick(k,i,j)*temp_area(j);
                end
                kinetic_E_v(L,j,k)=kinetic_E_v(L,j,k).*v_thick(k,i,j).*temp_area(j);
            end
            kinetic_E_u(1:L-1,M,k)=kinetic_E_u(1:L-1,M,k).*u_thick(k,i,M).*temp_area(M-1);
        end
        a_n_u = isnan(kinetic_E_u);
        kinetic_E_u(a_n_u)=0;
        KE_TS_u=sum(sum(sum(kinetic_E_u)));
        a_n_v = isnan(kinetic_E_v);
        kinetic_E_v(a_n_v)=0;
        KE_TS_v=sum(sum(sum(kinetic_E_v)));
        KE_TS(month)=KE_TS_u + KE_TS_v;
% %     %           Reference vector value setting
% %         u_rho(200:205,700:705) = 2;  %% for nwp
% %         v_rho(200:205,700:705) = 0.001;  %% for nwp
% % %         u_rho(300:305,870:875) = 1;  %% for ES
% % %         v_rho(300:305,870:875) = 0.001;  %% for ES

        
        % plot
% %         m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
% %         m_grid('fontsize',20);   %% for nwp = 25, for es = 20
% %         hold on;
% %         m_pcolor(lon,lat,slp_ind);
% %         shading interp;
%         m_gshhs_i('color','k')  
% %         if (fast==0)
%             m_gshhs_i('patch',[.8 .8 .8]);   % gray colored land (need much time)
% %         end
%         titlename = strcat('slope & uv (',num2str(tempyear),' year,',' ',char(calendar(tempmonth)), ')');
%         title(titlename,'fontsize',25);  %%title

        % contour
%         [C,h2]=m_contour(lon,lat,SST,-2:2:33,'k','linewidth',0.5);      %%for temp (East Sea)
%         clabel(C,h2,'FontSize',15,'Color','k','labelspacing',300,'Rotation',0,'fontweight','bold');
%         clabel(C,h2,'FontSize',18,'Color','k','labelspacing',100000,'Rotation',0,'fontweight','bold');

%         set(h2,'LineWidth',1.5);
        
        % set colorbar 
% %         h = colorbar;
% %         colormap(jet_mod);
% %         set(h,'fontsize',20);
% %         title(h,'( )','fontsize',20);
%         caxis([-2 33]);

%        quiver
% %         interval = 5;
% %         m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
% %         m_grid('fontsize',20);
% %         m_gshhs_i('color','k')
% %         m_gshhs_i('patch',[.8 .8 .8]);  
% %     uvplot=m_quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end),...
% %          u_rho(1:interval:end,1:interval:end),v_rho(1:interval:end,1:interval:end),6,'k');  %% for nwp
%     uvplot=m_quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end),...
%              u_rho(1:interval:end,1:interval:end),v_rho(1:interval:end,1:interval:end),9,'k');  %% for ES

%           Reference vector text setting
% %     m_text(123.9,43.9,'2 m/s','FontSize',25)  %% for nwp
%     m_text(129,50,'1 m/s','FontSize',25)  %% for ES

        % make jpg file
% %         xscale=lonlat(2)-lonlat(1);
% %         yscale=lonlat(4)-lonlat(3);
% %         halt = 1;
% %         while(halt)
% %             if (xscale > 1000 || yscale > 1000)
% %                 halt = 0;
% %             else
% %                 xscale = xscale * 1.2; yscale = yscale * 1.2;
% %             end
% %         end
% %         xscale = 980; yscale = 920; %% temporary scale

%         set(gcf,'Position',[200 100 xscale yscale])  % [hor_position ver_position width height]
%             jpgname=strcat(workdir,'figure\','nwp_1_20_SST_',num2str(tempyear),'_',num2str(tempmonth),'.jpg'); % year_month            
%         jpgname=strcat(outfile,num2str(month,'%04i'),'.jpg'); %% ~_month.jpg
%         saveas(gcf,jpgname,'jpg');

        disp(' ')
        disp([' cal of KE is completed.', month])
        disp(' ')
%         disp([' File path is : ',jpgname])
%         disp(' ')
%         disp(' ')

        close all;
        status=1; 
    end
    
plot(KE_TS);
jpgname=strcat(outfile,'_',num2str(month,'%04i'),'.jpg'); %% ~_month.jpg
saveas(gcf,jpgname,'jpg');
% return
