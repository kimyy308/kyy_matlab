clc; clear all; close all;

gridpath = 'F:\Toolkit\Backup\D\ROMS\data\MICT\NWP\';
gridname = 'roms_grid_nwp_1_10_test06.nc';
nc = netcdf([gridpath,gridname]);
lat_rho = nc{'lat_rho'}(:);
lon_rho = nc{'lon_rho'}(:);
maskr = nc{'mask_rho'}(:);
mon = {'Jan','Feb','Mar'};
[L,M] = size(lon_rho);
variable_u = 'Uwind';
variable_v = 'Vwind';
close(nc)

for year = 1993:1998
    
    yy = num2str(year);
    lat = lat_rho(:,1)'; lat = permute(lat, [2 1]);
    lon = lon_rho(1,:); lon = permute(lon, [2 1]);
    
    for month = 2
        mm = num2char(month,2);
        
        filepath_u = ['F:\Toolkit\Backup\D\ROMS\data\MICT\NWP\monthly\Uwind\'];
        filepath_v = ['F:\Toolkit\Backup\D\ROMS\data\MICT\NWP\monthly\Vwind\'];
        filename_u = ['Uwind_',yy,mm,'.nc'];
        filename_v = ['Vwind_',yy,mm,'.nc'];
        
        nc_u = netcdf([filepath_u,filename_u]);
        nc_v = netcdf([filepath_v,filename_v]);
        value_u = nc_u{variable_u}(:);
        value_v = nc_v{variable_v}(:);
        u = value_u .* maskr ./maskr;
        v = value_v .* maskr ./maskr;
        text_mon = mon(month);
        
        [curlZ]=ra_windstrcurl_lat(lat,lon,u,v);
        curlZ = curlZ * 10^7;
        curlZ = curlZ .* maskr ./ maskr;
        
        [taux tauy] = ra_windstr(u,v);
        taux = taux .* maskr ./ maskr;
        tauy = tauy .* maskr ./ maskr;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        domaxis = [127 142.5 34 52];
        domain = 'eastsea';
        text_lon = 129;
        text_lat = 50;
        
        figure('position',[400 100 800 600],'PaperUnits','inches','PaperPosition',[0 0 5.7 5.8]);
        set(gca,'Position',[0.15 0.06 0.8 0.8]);
        m_proj('mercator',...
            'lon',[domaxis(1) domaxis(2)],...
            'lat',[domaxis(3) domaxis(4)]);
        hold on
        m_gshhs_h('color','k');
        m_gshhs_h('patch',[.7 .7 .7]);
       m_grid('linewi',2,'linest','none','tickdir','out','fontsize',15);
        m_text(text_lon,text_lat,text_mon,'fontsize',17,'fontweight','bold')
        
        m_pcolor(lon_rho, lat_rho, curlZ); shading interp;
        h = colorbar;
%         caxis([-4 4])
        colormap(bluewhitered)
        title(h,'10^-7(N/m^3)'); ylabel(h,'Wind stress curl','fontsize',16);
        set(gca,'fontsize',15);
        
        limit = 10; interval = 5;
        vector_size = 1/10;
        
        H = m_quiver(lon_rho(1:interval:end,1:interval:end), lat_rho(1:interval:end,1:interval:end),...
            u(1:interval:end,1:interval:end)*vector_size, v(1:interval:end,1:interval:end)*vector_size,...
            'Color','k','AutoScale','off','maxheadsize',2);
        
        m_quiver(129,49,limit,0,vector_size,'color','b','AutoScale','off',...
            'linewidth',2,'maxheadsize',2);
        m_text(129,48.6,[num2str(limit),'m/s'])
        title(['Wind(',yy,')'],'fontsize',25)
        
        savepath = ['D:\ROMS\data\MICT\NWP\monthly\curl\'];
        savename = [domain,'_curl_',yy,mm];
        saveas(gcf,[savepath,savename],'jpg'); close;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         domaxis = [127 135 33.5 45];
%         domain = 'EKB';
%         text_lon = 127.7;
%         text_lat = 44;
%         
%         figure('position',[400 100 800 600],'PaperUnits','inches','PaperPosition',[0 0 5.7 5.8]);
%         set(gca,'Position',[0.15 0.06 0.8 0.8]);
%         m_proj('mercator',...
%             'lon',[domaxis(1) domaxis(2)],...
%             'lat',[domaxis(3) domaxis(4)]);
%         hold on
%         
%         m_gshhs_h('color','k');
%         m_gshhs_h('patch',[.7 .7 .7]);
%         m_grid('linestyle','none');
%         m_text(text_lon,text_lat,[yy,'.',mm],'fontsize',17,'fontweight','bold')
%         
%         m_pcolor(lon_rho, lat_rho, curlZ); shading interp;
%         h = colorbar;
%         colormap(bluewhitered)
%         title(h,'10^-7(N/m^3)'); ylabel(h,'Wind stress curl','fontsize',16);
%         set(gca,'fontsize',15);
%         
%         if month == 1
%             caxis([-3 3])
%         elseif month == 2
%             caxis([-3 3])
%         end
%           
%         limit = 5; interval = 5;
%         vector_size = 1/10;
%         
%         limit = 10; interval = 5;
%         vector_size = 1/10;
%         
%         H = m_quiver(lon_rho(1:interval:end,1:interval:end), lat_rho(1:interval:end,1:interval:end),...
%             u(1:interval:end,1:interval:end)*vector_size, v(1:interval:end,1:interval:end)*vector_size,...
%             'Color','k','AutoScale','off','maxheadsize',2);
%         
%         m_quiver(129,49,limit,0,vector_size,'color','b','AutoScale','off',...
%             'linewidth',2,'maxheadsize',2);
%         m_text(129,48.6,[num2str(limit),'m/s'])
%         title(['Wind(',yy,')'],'fontsize',25)
%         
%         savepath = ['D:\ROMS\data\MICT\NWP\monthly\curl\'];
%         savename = [domain,'_curl2_',yy,mm];
%         saveas(gcf,[savepath,savename],'jpg'); close;
        
    end
end