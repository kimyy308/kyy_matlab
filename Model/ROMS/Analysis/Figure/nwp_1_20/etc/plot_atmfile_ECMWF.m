clear all; clc; close all;
windows=1;
if (windows ==1)
    % % for windows
    dropboxpath='C:\Users\KYY\Dropbox';
    addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
    addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
    addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
elseif (linux==1)
    dropboxpath='/home01/kimyy/Dropbox';
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
    addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
end

load C:\Users\KYY\Dropbox\source\matlab\Common\Figure\jet_mod

% % read grid

filename = 'D:\MEPL\project\NWP\make_init\grid\roms_grid_combine2_test34.nc';
lon = ncread(filename,'lon_rho');
lat = ncread(filename,'lat_rho');
lon_u = ncread(filename,'lon_u');
lat_u = ncread(filename,'lat_u');
lon_v = ncread(filename,'lon_v');
lat_v = ncread(filename,'lat_v');
mask_rho = ncread(filename, 'mask_rho');

% lonlat=[115 164 15 52];  %% nwp
% lonlat = [115 130 30 42]  %% yellow sea
% lonlat=[127 144 33 52];  %% east sea
lonlat=[140 162 32 52];  %% kuroshio Extension

name = 'nwp_1_20_';
year = '2001_'
yday = '365'
workdir = 'D:\MEPL\project\SSH\figure\atm\';
inputdir = 'D:\MEPL\project\NWP\forcing_code\'
filename_suffix = '.nc';

% % read Tair
vname = 'Tair_'
filename=strcat(inputdir,name,year,vname,yday,filename_suffix);
Tair =  ncread(filename,'Tair');

% % read Uwind
vname = 'Uwind_'
filename=strcat(inputdir,name,year,vname,yday,filename_suffix);
Uwind =  ncread(filename,'Uwind');

% % read Vwind
vname = 'Vwind_'
filename=strcat(inputdir,name,year,vname,yday,filename_suffix);
Vwind =  ncread(filename,'Vwind');

name = 'Kuro_1_20_';  %% for jpgname


calendar=cell(1,12);
calendar{1} = 'Jan'; calendar{2} = 'Feb'; calendar{3} = 'Mar'; calendar{4} = 'Apr'; calendar{5} = 'May'; calendar{6} = 'Jun';
calendar{7} = 'Jul'; calendar{8} = 'Aug'; calendar{9} = 'Sep'; calendar{10} = 'Oct'; calendar{11} = 'Nov'; calendar{12} = 'Dec';



for month =2:2
    if strcmp(year,'2001_')
        day = [31 28 31 30 31 30 31 31 30 31 30 31];
    elseif strcmp(year,'2012_')
        day = [31 29 31 30 31 30 31 31 30 31 30 31];
    end

    % plot Tair
    m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
    
    if (strcmp(name,'ES_1_20_'))
        m_grid('fontsize',15);
    elseif (strcmp(name,'nwp_1_20_'))
        m_grid('fontsize',25);
    else
        m_grid('fontsize',20);
    end
    hold on;
    if (month==1)
        m_pcolor(lon,lat,mean(Tair(:,:,1:day(1)),3));
    else
        m_pcolor(lon,lat,mean(Tair(:,:,sum(day(1:month-1))+1:sum(day(1:month))),3));
    end
    shading interp;
    m_gshhs_i('color','k')  
    m_gshhs_i('patch',[.8 .8 .8]);   % gray colored land (need much time)

    titlename = strcat(year,' ECMWF Tair, ',char(calendar(month)));  %% for Tair
    title(titlename,'fontsize',25);

    % set colorbar 
    h = colorbar;
    colormap(jet_mod);
    set(h,'fontsize',25);
    title(h,'(^oC)','fontsize',25);  %% for temp
    if (strcmp(name,'Kuro_1_20_'))
        caxis([-20 15]);  %% for temp
    end
% %     contour
    hold on;
    if (month==1)
        [C,h2]=m_contour(lon,lat,mean(Tair(:,:,1:day(1)),3),'k','linewidth',0.5);
    else
        [C,h2]=m_contour(lon,lat,mean(Tair(:,:,sum(day(1:month-1))+1:sum(day(1:month))),3),'k','linewidth',0.5);
    end
%         [C,h2]=m_contour(lon,lat,squeeze(comb_modelSSS(:,:,month)),0:2:35,'k','linewidth',0.5);    %% for salt
%         [C,h2]=m_contour(lon,lat,squeeze(comb_modelSSS(:,:,month)),0:5:35,'k','linewidth',0.5);    %% for salt   
        clabel(C,h2,'FontSize',18,'Color','k','labelspacing',100000,'Rotation',0,'fontweight','bold');
%         clabel(C,h,'FontSize',15,'Color','k','labelspacing',300,'Rotation',0,'fontweight','bold');
        set(h2,'LineWidth',1.5);
    % make jpg file
    xscale=lonlat(2)-lonlat(1);
    yscale=lonlat(4)-lonlat(3);
    halt = 1;
    while(halt)
        if (xscale > 1000 || yscale > 1000)
            halt = 0;
        else
            xscale = xscale * 1.2; yscale = yscale * 1.2;
        end
    end
    xscale = 980; yscale = 920; %% temporary scale
    if (strcmp(name,'ES_1_20_'))
        xscale = 600; yscale =920;
    end
    set(gcf,'Position',[200 100 xscale yscale])  % [hor_position ver_position width height]
    jpgname=strcat(workdir,'\',name,'Tair_',year,num2str(month,'%04i'),'.jpg'); %% ~_month.jpg for temp
    saveas(gcf,jpgname,'jpg');

    disp(' ')
    disp([' Making Tair plot is completed.'])
    disp(' ')
    disp([' File path is : ',jpgname])
    disp(' ')
    disp(' ')
    hold off;
    close all;
    status=1; 

% % %         plot wind
    if (month==1)
        u(:,:,month)=mean(Uwind(:,:,1:day(1)),3);
        v(:,:,month)=mean(Vwind(:,:,1:day(1)),3);
    else
        u(:,:,month)=mean(Uwind(:,:,sum(day(1:month-1))+1:sum(day(1:month))),3);
        v(:,:,month)=mean(Vwind(:,:,sum(day(1:month-1))+1:sum(day(1:month))),3);
    end
%     u_rho=griddata(lon_u, lat_u, squeeze(u(:,:,month)), lon, lat);
%     v_rho=griddata(lon_v, lat_v, squeeze(v(:,:,month)), lon, lat);
    u_rho=u(:,:,month).*mask_rho;
    v_rho=v(:,:,month).*mask_rho;
    u_rho(200:205,700:705) = 10;  %% for nwp
    v_rho(200:205,700:705) = 0.001;  %% for nwp
%     u_rho(300:305,870:875) = 10;  %% for ES
%     v_rho(300:305,870:875) = 0.001;  %% for ES
    interval = 20;
    if (strcmp(name,'Kuro_1_20_'))
        interval = 40;
    end
    m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
    if (strcmp(name,'ES_1_20_'))
        m_grid('fontsize',15);
    else
        m_grid('fontsize',20);
    end
    m_gshhs_i('color','k')
    m_gshhs_i('patch',[.8 .8 .8]);
    hold on
    uvplot=m_quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end),...
         u_rho(1:interval:end,1:interval:end),v_rho(1:interval:end,1:interval:end),3,'k');  %% for nwp
%     uvplot=m_quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end),...
%              u_rho(1:interval:end,1:interval:end),v_rho(1:interval:end,1:interval:end),9,'k');  %% for ES
    titlename = strcat(year,' ECMWF Wind, ',char(calendar(month)));
    title(titlename,'fontsize',25);
    if (strcmp(name,'ES_1_20_'))
        %     m_text(129,50,'10 m/s','FontSize',25)  %% for ES
    else
         m_text(123.9,43.9,'10 m/s','FontSize',25)  %% for nwp
    end

    
    set(gcf,'Position',[200 100 xscale yscale])
%             jpgname=strcat(workdir,'figure\','nwp_1_20_SST_',num2str(tempyear),'_',num2str(tempmonth),'.jpg'); % year_month            
    jpgname=strcat(workdir,'\',name,'wind_',year,num2str(month,'%04i'),'.jpg'); %% ~_month.jpg
    saveas(gcf,jpgname,'jpg');
    hold off;

    disp(' ')
    disp([' Making UV plot is completed.'])
    disp(' ')
    disp([' File path is : ',jpgname])
    disp(' ')
    disp(' ')

    close all;
    status=1; 
end
