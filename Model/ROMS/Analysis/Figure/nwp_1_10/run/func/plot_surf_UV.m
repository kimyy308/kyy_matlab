function status=plot_surf_UV(outfile, workdir, lonlat, year, inputmonth, inputdir, regionflag)
% addpath(genpath('D:\MEPL\project\NWP\m_map'))

% % lonlat   : [lon_start lon_end lat_start lat_end]
% % year     : [year_start year_end]
% % filename : optional input

startmonth=(year(1)-1992)*12+inputmonth(1)
endmonth=(year(2)-1992)*12+inputmonth(2)
calendar=cell(1,12);
calendar{1} = ' January'; calendar{2} = ' February'; calendar{3} = ' March'; calendar{4} = ' April'; calendar{5} = ' May'; calendar{6} = ' June';
calendar{7} = ' July'; calendar{8} = ' August'; calendar{9} = ' September'; calendar{10} = ' October'; calendar{11} = ' November'; calendar{12} = ' December';

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
    filename = strcat(inputdir,name,num2str(month,'%04i'),filename_suffix);
    % read data
    if (exist('lon_u')==0)
        lon_u = ncread(filename,'lon_u');
        lat_u = ncread(filename,'lat_u');
        lon_v = ncread(filename,'lon_v');
        lat_v = ncread(filename,'lat_v');
        lon = ncread(filename,'lon_rho');
        lat = ncread(filename,'lat_rho');
    end
    u = ncread(filename,'u',[1 1 40 1], [length(lon_u(:,1)) length(lat_u(1,:)) 1 1]);
    v = ncread(filename,'v',[1 1 40 1], [length(lon_v(:,1)) length(lat_v(1,:)) 1 1]);
    u_rho=griddata(lon_u, lat_u, u, lon, lat);
    v_rho=griddata(lon_v, lat_v, v, lon, lat);

%           Reference vector value setting
if (regionflag==1)
elseif (regionflag==2)
end
if (regionflag==1)
    u_rho(200:205,700:705) = 2;  %% for nwp
    v_rho(200:205,700:705) = 0.001;  %% for nwp
elseif (regionflag==2)
    u_rho(300:305,870:875) = 1;  %% for ES
    v_rho(300:305,870:875) = 0.001;  %% for ES
end


%            
% plot

m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
m_gshhs_i('color','k')
m_gshhs_i('patch',[.8 .8 .8]);  
hold on

if (regionflag==1)
        m_grid('fontsize',20, 'box', 'fancy');  %% for nwp
        interval = 10; %% for nwp
        uvplot=m_quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end),...
               u_rho(1:interval:end,1:interval:end),v_rho(1:interval:end,1:interval:end),6,'k');  %% for nwp
elseif (regionflag==2)
        m_grid('fontsize',15, 'box', 'fancy');  %% for ES
        interval = 5;  %% for ES
            uvplot=m_quiver(lon(1:interval:end,1:interval:end),lat(1:interval:end,1:interval:end),...
       u_rho(1:interval:end,1:interval:end),v_rho(1:interval:end,1:interval:end),9,'k');  %% for ES
end


%             geoshow('landareas.shp', 'FaceColor', [0.8 0.8 0.8]);



    titlename = strcat('Surface vector(',num2str(tempyear),' year,',' ',char(calendar(tempmonth)), ')');
    title(titlename,'fontsize',25);

%           Reference vector text setting
if (regionflag==1)
        m_text(123.9,43.9,'2 m/s','FontSize',25)  %% for nwp
elseif (regionflag==2)
        m_text(129,50,'1 m/s','FontSize',25)  %% for ES
end


    

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
    xscale = 800; yscale = 920; %% temporary scale
    set(gcf,'Position',[200 100 xscale yscale])
%             jpgname=strcat(workdir,'figure\','nwp_1_20_SST_',num2str(tempyear),'_',num2str(tempmonth),'.jpg'); % year_month            
    jpgname=strcat(outfile,'nwp_1_20_UV_',num2str(month,'%04i'),'.jpg'); %% ~_month.jpg
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
return