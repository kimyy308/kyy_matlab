% function status=plot_SSS(sssfile, workdir, lonlat, year, inputmonth, inputdir)
% addpath(genpath('D:\MEPL\project\NWP\m_map'))

% % lonlat   : [lon_start lon_end lat_start lat_end]
% % year     : [year_start year_end]
% % filename : optional input

sssfile='D:\MEPL\project\NWP\Roms_tools\WOA1998\woa2013_salt.nc'
sstfile='D:\MEPL\project\NWP\Roms_tools\WOA1998\woa2013_temp.nc'

load C:\Users\KYY\Dropbox\source\matlab\Common\Figure\jet_mod
lonlat=[128 144 34 52] %% EJS
% lonlat=[115 164 15 52] %% NWP
calendar=cell(1,12);
calendar{1} = 'January'; calendar{2} = 'February'; calendar{3} = 'March'; calendar{4} = 'April'; calendar{5} = 'May'; calendar{6} = 'June';
calendar{7} = 'July'; calendar{8} = 'August'; calendar{9} = 'September'; calendar{10} = 'October'; calendar{11} = 'November'; calendar{12} = 'December';
lon_salt = ncread(sssfile,'LON');
lat_salt = ncread(sssfile,'LAT');
SSS = ncread(sssfile,'S_AN',[1 1 1 1], [length(lon_salt(:,1)) length(lat_salt(:,1)) 1 12]);
lon_temp = ncread(sstfile,'LON');
lat_temp = ncread(sstfile,'LAT');
SST = ncread(sstfile,'T_AN',[1 1 1 1], [length(lon_temp(:,1)) length(lat_temp(:,1)) 1 12]);
size(SSS)
for month = 1:12
    % plot
    m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
    m_grid('fontsize',20, 'box', 'fancy');
    hold on;
    m_pcolor(lon_salt,lat_salt,squeeze(SSS(:,:,1,month))');
    shading interp;
    m_gshhs_i('color','k')  
    m_gshhs_i('patch',[.8 .8 .8]);   % gray colored land (need much time)
    titlename = strcat('WOA 2013 SSS (',char(calendar(month)), ')');
    title(titlename,'fontsize',25);

    % set colorbar 
    h = colorbar;
    colormap(jet_mod);
    set(h,'fontsize',20);
    title(h,'','fontsize',15);
    caxis([33 35]);
%             colorbar('FontSize',11,'YTick',30:36,'YTickLabel',30:36);
% % log scale
%             caxis(log([1 36]));
%             colorbar('FontSize',11,'YTick',log(1:36),'YTickLabel',1:36);

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
    set(gcf,'Position',[200 100 xscale yscale])  % [hor_position ver_position width height]
%             jpgname=strcat(workdir,'figure\','nwp_1_20_SSS_',num2str(tempyear),'_',num2str(tempmonth),'.jpg'); % year_month       

    jpgname=strcat('D:\MEPL\project\SSH\1st_year\figure\WOA\2013\SSS_EJS_',num2str(month,'%04i'),'.jpg'); %% ~_month.jpg
    saveas(gcf,jpgname,'jpg');

    disp(' ')
    disp([' Making SSS plot is completed.'])
    disp(' ')
    disp([' File path is : ',jpgname])
    disp(' ')
    disp(' ')

    close all;
    
    
    % plot
    m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
    m_grid('fontsize',20, 'box', 'fancy');
    hold on;
    m_pcolor(lon_temp,lat_temp,squeeze(SST(:,:,1,month))');
    shading interp;
    m_gshhs_i('color','k')  
    m_gshhs_i('patch',[.8 .8 .8]);   % gray colored land (need much time)
    titlename = strcat('WOA 2013 SST (',char(calendar(month)), ')');
    title(titlename,'fontsize',25);

    % set colorbar 
    h = colorbar;
    colormap(jet_mod);
    set(h,'fontsize',20);
    title(h,'','fontsize',15);
    caxis([-2 33]);
%             colorbar('FontSize',11,'YTick',30:36,'YTickLabel',30:36);
% % log scale
%             caxis(log([1 36]));
%             colorbar('FontSize',11,'YTick',log(1:36),'YTickLabel',1:36);

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
    set(gcf,'Position',[200 100 xscale yscale])  % [hor_position ver_position width height]
%             jpgname=strcat(workdir,'figure\','nwp_1_20_SSS_',num2str(tempyear),'_',num2str(tempmonth),'.jpg'); % year_month       

    jpgname=strcat('D:\MEPL\project\SSH\1st_year\figure\WOA\2013\SST_EJS_',num2str(month,'%04i'),'.jpg'); %% ~_month.jpg
    saveas(gcf,jpgname,'jpg');

    disp(' ')
    disp([' Making SSS plot is completed.'])
    disp(' ')
    disp([' File path is : ',jpgname])
    disp(' ')
    disp(' ')

    close all;
    status=1; 
end
%         end
% status= 1;
% return