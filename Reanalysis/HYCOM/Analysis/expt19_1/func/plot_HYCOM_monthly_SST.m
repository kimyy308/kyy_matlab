function status=plot_ROMS_monthly_SST(testname, outfile, filedir, lonlat, tempyear, inputmonth, shadlev, conlev)
load C:\Users\KYY\Dropbox\source\matlab\Common\Figure\jet_mod
    
    calendar=cell(1,12);
    calendar{1} = ' January'; calendar{2} = ' February'; calendar{3} = ' March'; calendar{4} = ' April'; calendar{5} = ' May'; calendar{6} = ' June';
    calendar{7} = ' July'; calendar{8} = ' August'; calendar{9} = ' September'; calendar{10} = ' October'; calendar{11} = ' November'; calendar{12} = ' December';
    for monthij = 1:length(inputmonth)
        tempmonth = inputmonth(monthij)
        % ex : rootdir\test37\data\2001\test37_monthly_2001_01.nc
        filename = strcat(filedir, num2str(tempyear,'%04i'), '\', ...
                testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');

        % read data
        lon = ncread(filename,'lon_rho');
        lat = ncread(filename,'lat_rho');
        SST = ncread(filename,'temp',[1 1 40 1], [length(lon(:,1)) length(lat(1,:)) 1 1]);
%         size(SST);

        % plot
        m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        m_grid('fontsize',20);   %% for nwp = 25, for es = 20
        hold on;
        m_pcolor(lon,lat,SST);
        shading interp;
        m_gshhs_i('color','k')  
        m_gshhs_i('patch',[.8 .8 .8]);   % gray colored land
        titlename = strcat('SST (',char(calendar(tempmonth)), ', ', num2str(tempyear),')');
        title(titlename,'fontsize',20);  %%title

        % contour
        [C,h2]=m_contour(lon,lat,SST,conlev,'k','linewidth',0.5);      %%for temp
        clabel(C,h2,'FontSize',18,'Color','k','labelspacing',100000,'Rotation',0,'fontweight','bold');

        set(h2,'LineWidth',1.5);
        
        % set colorbar 
        h = colorbar;
        colormap(jet_mod);
        set(h,'fontsize',20);
        title(h,'(^oC)','fontsize',20);
        caxis(shadlev);

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
%         xscale = 800; yscale = 920; %% temporary scale

        set(gcf,'Position',[200 100 xscale yscale])  % [hor_position ver_position width height]
        jpgname=strcat(outfile, '_', testname, '_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.jpg'); %% ~_year_month.jpg
        saveas(gcf,jpgname,'jpg');

        disp(' ')
        disp([' Making SST plot is completed.'])
        disp(' ')
        disp([' File path is : ',jpgname])
        disp(' ')
        disp(' ')

        close all;
        status=1; 
    end
return