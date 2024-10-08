function status=plot_SSS(sssfile, workdir, lonlat, year, inputmonth, inputdir)
% addpath(genpath('D:\MEPL\project\NWP\m_map'))

% % lonlat   : [lon_start lon_end lat_start lat_end]
% % year     : [year_start year_end]
% % filename : optional input

load C:\Users\KYY\Dropbox\source\matlab\Common\Figure\jet_mod

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
            filename = strcat(workdir,name,num2str(month,'%04i'),filename_suffix);

            % read data
            lon = ncread(filename,'lon_rho');
            lat = ncread(filename,'lat_rho');
            SSS = ncread(filename,'salt',[1 1 40 1], [length(lon(:,1)) length(lat(1,:)) 1 1]);
% %  log scale
%             SSS = log(abs(ncread(filename,'salt',[1 1 40 1], [length(lon(:,1)) length(lat(1,:)) 1 1])));
            
            size(SSS);
            
            % plot
            m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            m_grid('fontsize',20, 'box', 'fancy');
            hold on;
            m_pcolor(lon,lat,SSS);
            shading interp;
            m_gshhs_i('color','k')  
            m_gshhs_i('patch',[.8 .8 .8]);   % gray colored land (need much time)
            titlename = strcat('SSS (',num2str(tempyear),' year,',' ',char(calendar(tempmonth)), ')');
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
            jpgname=strcat(sssfile,num2str(month,'%04i'),'.jpg'); %% ~_month.jpg
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
status= 1;
return