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



% % read model domain

name = 'monthly_spinup_';
filename_suffix = '.nc';
% ex : ~workdir\monthly_spinup_0001.nc
filename = strcat('D:\MEPL\project\SSH\중간보고\smooth13_vtvs\',name,num2str(0001,'%04i'),filename_suffix);
workdir = 'D:\MEPL\project\SSH\figure\avhrr';

lon = ncread(filename,'lon_rho');
lat = ncread(filename,'lat_rho');
% SST = ncread(filename,'temp',[1 1 40 1], [length(lon(:,1)) length(lat(1,:)) 1 1]);
% size(SST);
% lonlat=[115 164 15 52];  %% nwp
lonlat=[127 144 33 52];  %% east sea
 
calendar=cell(1,12);
calendar{1} = 'Jan'; calendar{2} = 'Feb'; calendar{3} = 'Mar'; calendar{4} = 'Apr'; calendar{5} = 'May'; calendar{6} = 'Jun';
calendar{7} = 'Jul'; calendar{8} = 'Aug'; calendar{9} = 'Sep'; calendar{10} = 'Oct'; calendar{11} = 'Nov'; calendar{12} = 'Dec';

% % for climate
load 'D:\MEPL\project\SSH\data\OISST_monthly\avhrr_n_model_SST.mat'


for month =1:12

    oisstlon = ncread('D:\MEPL\project\SSH\data\OISST_monthly\avhrr-monthly_201312.nc','lon');
    oisstlat = ncread('D:\MEPL\project\SSH\data\OISST_monthly\avhrr-monthly_201312.nc','lat');

    % % lonlat   : [lon_start lon_end lat_start lat_end]
    % % year     : [year_start year_end]
    % % month    : [first month of the year_start, last month of the year_end]
    % % filename : optional input

    load D:\MEPL\project\SSH\중간보고\smooth13_vtvs\kyy_plot_subroutine\jet_mod
    
%     startmonth=(year(1)-1992)*12+inputmonth(1)
%     endmonth=(year(2)-1992)*12+inputmonth(2)


%     for month = startmonth:endmonth
%         tempyear = int32(fix(month/12) +1);
%         tempmonth = mod(month,12);
%         if (tempmonth==0) 
%             tempmonth=12;
%             tempyear=tempyear-1;
%         end

        % plot
        m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        m_grid('fontsize',20);  %% for ES = 20 for NWP = 25
        hold on;
%         m_pcolor(lon,lat,SST);
        m_pcolor(lon,lat,squeeze(comb_OISST(:,:,month)));
        shading interp;
        m_gshhs_i('color','k')  
%         if (fast==0)
            m_gshhs_i('patch',[.8 .8 .8]);   % gray colored land (need much time)
%         end
%         tempyear = 5;
%         titlename = strcat('SST (',num2str(tempyear),' year,',' ',char(calendar(tempmonth)), ')');
%         titlename = strcat('AVHRR SST (',char(calendar(month)), ', 2012)');
        titlename = strcat('AVHRR SST (',char(calendar(month)), ',(1982-2011))');
%         title(titlename,'fontsize',25);  %% for nwp
        title(titlename,'fontsize',20);  %% for ES

        % set colorbar 
        h = colorbar;
        colormap(jet_mod);
        set(h,'fontsize',20);
        title(h,'(^oC)','fontsize',20);
        caxis([-2 33]);

% %     contour
        hold on;
        [C,h2]=m_contour(lon,lat,squeeze(comb_OISST(:,:,month)),-2:2:33,'k','linewidth',0.5);            
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

        set(gcf,'Position',[200 100 xscale yscale])  % [hor_position ver_position width height]
%             jpgname=strcat(workdir,'figure\','nwp_1_20_SST_',num2str(tempyear),'_',num2str(tempmonth),'.jpg'); % year_month            
        jpgname=strcat(workdir,'\','OISST_ES_1_20_SST_',num2str(month,'%04i'),'.jpg'); %% ~_month.jpg
        saveas(gcf,jpgname,'jpg');

        disp(' ')
        disp([' Making SST plot is completed.'])
        disp(' ')
        disp([' File path is : ',jpgname])
        disp(' ')
        disp(' ')
        hold off;
        close all;
        status=1; 

%     end
end
