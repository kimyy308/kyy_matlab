close all; clc; clear all;


grdname='E:\Data\Model\ROMS\ysecs\input\roms_grd.nc';
mask_rho=ncread(grdname, 'mask_rho');
lon_rho=ncread(grdname, 'lon_rho');
lat_rho=ncread(grdname, 'lat_rho');
lon_u=ncread(grdname, 'lon_u');
lat_u=ncread(grdname, 'lat_u');
lon_v=ncread(grdname, 'lon_v');
lat_v=ncread(grdname, 'lat_v');
h=ncread(grdname, 'h');
mask_rho(mask_rho==0)=NaN;

Vtransform=1;
Vstretching=1;
theta_s=5;
theta_b=0.4;
N=20;
hc=4;
var={'salt'};

calendar{1} = ' Jan'; calendar{2} = ' Feb'; calendar{3} = ' Mar'; 
calendar{4} = ' Apr'; calendar{5} = ' May'; calendar{6} = ' Jun';
calendar{7} = ' Jul'; calendar{8} = ' Aug'; calendar{9} = ' Sep'; 
calendar{10} = ' Oct'; calendar{11} = ' Nov'; calendar{12} = ' Dec';
for month = 1:12
    for varind=1:length(var)
%         tempyear=year;
        tempmonth=month;
        bryname=['E:\Data\Model\ROMS\ysecs\input\roms_bndy', '.nc'];
        bryinfo=ncinfo(bryname);
        % Vtransform=ncread(bryname,'Vtransform');
        % Vstretching=ncread(bryname,'Vstretching');
        % hc=ncread(bryname,'hc');
        varname=var{varind};
        switch varname %%read data, set name and unit
            case 'temp'
                unit = '^oC';
               data_south=ncread(bryname,'temp_south');
            case 'salt'
                unit = ' ';
               data_south=ncread(bryname,'salt_south');
            case 'v'
               unit = 'm/s';
               databar_south=ncread(bryname,'vbar_south');
               data3_south=ncread(bryname,'v_south');
               for k=1:20
                 databar3_south(:,k,:)=databar_south(:,:);
               end
               data_south=data3_south + databar3_south;
        end

        zeta_south=ncread(bryname,'zeta_south');
        z_south=zlevs(Vtransform, Vstretching, h(:,1), zeta_south(:,tempmonth), theta_s, theta_b, hc, N,'r');
        lon_rho_vert=repmat(lon_rho(:,1), [1 20]);
        mask_rho_vert=repmat(mask_rho(:,1), [1 20]);

        pcolor(lon_rho_vert,z_south',squeeze(data_south(:,:,tempmonth)).*mask_rho_vert)
        shading flat
        colorbar
        colormap(jet)
        xlabel('lon(^o)','color','k','FontSize',17,'fontweight','bold')
        ylabel('Depth(m)','color','k','FontSize',17,'fontweight','bold')
        title([varname, ',', 'clim', ',', calendar{tempmonth}]);

        switch varname %%read data, set name and unit
            case 'temp'
               caxis([0 35])
            case 'salt'
               caxis([32 35])
            case 'v'
               caxis([-0.1 0.1])
        end
        saveas(gcf, ['D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\4th_year\Figure\YSECS\south_', ...
            varname,'_','clim','_',num2str(tempmonth,'%02i')], 'tif');

        ylim([-200 0])   
        saveas(gcf, ['D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\4th_year\Figure\YSECS\south_', ...
            varname,'_200_','clim','_',num2str(tempmonth,'%02i')], 'tif');
        close all;
    end
end
