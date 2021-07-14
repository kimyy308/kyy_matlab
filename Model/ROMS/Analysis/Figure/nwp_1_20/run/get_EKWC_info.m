% function status=get_EKWC_info(outfile, year, inputdir)

% horsection = [128 135 35.5 44 -235 -235];
% horfile =[figdir,'\EJS\NKCC_',num2str(abs(horsection(5))),'m_temp_uv_']; 

% clear all;close all;

%==========================================================================
% % This function needs 
% % parameter : vert_param.m (Vtransform, Vstretching, theta_s, theta_b, hc(Tcline), N(the number of vertical level);
% % function : 'read_grid.m(similar to grd.m)', 'stretching.m', 'zlevs.m'.
% % library : 'netcdf_old'

% % outfile       : [path] outfilename
% % section       : [lon_start, lon_end, lat_start, lat_end, depth_deep(negative), depth_shallow(negative or zero)]
% % year          : [year_start, year_end]
% % inputmonth    : [first month of the year_start, last month of the year_end]
% % fast          : [switch] option for 'm_gshhs_i' (gray colored land). 0 or 1(use)
% % var           : [var] kind of variables. 1=temp, 2=salt, 3=u, 4=v, 5=rho
% % var_lim       : [lowest value, highest value] color scale 
% % inputdir      : [path] inputfile path. 
% % full filename ex) inputdir\monthly_spinup_0001.nc
warning off;
vert_param;
% addpath(genpath('D:\MEPL\project\NWP\netcdf_old'))

% max_level= 40;
% var = 2;  %  1 = temperature , 2 = salinity ,3 = density, 4= meridional vel

% current   = 1;   skip_v = 3;     size_v = 2;    color_v = 'k';

% plot_pcolorjw = 1;    temp_lim = [5 20];    salt_lim = [31.0 35]; den_lim = [25 27.5]; v_lim = [-10 10];
 
% plot_contour  = 1;    color_c  ='-k' ;      temp_c  = [5:1:20];  salt_c  =[31:1:35]; v_c =[-10:2:10];

% plot_geoshow  = 0;    color_g = [.7 .7 .7];'black';

% switch_save   = 1;    out_type = 'png';

% section = 1; % 0 => whole, 1=> yellowsea, 2=> eastsea

% grdfile       = 'D:\MEPL\project\NWP\make_init\grid\roms_grid_combine2_test17.nc';
% Vtransform = 2;  Vstretching = 4; theta_s = 7; theta_b = 2; hc = 250; N = max_level;
% yy = 1999;
% start_mm=1;
% end_mm=1;


time_step=1;

% file_dir='D:\MEPL\project\SSH\占쌩곤옙占쏙옙占쏙옙\smooth13_vtvs\';
% mm=start_mm;
%==========================================================================

%% MAX v velocity at 102
%% width of EKWC (102에서 max vel 기준 양옆으로 efolding vel)
%% 102 단면에서의 temp & geostrophic velocity (500m 기준)
%% 모델에서의 유속과 geostrophic velocity와의 차이
%% 수평 vorticity, pycnocline depth 분포

%% 유속이 증가하는 위치가 있는지? 증가한다면 냉수가 있는 위치인가?
%% separation latitude (102 기준 e folding 감소 vel latitude)

%% 102 : (129.7967, 36.0767) , (130.9183, 36.0767)
%% year를 받아서, 한 text파일에 모든 정보 기록하도록.


% . text file open
% for
%     read data(section)
%     0. 102 라인 단면 유속(v), 수온 정보 읽어오기.
%     1. 102라인에서의 maximum velocity 위치.
%     2. efolding vel 구함
%     3. vel - efolding -> min 위치 찾음 max 동쪽, max 서쪽에서.
%     4. min 위치간의 거리 구함.
%     
%     . standard depth로 interpolation
%     . sea level 
%     . gvel, maximum surface gvel
%     . temp plot
%     . rho plot
%     . salt plot
%     . gvel plot
%     . eta plot
%     . thermocline, pycnocline depth plot
%     . uwind, vwind plot
%     . gvel 과 model vel 차이 plot
%     
%     .surface hor temp. salt, u, v read
%     .vorticity
%     .separation latitude
% end

startmonth=(year(1)-1992)*12+1
endmonth=(year(2)-1992)*12+12
mid=num2str(startmonth,'%04i');
file = [inputdir,'monthly_spinup_',mid,'.nc'];

% if (exist('gd.lon_rho')==0)
    gd = read_grid(file,Vtransform,Vstretching,theta_s,theta_b,hc,N);
% end

lon_rho  = gd.lon_rho;
lat_rho  = gd.lat_rho; 
mask_rho = gd.mask_rho;
h=gd.h;
N = gd.N;
% depth=zlevs(h,gd.zeta,gd.theta_s,gd.theta_b,gd.hc,N,'r',1);
% depth=zlevs(h,gd.zeta,gd.theta_s,gd.theta_b,gd.hc,N,'r');
% angle    = gd{'angle'}(:);
% mask_u = gd{'mask_u'}(:);
% mask_v = gd{'mask_v'}(:);
depth=gd.z_r;


vname = {'temp','salt', 'u', 'v', 'rho'};

calendar{1} = ' Jan'; calendar{2} = ' Feb'; calendar{3} = ' Mar'; 
calendar{4} = ' Apr'; calendar{5} = ' May'; calendar{6} = ' Jun';
calendar{7} = ' Jul'; calendar{8} = ' Aug'; calendar{9} = ' Sep'; 
calendar{10} = ' Oct'; calendar{11} = ' Nov'; calendar{12} = ' Dec';

for im=startmonth:time_step:endmonth
    mid=num2str(im,'%04i');   %% ex) '0036'
    file = [inputdir,'monthly_spinup_',mid,'.nc'];   %% ex) 'monthly_spinup_0036.nc'
    disp(['read  ' file])
    nc=netcdf(file);
    tempyear = int32(fix(im/12) +1);
    tempmonth = mod(im,12);
    if (tempmonth==0) 
        tempmonth=12;
        tempyear=tempyear-1;
    end
    date=[num2str(tempyear),' year, ',char(calendar(tempmonth))]; %% ex) 2 year, December

% %     read data (temp, salt, u, v, rho)
    temp=nc{'temp'}(:);
    temp(find(temp<-1000))=NaN; 
    temp_name='Temperature';
    temp_unit = '^oC';
    
    salt=nc{'salt'}(:);
    salt(find(salt<-1000))=NaN; 
    salt_name='Salinity';
    salt_unit = ' ';
    
    u=nc{'u'}(:);
    u(find(u<-1000))=NaN; 
    u_name='Zonal velocity';
    u_unit = 'm/s';
    
    v=nc{'v'}(:);
    v(find(v<-1000))=NaN; 
    v_name='Zonal velocity';
    v_unit = 'm/s';
    
    rho=nc{'rho'}(:);
    rho(find(rho<-1000))=NaN; 
    rho_name='Zonal velocity';
    rho_unit = 'm/s';
    
    section = [129.7967, 130.9183, 36.0767, 36.0767 -500 0]; %% NIFS data 102-06 ~ 102-10
    
    dist=sqrt((lon_rho-section(1)).^2+(lat_rho-section(3)).^2); %% get distance from station 1
    min_dist=min(min(dist));
    dist2=sqrt((lon_rho-section(2)).^2+(lat_rho-section(4)).^2);  %% get distance from station 2
    min_dist2=min(min(dist2));                
    [x1,y1]=find(dist==min_dist);  %% find closest point from station 1. [lat lon]
    [x2,y2]=find(dist2==min_dist2); %% find closest point from station 2. [lat lon]

    lat1=lat_rho(x1(1),y1(1));  lon1=lon_rho(x1(1),y1(1));
    lat2=lat_rho(x2(1),y2(1));  lon2=lon_rho(x2(1),y2(1));
%     if (lon2-lon1) >= (lat2-lat1)
        lon_line = lon1:mean(gradient(lon_rho(1,:))):lon2;  %% for 1/20^o horizontal resolution
        lat_line = (lon_line-lon1)/((lon2-lon1)/(lat2-lat1))+lat1; %% weight for latitude index (lat1 : 0.05/((lon2-lon1) * (lat2-lat1) : lat2 
        x=repmat(lon_line,gd.N,1);  %% copy lon_line (gd.N times) to make matrix 
        x_label='Longitude(^oE)';
%                 domaxis=[domaxis(3) domaxis(4) domaxis(5) domaxis(6) domaxis(5) domaxis(6)];
        temp_102=zeros(gd.N,length(lon_line)); %% initialize temp matrix that size is same with x
        salt_102=zeros(gd.N,length(lon_line));
        u_102=zeros(gd.N,length(lon_line));
        v_102=zeros(gd.N,length(lon_line));
        rho_102=zeros(gd.N,length(lon_line));
%     else
%         lat_line=[min(lat1,lat2):mean(gradient(lat_rho(:,1))):max(lat1,lat2)];
%         lon_line=(lat_line-lat1)*((lon2-lon1)/(lat2-lat1))+lon1;
%         x=repmat(lat_line,gd.N,1);
%         x_label='Latitude(^oN)';
%         temp_102=zeros(gd.N,length(lat_line)); %% initialize temp matrix that size is same with x
%         salt_102=zeros(gd.N,length(lat_line));
%         u_102=zeros(gd.N,length(lat_line));
%         v_102=zeros(gd.N,length(lat_line));
%         rho_102=zeros(gd.N,length(lat_line));
%     end
            
            
%     if (x1(1)==x2(1) || y1(1)==y2(1)) %% fixed lon or lat
        temp_102(:,:) = squeeze(temp(:,min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1))));
        salt_102(:,:) = squeeze(salt(:,min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1))));
        u_102(:,:) = squeeze(u(:,min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1))));
        v_102(:,:) = squeeze(v(:,min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1))));
        rho_102(:,:) = squeeze(rho(:,min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1))));
        
        Yi(:,:)= squeeze(depth(:,min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1))));
%     else
%         for k=1:1:gd.N
%             lon_range=lon_rho(min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1)));
%             lat_range=lat_rho(min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1)));
%             data_range=squeeze(data(k,min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1))));  %%data(zlevel, xmin : xmax, ymin : ymax)
%             depth_range=squeeze(depth(k,min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1)))); %%depth(zlevel, xmin : xmax, ymin : ymax)
%             Temp(k,:)=griddata(lon_range,lat_range,data_range,lon_line,lat_line); %%get interpolated data section
%             Yi(k,:)=griddata(lon_range,lat_range,depth_range,lon_line,lat_line); %%get interpolated depth section
%         end
%     end
            
%     data=Temp;
%             left bottom width height

    for var=1:5
        figure('position',[400 200 1500 550],'PaperUnits','inches','PaperPosition',[0 0 9.5 5]); %%   %%figure window, figure file 
        set(gca,'Position',[0.2 0.25 0.65 0.5]);  %% figure
        hold on
        
%         if (lon2-lon1) >= (lat2-lat1)        
            axis([section(1) section(2) section(5) section(6)]);
%         else
%             axis([section(3) section(4) section(5) section(6)]);
%         end
        if (var ==1)
            pcolor(x,Yi,temp_102)
            titlename = strcat('Temp (',num2str(tempyear),'th year,',' ',char(calendar(tempmonth)), ')');
            outfile =[figdir,'\EJS\EKWC\model_102_',num2str(abs(section(5))),'m_temp_']; 
            caxis([-2 33])
            level_c =5:5:30;
            [C,h]=contour(x,Yi,temp_102,level_c,'k','linewidth',1); 
            bar = colorbar('fontsize',17,'fontweight','bold');
            set(get(bar,'title'),'string',temp_unit,'FontSize',17,'fontweight','bold')
        elseif (var ==2)
            pcolor(x,Yi,salt_102)
            titlename = strcat('Salt (',num2str(tempyear),'th year,',' ',char(calendar(tempmonth)), ')');
            outfile =[figdir,'\EJS\EKWC\model_102_',num2str(abs(section(5))),'m_salt_']; 
            caxis([33.95 34.15])
            level_c =31:0.5:35;
            [C,h]=contour(x,Yi,salt_102,level_c,'k','linewidth',1); 
            bar = colorbar('fontsize',17,'fontweight','bold');
            set(get(bar,'title'),'string',salt_unit,'FontSize',17,'fontweight','bold')
        elseif (var ==3)
            pcolor(x,Yi,u_102)
            titlename = strcat('U (',num2str(tempyear),'th year,',' ',char(calendar(tempmonth)), ')');
            outfile =[figdir,'\EJS\EKWC\model_102_',num2str(abs(section(5))),'m_u_']; 
            caxis([-0.5 0.5])
            level_c =-0.5:0.1:0.5;
            [C,h]=contour(x,Yi,u_102,level_c,'k','linewidth',1); 
            bar = colorbar('fontsize',17,'fontweight','bold');
            set(get(bar,'title'),'string',u_unit,'FontSize',17,'fontweight','bold')
        elseif (var ==4)
            pcolor(x,Yi,v_102)
            titlename = strcat('V (',num2str(tempyear),'th year,',' ',char(calendar(tempmonth)), ')');
            outfile =[figdir,'\EJS\EKWC\model_102_',num2str(abs(section(5))),'m_v_']; 
            caxis([-0.5 0.5])
            level_c =-0.5:0.1:0.5;
            [C,h]=contour(x,Yi,v_102,level_c,'k','linewidth',1); 
            bar = colorbar('fontsize',17,'fontweight','bold');
            set(get(bar,'title'),'string',v_unit,'FontSize',17,'fontweight','bold')
        elseif (var ==5)
            pcolor(x,Yi,rho_102)
            titlename = strcat('Rho (',num2str(tempyear),'th year,',' ',char(calendar(tempmonth)), ')');
            outfile =[figdir,'\EJS\EKWC\model_102_',num2str(abs(section(5))),'m_rho_']; 
            caxis([0 35])
            level_c =20:1:35;
            [C,h]=contour(x,Yi,rho_102,level_c,'k','linewidth',1); 
            bar = colorbar('fontsize',17,'fontweight','bold');
            set(get(bar,'title'),'string',rho_unit,'FontSize',17,'fontweight','bold')
        end
        shading flat;
        
        set(gca,'box','on','linewidth',1.5,'fontsize',17)
        xlabel(x_label,'color','k','FontSize',17,'fontweight','bold')
        ylabel('Depth(m)','color','k','FontSize',17,'fontweight','bold')
        title(titlename,'fontsize',17); % stlee
        hold on

        
        clabel(C,h,'FontSize',15,'Color','k','labelspacing',100000,'Rotation',0,'fontweight','bold');
        
        load C:\Users\KYY\Dropbox\source\matlab\Common\Figure\jet_mod
        colormap(jet_mod);
    %             colormap(jet);
        
        
        saveas(gcf,[outfile,mid,'.tif'],'tiff');
        close all
    end
end

status = 1;
% end