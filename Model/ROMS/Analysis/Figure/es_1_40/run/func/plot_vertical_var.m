function status=plot_vertical_var(outfile, section, year, inputmonth, var, var_lim, inputdir)
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

% vert_param;
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

% file_dir='D:\MEPL\project\SSH\�߰�����\smooth13_vtvs\';
% mm=start_mm;
%==========================================================================

startmonth=(year(1)-1992)*12+inputmonth(1)
endmonth=(year(2)-1992)*12+inputmonth(2)
mid=num2str(startmonth,'%04i');
file = [inputdir,'monthly_spinup_',mid,'.nc'];

% gd = read_grid(grdfile,Vtransform,Vstretching,theta_s,theta_b,hc,N);
gd = read_grid(file,Vtransform,Vstretching,theta_s,theta_b,hc,N);
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
warning off
% mask_rho = mask_rho./mask_rho;
% mask_u = mask_u./mask_u;
% mask_v = mask_v./mask_v;
% warning on
% vname = {'temp','salt'};%,'zeta','ubar','vbar','u','v','omega'};
% vname = {'temp','salt', 'rho','v'};%,'zeta','ubar','vbar','u','v','omega'};
vname = {'temp','salt', 'u', 'v', 'rho'};
% calendar{1} = ' January'; calendar{2} = ' February'; calendar{3} = ' March'; 
% calendar{4} = ' April'; calendar{5} = ' May'; calendar{6} = ' June';
% calendar{7} = ' July'; calendar{8} = ' August'; calendar{9} = ' September'; 
% calendar{10} = ' October'; calendar{11} = ' November'; calendar{12} = ' December';

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
        
        switch var %%read data, set name and unit
            case 1
                value=nc{char(vname(var))}(:);
                value(find(value<-1000))=NaN; 
                val_name='Temperature';
                unit = '^oC';
%                     out_name_1=['vertical_NS_Temp-',num2str(yy),'-'];
%                     val_caxis=temp_lim;
%                     level_c=temp_c;
                data=value;
               clear value;         
            case 2
                value=nc{char(vname(var))}(:);
                val_name='Salinity';
%                 unit = 'psu';
                unit = ' ';
%                     out_name_1=['vertical_NS_Salt-',num2str(yy),'-'];
%                     val_caxis=salt_lim;
%                     level_c=salt_c;
                data=value;
               clear value;     
           case 3
                value=nc{char(vname(var))}(:);
                value=value*100;
                val_name='u';
                unit = 'cm/s';
%                     out_name_1=['vertical_meridional_v-',num2str(yy),'-'];
%                     val_caxis=v_lim;
%                     level_c=v_c;
                data=value;
               clear value;
           case 4
                value=nc{char(vname(var))}(:);
                value=value*100;
                val_name='v';
                unit = 'cm/s';
%                     out_name_1=['vertical_meridional_v-',num2str(yy),'-'];
%                     val_caxis=v_lim;
%                     level_c=v_c;
                data=value;
               clear value;
        end
   

%         if (plot_pcolorjw)
%             for i=1:1:length(data(:,1))
%                 for j=1:1:length(data(1,:))
%                     if data(i,j) > 10000
%                         data(i,j) = NaN;
%                     end
%                 end
%             end
%             switch section
%                 case 1
%                 domaxis=[38 39.5 120 124 -100 0]; % ����
%                 domaxis=[35.28 33.42 128.88 130.09 -100 0]; % �������� along_stlee  
%                 domaxis=[35.28 33.42 128.88 130.09 -180 0]; % �������� (NWP����) along_stlee  
%                 domaxis=[41 42 140.8 141 -50 0]; % ������ stlee
%                 domaxis=[41 42 140.5 140.7 -190 0]; % ������ stlee
%                   domaxis=[45.1872096987 46.4790236999 142.043821201 142.151473264  -100 0]; % �Ҿ� stlee
%                   domaxis=[36 36.05 121 126  -90 0]; % Ȳ�� kimyy
%                     domaxis=[33 36.5 124.4 124.45  -100 0]; % tak_124,4 kimyy
%                        domaxis=[35 35.05 119 127  -100 0]; % tak_124,4 kimyy
%                 domaxis=[34.85 33.2265 128.6 129.5 -100 0]; % �������� along_stlee2
%                 domaxis=[128.88 35.28 130.09 33.42] % �������� along_stlee  
%                 domaxis=[30 36 122 125 -200 0]; % ����
%                 domaxis=[23 30 118.3 124 -100 0]; % �븸���� along
%                 domaxis=[25 24 118.5 121 -100 0]; % �븸���� cross
            dist=sqrt((lon_rho-section(1)).^2+(lat_rho-section(3)).^2); %% get distance from station 1
            min_dist=min(min(dist));
            dist2=sqrt((lon_rho-section(2)).^2+(lat_rho-section(4)).^2);  %% get distance from station 2
            min_dist2=min(min(dist2));                
            [x1,y1]=find(dist==min_dist);  %% find closest point from station 1. [lat lon]
            [x2,y2]=find(dist2==min_dist2); %% find closest point from station 2. [lat lon]
            
            lat1=lat_rho(x1(1),y1(1));  lon1=lon_rho(x1(1),y1(1));
            lat2=lat_rho(x2(1),y2(1));  lon2=lon_rho(x2(1),y2(1));
            if (lon2-lon1) >= (lat2-lat1)
                lon_line = lon1:mean(gradient(lon_rho(1,:))):lon2;  %% for 1/20^o horizontal resolution
                lat_line = (lon_line-lon1)/((lon2-lon1)/(lat2-lat1))+lat1; %% weight for latitude index (lat1 : 0.05/((lon2-lon1) * (lat2-lat1) : lat2 
                x=repmat(lon_line,gd.N,1);  %% copy lon_line (gd.N times) to make matrix 
                x_label='Longitude(^oE)';
%                 domaxis=[domaxis(3) domaxis(4) domaxis(5) domaxis(6) domaxis(5) domaxis(6)];
                Temp=zeros(gd.N,length(lon_line)); %% initialize temp matrix that size is same with x
            else
                lat_line=[min(lat1,lat2):mean(gradient(lat_rho(:,1))):max(lat1,lat2)];
                lon_line=(lat_line-lat1)*((lon2-lon1)/(lat2-lat1))+lon1;
                x=repmat(lat_line,gd.N,1);
                x_label='Latitude(^oN)';
                Temp=zeros(gd.N,length(lat_line)); %% initialize temp matrix that size is same with x
            end
            
            
            if (x1(1)==x2(1) || y1(1)==y2(1)) %% fixed lon or lat
                Temp(:,:) = squeeze(data(:,min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1))));
                Yi(:,:)= squeeze(depth(:,min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1))));
            else
                for k=1:1:gd.N
                    lon_range=lon_rho(min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1)));
                    lat_range=lat_rho(min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1)));
                    data_range=squeeze(data(k,min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1))));  %%data(zlevel, xmin : xmax, ymin : ymax)
                    depth_range=squeeze(depth(k,min(x1(1),x2(1)):max(x1(1),x2(1)),min(y1(1),y2(1)):max(y1(1),y2(1)))); %%depth(zlevel, xmin : xmax, ymin : ymax)
                    Temp(k,:)=griddata(lon_range,lat_range,data_range,lon_line,lat_line); %%get interpolated data section
                    Yi(k,:)=griddata(lon_range,lat_range,depth_range,lon_line,lat_line); %%get interpolated depth section
                end
            end
            
            data=Temp;
%             end
%         end    
%             left bottom width height
        figure('position',[400 200 1500 550],'PaperUnits','inches','PaperPosition',[0 0 9.5 5]); %%   %%figure window, figure file 
        set(gca,'Position',[0.2 0.25 0.65 0.5]);  %% figure
% %         text_posi_x=(domaxis(2)-domaxis(1))/20+domaxis(1);
% %         text_posi_y1=(domaxis(6)-domaxis(5))/20+domaxis(5);
% %         text_posi_y2=2*(domaxis(6)-domaxis(5))/20+domaxis(5);
% %         text_posi_y3=3*(domaxis(6)-domaxis(5))/20+domaxis(5);
        text_posi_x=(section(4)-section(3))/20+section(3);
        text_posi_y1=(section(6)-section(5))/20+section(5);
        text_posi_y2=2*(section(6)-section(5))/20+section(5);
        text_posi_y3=3*(section(6)-section(5))/20+section(5);
%             switch section
%                 case 1
                hold on
                pcolor(x,Yi,data)
                if (lon2-lon1) >= (lat2-lat1)        
                    axis([section(1) section(2) section(5) section(6)]);
                else
                    axis([section(3) section(4) section(5) section(6)]);
                end
                shading flat;
%                 caxis(val_caxis)
                caxis(var_lim)
                set(gca,'box','on','linewidth',1.5,'fontsize',17)
                xlabel(x_label,'color','k','FontSize',17,'fontweight','bold')
                ylabel('Depth(m)','color','k','FontSize',17,'fontweight','bold')

% %                 text(text_posi_x,text_posi_y2,val_name,'color','k','FontSize',17,'fontweight','bold')
% %                 text(text_posi_x,text_posi_y3,date,'color','k','FontSize',17,'fontweight','bold')


            if (var ==1)
                titlename = strcat('Temp (',num2str(tempyear),'th year,',' ',char(calendar(tempmonth)), ')');
            elseif (var ==2)
                titlename = strcat('Salt (',num2str(tempyear),'th year,',' ',char(calendar(tempmonth)), ')');
            elseif (var ==3)
                titlename = strcat('U (',num2str(tempyear),'th year,',' ',char(calendar(tempmonth)), ')');
            elseif (var ==4)
                titlename = strcat('V (',num2str(tempyear),'th year,',' ',char(calendar(tempmonth)), ')');
            end
%                 title('Vertical Temperature','fontsize',17); % stlee
                title(titlename,'fontsize',17); % stlee
%                 out_name_1=['YellowSea',out_name_1];
%                 if (plot_contour)
                  hold on
                  level_c =5:5:20;
%                   level_c =31:0.5:35;
                  [C,h]=contour(x,Yi,data,level_c,'k','linewidth',1);
% %                   [C2,h2]=contour(x,Yi,data,level_c-0.5,color_c,'linewidth',1);
% %                   [C2,h2]=contour(x,Yi,data,[-1:2:1],'-w','linewidth',1);                  
                  clabel(C,h,'FontSize',15,'Color','k','labelspacing',100000,'Rotation',0,'fontweight','bold');
%                 end
%             end
%             caxis(val_caxis);
            bar = colorbar('fontsize',17,'fontweight','bold');
%             colormap(jet);
            load C:\Users\KYY\Dropbox\source\matlab\Common\Figure\jet_mod
%             colormap(jet_mod);
            colormap(jet);
            set(get(bar,'title'),'string',unit,'FontSize',17,'fontweight','bold')
            
%             out_name=['09_13YellowSea_salt_',num2str(mid,'%02i')];
%             out_name=['YellowSea_mer34_v_tak_',num2str(mid,'%04i')];

%             if (switch_save)
%                 saveas(gcf,out_name,out_type);
                  saveas(gcf,[outfile,mid,'.tif'],'tiff');
%             end
            close all
end
%      close all
% % % end
status = 1;
end