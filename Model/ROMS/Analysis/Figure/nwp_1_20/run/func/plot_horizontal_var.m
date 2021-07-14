function status=plot_horizontal_var(outfile, section, year, inputmonth, var, var_lim, inputdir)
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

vert_param;

time_step=1;

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

    dist=sqrt((lon_rho-section(1)).^2+(lat_rho-section(3)).^2); %% get distance from station 1
    min_dist=min(min(dist));
    dist2=sqrt((lon_rho-section(2)).^2+(lat_rho-section(4)).^2);  %% get distance from station 2
    min_dist2=min(min(dist2));                
    [x1,y1]=find(dist==min_dist);  %% find closest point from station 1. [lat lon]
    [x2,y2]=find(dist2==min_dist2); %% find closest point from station 2. [lat lon]

    cut_data=data(:,x1(1):x2(1),y1(1):y2(1));
    cut_depth=depth(:,x1(1):x2(1),y1(1):y2(1));
    cut_lon=lon_rho(x1(1):x2(1),y1(1):y2(1));
    cut_lat=lat_rho(x1(1):x2(1),y1(1):y2(1));
    cut_mask_rho=mask_rho(x1(1):x2(1),y1(1):y2(1));
    refdepth=(section(5)+section(6))/2.0
    vert_dist=abs(cut_depth-refdepth);
   
    for i=1:length(vert_dist(1,:,1))
%         tic;
        for j=1:length(vert_dist(1,1,:))
            if (cut_mask_rho(i,j)==0)
                cut_data3(i,j)=NaN;
            else
                if (cut_depth(1,i,j)>refdepth)
                    cut_data3(i,j)=NaN;
                else
                    z_ind=find(squeeze(vert_dist(:,i,j))==min(squeeze(vert_dist(:,i,j))));
                    if (z_ind==1)
                        cut_depth2(1:3)=cut_depth(z_ind:z_ind+2,i,j);
                        cut_data2(1:3)=cut_data(z_ind:z_ind+2,i,j);
                    else
                        cut_depth2(1:3)=cut_depth(z_ind-1:z_ind+1,i,j);
                        cut_data2(1:3)=cut_data(z_ind-1:z_ind+1,i,j);
                    end
                    cut_data3(i,j)=interpn(cut_depth2,cut_data2,refdepth,'cubic');
                end
            end
        end
%        toc;
    end


    % make jpg file
    xscale=section(2)-section(1);
    yscale=section(4)-section(3);
    halt = 1;
    while(halt)
        if (xscale > 1000 || yscale > 1000)
            halt = 0;
        else
            xscale = xscale * 1.2; yscale = yscale * 1.2;
        end
    end
%         xscale = 800; yscale = 920; %% temporary scale
    set(gcf,'Position',[200 100 xscale yscale])
    
    hold on
%     pcolor(x,Yi,data)
    m_proj('mercator','lon',[section(1) section(2)],'lat',[section(3) section(4)]);
    m_grid('fontsize',17, 'box', 'fancy');
    hold on;
    m_pcolor(cut_lon,cut_lat,cut_data3)
%     if (lon2-lon1) >= (lat2-lat1)        
%         axis([section(1) section(2) section(5) section(6)]);
%     else
%         axis([section(3) section(4) section(5) section(6)]);
%     end
    shading interp;

    caxis(var_lim)
    set(gca,'box','on','linewidth',1.5,'fontsize',17)
    xlabel('Longitude(^oE)','color','k','FontSize',17,'fontweight','bold')
    ylabel('Latitude(^oN)','color','k','FontSize',17,'fontweight','bold')

    if (var ==1)
        titlename = strcat(num2str(abs(refdepth)),'m Temp (',num2str(tempyear),'th year,',' ',char(calendar(tempmonth)), ')');
    elseif (var ==2)
        titlename = strcat('Salt (',num2str(tempyear),'th year,',' ',char(calendar(tempmonth)), ')');
    elseif (var ==3)
        titlename = strcat('U (',num2str(tempyear),'th year,',' ',char(calendar(tempmonth)), ')');
    elseif (var ==4)
        titlename = strcat('V (',num2str(tempyear),'th year,',' ',char(calendar(tempmonth)), ')');
    end
    title(titlename,'fontsize',17); 
      
    level_c =4:3:34;
    [C,h]=m_contour(cut_lon,cut_lat,cut_data3,level_c,'k','linewidth',1);             
    clabel(C,h,'FontSize',15,'Color','k','labelspacing',100000,'Rotation',0,'fontweight','bold');
    m_gshhs_i('color','k')  
    m_gshhs_i('patch',[.8 .8 .8]);   % gray colored land (need much time)
    bar = colorbar('fontsize',17,'fontweight','bold');

    load C:\Users\KYY\Dropbox\source\matlab\Common\Figure\jet_mod
    colormap(jet_mod);
    set(get(bar,'title'),'string',unit,'FontSize',17,'fontweight','bold')

    saveas(gcf,[outfile,mid,'.tif'],'tiff');
    close all
end
status = 1;
end