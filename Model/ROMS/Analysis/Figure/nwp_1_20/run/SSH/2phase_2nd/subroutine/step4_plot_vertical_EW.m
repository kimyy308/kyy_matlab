%--------------------------------------------------------------------------
%  using this code when z coordinate model
%
%                        Plotting ROMS OUT ( Climate monthly )
%
%                                             date : 2021. 05. 21
%                                             edited by Y.K.Kang
%
%--------------------------------------------------------------------------
clear all;close all;
modelname = 'CMIP6';
testnames = {'test2102','test2103','test2104','test2105','test2106'};
% testnames = {'test2101'};
st_line  = 37.5  %latitude of stations line
fig_salt = 0;

syear = 2014; ss = num2str(syear);
domaxis=[128 136 -1000 0]; % [min_lat max_lat max_depth min_depth]

polypath = 'H:\kisti_paper\data\data_mat\';
polyname = 'location_48ES.mat';
load([polypath, polyname]);

gridpath = ['I:\CMIP6\input\'];
gridname = ['ES_grid.nc'];
nc = netcdf([gridpath, gridname]);
lon_rho = nc{'lon_rho'}(:); lon = lon_rho(1,:);
lat_rho = nc{'lat_rho'}(:); lat = lat_rho(:,1);
mask_rho = nc{'mask_rho'}(:);
angle = nc{'angle'}(:);
h = nc{'h'}(:);
close(nc);

gpath = gridpath;
gname = ['ES_grid_info.nc'];
nc = netcdf([gpath,gname]);
hc = nc{'hc'}(:);
Vstretching = nc{'Vstretching'}(:);
Vtransform = nc{'Vtransform'}(:);
theta_s = nc{'theta_s'}(:);
theta_b = nc{'theta_b'}(:);
close(nc);

for i = 1 : length(testnames)
    testname = testnames{i};
    
    fpath = ['I:\',modelname,'\output\',testname,'\monthly\'];
    flist = dir(fpath);
    ee = flist(end).name; eyear = str2num(ee);
    
    dist = abs(lat - st_line);
    min_dist = min(dist(:));
    idx_lon = find(lon>=domaxis(1) & lon<=domaxis(2));
    idx_lat = find(dist == min_dist);
    fig_lon = lon(idx_lon);
    warning off
    
    for year = syear : eyear
        yy = num2str(year);
        filepath = [fpath, yy,'\stddepth\'];
        
        out_path = ['I:\',modelname,'\output\',testname,'\monthly\',yy,'\fig\vertical\EW\'];
        if exist(out_path) == 0
            mkdir(out_path)
        end
        
        for month  = 1 : 12
            mm = num2char(month,2);
            disp([testname,':',yy,'.',mm])
            
            file = [filepath,testname,'_monthly_std_',yy,'_',mm];
            load(file)
            [N,L,M] = size(temp);
            
            mask_Wh = inpolygon(lon_rho,lat_rho,polygon_Wh(:,1),polygon_Wh(:,2));
            
            for i = 1 : N
                c_temp(i,:,:) = squeeze(temp(i,:,:)) .*mask_Wh./mask_Wh;
                c_salt(i,:,:) = squeeze(salt(i,:,:)) .*mask_Wh./mask_Wh;
            end
            
            m_temp = squeeze(c_temp(:,idx_lat,idx_lon));
            m_salt = squeeze(c_salt(:,idx_lat,idx_lon));
            
            % fig temperature
            figure('position',[400 100 550 550],'PaperUnits','inches','PaperPosition',[0 0 7 5 ],...
                'visible','off');
            set(gca,'Position',[0.15 0.15 0.80 0.75]);
            text_posi_x=(domaxis(2)-domaxis(1))/20+domaxis(1);
            text_posi_y1=(domaxis(4)-domaxis(3))/20+domaxis(3);
            text_posi_y2=2*(domaxis(4)-domaxis(3))/20+domaxis(3);
            text_posi_y3=3*(domaxis(4)-domaxis(3))/20+domaxis(3);
            text_posi_y4=(domaxis(4)-domaxis(3))/20+domaxis(4);
            hold on
            pcolor(fig_lon,depth,m_temp); shading flat;
            caxis([0 15]); colormap jet
            
            set(gca,'box','on','linewidth',1.5,'fontsize',15)
            axis([domaxis(1) domaxis(2) domaxis(3) domaxis(4)])
            lab_lat=['Latitude ',num2str(st_line),'(^oE)'];
            xlabel('Longitude(^oE)','color','k','FontSize',17,'fontweight','bold')
            ylabel('Depth(m)','color','k','FontSize',17,'fontweight','bold')
            text(text_posi_x,text_posi_y1,lab_lat,'color','k','FontSize',17,'fontweight','bold')
            text(text_posi_x,text_posi_y2,'Temperature (\circC)','color','k','FontSize',17,'fontweight','bold')
            text(text_posi_x,text_posi_y3,[yy,'.',mm],'color','k','FontSize',17,'fontweight','bold')
            out_name_1=[modelname,'_',testname,'_monthly_EW_temp_',yy,mm];
            
            [C,h]=contour(fig_lon,depth,m_temp,[1:1:30],'k','linewidth',1);
            [C2,h2]=contour(fig_lon,depth,m_temp,[-1 -1],'-w','linewidth',2);
%             [C3,h3]=contour(fig_lon,depth,m_temp,[0.6 0.6],'-w','linewidth',2);
            clabel(C,h,'FontSize',15,'Color','k','labelspacing',100000,'fontweight','bold');
            clabel(C2,h2,'FontSize',15,'Color','w','labelspacing',100000,'fontweight','bold');
%             clabel(C3,h3,'FontSize',15,'Color','w','labelspacing',100000,'fontweight','bold');
            
            bar = colorbar('fontsize',15);
            set(get(bar,'title'),'string','\circC','FontSize',15,'fontweight','bold')
            out_name=[out_path,out_name_1];
            saveas(gcf,out_name,'jpg');
            close all
            
            if fig_salt == 1
                % fig salinity
                figure('position',[400 100 550 550],'PaperUnits','inches','PaperPosition',[0 0 7 5 ],...
                    'visible','off');
                set(gca,'Position',[0.15 0.15 0.80 0.75]);
                text_posi_x=(domaxis(2)-domaxis(1))/20+domaxis(1);
                text_posi_y1=(domaxis(4)-domaxis(3))/20+domaxis(3);
                text_posi_y2=2*(domaxis(4)-domaxis(3))/20+domaxis(3);
                text_posi_y3=3*(domaxis(4)-domaxis(3))/20+domaxis(3);
                text_posi_y4=(domaxis(4)-domaxis(3))/20+domaxis(4);
                hold on
                pcolor(fig_lon,depth,m_salt); shading flat;
                caxis([34 34.3]); colormap jet
                
                set(gca,'box','on','linewidth',1.5,'fontsize',15)
                axis([domaxis(1) domaxis(2) domaxis(3) domaxis(4)])
                lab_lat=['Latitude ',num2str(st_line),'(^oE)'];
                xlabel('Longitude(^oE)','color','k','FontSize',17,'fontweight','bold')
                ylabel('Depth(m)','color','k','FontSize',17,'fontweight','bold')
                text(text_posi_x,text_posi_y1,lab_lat,'color','k','FontSize',17,'fontweight','bold')
                text(text_posi_x,text_posi_y2,'Salinity','color','k','FontSize',17,'fontweight','bold')
                text(text_posi_x,text_posi_y3,[yy,'.',mm],'color','k','FontSize',17,'fontweight','bold')
                out_name_2=[modelname,'_',testname,'_monthly_EW_salt_',yy,mm];
                
                [C,h]=contour(fig_lon,depth,m_salt,[30:0.1:40],'k','linewidth',1);
                [C2,h2]=contour(fig_lon,depth,m_salt,[34.06 34.06],'-w','linewidth',1);
                clabel(C,h,'FontSize',15,'Color','k','labelspacing',100000,'fontweight','bold');
                clabel(C2,h2,'FontSize',15,'Color','w','labelspacing',100000,'fontweight','bold');
                
                bar = colorbar('fontsize',15);
                out_name=[out_path,out_name_2];
                saveas(gcf,out_name,'jpg');
                close all
                
            end % fig salt
        end % month
    end %year
end % testname
