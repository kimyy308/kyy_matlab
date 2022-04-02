%--------------------------------------------------------------------------
%  using this code when model grid is tilted.
%
%                        Plotting ROMS EW section
%
%                                             date : 2019. 05. 09
%                                             edited by Y.K.
%
%--------------------------------------------------------------------------
clear all;close all;
%==========================================================================
st_line  = 37.5;  %latitude of stations line
max_level= 40;
Vtransform=2;
Vstretching=4;

modelname = 'NWP';
test = [16, 19];
syear = 1993; eyear = 2018;

for tt = 1 : 2
    testname = ['test',num2str(test(tt))];
    progress = 'run';
    
    domaxis=[128 139 -1000 0]; % [min_lon max_lon max_depth min_depth]
    % domaxis=[129.5 132.5 -500 0]; %KODC 105 line
    
    variation =1;  %  1 = temperature , 2 = salinity ,3 = density
    current   = 1;   skip_v = 3;     size_v = 2;    color_v = 'k';
    plot_pcolorjw = 1;  temp_lim = [0 15];    salt_lim = [33 34.5];
    plot_contour  = 1;  color_c  ='-k' ;  level_temp =[2:1:15]; level_c  = [33:0.1:34.5];
    switch_save   = 1;    out_type = 'jpg';
    
    for year = syear : eyear
        yy = num2str(year);
        
        smonth=1;
        emonth=12;
        
        filepath=['H:\new_deep_eastsea\',modelname,'\',testname,'\output\',progress,'\',yy,'\'];
        fname = [testname,'_monthly_',yy,'_01.nc'];
        nw = netcdf([filepath, fname]);
        lon_rho  = nw{'lon_rho'}(:);
        lat_rho  = nw{'lat_rho'}(:);
        mask_rho = nw{'mask_rho'}(:);
        h=nw{'h'}(:); hc = nw{'hc'}(:);
        temp = nw{'temp'}(:);
        theta_s = nw{'theta_s'}(:); theta_b = nw{'theta_b'}(:);
        dist = abs(st_line-lat_rho);
        min_dist = min(dist);
        [N,L,M] = size(temp);
        
        depth = zlevs(h,0,theta_s,theta_b,hc,N,Vstretching,Vtransform,'r');
        warning off
        mask_rho = mask_rho./mask_rho;
        close(nw);
        
        outpath=[filepath,'fig\vertical\'];
        if exist(outpath) == 0
            mkdir(outpath)
        end
        
        for ll = 1 : M
            num_stline(ll) = find(dist(:,ll) == min_dist(ll));
            find_lat(ll) = lat_rho(num_stline(ll),ll);
            find_lon(ll) = lon_rho(num_stline(ll),ll);
        end
        
        idx_nan = abs(st_line - find_lat);
        find_nan = find(idx_nan > 1);
        
        warning on
        vname = {'temp','salt'};%,'zeta','ubar','vbar','u','v','omega'};
        
        for month = smonth : emonth
            mm = num2char(month,2);
            file = [filepath,testname,'_monthly_',yy,'_',mm,'.nc'];
            disp([file,' : ', mm])
            nc=netcdf(file);
            
            switch variation
                case 1
                    value = nc{char(vname(variation))}(:);
                    val_name = 'Temperature';
                    unit = '^oC';
                    var = 'temp';
                    val_caxis = temp_lim;
                    
                    for kk = 1:M
                        data(:,kk) = squeeze(value(:,num_stline(kk),kk));
                        Yi(:,kk) = squeeze(depth(:,num_stline(kk),kk));
                    end
                    level_c=level_temp;
                    clear value;
                    
                case 2
                    value=nc{char(vname(variation))}(:);
                    val_name='Salinity';
                    unit = 'psu';
                    var = 'salt';
                    val_caxis=salt_lim;
                    
                    for kk = 1:M
                        data(:,kk) = squeeze(value(:,num_stline(kk),kk));
                        Yi(:,kk) = squeeze(depth(:,num_stline(kk),kk));
                    end
                    clear value;
            end
            
            if (plot_pcolorjw)
                for i=1:1:length(data(:,1))
                    for j=1:1:length(data(1,:))
                        if data(i,j) > 10000
                            data(i,j) = NaN;
                        end
                    end
                end
            end
            
            data(:,find_nan) = NaN;
            x_1 = find(find_lon >= domaxis(1));
            x_2 = find(find_lon >= domaxis(2));
            x = find_lon(x_1(1):x_2(1));
            x = repmat(x,40,1);
            data = data(:,x_1(1):x_2(1));
            Yi = Yi(:,x_1(1):x_2(1));
            
            text_posi_x=(domaxis(2)-domaxis(1))/20+domaxis(1);
            text_posi_y1=(domaxis(4)-domaxis(3))/20+domaxis(3);
            text_posi_y2=2*(domaxis(4)-domaxis(3))/20+domaxis(3);
            text_posi_y3=3*(domaxis(4)-domaxis(3))/20+domaxis(3);
            text_posi_y4=4*(domaxis(4)-domaxis(3))/20+domaxis(3);
            data(data>1000) = NaN;
            
            figure('position',[400 100 550 550],'PaperUnits','inches','PaperPosition',[0 0 7 5 ]);
            set(gca,'Position',[0.15 0.15 0.80 0.75]); hold on;
            pcolor(x,Yi,data); shading flat;caxis(val_caxis); colormap jet;
            set(gca,'box','on','linewidth',1.5,'fontsize',17);
            axis([domaxis(1) domaxis(2) domaxis(3) domaxis(4)]);
            
            lab_lat=['Latitude',num2str(st_line),'(^oN)'];
            xlabel('Longitude(^oE)','color',color_v,'FontSize',17,'fontweight','bold')
            ylabel('Depth(m)','color',color_v,'FontSize',17,'fontweight','bold')
            text(text_posi_x,text_posi_y1,lab_lat,'color',color_v,'FontSize',17,'fontweight','bold')
            text(text_posi_x,text_posi_y2,val_name,'color',color_v,'FontSize',17,'fontweight','bold')
            switch progress
                case 'spinup'
                    date=[num2str(syear),'. ',mm];
                    text(text_posi_x,text_posi_y3,[date,'(spinup ',num2str(spin),'yrs)'],'color',color_v,'FontSize',17,'fontweight','bold')
                case 'run'
                    date=[yy,'. ',mm];
                    text(text_posi_x,text_posi_y3,[date],'color',color_v,'FontSize',17,'fontweight','bold')
            end
            %             text(text_posi_x,text_posi_y4,testname,'color',color_v,'FontSize',17,'fontweight','bold')
            outname=['EW_',testname,'_monthly_',var,'_',yy,mm];
            
            if (plot_contour)
                hold on
                [C,h]=contour(x,Yi,data,level_c,color_c,'linewidth',1);
                [C2,h2]=contour(x,Yi,data,level_c-1,color_c,'linewidth',1);
                [C2,h2]=contour(x,Yi,data,[-1:2:1],'-w','linewidth',1);
                [C2,h2]=contour(x,Yi,data,[30:4.06:34.06],'-w','linewidth',2);
                clabel(C,h,'FontSize',15,'Color','k','labelspacing',100000,'fontweight','bold');
            end
            
            bar = colorbar('fontsize',17);
            set(get(bar,'title'),'string',unit,'FontSize',17,'fontweight','bold')
            
            if (switch_save)
                saveas(gcf,[outpath, outname],out_type);
            end
            close all
        end
    end
end
