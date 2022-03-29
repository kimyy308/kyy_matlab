%--------------------------------------------------------------------------
%  using this code when model grid is tilted.
%
%                        Plotting ROMS OUT ( 365daily )
%
%                                             date : 2019. 05. 09
%                                             edited by Y.K.
%
%--------------------------------------------------------------------------
clear all;close all;
%==========================================================================
st_line  = 132;  %latitude of stations line
max_level= 40;
Vtransform=2;
Vstretching=4;

modelname = 'new_deep_eastsea';
testname = 'v08';
domaxis=[35 43 -1000 0]; % [min_lat max_lat max_depth min_depth]

variation =1;  %  1 = temperature , 2 = salinity ,3 = density
current   = 1;   skip_v = 3;     size_v = 2;    color_v = 'k';
plot_pcolorjw = 1;  temp_lim = [0 15];    salt_lim = [33 34.5];
plot_contour  = 1;  color_c  ='-k' ;  level_temp =[2:1:20]; level_c  = [33:0.1:34.5];
switch_save   = 1;    out_type = 'jpg';

for yy = 1992: 1992
    year = num2str(yy);
    start_mm=1;
    end_mm=12;
    time_step=1;
    
    file_dir=['H:\',modelname,'\output\',testname,'\spinup\',year,'\'];
    file_dir_out=[file_dir,'fig\'];
    if exist(file_dir_out) == 0
        mkdir(file_dir_out)
    end
    
    mm = start_mm;
    
    gd = grd([modelname]);
    lon_rho  = gd.lon_rho;
    lat_rho  = gd.lat_rho;
    mask_rho = gd.mask_rho;
    h=gd.h;
    N = gd.N;
    dist = abs(st_line-lon_rho);
    min_dist = min(dist,[],2);
    [L,M] = size(lat_rho);
    
    for ll = 1 : L
        num_stline(ll,1) = find(dist(ll,:) == min_dist(ll));
        find_lat(ll,1) = lat_rho(ll,num_stline(ll));
        find_lon(ll,1) = lon_rho(ll,num_stline(ll));
    end
    
    idx_nan = abs(st_line-find_lon);
    find_nan = find(idx_nan > 1);
    
    depth = zlevs(h,0,gd.theta_s,gd.theta_b,gd.hc,N,Vstretching,Vtransform,'r');
    warning off
    mask_rho = mask_rho./mask_rho;
    
    warning on
    vname = {'temp','salt'};%,'zeta','ubar','vbar','u','v','omega'};
    
    for im = start_mm : time_step : end_mm
        mid=[num2char(im,2)];
        file = [file_dir,testname,'_monthly_',year,'_',mid,'.nc'];
        disp([file,' : ', num2char(im,4)])
        nc=netcdf(file);
        date=[num2str(yy),'. ',num2str(im)];
        
        switch variation
            case 1
                value = nc{char(vname(variation))}(:);
                val_name = 'Temperature';
                unit = '^oC';
                out_name_1 = [testname,'_monthly_temp_',num2str(yy),mid];
                val_caxis = temp_lim;
                
                for kk = 1:L
                    data(:,kk) = squeeze(value(:,kk,num_stline(kk)));
                    Yi(:,kk) = squeeze(depth(:,kk,num_stline(kk)));
                end
                level_c=level_temp;
                clear value;
                
            case 2
                value=nc{char(vname(variation))}(:);
                val_name='Salinity';
                unit = 'psu';
                out_name_1 = [testname,'_monthly_salt_',num2str(yy),mid];
                val_caxis=salt_lim;
                
                for kk = 1:L
                    data(:,kk) = squeeze(value(:,kk,num_stline(kk)));
                    Yi(:,kk) = squeeze(depth(:,kk,num_stline(kk)));
                end
                clear value;
                
                if (plot_pcolorjw)
                    for i=1:1:length(data(:,1))
                        for j=1:1:length(data(1,:))
                            if data(i,j) > 10000
                                data(i,j) = NaN;
                            end
                        end
                    end
                end
        end
        
        data(:,find_nan) = NaN;
        x_1 = find(find_lat >= domaxis(1));
        x_2 = find(find_lat >= domaxis(2));
        x = find_lat(x_1(1):x_2(1));
        x = repmat(x',40,1);
        data = data(:,x_1(1):x_2(1));
        Yi = Yi(:,x_1(1):x_2(1));
        
        figure('position',[400 100 550 550],'PaperUnits','inches','PaperPosition',[0 0 7 5 ]);
        set(gca,'Position',[0.15 0.15 0.80 0.75]);
        text_posi_x=(domaxis(2)-domaxis(1))/20+domaxis(1);
        text_posi_y1=(domaxis(4)-domaxis(3))/20+domaxis(3);
        text_posi_y2=2*(domaxis(4)-domaxis(3))/20+domaxis(3);
        text_posi_y3=3*(domaxis(4)-domaxis(3))/20+domaxis(3);
        text_posi_y4=(domaxis(4)-domaxis(3))/20+domaxis(4);
        
        hold on
        pcolor(x,Yi,data)
        shading flat;caxis(val_caxis)
        colormap jet
        set(gca,'box','on','linewidth',1.5,'fontsize',17)
        axis([domaxis(1) domaxis(2) domaxis(3) domaxis(4)])
        lab_lat=['Longitude ',num2str(st_line),'(^oE)'];
        xlabel('Latitude(^oN)','color',color_v,'FontSize',17,'fontweight','bold')
        ylabel('Depth(m)','color',color_v,'FontSize',17,'fontweight','bold')
        text(text_posi_x,text_posi_y1,lab_lat,'color',color_v,'FontSize',17,'fontweight','bold')
        text(text_posi_x,text_posi_y2,val_name,'color',color_v,'FontSize',17,'fontweight','bold')
        text(text_posi_x,text_posi_y3,date,'color',color_v,'FontSize',17,'fontweight','bold')
        %                 plot([40.5 40.5],[0 -3500],'--r','linewidth',1.5);
        out_name_1=['NS_MICT_',out_name_1];
        
        if (plot_contour)
            hold on
            [C,h]=contour(x,Yi,data,level_c,color_c,'linewidth',1);
            [C2,h2]=contour(x,Yi,data,level_c-1,color_c,'linewidth',1);
            [C2,h2]=contour(x,Yi,data,[-1:2:1],'-w','linewidth',1);
            [C2,h2]=contour(x,Yi,data,[30:4.06:34.06],'-w','linewidth',2);
            clabel(C,h,'FontSize',15,'Color','k','labelspacing',100000,'fontweight','bold');
        end
        %             caxis(val_caxis);
        bar = colorbar('fontsize',17);
        set(get(bar,'title'),'string',unit,'FontSize',17,'fontweight','bold')
        
        out_name=[file_dir_out,out_name_1];
        
        if (switch_save)
            saveas(gcf,out_name,out_type);
        end
        close all
    end
    %      close all
    % % % end
end