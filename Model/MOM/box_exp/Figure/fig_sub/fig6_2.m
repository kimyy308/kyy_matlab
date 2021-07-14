[lon44,ilevel44] = meshgrid(ilevel2(1:23),lon2(1:30));

gifname = sprintf('Meridional Vector & Vertical Temperature along 129^oE, 10th year averaged value');
    sbp1=subplot(2,1,1);
    pcolor(lat2(31:50),ilevel2(1:20)/20,(reshape(u2mean(51,31:50,1:20,1),20,20))');
    hold on;
    shading interp;
    set(gca, 'color', [0.8,0.8,0.8]);
    caxis([-0.35,0.35]);
    c=colorbar;
    c.Label.String= 'Zonal velocity (m/s)';
    c.Label.FontSize= 12;
    colormap(bwrmap);
    [C,h]=contour(lat2(31:50),ilevel2(1:20)/20,(reshape(u2mean(51,31:50,1:20,1),20,20))','ShowText','on');
    h.LineColor='black';
    %     title(gifname,'fontsize',27);
    xlabel('Latitude (^oN)','fontsize',15);
    ylabel('Depth (M)','fontsize',15);
    set(gcf,'PaperPosition',[0 0 20 15]);  
    set(gca,'fontsize',12);
    
    ar=annotation(gcf,'rectangle',...
    [0.308958333333333 0.629692440844922 0.410979166666667 0.293117647058822],...
    'Color',[1 0 0],...
    'LineWidth',1.5);
    tex=text(37.4,-1.5,'A');
    set(tex,'fontsize',23);
    set(tex,'Color',[0 0 0]);    
    hold off;
    
    set(gca,'yTickLabel',{-160:40:0})

 gifname = sprintf('Meridional Vector & Vertical Temperature along 37^oN, 10th year averaged value');
    sbp2=subplot(2,1,2);
    pcolor(ilevel44,lon44,reshape(v2mean(1:30,60,1:23,1),30,23));
    hold on;
    shading interp;
    set(gca, 'color', [0.8,0.8,0.8]);
    caxis([-0.35,0.35]);
    c=colorbar;
    c.Label.String= 'Meridional velocity (m/s)';
    c.Label.FontSize= 12;
    colormap(bwrmap);
%     title(gifname,'fontsize',27);
    [C,h]=contour(ilevel44,lon44,reshape(v2mean(1:30,60,1:23,1),30,23),'ShowText','on');
    h.LineColor='black';
    xlabel('Longitude (^oE)','fontsize',15);
    ylabel('Depth (M)','fontsize',15);
    set(gcf,'PaperPosition',[0 0 20 15]);  
%     set(gcf,'PaperPositionMode','auto'); 
%     set(sbp2,'Position',[200 100 600 700]);  
    set(gca,'fontsize',12);
%    [C,h]=contour(ilevel4,lon4,reshape(v2(1:25,30,1:15,1),25,15),'ShowText','on');
%    set(gca,'YTick',[-300 -200 -100 0])
%    set(gca,'yTickLabel',{-300:100:0})
%     h.LineColor=[0.5 0.5 0.5];
%     clabel(C,h,'Color',[0.5 0.5 0.5]);
    br=annotation(gcf,'rectangle',...
    [0.200770833333333 0.229067930489731 0.295354166666667 0.220589257503949],...
    'Color',[1 0 0],...
    'LineWidth',1.5);    
    tex=text(128.15,-40,'B');
    set(tex,'fontsize',23);
    set(tex,'Color',[0 0 0]);    
    hold off;
    fig=gcf;
    fig.InvertHardcopy='off';
    saveas(gcf,filename,'tiff');
    delete(sbp1);
    delete(sbp2);
    delete(ar);
    delete(br); 