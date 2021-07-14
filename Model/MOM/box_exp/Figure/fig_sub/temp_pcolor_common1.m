if (pcolindex==1)
    hold on;
    shading interp;
    set(gca, 'color', [0.8,0.8,0.8]);
    caxis([0.2,20]);
%     c=colorbar('location','north');
%      c=colorbar
%     c.Label.String= 'Temp (^oC)';
%     c.Label.FontSize= 12;
    colormap jet;
    set(gcf,'PaperPosition',[0 0 20 15]);
    axis equal;
    set(gca,'fontsize',13);
else
    hold on;
    shading interp;
    set(gca, 'color', [0.8,0.8,0.8]);
    caxis([0.2,20]);
    c=colorbar;
    c.Label.String= 'Temp (^oC)';
    c.Label.FontSize= 12;
    colormap jet;
    xlabel('Longitude (^oE)','fontsize',27);
    ylabel('Latitude (^oN)','fontsize',27);
    set(gcf,'PaperPosition',[0 0 20 15]);
    axis equal;
    set(gca,'fontsize',15);
end