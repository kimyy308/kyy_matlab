%     hold on;
    shading interp;
    xlim([min(lat2(1:latlimind)) max(lat2(1:latlimind))]);
%     set(gca,'yTickLabel',{'-160m', '-120m', '-80m', '-40m', '0m'}, ...
%         'YTick', [-160, -120, -80, -40, 0])
    set(gca,'yTickLabel',{'-160m', '-120m', '-80m', '-40m', '0m'}, ...
        'YTick', [-8, -6, -4, -2, 0])
    set(gca,'xTickLabel',{'35^oN', '36^oN', '37^oN', '38^oN', '39^oN', '40^oN', '41^oN', '42^oN'}, ...
        'XTick', [35, 36, 37, 38, 39, 40, 41, 42])
    set(gca, 'color', [0.8,0.8,0.8]);
    colormap(gca,jet);
    set(gcf,'PaperPosition',[0 0 30 12]);
    set(gca, 'TickDir','out');
%     axis equal;
    set(gca,'fontsize',20);