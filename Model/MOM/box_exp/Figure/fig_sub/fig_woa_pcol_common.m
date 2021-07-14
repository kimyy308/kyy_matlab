%     hold on;
    shading interp;
    xlim([min(lat2_woa(latminind:latmaxind)) max(lat2_woa(latminind:latmaxind))]);
%     set(gca,'yTickLabel',{'-160m', '-120m', '-80m', '-40m', '0m'}, ...
%         'YTick', [-160, -120, -80, -40, 0])
    set(gca,'yTickLabel',{'-160m', '-120m', '-80m', '-40m', '0m'}, ...
        'YTick', [-8, -6, -4, -2, 0])
%     set(gca,'xTickLabel',{'128.3^oE', '128.8^oE', '129.3^oE', '129.8^oE', '130.3^oE'}, ...
    set(gca,'xTickLabel',{'35^oN','36^oN','37^oN','38^oN'}, ...
        'XTick', (35:1:39))
    set(gca, 'color', [0.8,0.8,0.8]);
%     load jet_mod;
%     colormap(gca,jet_mod);
    colormap(gca,jet);
    set(gcf,'PaperPosition',[0 0 16 12]);
%     axis equal;
    set(gca,'fontsize',20);