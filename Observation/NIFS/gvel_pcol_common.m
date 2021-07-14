%     hold on;
    shading interp;
    xlim([min(lon2(lonminind:5)) max(lon2(lonminind:5))]);
%     set(gca,'yTickLabel',{'-160m', '-120m', '-80m', '-40m', '0m'}, ...
%         'YTick', [-160, -120, -80, -40, 0])
    set(gca,'yTickLabel',{'-160m', '-120m', '-80m', '-40m', '0m'}, ...
        'YTick', [-8, -6, -4, -2, 0])
%     set(gca,'xTickLabel',{'128.3^oE', '128.8^oE', '129.3^oE', '129.8^oE', '130.3^oE'}, ...
    set(gca,'xTickLabel',{'129.8^oE','130.0^oE','130.2^oE','130.4^oE','130.6^oE','130.8^oE'}, ...
        'XTick', (129.8:0.2:130.8))
    set(gca, 'color', [0.8,0.8,0.8]);
%     load jet_mod;
%     colormap(gca,jet_mod);
bwr_map;
    colormap(gca,bwrmap);
    set(gcf,'PaperPosition',[0 0 28 12]);
%     axis equal;
    set(gca,'fontsize',20);
    set(gca, 'TickDir','out');