    hold on;
    shading interp;
    set(gca, 'color', [0.8,0.8,0.8]);
%     caxis([0.2,20]);
% c=colorbar;
% c.Label.String= 'Upper layer thickness(m)';
% c.Label.FontSize= 12;
    colormap(gca,jet);
%     set(gcf,'PaperPosition',[0 0 20 15]);
%     axis equal;
    set(gca,'fontsize',13);
   
    axis equal;
    axis tight;