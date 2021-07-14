
zlim([-201 0]);
grid on
hold off
colormap jet;
colormap(flipud(colormap))
set(gca,'PlotBoxAspectRatio',[2.5 10 2])
set(gcf,'PaperPosition',[0 0 22 15]);
set(gca,'fontsize',15);
xlabel('Longitude (^oE)','fontsize',12);
ylabel('Latitude (^oN)','fontsize',12);
zlabel('Depth (M)','fontsize',12);
fig=gcf;
fig.InvertHardcopy='off';
c=colorbar;
c.Label.String= 'Depth (m)';
c.Label.FontSize= 12;