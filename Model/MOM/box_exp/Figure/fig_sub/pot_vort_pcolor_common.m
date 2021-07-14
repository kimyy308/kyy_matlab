pcolor(lon5,lat5,(pot_vort5mean(:,:,1))');
hold on;
shading interp;
set(gca, 'color', [0.8,0.8,0.8]);
 caxis([4.1e-7,4.3e-6]);
c=colorbar;
c.Label.String= 'Potential Vorticity';
c.Label.FontSize= 12;
xlabel('Longitude (^oE)','fontsize',27);
ylabel('Latitude (^oN)','fontsize',27);
set(gca,'fontsize',15); 
colormap(jet); 
set(gcf,'PaperPosition',[0 0 20 15]);
axis equal;
uu2mean(:,:)=u2mean(:,:,1,1);
vv2mean(:,:)=v2mean(:,:,1,1);
lonn=128.05:0.1:139.95;
latt=34.05:0.1:44.95;
cc=streamslice(lonn,latt,uu2mean',vv2mean');
set(cc,'color','k')
hold off;
fig.InvertHardcopy='off';
set(gca, 'color', [0.8,0.8,0.8]);
saveas(gcf,filename,'tiff');