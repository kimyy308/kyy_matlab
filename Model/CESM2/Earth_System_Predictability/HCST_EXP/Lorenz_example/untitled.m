two axis example
figH = figure;
axLH = gca;
axRH = axes('color','none');
mslplot{1}=plot(inputyear,Model.amp_S2(3,:), 'b','parent',axLH);
mslplot{2}=plot(inputyear,Obs.amp_S2(3,:), 'k','parent',axRH);
ylabel(axLH,'Model S2 Tidal amplitude (cm)')
ylabel(axRH,'Obs S2 Tidal amplitude (cm)')
ax_pos = get(axLH,'position');
set(axLH,'yaxislocation','left','position',ax_pos+[0 0.02 -0.01 -0.02]);
set(axRH,'color','none','yaxislocation','right','xtick', inputyear, 'position', ax_pos+[0 0.02 -0.01 -0.02]);
%             set(axRH,'color','none','yaxislocation','right');

%             set(axRH,'xticklabel',[], 'xtick', get(axLH,'xtick'));
set(axLH,'xlim', get(axRH, 'xlim'), 'xtick', get(axRH,'xtick'),'xticklabel',[], 'xaxislocation','top');
set(axLH,'ycolor','b', 'box', 'off', 'FontSize',15);
set(axRH,'ycolor','k', 'box', 'off', 'FontSize',15);
xlabel(axRH, 'Year');

title([num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i')])
% datetick('x','yyyy','keepticks')
axis tight;
% ylim(meanplotlev2)
set(mslplot{1},'LineWidth',2);
set(mslplot{2},'LineWidth',2);
grid on

lgd=legend([mslplot{1} mslplot{2}], 'Model','TG-UST');

%             lgd=legend('Model','TG-UST');
set(lgd,'FontSize',15);
set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
set(lgd,'Orientation','horizontal');

set(gcf,'PaperPosition', [0 0 36 12]) 
saveas(gcf,pngname,'tif');
grid off
close all;