[lon44,ilevel44] = meshgrid(ilevel2(1:23),lon2(1:30));
pcolindex=1;
% % % 
% % % Exp1, temp(surf) + vec(surf)
% % % 

sbp1=subplot(2,2,1);
pcolor(lon2,lat2,(temp2mean_1(:,:,colordepthindex,1)'));
temp_pcolor_common1;
vector_common1;
hold on;
quiver(lon(1:3:110,1:3:120),lat(1:3:110,1:3:120), ...
  squeeze(u2mean_1(1:3:120,1:3:110,quivdepthindex,1))'*2,squeeze(v2mean_1(1:3:120,1:3:110,quivdepthindex,1))'*2, ...
  'k','AutoScale','off');
% ar_1=annotation(gcf,'rectangle',...
% [0.465 0.3815 0.00353378378378377 0.07],...
% 'FaceColor','red');
% tex=text(132.815,37.1,'A');
% set(tex,'fontsize',20);
% br_1=annotation(gcf,'rectangle',...
% [0.225 0.55 0.07010055865921788 0.00353378378378377 ],...    
%  'FaceColor','red');
% tex=text(127.8,40.05,'B');
% set(tex,'fontsize',20);
vector_common2;
hold off;

    

sbp2=subplot(2,2,2);
pcolor(lon2,lat2,(temp2mean_2(:,:,colordepthindex,1)'));
temp_pcolor_common1;
vector_common1;
hold on;
quiver(lon(1:3:110,1:3:120),lat(1:3:110,1:3:120), ...
  squeeze(u2mean_2(1:3:120,1:3:110,quivdepthindex,1))'*2,squeeze(v2mean_2(1:3:120,1:3:110,quivdepthindex,1))'*2, ...
  'k','AutoScale','off');
% ar_2=annotation(gcf,'rectangle',...
% [0.465 0.3815 0.00353378378378377 0.07],...
% 'FaceColor','red');
% tex=text(132.815,37.1,'A');
% set(tex,'fontsize',20);
% br_2=annotation(gcf,'rectangle',...
% [0.225 0.55 0.07010055865921788 0.00353378378378377 ],...    
%  'FaceColor','red');
% tex=text(127.8,40.05,'B');
% set(tex,'fontsize',20);
vector_common2;
hold off;



sbp3=subplot(2,2,3);

vector_common1;
hold on;
pcolor(lon2,lat2,(temp2mean_3(:,:,colordepthindex,1)'));
temp_pcolor_common1;
quiver(lon(1:3:110,1:3:120),lat(1:3:110,1:3:120), ...
  squeeze(u2mean_3(1:3:120,1:3:110,quivdepthindex,1))'*2,squeeze(v2mean_3(1:3:120,1:3:110,quivdepthindex,1))'*2, ...
  'k','AutoScale','off');
% ar_3=annotation(gcf,'rectangle',...
% [0.465 0.3815 0.00353378378378377 0.07],...
% 'FaceColor','red');
% tex=text(132.815,37.1,'A');
% set(tex,'fontsize',20);
% br_3=annotation(gcf,'rectangle',...
% [0.225 0.55 0.07010055865921788 0.00353378378378377 ],...    
%  'FaceColor','red');
% tex=text(127.8,40.05,'B');
% set(tex,'fontsize',20);
vector_common2;
hold off;


sbp3=subplot(2,2,4);
vector_common1;
hold on;
pcolor(lon2,lat2,(temp2mean_4(:,:,colordepthindex,1)'));
temp_pcolor_common1;
quiver(lon(1:3:110,1:3:120),lat(1:3:110,1:3:120), ...
  squeeze(u2mean_4(1:3:120,1:3:110,quivdepthindex,1))'*2,squeeze(v2mean_4(1:3:120,1:3:110,quivdepthindex,1))'*2, ...
  'k','AutoScale','off');
% ar_4=annotation(gcf,'rectangle',...
% [0.465 0.3815 0.00353378378378377 0.07],...
% 'FaceColor','red');
% tex=text(132.815,37.1,'A');
% set(tex,'fontsize',20);
% br_4=annotation(gcf,'rectangle',...
% [0.225 0.55 0.07010055865921788 0.00353378378378377 ],...    
%  'FaceColor','red');
% tex=text(127.8,40.05,'B');
% set(tex,'fontsize',20);
h=colorbar;
set(h, 'Position', [.9150 .11 .0181 .8150])  % right, up, width, height
set(h,...
    'YTickLabel',{'2^oC','4^oC','6^oC','8^oC','10^oC','12^oC','14^oC','16^oC','18^oC','20^oC'},...
        'YTick',[2 4 6 8 10 12 14 16 18 20]);
h.Label.FontSize=7;
% set(h,'Visible','off');
vector_common2;
hold off;
    

% delete(ar_1);
% delete(br_1);
% delete(ar_2);
% delete(br_2);     
% delete(ar_3);
% delete(br_3);     
% delete(ar_4);
% delete(br_4);     
    
    
    
    
    
    
    
  