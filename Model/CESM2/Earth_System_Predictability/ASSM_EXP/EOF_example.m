clc; clear all; close all;
%%
%4�� ����

 %read ssh data
X = textread('Eta_bp.dat');
X = X(:,2:end);
X = X' ; %transpose : (t,x) -> (x,t)
[M N] = size(X);

 %read station location data
[lon lat station] = textread('station.dat','%f %f %s');

 %solve eigen value problem
R = 1/N * X * X';
[V, D] = eig(R); %eigen vector , eigen value

%sort eigen values 
L = diag(D);L = sort(L,'descend'); %sorted eigen value
for i = 1:length(L);
    ind = find(L == L(i));
    E(:,i) = V(:,ind); %E : sorted eigen vector
end;

 % find PC 
Z = E' * X; %Z: PC

 %find variance / total variance
total_variance = sum(L) ;
mode_variance = abs(L)/total_variance*100; %�� mode�� ��������

%find K number 
K = 5;

%%%plot %%%%
%3) plot L 
plot(mode_variance,'*','LineWidth',3);set(gca,'XTick',[1:1:M]);
title('Explained variance of Each Mode','fontsize',30);
xlabel('Mode number','fontsize', 20);
ylabel('variance/total variance (%)','fontsize',20);set(gca,'FontSize',18);
saveas(gcf,'ex4_varialbes.png','png');
close;

%1,2) plot E,Z

% Station positon
m_proj('mercator','lon',[125 145],'lat',[30 50]);
m_grid('box','fancy','tickdir','in');
m_gshhs_i('color','k');
m_gshhs_i('patch',[.6 .6 .6]);
hold on
m_plot(lon,lat,'bo'); hold on
title('Station','fontsize',30);set(gca,'FontSize',18);
saveas(gcf,'Station.png','png');
close;

%1) plot E 
 %set interpolation grid
[xx, yy] = meshgrid([125:0.1:145],[30:0.1:50]);
for i = 1 : K
    %interpolation
    e = griddata(lon,lat, E(:,i), xx, yy);
    %set title name
    title_name = strcat('Eigen Vector-','MODE',num2str(i));
    %plot 
    figure; hold on;
    m_proj('mercator','lon',[125 145],'lat',[30 50]);
    m_contourf(xx,yy,e,20); m_contour(xx,yy,e,20);
    m_grid('box','fancy','tickdir','in'); 
    m_gshhs_i('color','k'); m_gshhs_i('patch',[.6 .6 .6]);
    title(title_name,'fontsize',30);set(gca,'FontSize',18);
    xlabel('Lon','fontsize', 20);ylabel('Lat(dimension)','fontsize',20);
    colorbar;
    saveas(gcf,strcat('ex4',title_name,'.png'),'png');
    close;
end

%2) plot Z
for i = 1:K
    plot(Z(i,:)); 
    title_name = strcat('PC- ','MODE',num2str(i));
    title(title_name,'fontsize',30);set(gca,'FontSize',18);
    xlabel('time','fontsize', 20);ylabel('value(dimension)','fontsize',20);
    saveas(gcf,strcat('ex4',title_name,'.png'),'png');
    close;
end

close;