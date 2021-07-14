x=[1,2,4,6,8]';
y=[100,140,160,170,175].';
g = fittype('a-b*exp(-c*x)');
f0 = fit(x,y,g,'StartPoint',[[ones(size(x)), -exp(-x)]\y; 1]);
xx = linspace(1,8,50);
plot(x,y,'o',xx,f0(xx),'r-');



inputyear=[2006 : 2100]
testname='test66'
regionname='AKP4'
filedir = strcat('D:\', 'Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
interpedfilename = strcat(filedir, testname,'_',regionname, 'cmems_interped_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
interped_ssh=ncread(interpedfilename, 'interped_ssh');
reshap_interped_ssh=reshape(interped_ssh,[104 90 12 1140/12]);
yearly_interped_ssh=mean(reshap_interped_ssh,3);

x=squeeze(inputyear)'-2006;
y=squeeze(yearly_interped_ssh(28,22,:));
% y = y + 0.1*rand(size(y));
% x = x(3:2:end-5);
% y = y(3:2:end-5);
% y(end)=y(end)+5;
f = fit(x,y,'exp1');
plot(f,x,y)
SLR1=[num2str(round((f(x(end))-f(x(1)))*100,2)), ' cm']
hold on
[f2, gof2] = fit(x,y,'poly2');
plot(f2, '-')
SLR2=[num2str(round((f2(x(end))-f2(x(1)))*100,2)), ' cm']
f3 = fit(x,y,'poly1');
plot(f3, '--')
SLR3=[num2str(round((f3(x(end))-f3(x(1)))*100,2)), ' cm']
% a=polyfit(x,y,1)
% y=a(1).*x+a(2) +0.3465-0.2748
% plot(x,y, '--')
hold off
legend('ssh', ['exp, ', SLR1], ['poly2, ', SLR2], ['poly1, ', SLR3])
set(gca, 'Fontsize', 20)

% % % 
% % % 
% % % % g = fittype('a-b*exp(-c*x)');
% % % % fo = fitoptions('Method','NonlinearLeastSquares',...
% % % %                'Lower',[0,0],...
% % % %                'Upper',[Inf,max(x)],...
% % % %                'StartPoint',[1 1]);
% % % % g = fittype('a*(x-b)^n','problem','n','options',fo);
% % % % f0 = fit(x,y,g,'problem',3);
% % % 
% % % f0 = fit(x,y,'poly2');
% % % 
% % % xx = x;
% % % plot(x,y,'o',xx,f0(xx),'r-');



x=squeeze(inputyear)';
for i=1:size(yearly_interped_ssh,1)
    for j=1:size(yearly_interped_ssh,2)
        if isfinite(yearly_interped_ssh(i,j,1))
            y=squeeze(yearly_interped_ssh(i,j,:));
            f = fit(x,y,'exp1');
            comb_SLR1(i,j)=f(x(end))-f(x(1));
            f = fit(x,y,'poly2');
            comb_SLR2(i,j)=f(x(end))-f(x(1));
            f = fit(x,y,'poly1');
            comb_SLR3(i,j)=f(x(end))-f(x(1));
        else
            comb_SLR1(i,j)=NaN;
            comb_SLR2(i,j)=NaN;
            comb_SLR3(i,j)=NaN;
        end
    end
end
mean(comb_SLR1(:),'omitnan')
mean(comb_SLR2(:),'omitnan')
mean(comb_SLR3(:),'omitnan')

pcolor(comb_SLR1' *100)
colormap(jet)
colorbar
caxis([50 110])
shading flat

pcolor(comb_SLR2' *100)
colormap(jet)
colorbar
caxis([50 110])
shading flat

pcolor(comb_SLR3' *100)
colormap(jet)
colorbar
caxis([50 110])
shading flat

