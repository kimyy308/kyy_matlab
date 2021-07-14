[r,c,v] = ind2sub(size(comb_MyOcean_data),find(comb_MyOcean_data >= 100));
%[lon lat 

%aa=comb_MyOcean_data(r(1),c(1),v(1),:);
% [q w] = find(comb_MyOcean_data >= 200);
% comb_MyOcean_data(q,w)
      


[a b]=meshgrid(MyOcean_lon,MyOcean_lat);

figure
pcolor(a,b,comb_MyOcean_data(:,:,v(1))')
shading flat
colorbar


yearij=5
monthji=8
aa=166
bb=43
plot(point_temp2(1:23),MyOcean_st(1:23))
set(gca,'ydir','reverse')


%%

aaa=squeeze(comb_MyOcean_data(aa,bb,:));
aa=8; bb=79


%% 
> apink=squeeze(nanmean(comb_MyOcean_data,3));
pcolor(aaa')
shading flat
colorbar
%% regression


tt3 = datenum(tt1,tt2,15,0,0,0);

MyOcean_temp_divided=zeros(360,240);
error=[];
for aa=1:length(MyOcean_lon)
    for bb=1:length(MyOcean_lat)
        data2 = comb_MyOcean_data(aa,bb,:); %¿Âµµ
        data = squeeze(data2);
        [r,m,d]=regression(double(tt3'),double(data'));%tt¿Ídata°¡ ¼¼·ÎÇàÇüÅÂ·Î µé¾î°¡¾ßÇÔ
        
        notnan=find(isnan(data)==0);%data°¡ nanÀÌ ¾Æ´Ñ°Í
        if length(notnan)==0
            m=NaN;
        elseif length(notnan)>0&&length(notnan)<5
           error= [error;aa bb] ;
        else
        end
        
        MyOcean_temp_divided(aa,bb)=m*365; % Æ¯Á¤À§°æµµ¿ùÀÇ ¼ö¿ÂÆ®·»µå
    end
end

%% figure regression
myocean_temp_data=MyOcean_temp_divided;
[x,y]=find(myocean_temp_data==0);
for i=1:length(x)
    myocean_temp_data(x(i),y(i))=NaN;
end
pcolor(myocean_temp_data')
shading flat
colorbar

%% steric sea level


% 35g/kg,10¡ÆC,1000bar
% rho(t) - rho(10)/rho(10) * z
rho_10 = swp('rho',10,35000) % kg/m3
t = myocean_temp_data +10; %¡ÆC 

rho_t=NaN(length(MyOcean_lon),length(MyOcean_lat));
for ii=1:length(MyOcean_lon)
    for jj=1:length(MyOcean_lat)
       rho_t(ii,jj)=swp('rho',t(ii,jj),35000);
    end
end

steric_sea_lv=zeros(306,240);
for kk=1:length(MyOcean_lon)
    for jj=1:length(MyOcean_lat)
        if isnan(rho_t(kk,jj))==0%nanÀÌ¾Æ´Ñ°Í
            rho_cal=(rho_t(kk,jj) - rho_10)/rho_10;
            steric_sea_lv(kk,jj) = rho_cal*depth(kk,jj);
        else
            steric_sea_lv(kk,jj) = NaN;
        end
    end
end

pcolor(steric_sea_lv')
shading flat; colorbar



%%

apink= nanmean(v5,3);

lonlat=[115 145 30 52];
% val_min=-0.2; val_max=0.2; val_con=0.05; % color
m_proj('miller','lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
hold on;
m_gshhs_c('patch',[.7 .7 .7],'edgecolor','none');
m_grid('linestyle','none','Tickdir','out', 'xtick',lonlat(1):5:lonlat(2),...
'ytick',lonlat(3):3:lonlat(4), 'fontsize', 8, 'XaxisLocation', 'bottom');
[c,h]=m_contourf(MyOcean_lon,MyOcean_lat,apink'); %contour style
set(h,'linestyle','none');
% m_contour(soda_lon,soda_lat,soda_data',[val_min:val_con:val_max],'k');
% colormap('jet'); caxis([val_min val_max]); %color qjadnl wlwjd
h = colorbar;



