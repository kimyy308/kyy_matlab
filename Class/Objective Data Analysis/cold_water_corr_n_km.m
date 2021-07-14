% % % 
% % % exp 1.
% % % 
clc; clear all; close all;

addpath(genpath('D:\MEPL\project\NWP\m_map'))

% rname='D:\MEPL\Ph_D_course\term_paper\ES_clim2.nc';
% temp1=ncread(rname,'CLIM_TEMP');
% u1=ncread(rname,'CLIM_U');
% v1=ncread(rname,'CLIM_V');

%%% only temperature
rname='D:\MEPL\Ph_D_course\term_paper\temp.nc';
temp2=ncread(rname,'ES_TEMP');
lon2=ncread(rname,'XT_OCEAN293_442');
lat2=ncread(rname,'YT_OCEAN281_470');
level2=ncread(rname,'ST_OCEAN1_33');
ilevel2 = -level2;

numx=length(lon2);
numy=length(lat2);
numz=length(ilevel2);
numt=length(temp2(1,1,1,:));

% % % monthly transport data
rname2='D:\MEPL\Ph_D_course\term_paper\transport2.nc';
tr=ncread(rname2,'TR');
% y_tr=ncread(rname2,'Y_TR');


% % masking land & ocean and change missing value to NaN
[Plg,Plt]=meshgrid(lon2,lat2);
for i=1:numx
    for j=1:numy
        for k=1:numz
            for l=1:numt
                if (temp2(i,j,k,l)< -1.0e+5)
                    mask(i,j)=0;
                    temp2(i,j,k,l)=NaN;
                else
                    mask(i,j)=1;
                end
            end
        end
    end
end
% m_proj('UTM','long',[128 143],'lat',[33 52]);
% m_proj('mercator','lon',[128 143],'lat',[33 52]);
% m_pcolor(Plg,Plt,squeeze(temp2(1:150,1:190,1,335))');



% % masking land & ocean.
% % land = 0, ocean = 1
for i=1:numx
    for j=1:numy
      mask(i,j)=isfinite(temp2(i,j,1,1));
%         if(isnan(temp2(i,j,1,1))) 
%             mask(i,j)=0;
%         else
%             mask(i,j)=1;
%         end
    end
end

mask_raw=mask';

% % % searching coastline index
ni=1;
startj=23;
startx=find(mask(:,startj)==1)
tempi(ni)=startx(1);
tempj(ni)=startj;
west=0; east=0; south=0; north=0;
if(mask(tempi(ni)-1,tempj(ni))) west=1; end %%west
if(mask(tempi(ni)+1,tempj(ni))) east=1; end %%east
if(mask(tempi(ni),tempj(ni)-1)) south=1; end %%south
if(mask(tempi(ni),tempj(ni)+1)) north=1; end %%north

if (south ==1 && east==1)  %%% past i, j --> west point
    pasti(ni)= tempi(ni)-1;
    pastj(ni)= tempj(ni);
    south=-1;
end
a=1;
tmask=mask_raw;
tmask(:,1)=0;  %%% -> close the open boundary

while(a<400)
    if(west==-1) %%% from west, curve to north
        if(tmask(tempj(ni)+1,tempi(ni))==1)
            ni=ni+1;
            tempi(ni)=tempi(ni-1);
            tempj(ni)=tempj(ni-1)+1;
            west=1;
            east=1;
            south=-1;
            north=1;
            a=a+1;
            display('from west, curve to north')
            tempi(ni-1)
            tempj(ni-1)
            tempi(ni)
            tempj(ni)
        else 
            if(tmask(tempj(ni),tempi(ni)+1)==0) %%% from west, curve to south
                ni=ni+1;
                tempi(ni)=tempi(ni-1);
                tempj(ni)=tempj(ni-1)-1;
                west=1;
                east=0;
                south=1;
                north=-1;
                a=a+1;
                display('from west, curve to south')
                tempi(ni-1)
                tempj(ni-1)
                tempi(ni)
                tempj(ni)
            else  %%% from west, go straight to east
                ni=ni+1;
                tempi(ni)=tempi(ni-1)+1;
                tempj(ni)=tempj(ni-1);
                west=-1;
                east=1;
                south=1;
                north=0;
                a=a+1;
                display('from west, go straight to east')
                tempi(ni-1)
                tempj(ni-1)
                tempi(ni)
                tempj(ni)
            end
        end
    elseif(south==-1) 
        if(tmask(tempj(ni),tempi(ni)-1)==1) %%% from south, curve to west
            ni=ni+1;
            tempi(ni)=tempi(ni-1)-1;
            tempj(ni)=tempj(ni-1);
            west=1;
            east=-1;
            south=1;
            north=1;
            a=a+1;
            display('from south, curve to west')
            tempi(ni-1)
            tempj(ni-1)
            tempi(ni)
            tempj(ni)
        else 
            if(tmask(tempj(ni)+1,tempi(ni))==0) %%% from south, curve to east
                ni=ni+1;
                tempi(ni)=tempi(ni-1)+1;
                tempj(ni)=tempj(ni-1);
                west=-1;
                east=1;
                south=1;
                north=0;
                a=a+1;
                display('from south, curve to east')
                tempi(ni-1)
                tempj(ni-1)
                tempi(ni)
                tempj(ni)
            else  %%% from south, go straight to north
                ni=ni+1;
                tempi(ni)=tempi(ni-1);
                tempj(ni)=tempj(ni-1)+1;
                west=0;
                east=1;
                south=-1;
                north=1;
                a=a+1;
                display('from south, go straight to north')
                tempi(ni-1)
                tempj(ni-1)
                tempi(ni)
                tempj(ni)
            end
        end
    elseif(east==-1) 
        if(tmask(tempj(ni)-1,tempi(ni))==1) %%% from east, curve to south
            ni=ni+1;
            tempi(ni)=tempi(ni-1);
            tempj(ni)=tempj(ni-1)-1;
            west=1;
            east=1;
            south=1;
            north=-1;
            a=a+1;
            display('from east, curve to south')
            tempi(ni-1)
            tempj(ni-1)
            tempi(ni)
            tempj(ni)
        else 
            if(tmask(tempj(ni),tempi(ni)-1)==0) %%% from east, curve to north 
                ni=ni+1;
                tempi(ni)=tempi(ni-1);
                tempj(ni)=tempj(ni-1)+1;
                west=0;
                east=1;
                south=-1;
                north=1;
                a=a+1;
                display('from east, curve to north')
                tempi(ni-1)
                tempj(ni-1)
                tempi(ni)
                tempj(ni)
            else  %%% from east, go straight to west
                ni=ni+1;
                tempi(ni)=tempi(ni-1)-1;
                tempj(ni)=tempj(ni-1);
                west=1;
                east=-1;
                south=0;
                north=1;
                a=a+1;
                display('from east, go straight to west')
                tempi(ni-1)
                tempj(ni-1)
                tempi(ni)
                tempj(ni)
            end
        end
    elseif(north==-1) 
        if(tmask(tempj(ni),tempi(ni)+1)==1) %%% from north, curve to east
            ni=ni+1;
            tempi(ni)=tempi(ni-1)+1;
            tempj(ni)=tempj(ni-1);
            west=-1;
            east=1;
            south=1;
            north=1;
            a=a+1;
            display('from north, curve to east')
            tempi(ni-1)
            tempj(ni-1)
            tempi(ni)
            tempj(ni)
        else 
            if(tmask(tempj(ni),tempi(ni)-1)==0) %%% from north, curve to west
                ni=ni+1;
                tempi(ni)=tempi(ni-1)-1;
                tempj(ni)=tempj(ni-1);
                west=1;
                east=-1;
                south=0;
                north=1;
                a=a+1;
                display('from north, curve to west')
                tempi(ni-1)
                tempj(ni-1)
                tempi(ni)
                tempj(ni)
            else  %%% from north, go straight to south
                ni=ni+1;
                tempi(ni)=tempi(ni-1);
                tempj(ni)=tempj(ni-1)-1;
                west=1;
                east=0;
                south=1;
                north=-1;
                a=a+1;
                display('from north, go straight to south')
                tempi(ni-1)
                tempj(ni-1)
                tempi(ni)
                tempj(ni)
            end
        end
    end
    if (tempi(ni)==numx ||tempj(ni)==numy)
        a=401;
        display('search complete')
    end
end

% % % draw circle and mask grid which distance from center of circle is
% % % closer than 500m
% % % this uses haversine formula to calculate the circle distance between
% % % two points (shortest distance over the earth's surface)
R = 6373.0 *1000  %%%% earth radius(m)
maxdist=70000;
tmask2=tmask;
for k=1:length(tempi)
    for i=tempi(k)-8 : tempi(k)+8
        for j=tempj(k)-8 : tempj(k)+8
            if (i<1 || i>numx || j<1 || j>numy)
            else
                clat=deg2rad(lat2(tempj(ni)));
                clon=deg2rad(lon2(tempi(ni)));
                tlat=deg2rad(lat2(j));
                tlon=deg2rad(lon2(i));
                dlon=tlon-clon;
                dlat=tlat-clat;
                dlon=deg2rad(lon2(i)-lon2(tempi(k)));
                dlat=deg2rad(lat2(j)-lat2(tempj(k)));
                a=sin(dlat / 2)^2 + cos(clat) * cos(tlat) * sin(dlon / 2)^2;
                c = 2 * atan2(sqrt(a), sqrt(1 - a));
                dist= R * c;
                if (abs(dist)<=maxdist)    %%%% distance between station and coastline is 40km
                    tmask2(j,i)=0;
                end
            end
        end
    end
end

% pcolor(tmask);
% figure;
% pcolor(tmask2);



% % % searching fixed line index
ni=1;
startj=23;
startx=find(tmask2(startj,:)==1)
tempi2(ni)=startx(1);
tempj2(ni)=startj;
% west=0; east=0; south=0; north=0;
% if(mask(tempi2(ni)-1,tempj2(ni))) west=1; end %%west
% if(mask(tempi2(ni)+1,tempj2(ni))) east=1; end %%east
% if(mask(tempi2(ni),tempj2(ni)-1)) south=1; end %%south
% if(mask(tempi2(ni),tempj2(ni)+1)) north=1; end %%north
% 
% if (south ==1 && east==1)  %%% past i, j --> west point
%     pasti(ni)= tempi2(ni)-1;
%     pastj(ni)= tempj2(ni);
%     south=-1;
% end
south=-1;
a=1;

while(a<400)
    if(west==-1) %%% from west, curve to north
        if(tmask2(tempj2(ni)+1,tempi2(ni))==1)
            ni=ni+1;
            tempi2(ni)=tempi2(ni-1);
            tempj2(ni)=tempj2(ni-1)+1;
            west=1;
            east=1;
            south=-1;
            north=1;
            a=a+1;
            display('from west, curve to north')
            tempi2(ni-1)
            tempj2(ni-1)
            tempi2(ni)
            tempj2(ni)
        else 
            if(tmask2(tempj2(ni),tempi2(ni)+1)==0) %%% from west, curve to south
                ni=ni+1;
                tempi2(ni)=tempi2(ni-1);
                tempj2(ni)=tempj2(ni-1)-1;
                west=1;
                east=0;
                south=1;
                north=-1;
                a=a+1;
                display('from west, curve to south')
                tempi2(ni-1)
                tempj2(ni-1)
                tempi2(ni)
                tempj2(ni)
            else  %%% from west, go straight to east
                ni=ni+1;
                tempi2(ni)=tempi2(ni-1)+1;
                tempj2(ni)=tempj2(ni-1);
                west=-1;
                east=1;
                south=1;
                north=0;
                a=a+1;
                display('from west, go straight to east')
                tempi2(ni-1)
                tempj2(ni-1)
                tempi2(ni)
                tempj2(ni)
            end
        end
    elseif(south==-1) 
        if(tmask2(tempj2(ni),tempi2(ni)-1)==1) %%% from south, curve to west
            ni=ni+1;
            tempi2(ni)=tempi2(ni-1)-1;
            tempj2(ni)=tempj2(ni-1);
            west=1;
            east=-1;
            south=1;
            north=1;
            a=a+1;
            display('from south, curve to west')
            tempi2(ni-1)
            tempj2(ni-1)
            tempi2(ni)
            tempj2(ni)
        else 
            if(tmask2(tempj2(ni)+1,tempi2(ni))==0) %%% from south, curve to east
                ni=ni+1;
                tempi2(ni)=tempi2(ni-1)+1;
                tempj2(ni)=tempj2(ni-1);
                west=-1;
                east=1;
                south=1;
                north=0;
                a=a+1;
                display('from south, curve to east')
                tempi2(ni-1)
                tempj2(ni-1)
                tempi2(ni)
                tempj2(ni)
            else  %%% from south, go straight to north
                ni=ni+1;
                tempi2(ni)=tempi2(ni-1);
                tempj2(ni)=tempj2(ni-1)+1;
                west=0;
                east=1;
                south=-1;
                north=1;
                a=a+1;
                display('from south, go straight to north')
                tempi2(ni-1)
                tempj2(ni-1)
                tempi2(ni)
                tempj2(ni)
            end
        end
    elseif(east==-1) 
        if(tmask2(tempj2(ni)-1,tempi2(ni))==1) %%% from east, curve to south
            ni=ni+1;
            tempi2(ni)=tempi2(ni-1);
            tempj2(ni)=tempj2(ni-1)-1;
            west=1;
            east=1;
            south=1;
            north=-1;
            a=a+1;
            display('from east, curve to south')
            tempi2(ni-1)
            tempj2(ni-1)
            tempi2(ni)
            tempj2(ni)
        else 
            if(tmask2(tempj2(ni),tempi2(ni)-1)==0) %%% from east, curve to north 
                ni=ni+1;
                tempi2(ni)=tempi2(ni-1);
                tempj2(ni)=tempj2(ni-1)+1;
                west=0;
                east=1;
                south=-1;
                north=1;
                a=a+1;
                display('from east, curve to north')
                tempi2(ni-1)
                tempj2(ni-1)
                tempi2(ni)
                tempj2(ni)
            else  %%% from east, go straight to west
                ni=ni+1;
                tempi2(ni)=tempi2(ni-1)-1;
                tempj2(ni)=tempj2(ni-1);
                west=1;
                east=-1;
                south=0;
                north=1;
                a=a+1;
                display('from east, go straight to west')
                tempi2(ni-1)
                tempj2(ni-1)
                tempi2(ni)
                tempj2(ni)
            end
        end
    elseif(north==-1) 
        if(tmask2(tempj2(ni),tempi2(ni)+1)==1) %%% from north, curve to east
            ni=ni+1;
            tempi2(ni)=tempi2(ni-1)+1;
            tempj2(ni)=tempj2(ni-1);
            west=-1;
            east=1;
            south=1;
            north=1;
            a=a+1;
            display('from north, curve to east')
            tempi2(ni-1)
            tempj2(ni-1)
            tempi2(ni)
            tempj2(ni)
        else 
            if(tmask2(tempj2(ni),tempi2(ni)-1)==0) %%% from north, curve to west
                ni=ni+1;
                tempi2(ni)=tempi2(ni-1)-1;
                tempj2(ni)=tempj2(ni-1);
                west=1;
                east=-1;
                south=0;
                north=1;
                a=a+1;
                display('from north, curve to west')
                tempi2(ni-1)
                tempj2(ni-1)
                tempi2(ni)
                tempj2(ni)
            else  %%% from north, go straight to south
                ni=ni+1;
                tempi2(ni)=tempi2(ni-1);
                tempj2(ni)=tempj2(ni-1)-1;
                west=1;
                east=0;
                south=1;
                north=-1;
                a=a+1;
                display('from north, go straight to south')
                tempi2(ni-1)
                tempj2(ni-1)
                tempi2(ni)
                tempj2(ni)
            end
        end
    end
    if (tempi2(ni)>130 ||tempj2(ni)==numy)
        a=401;
        display('search station that distance from coastline is uniform complete')
    end
end


% % % make new temperature variable
ksbcw=squeeze(temp2(17,26,16,:));
m_ksbcw=ksbcw-mean(ksbcw);
for k=1:length(tempi2)
    coasttemp(k,:,:)=temp2(tempi2(k),tempj2(k),:,:);
end
m_coasttemp=coasttemp-mean(coasttemp,3);
for j=1:length(squeeze(coasttemp(:,1,1)))
    for k=1:length(squeeze(coasttemp(1,:,1)))
%         corcoef_coast(j,k)=corr(m_ksbcw,squeeze(m_coasttemp(j,k,:)));   
        corcoef_coast2(j,k)=corr(ksbcw,squeeze(coasttemp(j,k,:)));   
    end
end

% % % coastline latitude plot
plot(1:length(tempj2),lat2(tempj2))
filename=strcat('D:\MEPL\Ph_D_course\term_paper\coastlat_',num2str(maxdist/1000),'km.tif')
saveas(gcf,filename,'tiff');

% % % coastline
tmask3=tmask;
for k=1:length(tempi2)
    tmask3(tempj2(k),tempi2(k))=10;
end

pcolor(lon2,lat2,tmask3);
% shading interp; 
% colormap jet; 
set(gca, 'color', [0.8,0.8,0.8]);
% c=colorbar;
% c.Label.String= 'Temp (^oC)';
xlabel('Longitude (^oE)','fontsize',15);
ylabel('Latitude (^oN)','fontsize',15);
set(gca,'fontsize',10);
shading flat
hold on;
% [C hT]= contour(lon2(5:65),lat2(15:70),squeeze(temp2(5:65,15:70,16,335))', '-k','ShowText','on');
filename=strcat('D:\MEPL\Ph_D_course\term_paper\coastline_',num2str(maxdist/1000),'km.tif')
hold off;
fig=gcf;
fig.InvertHardcopy='off';
set(gca, 'color', [0.8,0.8,0.8]);
axis equal;
axis tight;
saveas(gcf,filename,'tiff');






% % % ksbcw & vertical section along coastline corr
pcolor(1:length(tempi2),ilevel2(1:33),corcoef_coast2');
shading interp; colormap jet; set(gca, 'color', [0.8,0.8,0.8]);
c=colorbar;
c.Label.String= 'R';
caxis([-1 1])
% xlabel('Latitude (^oN)','fontsize',10);
ylabel('depth(m)','fontsize',10);
set(gca,'fontsize',15);
filename=strcat('D:\MEPL\Ph_D_course\term_paper\coastlinecorr_',num2str(maxdist/1000),'km.tif')
hold off;
fig=gcf;
fig.InvertHardcopy='off';
set(gca, 'color', [0.8,0.8,0.8]);
% axis equal;
axis tight;
saveas(gcf,filename,'tiff');


% % % % ksbcw & vertical section along coastline corr
% pcolor(1:length(tempi2),ilevel2(1:33),corcoef_coast');
% shading interp; colormap jet; set(gca, 'color', [0.8,0.8,0.8]);
% c=colorbar;
% c.Label.String= 'R';
% caxis([-1 1])
% % xlabel('Latitude (^oN)','fontsize',10);
% ylabel('depth(m)','fontsize',10);
% set(gca,'fontsize',15);
% filename='D:\MEPL\Ph_D_course\term_paper\m_coastlinecorr.tif'
% hold off;
% fig=gcf;
% fig.InvertHardcopy='off';
% set(gca, 'color', [0.8,0.8,0.8]);
% % axis equal;
% axis tight;
% saveas(gcf,filename,'tiff');




