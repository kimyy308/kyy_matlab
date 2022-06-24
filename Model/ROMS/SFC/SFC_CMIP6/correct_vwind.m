clc; close all; clear all;


rootdir='/data1/RCM/CMIP6/input_hdd/nwp_1_20/input/SBC/';
scenname='ssp585';
modelname ='EC-Earth3-Veg';
expname ='nwp_1_20';
filedir=[rootdir, scenname, '/', modelname, '/'];
year = 2036;
filename = [modelname, '_', expname, '_', num2str(year), '_', 'Vwind.nc.org'];
writefilename = [modelname, '_', expname, '_', num2str(year), '_', 'Vwind.nc'];

tgtday = 337;
for ti=1:7
   Vwind(:,:,ti)= ncread([filedir, filename], 'Vwind', [1, 1, tgtday-4+ti], [inf, inf, 1]);
end

% pcolor(Vwind(:,:,1)'); shading flat; colorbar;


% for ti=1:3
%     Vwind_mod(:,:,ti) = Vwind(:,:,ti+2);
%     for xi= 50: 150
%         for yi= 500 : 600
%             Vwind_mod(xi,yi,ti)= mean(Vwind(xi,yi,ti:ti+4));
%         end
%     end
% end
% pcolor(Vwind_mod(:,:,1)'); shading flat; colorbar;


[xf, yf] = find(Vwind(:,:,4)<-15);
for ti=1:5
    Vwind_mod2(:,:,ti) = Vwind(:,:,ti+2);
    for fi= 1:length(xf)
        Vwind_mod2(xf(fi),yf(fi),ti)= mean(Vwind(xf(fi),yf(fi),ti:ti+2));
    end
end
pcolor(Vwind_mod2(:,:,2)'); shading flat; colorbar;

Vwind_mod3 = Vwind_mod2;
hor_filt_range=10;
for ti=1:5
    for xfilti = 50 : 200
        for yfilti = 450 : 600
            Vwind_mod3(xfilti, yfilti,ti)= ...
                mean(Vwind_mod2(xfilti-hor_filt_range:xfilti+hor_filt_range, ...
                yfilti-hor_filt_range : yfilti + hor_filt_range,ti),'all');
        end
    end
end
pcolor(Vwind_mod3(:,:,2)'); shading flat; colorbar;



for ti=1:7
    Vwind_mod4(:,:,ti) = Vwind(:,:,ti);
    [xf, yf] = find(Vwind(:,:,ti)<-15);
    for fi= 1:length(xf)
        Vwind_mod4(xf(fi),yf(fi),ti)= -15;
    end
end
pcolor(Vwind_mod4(:,:,4)'); shading flat; colorbar;

pcolor(Vwind(:,:,4)'); shading flat; colorbar;


% %  write
for ti= 1:7
    ncwrite([filedir, writefilename], 'Vwind', Vwind_mod4(:,:,ti), [1, 1, tgtday-4+ti]);
end


% filt_ex=imgaussfilt(Vwind_mod2(:,:,2), 'FilterSize', 9);