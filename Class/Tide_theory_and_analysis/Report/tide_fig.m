clear all; close all;
% 
% Netcdf_old 
% 
% addpath(genpath('C:\Users\kyy\Desktop\nwp\netcdf_old'))

% % http://volkov.oce.orst.edu/tides/tpxo8_atlas.html
% % TPX08 + atlas reanalysis data

% 
% m_map
% 
% addpath(genpath('C:\Users\kyy\Desktop\nwp\m_map'))

% 
% t_tide
% 
addpath(genpath('D:\class\PhD\Tide_theory_and_analysis\harmonic_analysis\t_tide_v1.3beta'))

% UH# GLOSS#    Version Location    Country Latitude    Longitude   Start       End
% A : 021	176	a	Juan Fernandez	Chile	-33.61700	-78.83300	1977-01-26	1984-05-03
% B : 021	176	b	Juan Fernandez	Chile	-33.61700	-78.83300	1985-09-04	2014-12-31
% Units        : millimeters  (change to cm in the code)

% 
% Time ref : GMT (region time : GMT - 4)

% previously within
%                              separate series (022B and 022C) were 
%                              replaced in March 2004 with the originator's 
%                              updated series, now referred to as 022B.

% The mean level in 2005 is low relative to
%                              satellite altimetry and comparison with 
%                              Antofagasta by about 20 cm.

time=ncread('D:\class\PhD\Tide_theory_and_analysis\harmonic_analysis\OS_UH-RQH021B_20160323_R.nc','time');
% time=time+4;
depth=ncread('D:\class\PhD\Tide_theory_and_analysis\harmonic_analysis\OS_UH-RQH021B_20160323_R.nc','depth');
lat=ncread('D:\class\PhD\Tide_theory_and_analysis\harmonic_analysis\OS_UH-RQH021B_20160323_R.nc','latitude');
lon=ncread('D:\class\PhD\Tide_theory_and_analysis\harmonic_analysis\OS_UH-RQH021B_20160323_R.nc','longitude');
ssh=ncread('D:\class\PhD\Tide_theory_and_analysis\harmonic_analysis\OS_UH-RQH021B_20160323_R.nc','sea_surface_height_above_reference_level');
sensor_type=ncread('D:\class\PhD\Tide_theory_and_analysis\harmonic_analysis\OS_UH-RQH021B_20160323_R.nc','sensor_type_code');
ssh=squeeze(ssh);
ssh=ssh/10; %%% unit (mm -> cm)
sensor_type=squeeze(sensor_type);
ncdisp('D:\class\PhD\Tide_theory_and_analysis\harmonic_analysis\OS_UH-RQH021B_20160323_R.nc');

ind=1;
checkssh(1)=1;
startssh(1)=0;
validssh=isnan(ssh);
for i=2:(size(ssh)-1)
    if (validssh(i)==0)
        if(validssh(i-1)==1)
            startssh(ind)=i;
        elseif(validssh(i+1)==1)
            ind=ind+1;
            checkssh(ind)=1;
        else
            checkssh(ind)=checkssh(ind)+1;
        end
    end
end
clear validssh ind
maxperiod=max(checkssh);
startperiod=startssh(find(checkssh==max(checkssh)));
clear startssh checkssh
anassh=ssh(startperiod:startperiod+maxperiod-1);
anatime=time(startperiod:startperiod+maxperiod-1);
% t_tide(XIN,INTERVAL,START_TIME,LATITUDE,RAYLEIGH)
double(time) + datenum('1700-01-01 00:00:00');
% disp('starting time :');  % 152490, 25-May-2002 17:00:00
% datestr(ans(startperiod))
% [tname,tfreq,tcon,tout] =t_tide(anassh(1:8760), ...
%                                 'interval', 1, ...
%                                 'start time', [2002,5,25,17,0,0], ...
%                                 'latitude', lat, ...
%                                 'rayleigh', 1);
% plot ((double(anatime(1:8760))+datenum('1700-01-01 00:00:00')),anassh(1:8760));
% datetick('x','mmmm')
disp('starting time :'); % 157777, 01-Jan-2003 00:00:00
datestr(time(157777)+datenum('1700-01-01 00:00:00'))
[tname,tfreq,tcon,tout] =t_tide(anassh(5289:14048), ...
                                'interval', 1, ...
                                'start time', [2003,1,1,0,0,0], ...
                                'latitude', lat, ...
                                'rayleigh', 1);
                            
tsnr=(tcon(:,1)./tcon(:,2)).^2;  % signal to noise ratio
% as long as the SNR > 10; and is
% probably not bad for SNR as low as 2 or 3. The nonlinear
% procedure gives similar results to the linearized
% procedure at high SNR, and is more accurate at low
% SNR.        
         

% % %  figure 2, 1year(2003)
plot ((double(anatime(5289:14044))+datenum('1700-01-01 00:00:00')),anassh(5289:14044));
datetick('x','mmmm')
title('Sea level at the Juan Fernandez-B (2003)');
ylabel('Sea level (cm)');
xlabel('Month');

% % %  figure 2_2, 1month(Jan, 2003)
plot ((double(anatime(5289:6032))+datenum('1700-01-01 00:00:00')),anassh(5289:6032));
datetick('x','dd')
title('Sea level at the Juan Fernandez-B (January, 2003)');
ylabel('Sea level (cm)');
xlabel('Day');

% % % Figure 3
fsig=tcon(:,1)>tcon(:,2); % Significant peaks
% fsig=tsnr>2; % Significant peaks
semilogy([tfreq(~fsig),tfreq(~fsig)]',[.0005*ones(sum(~fsig),1),tcon(~fsig,1)]','.-r');
line([tfreq(fsig),tfreq(fsig)]',[.0005*ones(sum(fsig),1),tcon(fsig,1)]','marker','.','color','b');
line(tfreq,tcon(:,2),'linestyle',':','color',[0 .5 0]);
% set(gca,'ylim',[.0005 1],'xlim',[0 .5]);
xlabel('frequency (cph)');
text(tfreq,tcon(:,1),tname,'rotation',45,'vertical','base');
ylabel('Amplitude (cm)');
text(.2,50,'Analyzed lines with 95% significance level', 'Fontsize', 15);
text(.2,20,'Significant Constituents','color','b', 'Fontsize', 15);
text(.2,10,'Insignificant Constituents','color','r', 'Fontsize', 15);
text(.2,5,'95% Significance Level','color',[0 .5 0], 'Fontsize', 15);
title('Amplitude of tidal constituents');

% % % FIgure 4
errorbar(tfreq(~fsig),tcon(~fsig,3),tcon(~fsig,4),'.r');
hold on;
errorbar(tfreq(fsig),tcon(fsig,3),tcon(fsig,4),'o');
% text(tfreq(find(fsig),1),tcon(find(fsig),1),tname(find(fsig),:,1),'rotation',45);
hold off;
set(gca,'ylim',[-45 360+45],'xlim',[0 .5],'ytick',[0:90:360]);
xlabel('frequency (cph)');
ylabel('Greenwich Phase (deg)');
text(.34,330,'Analyzed Phase angles with 95% CI', 'Fontsize', 15);
text(.34,290,'Significant Constituents','color','b', 'Fontsize', 15);
text(.34,250,'Insignificant Constituents','color','r', 'Fontsize', 15);
title('Phase of tidal constituents');

disp('maximum value')
max(anassh(5289:14044))
disp('minimum value')
min(anassh(5289:14044))