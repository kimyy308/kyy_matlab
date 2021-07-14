
%
% mean seasoanl cycle of river flow from ChangJiang River
% from 1976 to 1979. station is Datong, China (117.61E, 30.76N).
% The seasonal cyle is silmar to Fig. 5 of 
% (Dai and Trenberth, J. of Hydro., 2002)
%


%  River.flag can have any of the following values:
%             = 0,  All Tracer source/sink are off.
%             = 1,  Only temperature is on.
%             = 2,  Only salinity is on.
%             = 3,  Both temperature and salinity are on. 


season_cycle=[9650 9257.5 11632.5 21150 32450 41425 ...
	      47300 36925 33025 27000 19350 11215];
  
y1992=[15877.0 12048.0 14160.0 39917.0 28541.0 34698.0 ...
        41284.0 27360.0 24472.0 26016.0 19410.0 13372.0];
y1993=[14951.0 17125.0 22780.0 27089.0 25381.0 37336.0 ...
        39304.0 49309.0 37186.0 26415.0 22558.0 23235.0];
y1994=[12478.0 14863.0 17800.0 18786.0 32471.0 26562.0 ...
        40010.0 41306.0 33391.0 32542.0 23916.0 20490.0];
y1995=[16897.0 17446.0 16498.0 17797.0 35981.0 36067.0 ...
        53968.0 47684.0 38738.0 20251.0 32647.0 12401.0];
y1996=[12748.0 19299.0 11732.0 27175.0 20677.0 31979.0 ...
        53573.0 71507.0 42672.0 24483.0 25156.0 21083.0];
y1997=[10487.0 16134.0 19884.0 17368.0 25595.0 30261.0 ...
        37827.0 43122.0 22475.0 25982.0 19069.0 23962.0];
y1998=[19697.0 22207.0 15089.0 26035.0 37194.0 40313.0 ...
        55575.0 69106.0 42974.0 23043.0 23372.0 10601.0];
y1999=[13780.0 13638.0 10761.0 21968.0 43161.0 43064.0 ...
        61477.0 41113.0 51229.0 30299.0 27379.0 17498.0];
y2000=[12380.0 22674.0 15064.0 18525.0 26115.0 30300.0 ...
        39215.0 33924.0 37024.0 33622.0 29034.0 19227.0];
y2001=[14035.0 20892.0 14925.0 17833.0 33148.0 27650.0 ...
        43563.0 26776.0 34081.0 19971.0 27059.0 16781.0];
y2002=[18154.0 15298.0 17473.0 25617.0 47192.0 45848.0 ...
        54038.0 47352.0 42957.0 20511.0 22123.0 22616.0];

avgflow=mean(season_cycle); %25031.7
season_cycle=season_cycle./avgflow;

r = 1;
Name = 'ChangJiang+Huai';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 120.0;
River.lat(r) = 31.9;
River.flag(r) = 2;
filen='RivDis';
River.file(r,1:length(filen))=filen;
River.area(r)=1794+244;
River.pfact(r)=1;
River.flow_mean(r)=mean(y1992);
River.trans(r,1:12)=y1992;

r = 2;
Name = 'HuangHe';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 118.5;
River.lat(r) = 37.5;
River.flag(r) = 2;
filen='RivDis';
River.file(r,1:length(filen))=filen;
River.area(r)=894;
River.pfact(r)=1;
River.flow_mean(r)=1490.0;
River.trans(r,1:12)=season_cycle.*River.flow_mean(r);

save eas_rivers_discharge River
