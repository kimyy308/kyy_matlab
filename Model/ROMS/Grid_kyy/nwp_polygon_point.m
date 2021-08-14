% % test polygon
testpolygon = ...
    [150.0, 20.0;
     150.0, 20.8;
     150.8, 20.8;
     150.8, 20.0];



% % North Western Pacific polygon
nwppolygon = ...
    [115.0, 15.0;
     115.0, 52.0;
     164.0, 52.0;
     164.0, 15.0];


% % East Sea polygon for masking
espolygon = ...
    [130.2, 33.3;
     129.4, 35.5;
     129, 35.5;
     127, 39;
     127, 42;
     132, 44;
     140, 52;
     142.5, 52;
     142.3, 47;
     142, 46.5;
     142, 45;
     142, 43;
     141, 43;
     141, 42.8;
     140.2, 42.6;
     140.2, 42.2;
     140.4, 41.8;
     140.5, 41;
     140.5, 38;
     137, 36;
     136, 35;
     133, 35;
     131, 34];
 
 espolygon2 = ...
     [129.0, 35.2;
     129.3, 34.9;
     132.3, 37.1;
     132.3, 38.6;
     128.2, 38.6];
 es_khoapolygon = ...
    [130.5, 33.8;
     129.18, 35.16;
     129, 35.5;
     127, 39;
     127, 42;
     132, 44;
     140, 52;
     142.5, 52;
     142.3, 47;
     142, 46.5;
     142, 45;
     142, 43;
     141, 43;
     141, 42.8;
     140.2, 42.6;
     140.2, 42.2;
     140.4, 41.8;
     140.5, 41;
     140.5, 38;
     137, 36;
     136, 35;
     133, 35;
     131, 34];
 
 sespolygon = ...
    [130.2, 33.3;
     129.4, 35.5;
     129, 35.5;
     127, 39;
     127, 40;
     140.5, 40;
     140.5, 38;
     137, 36;
     136, 35;
     133, 35;
     131, 34];
 
nespolygon = ...
    [127, 40;
     132, 44;
     140, 52;
     142.5, 52;
     142.3, 47;
     142, 46.5;
     142, 45;
     142, 43;
     141, 43;
     141, 42.8;
     140.2, 42.6;
     140.2, 42.2;
     140.4, 41.8;
     140.5, 40];
 
 % % South Sea polygon for masking
 sspolygon = ...
     [130.2, 33.3;
      129.4, 35.5;
      125.5, 34.5;
      125.5, 32.6;
      128.6, 32.6];
  
  sspolygon2 = ...
     [125.0, 34.3;
      126.6, 34.75;
      128.0, 35.6;
      129.0, 35.2;
      129.3, 34.9;
      128.6, 34.1;
      126.6, 32.5;
      125.0, 33.3];
  ss_khoapolygon = ...
     [130.5, 33.8;
      129.18, 35.16;
      126.52, 34.29;
      124.0, 31.75];
  
  
  % % Yellow Sea polygon for masking
  yspolygon = ...
      [125.5, 32.6;
       125.5, 34.5;
       126.6, 34.5;
       127, 35;
       127.5, 36;
       127.5, 38;
       126, 40;
       122, 41;
       116, 41;
       116, 38;
       116, 32.6];

   yspolygon2 = ...
      [125.0, 34.3;
       126.6, 34.75;
       127, 35;
       127.5, 36;
       127.5, 38;
       126.7, 37.8;
       126.6, 37.8;
       126.4, 37.8;
       126.1, 37.8;
       126.2, 37.7;
       125.7, 37.7;
       125.5, 37.7];
   
   ys_khoapolygon = ...
      [124.0, 31.75;
       126.52, 34.29;
       127.5, 36;
       127.5, 38;
       126, 40;
       122, 41;
       116, 41;
       116, 38;
       119.0, 34.8;
       120.7, 32.6;
       120.75, 32.45;
       122, 31.75];
   
   ysecs_khoapolygon = ...    % land mask for changjiang river
      [117.0, 30;
       128.4, 30;
       131.0, 32;
       131.0, 33;
       131.0, 34;
       129.18, 35.16;
       129, 35.5;
       126.52, 34.29;
       127.5, 36;
       127.5, 38;
       126, 40;
       122, 41;
       116, 41;
       116, 38;
       119.0, 34.8;
       120.7, 32.6;
       120.75, 32.45;
       122, 31.72;
       121.9, 30.91;
       120, 30.91;
       120, 30;
       116, 30];
   
 % % East China Sea polygon for masking
   ecspolygon = ...
       [116, 32.6;
        126, 32.6
        128.6, 32.6;
        121.5, 25;
        116, 25];
  
 % % East China Sea polygon include south sea for masking
   ecs2polygon = ...
       [126.6, 34.75;
        127, 35;
        129.4, 35.5;
        132, 34;
        131, 34;
        131, 33;
        131.0, 32.0;
        128.4, 30.0; 
        116, 30.0;
        116, 32.6;
        125.5, 32.6;
        125.5, 34.5;
        126.6, 34.5];   
  
% %    Yellow Sea and East China Sea for masking
ysecspolygon = ...
       [127, 35;
        129.4, 35.5;
        132, 34;
        131, 34;
        131, 33;
        131.0, 32.0;
        128.4, 30.0; 
        116, 30.0;
        116, 32.6;
        116, 38;
        116, 41;
        122, 41;
        126, 40;
        127.5, 38;
        127.5, 36]; 
    
    yspolygon = ...
      [125.5, 32.6;
       125.5, 34.5;
       126.6, 34.5;
       127, 35;
       127.5, 36;
       127.5, 38;
       126, 40;
       122, 41;
       116, 41;
       116, 38;
       116, 32.6];
  % % East China Sea Changjiang river discharges path polygon for masking
   ecs_crdpolygon = ...
       [116, 41;
        128.6, 41;
        121.5, 28;
        116, 28];

% % Kuro polygon for masking (for DA)
   kuropolygon = ...
       [135.0, 32.0;
        135.0, 48.0
        150.0, 48.0;
        150.0, 32.0];

% % EKWC for masking (for DA)
   ekwcpolygon = ...
       [128.0, 40.0;
        133.0, 40.0
        133.0, 35.0;
        128.0, 35.0];
 
 % % North Korean coast polygon for walleye pollack
   nkcpolygon = ...
       [128.5, 40.0;
        128.7, 40.0;
	129.9, 41.0;
	128.5, 41.0];

 % % South Korean coast polygon for walleye pollack
   skcpolygon = ...
       [128.5, 38.5;
        128.7, 38.5;
	129.9, 37.5;
	128.5, 37.5];

 % % Wonsan Bay polygon for walleye pollack
   wbpolygon = ...
       [128.4, 38.5;
        128.4, 40.0;
	127.0, 40.0;
	127.0, 38.5];

% % East Korea Bay polygon for walleye pollack
   ekbpolygon = ...
       [130, 38.5;
        130, 41.0;
	127.0, 41.0;
	127.0, 38.5];



% %  Around Korea Peninsula polygon for SSH
   akppolygon = ...
       [117.0, 32.0;
        117.0, 52.0;
        140.0, 52.0;
        140.0, 32.0];
   akp2polygon = ...
       [117.0, 32.0;
        117.0, 52.0;
        142.5, 52;
        142.3, 47;
        142, 46.5;
        142, 45;
        142, 43;
        141, 43;
        141, 42.8;
        140.2, 42.6;
        140.2, 42.2;
        140.4, 41.8;
        140.5, 41;
        140.5, 38;
        137, 36;
        136, 35;
        133, 35;
        132, 34;
        131, 34;
        131, 33;
        131.0, 32.0]; 
    
    akp3polygon = ...
       [117.0, 25.0;
        117.0, 52.0;
        142.5, 52;
        142.3, 47;
        142, 46.5;
        142, 45;
        142, 43;
        141, 43;
        141, 42.8;
        140.2, 42.6;
        140.2, 42.2;
        140.4, 41.8;
        140.5, 41;
        140.5, 38;
        137, 36;
        136, 35;
        133, 35;
        132, 34;
        131, 34;
        131, 33;
        131.0, 32.0
        122.0, 25.0]; 
    
    
    akp4polygon = ...
       [117.0, 30.0;
        117.0, 52.0;
        142.5, 52;
        142.3, 47;
        142, 46.5;
        142, 45;
        142, 43;
        141, 43;
        141, 42.8;
        140.2, 42.6;
        140.2, 42.2;
        140.4, 41.8;
        140.5, 41;
        140.5, 38;
        137, 36;
        136, 35;
        133, 35;
        132, 34;
        131, 34;
        131, 33;
        131.0, 32.0;
        128.4, 30.0]; 
    
% %   Coastal stations polygon for SSH
    capolygon = ...
        [125.0, 38.0;
         125.0, 33.0;
         128.0, 33.0;
         130.5, 36.0;
         130.0, 38.0];

% %   Bohai bay polygon for SSH
    bohpolygon = ...
        [116.0, 41.0;
         123.0, 41.0;
         122.0, 39.3;
         121.0, 37.5;
         120.0, 37.0;
         116.0, 37.0];
     
% %      pollock_egg polygon
pollock_eggpolygon = ...
    [127.0, 36.0;
    127.0, 44.0;
    133.0, 44.0;
    133.0, 36.0];
    
% %      pollock_egg polygon
pollock_egg2polygon = ...
    [126.0, 34.0;
    127.0, 44.0;
    133.0, 44.0;
    133.0, 34.0];

% %      pollock_egg polygon
pollock_egg3polygon = ...
    [127.0, 36.0;
    127.0, 41.0;
    133.0, 41.0;
    133.0, 36.0];

% %      South Korean EEZ (East Sea, real fishing)
SK_EEZ_polygon = ...
    [128.5, 34.0;
    128.8, 34.4;
    129, 34.5;
    129, 34.7;
    129.5, 35;
    129.7, 35.2;
    130.25, 35.2;
    130.5, 35.5;
    131, 36;
    131.25, 36.2;
    131.42, 36.0;
    131.75, 35.6;
    132.25, 36;
    132.25, 36.3;
    132.9, 36.9;
    135.5, 36.9;
    135.5, 38.62;
    134.2, 39.8;
    133.0, 38.62;
    128.3, 38.62];

% %      South Korean EEZ (claimed)
SK_EEZ2_polygon = ...
    [128.36, 38.62;
    133.18, 39.85;
    133.75, 39.73;
    133.81, 39.09;
    132.1, 36.96;
    131.31, 36.21;
    130.57, 35.57;
    130.27, 35.12;
    129.68, 35.07;
    129.37, 34.96;
    129.21, 34.84;
    129.05, 34.67;
    129.01, 34.48;
    128.88, 34.31;
    128.69, 34.14;
    128.43, 33.83;
    128.42, 33.79;
    128.35, 33.75;
    127.87, 33.27;
    127.87, 33.23;
    127.7, 32.95];


% %      North Korean EEZ
NK_EEZ_polygon = ...
    [132.6, 38.7;
    133.0, 39.0;
    133.5, 40.0;
    130.0, 43.0;
    127.0, 40.0;
    127.0, 38.7];

% %      Russian Primorsky EEZ
RU_PR_EEZ_polygon = ...
    [130.0, 43.0;
    133.5, 40.0;
    133.3, 39.7;
    136.0, 40.5;
    136.3, 41.3;
    139.0, 44.5;
    140.0, 46.0;
    141.0, 45.9;
    141.7, 52.0;
    130.0, 52.0];

% %      Russian Sakhalin EEZ
RU_SH_EEZ_polygon = ...
    [141.0, 45.9;
    141.7, 52.0;
    142.5, 52.0;
    142.0, 45.8];


% %      Japan EEZ
JP_EEZ_polygon = ...
    [133.3, 39.7;
    136.0, 40.5;
    136.3, 41.3;
    139.0, 44.5;
    140.0, 46.0;
    141.0, 45.9;
    142.0, 45.8;
    142.0, 43.0;
    141.0, 43.0;
    140.0, 42.5;
    140.0, 36.0;
    137.0, 36.0;
    136.0, 35.0;
    130.0, 35.0;
    133.0, 38.0;
    132.6, 38.7;
    133.0, 39.0];