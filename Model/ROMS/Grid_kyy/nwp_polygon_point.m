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

% % southern Northwest Pacific polygon
snwppolygon = ...
    [120.0, 15.0;
     120.0, 20.0;
     164.0, 20.0;
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
 
 nes2polygon = ...
    [127, 39;
     132, 44;
     140, 45;
     142.5, 45;
     142.3, 45;
     142, 45
     142, 45;
     142, 43;
     141, 43;
     141, 42.8;
     140.2, 42.6;
     140.2, 42.2;
     140.4, 41.8;
     140.5, 39];
 
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
    kuro2polygon = ...
       [142.5, 32.0;
        142.5, 42.0
        147.5, 42.0;
        147.5, 32.0];

% % EKWC for masking (for DA)
    ekwcpolygon = ...
       [128.0, 40.0;
        133.0, 40.0
        133.0, 35.0;
        128.0, 35.0];
    ekwc2polygon = ...
       [127.0, 42.0;
        133.0, 42.0
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

    ekb2polygon = ...
       [129.5, 37.5;
        129.5, 41.0;
	127.0, 41.0;
	127.0, 37.5];

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

% %      pollock_egg polygon
pollock_egg4polygon = ...
    [126.1, 34.0;
    126.1, 41.0;
    133.0, 41.0;
    133.0, 34.0];

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
    131.65, 38.62;
    128.35, 38.62];

% SK_EEZ_polygon = ...
%     [128.5, 34.0;
%     128.8, 34.4;
%     129, 34.5;
%     129, 34.7;
%     129.5, 35;
%     129.7, 35.2;
%     130.25, 35.2;
%     130.5, 35.5;
%     131, 36;
%     131.25, 36.2;
%     131.42, 36.0;
%     131.75, 35.6;
%     132.25, 36;
%     132.25, 36.3;
%     132.9, 36.9;
%     135.5, 36.9;
%     135.5, 38.62;
%     134.2, 39.8;
%     131.65, 38.62;
%     131.65, 38.25;
%     128.6, 38.25]; % exclude prohibited area 


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


% %  SK_coastal (Most of islands are excluded, 
% %  but islands connected by bridge are not excluded.)
% % 'Gyeonggi-do', 'Chungcheongnam-do', 'Jeollabuk-do', 'Jeollanam-do', ...
% % 'Gyeongsangnam-do', 'Gyeongsangbuk-do', 'Gangwon-do', 'Jeju-do', ...
% % 'Ulleung-do', 'Dok-do', 'Ulleng-Dok'
SK_coastal_polygon = ...
    [126.34, 37.84; %%% Gyeonggi-do (northwestern edge)
    126.12, 37.7;
    126.36, 37.62;
    126.36, 37.46;
    126.40, 36.37;
    126.42, 37.23;
    126.35, 37.11;
    126.36, 37.08; %%% Gyeonggi-do (southwestern edge),'Chungcheongnam-do' (northwestern edge)
    126.00, 36.68;
    126.40, 36.08; 
    126.50, 36.08; 
    126.58, 36.06;
    126.60, 35.99; %%% 'Chungcheongnam-do' (southwestern edge), 'Jeollabuk-do' (northwestern edge)
    126.65, 35.99;
    126.45, 35.98;
    126.24, 35.86;
    126.38, 35.60;
    126.37, 35.45; %%% 'Jeollabuk-do' (southwestern edge), 'Jeollanam-do' (northwestern edge)
    126.22, 35.23;
    126.25, 35.20;
    125.98, 35.12;
    125.55, 34.63;
    125.98, 33.99; %%% 'Jeollanam-do' (southwestern edge)
    126.09, 33.26; %%% 'Jeju-do', (western edge)
    126.57, 32.11; %%% 'Jeju-do', (southern edge)
    127.10, 33.44; %%% 'Jeju-do', (eastern edge)
    127.10, 34.26; %%% 'Jeollanam-do' (southern edge)
    127.87, 34.37; %%% 'Jeollanam-do' (southeastern edge), 'Gyeongsangnam-do' (southwestern edge)
    128.88, 34.75;
    129.65, 35.59;
    129.49, 35.63; %%% 'Gyeongsangnam-do' (northeastern edge), 'Gyeongsangbuk-do' (southeastern edge)
    129.69, 36.13;
    129.50, 37.18; %%% 'Gyeongsangbuk-do' (northeastern edge) 'Gangwon-do' (southeastern edge)
    130.93, 36.02; %%% 'Ulleung-do' (southern edge)
    133.04, 37.63;
    128.40, 38.64; %%% 'Gangwon-do' (northwestern edge)
    126.68, 37.85; 
    126.61, 37.77];

adm_div_all_polygon = ...
    [126.66, 37.79; %%% Han river estuary
    126.34, 37.84; 
    126.19, 37.81;
    126.19, 37.73;
    126.12, 37.7;
    125.70, 37.70;
    125.35, 37.70;
    125.15, 37.60;
    124.95, 37.65;
    124.75, 37.80;
    124.75, 38.00;
    124.50, 38.15;
    123.35, 38.15;  %%% Gyeonggi-do (northwestern edge)
    123.50, 37.70;  
    123.80, 37.08;  %%% Gyeonggi-do (southwestern edge), 'Chungcheongnam-do' (northwestern edge)
    123.90, 36.70;
    123.75, 35.99;  %%% 'Chungcheongnam-do' (southwestern edge), 'Jeollabuk-do' (northwestern edge)
    123.73, 35.80;
    123.30, 35.45;  %%% 'Jeollabuk-do' (southwestern edge) 'Jeollanam-do' (northwestern edge)
    122.55, 34.50;
    122.70, 33.75;
    122.85, 33.40;
    123.60, 32.65;
    123.95, 32.45;
    126.05, 30.35; %%% 'Jeju-do' (southern edge)
    126.80, 31.00;
    127.30, 32.35;
    127.80, 33.25;
    128.95, 34.30;
    129.10, 34.80;
    129.49, 35.00; %%% 'Gyeongsangnam-do' (southeastern edge),
    130.25, 35.15;
    130.55, 35.63; %%% 'Gyeongsangnam-do' (northeastern edge), 'Gyeongsangbuk-do' (southeastern edge)
    131.20, 36.10;
    131.40, 36.10;
    131.55, 36.05;
    132.50, 36.70;
    132.95, 37.18; %%% 'Gyeongsangbuk-do' (northeastern edge) 'Gangwon-do' (southeastern edge)
    134.25, 38.50; 
    134.20, 38.62; %%% 'Gangwon-do' (northwestern edge)
    128.36, 38.62; %%% 'Gangwon-do' (land, northwestern edge)
    128.17, 37.92];

adm_div_YS_polygon = ...
    [126.66, 37.79; %%% Han river estuary
    126.34, 37.84; 
    126.19, 37.81;
    126.19, 37.73;
    126.12, 37.7;
    125.70, 37.70;
    125.35, 37.70;
    125.15, 37.60;
    124.95, 37.65;
    124.75, 37.80;
    124.75, 38.00;
    124.50, 38.15;
    123.35, 38.15;  %%% Gyeonggi-do (northwestern edge)
    123.50, 37.70;  
    123.80, 37.08;  %%% Gyeonggi-do (southwestern edge), 'Chungcheongnam-do' (northwestern edge)
    123.90, 36.70;
    123.75, 35.99;  %%% 'Chungcheongnam-do' (southwestern edge), 'Jeollabuk-do' (northwestern edge)
    123.73, 35.80;
    123.30, 35.45;  %%% 'Jeollabuk-do' (southwestern edge) 'Jeollanam-do' (northwestern edge)
    122.55, 34.50;
    122.70, 33.75;
    122.85, 33.40;
    123.60, 32.65;
    123.95, 32.45;
    124.05, 32.25;
    126.52, 34.29; %%% Haenam tower
    128.36, 38.62]; %%% 'Gangwon-do' (land, northwestern edge)

adm_div_SS_polygon = ...
    [126.52, 34.29; %%% haenam tower
    124.05, 32.25;
    126.05, 30.35; %%% 'Jeju-do' (southern edge)
    126.80, 31.00;
    127.30, 32.35;
    127.80, 33.25;
    128.95, 34.30;
    129.10, 34.80;
    129.49, 35.00; %%% 'Jeollanam-do' (southeastern edge), 'Gyeongsangnam-do' (southeastern edge),
    129.17, 35.15; %%% Haeundae Moontan hill
    128.36, 38.62]; %%% 'Gangwon-do' (land, northwestern edge)

adm_div_ES_polygon = ...
    [129.17, 35.15; %%% Haeundae Moontan hill
    129.49, 35.00; %%% 'Jeollanam-do' (southeastern edge), 'Gyeongsangnam-do' (southeastern edge),
    130.25, 35.15;
    130.55, 35.63; %%% 'Gyeongsangnam-do' (northeastern edge), 'Gyeongsangbuk-do' (southeastern edge)
    131.20, 36.10;
    131.40, 36.10;
    131.55, 36.05;
    132.50, 36.70;
    132.95, 37.18; %%% 'Gyeongsangbuk-do' (northeastern edge) 'Gangwon-do' (southeastern edge)
    134.25, 38.50; 
    134.20, 38.62; %%% 'Gangwon-do' (northwestern edge)
    128.36, 38.62]; %%% 'Gangwon-do' (land, northwestern edge)
    
    
adm_div_GGD_polygon = ...
    [126.66, 37.79; %%% Han river estuary
    126.34, 37.84; 
    126.19, 37.81;
    126.19, 37.73;
    126.12, 37.7;
    125.70, 37.70;
    125.35, 37.70;
    125.15, 37.60;
    124.95, 37.65;
    124.75, 37.80;
    124.75, 38.00;
    124.50, 38.15;
    123.35, 38.15;  %%% Gyeonggi-do (northwestern edge)
    123.50, 37.70;  
    123.80, 37.08;  %%% Gyeonggi-do (southwestern edge), 'Chungcheongnam-do' (northwestern edge)
    127.90, 37.08];  %%% Choongju city (center)

adm_div_CCND_polygon = ...
    [127.90, 37.08;  %%% Choongju city (center)
    123.80, 37.08;  %%% Gyeonggi-do (southwestern edge), 'Chungcheongnam-do' (northwestern edge)
    123.90, 36.70;
    123.75, 35.99;  %%% 'Chungcheongnam-do' (southwestern edge), 'Jeollabuk-do' (northwestern edge)
    126.73, 35.99; 
    127.20, 36.17];
    
adm_div_JBD_polygon = ...
    [126.73, 35.99; 
    123.75, 35.99;  %%% 'Chungcheongnam-do' (southwestern edge), 'Jeollabuk-do' (northwestern edge)
    123.73, 35.80;
    123.30, 35.45;  %%% 'Jeollabuk-do' (southwestern edge) 'Jeollanam-do' (northwestern edge)
    127.08, 35.45];
    
% adm_div_JND_polygon = ...   % Chuja island is included
%     [127.08, 35.45;
%     123.30, 35.45;  %%% 'Jeollabuk-do' (southwestern edge) 'Jeollanam-do' (northwestern edge)
%     122.55, 34.50;
%     122.65, 34.15;
%     125.05, 34.15;
%     126.20, 33.90;
%     126.90, 34.05;
%     127.81, 34.05;
%     127.81, 35.45];

adm_div_JND_polygon = ...   % considering 30 km range
    [127.08, 35.45;
    123.30, 35.45;  %%% 'Jeollabuk-do' (southwestern edge) 'Jeollanam-do' (northwestern edge)
    122.55, 34.50;
    122.65, 33.69;
    125.05, 33.69;
    126.20, 33.69;
    126.90, 33.78;
    127.81, 33.78;
    127.81, 35.45];    

% adm_div_JJD_polygon = ...
%     [122.65, 34.15;
%     125.05, 34.15;
%     126.20, 33.90;
%     126.90, 34.05;
%     127.81, 34.05;
%     127.80, 33.25;
%     127.30, 32.35;
%     126.80, 31.00;    
%     126.05, 30.35; %%% 'Jeju-do' (southern edge)
%     123.95, 32.45;
%     123.60, 32.65;
%     122.85, 33.40;
%     122.70, 33.75];


adm_div_JJD_polygon = ...   % considering 30 km range
    [122.65, 33.69;
    125.05, 33.69;
    126.20, 33.69;
    126.90, 33.78;
    127.81, 33.78;
    127.80, 33.25;
    127.30, 32.35;
    126.80, 31.00;    
    126.05, 30.35; %%% 'Jeju-do' (southern edge)
    123.95, 32.45;
    123.60, 32.65;
    122.85, 33.40;
    122.70, 33.78];


adm_div_GSND_polygon = ...
    [127.81, 35.63;
    127.80, 33.25;
    128.95, 34.30;
    129.10, 34.80;
    129.49, 35.00; %%% 'Gyeongsangnam-do' (southeastern edge),
    130.25, 35.15;
    130.55, 35.63]; %%% 'Gyeongsangnam-do' (northeastern edge), 'Gyeongsangbuk-do' (southeastern edge)

adm_div_GSBD_polygon = ...
    [127.81, 35.63;
    130.55, 35.63;
    131.20, 36.10;
    131.40, 36.10;
    131.55, 36.05;
    132.50, 36.70;
    132.95, 37.18; %%% 'Gyeongsangbuk-do' (northeastern edge) 'Gangwon-do' (southeastern edge)
    127.81, 37.18]; %%% 'Gyeongsangnam-do' (northeastern edge), 'Gyeongsangbuk-do' (southeastern edge)

adm_div_GSBD_coastal_polygon = ...
    [127.81, 35.63;
    130.20, 35.63;
    130.20, 37.18;
    127.81, 37.18]; %%% 'Gyeongsangnam-do' (northeastern edge), 'Gyeongsangbuk-do' (southeastern edge)

adm_div_ULD_polygon = ...
    [130.20, 35.63;
    130.20, 38.62;
    134.25, 38.62;
    132.95, 37.18;
    132.50, 36.70;
    131.55, 36.05;
    131.40, 36.10;
    131.20, 36.10]; %%% 'Gyeongsangnam-do' (northeastern edge), 'Gyeongsangbuk-do' (southeastern edge)

adm_div_ULD_only_polygon = ...
    [130.20, 35.63;
    130.20, 38.62;
    131.40, 38.62;
    131.40, 36.10]; %%% 'Gyeongsangnam-do' (northeastern edge), 'Gyeongsangbuk-do' (southeastern edge)

adm_div_DD_only_polygon = ...
    [131.40, 36.10;
    131.55, 36.05;
    132.50, 36.70;
    132.95, 37.18; %%% 'Gyeongsangbuk-do' (northeastern edge) 'Gangwon-do' (southeastern edge)
    134.25, 38.50; 
    134.20, 38.62;
    131.40, 38.62]; %%% 'Gyeongsangnam-do' (northeastern edge), 'Gyeongsangbuk-do' (southeastern edge)

adm_div_GWD_polygon = ...
    [130.20, 37.18; %%% 'Gyeongsangbuk-do' (northeastern edge) 'Gangwon-do' (southeastern edge)
    130.20, 38.50; 
    130.20, 38.62; %%% 'Gangwon-do' (northwestern edge)
    128.36, 38.62; %%% 'Gangwon-do' (land, northwestern edge)
    128.17, 37.92;
    127.80, 37.92;
    127.80, 37.18]; 

% %    reference salinity area in Isobe et al. (2002)
ref_sal_polygon = ...
    [122.2, 23.8;
    122.2, 25.0;
    123.4, 25.0;
    123.4, 23.8];

%% stleedrifter
stlee_drifter_polygon = ...
    [121, 24;
    122, 24.5;
    124, 25;
    126, 26;
    129, 28;
    131, 34;
    122, 34;
    122, 30;
    120, 28;
    120, 27];
