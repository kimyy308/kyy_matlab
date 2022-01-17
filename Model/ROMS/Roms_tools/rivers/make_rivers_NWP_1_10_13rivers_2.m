 %=======================================================================
% Make river Changjiang & Yellow river
%  add 10 rivers (rivers in the bohai sea and the Korea coast)
%                            editted by C.S KIM (10 Jan. 2012)
%    Updated 11-May-2018 by Yong-Yub Kim (makes possible to running at R2017b ~)
%=======================================================================
%  River.flag can have any of the following values:
%             = 0,  All Tracer source/sink are off.
%             = 1,  Only temperature is on.
%             = 2,  Only salinity is on.
%             = 3,  Both temperature and salinity are on. 
% 
%  + tumen river transport
% 
clear all
close all
clc
run('./../Run/start_damo');
romstools_param
testname='test11';
system_name=computer;
if (strcmp(system_name,'PCWIN64'))
    % % for windows
    outputdir=['E:\Data\Model\ROMS\nwp_1_10\input\',testname,'\'];
elseif (strcmp(system_name,'GLNXA64'))
    % % for linux
    outputdir=['/data1/kimyy/Model/ROMS/roms_nwp/nwp_1_10/input/',testname,'/'];
end
% addpath e:\working\roms\roms_tool\data_analysis\
warning off
% vert_param;

% season_cycle=[9650 9257.5 11632.5 21150 32450 41425 ...
% 	      47300 36925 33025 27000 19350 11215];

% year=[1993:2017];  %% if year == 0, make climate data!
year=[1983:2018];  %% if year == 0, make climate data!

year2 = year;



 cycle=425.25;
 time=[15:30:405];
% cycle for 425.25

% riv = load('cal_discharge_all.txt');
% riv = load('cal_discharge_all9210.txt');
riv = load('datong_1961_2019.txt');


% riv2= riv(12:end-1,:);
riv2= riv ;
riv3_yr= reshape(riv2(:,1),12,length(riv2)/12);
riv3_riv= reshape(riv2(:,2),12,length(riv2)/12);

riv3_yr(13,:)=riv3_yr(1,:);
riv3_yr(14,:)=riv3_yr(12,:);
riv3_riv(13,:)=riv3_riv(1,:);
riv3_riv(14,:)=riv3_riv(12,:);


if (year==0) % % make climate data;
    riv3_riv =  mean(riv3_riv,2);
    riv3_yr = repmat(0,14,1);
end


 for i=1:length(year)
% for i=1:1
disp(' ')
disp(['Year is ',num2str(year(i))])
disp(' ')

% cycle=yeardays(year(i));
% time=[1/24:1/12:23/24]*cycle;

idx = find(floor(riv3_yr(1,:)) == year(i));
dis = riv3_riv(:,idx);

avgflow=mean(dis); %25031.7
season_cycle_CRD = dis./avgflow;
sizedis = length(dis);

r = 1;
Name = 'ChangJiang+Huai';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 120.0;
River.lat(r) = 31.9;  % raw
% River.lon(r) = 117.7;
% River.lat(r) = 30.9; %edyz3
River.lon(r) = 119.5;
River.lat(r) = 32; % nwps12
River.flag(r) = 3;
filen='RivDis';
River.file(r,1:length(filen))=filen;
River.area(r)=1794+244;
River.pfact(r)=1;
River.flow_mean(r)=mean(dis);
River.trans(r,1:14)=dis;

season_cycle=[ 378  296  299  209  216  271  907  1481  1345  1084  619  440  629  378  629];
% % % season_cycle_CRD= [226 207 355 438 638 791 1208 1227 1156 756 559 391 226 391];93-17
% % % season_cycle_CRD= [258 270 432 528 750 868 1266 1023 891 870 620 268 258 268];83-18

avgflow=mean(season_cycle); %25031.7
season_cycle=season_cycle./avgflow;

r = 2;
Name = 'HuangHe';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 118.5;
River.lat(r) = 37.5;
River.flag(r) = 3;
filen='RivDis';
River.file(r,1:length(filen))=filen;
River.area(r)=894;
River.pfact(r)=1;
River.flow_mean(r)=  avgflow ;
River.trans(r,1:14)= season_cycle_CRD.*River.flow_mean(r);
% River.trans(r,1:14)= season_cycle.*River.flow_mean(r);

% r3 ~r6 
% Seasonal variations of the Yellow River plume in the Bohai Sea: A model study
% WangQ_2008JGR [Qiang Wang,1 Xinyu Guo,1 and Hidetaka Takeoka1]


r = 3;
Name = 'HaiheRiv.';
season_cycle=[12  19  22  13  8  18  283  651 273  66  53  27 12 27] ; % 
avgflow=mean(season_cycle); season_cycle=season_cycle./avgflow;
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 117.75; %126.7;
River.lat(r) = 39;
River.flag(r) = 3;
filen='RivDis';
River.file(r,1:length(filen))=filen;
River.area(r)=30;
River.pfact(r)=1;
River.flow_mean(r)= avgflow ;
River.trans(r,1:14)=season_cycle.*River.flow_mean(r);

r = 4;
Name = 'LuanheRiv.';
season_cycle=[3 25 37 33 12 23 341 467 134 76 63 37 3 37] ; % 
avgflow=mean(season_cycle); season_cycle=season_cycle./avgflow;
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 119.2; %126.7;
River.lat(r) = 39.5;
River.flag(r) = 3;
filen='RivDis';
River.file(r,1:length(filen))=filen;
River.area(r)=30;
River.pfact(r)=1;
River.flow_mean(r)= avgflow ;
River.trans(r,1:14)=season_cycle.*River.flow_mean(r);
 
r = 5;
Name = 'LiaoheRiv.';
season_cycle=[36  40 100 133 134 193 559 1144 513 279 158 90 36 90] ; % 
avgflow=mean(season_cycle); season_cycle=season_cycle./avgflow;
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 121.9; %121.8;
River.lat(r) = 41; %40.9
River.flag(r) = 3;
filen='RivDis';
River.file(r,1:length(filen))=filen;
River.area(r)=30;
River.pfact(r)=1;
River.flow_mean(r)= avgflow ;
River.trans(r,1:14)=season_cycle.*River.flow_mean(r);

r = 6;
Name = 'YaluRiv.';
season_cycle=[ 656 657 709 693 695 828 1353 2117 1158 699 678 685 656 685] ; % 
avgflow=mean(season_cycle); season_cycle=season_cycle./avgflow;
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 124.4; %126.7;
River.lat(r) = 40.2; %40.0
River.flag(r) = 3;
filen='RivDis';
River.file(r,1:length(filen))=filen;
River.area(r)=30;
River.pfact(r)=1;
River.flow_mean(r)= avgflow ;
River.trans(r,1:14)=season_cycle.*River.flow_mean(r);


r = 7;
Name = 'TaedongRiv.';
season_cycle=[ 25.43 23.21 36.86 89.57 84.86 181.36  513.79  800.29  259.71  42.36  48.79 34.57 25.43 34.57] ; % clim
avgflow=mean(season_cycle); season_cycle=season_cycle./avgflow;
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 125.8; %126.7;
River.lat(r) = 38.8;
River.flag(r) = 3;
filen='RivDis';
River.file(r,1:length(filen))=filen;
River.area(r)=30;
River.pfact(r)=1;
River.flow_mean(r)= avgflow   ;
River.trans(r,1:14)=season_cycle.*River.flow_mean(r);


r = 8;
Name = 'HanRiv.';
% season_cycle=[129.9 133.9 133.2  204.6   389.1  531.6 ...
%                4705.6  916.9  217.2  144.0  153.9  150.3]*1.2;    % Pal Dang Dam 2007
season_cycle=[ 137.79  150.58  240.04  482.46  302.92  426.28 ...
               1738.92 1268.71 920.92  273.36  214.92 148.72 137.79 148.72] ; % clim
           season_cycle(7:8) = season_cycle(7:8) + 500 ;
avgflow=mean(season_cycle); season_cycle=season_cycle./avgflow;
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 127.4; %126.7;
River.lat(r) =  37.4;
River.flag(r) = 3;
filen='RivDis';
River.file(r,1:length(filen))=filen;
River.area(r)=30;
River.pfact(r)=1;
River.flow_mean(r)= avgflow   ;
River.trans(r,1:14)=season_cycle.*River.flow_mean(r);

r = 9;
Name = 'KeumRiv.';
season_cycle=[71 89 100.5 140.5 100.5 155 295 347.5 212.5 103 87.5 82.25 71 82.25] ; % clim
%            season_cycle(7:8) = season_cycle(7:8) + 500 ;
avgflow=mean(season_cycle); season_cycle=season_cycle./avgflow;
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 127.2;  %126.7;
River.lat(r) = 36.00;
River.flag(r) = 3;
filen='RivDis';
River.file(r,1:length(filen))=filen;
River.area(r)=30;
River.pfact(r)=1;
River.flow_mean(r)= avgflow   ;
River.trans(r,1:14)=season_cycle.*River.flow_mean(r);

r = 10;
Name = 'YeongsanRiv.';
season_cycle=[20  36.75  26.75  53.75  32  93  52.25  76.5  54.5  25.25  22.75  21.75 20 21.75] ; % clim
%            season_cycle(7:8) = season_cycle(7:8) + 500 ;
avgflow=mean(season_cycle); season_cycle=season_cycle./avgflow;
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 126.6;%126.6;
River.lat(r) = 34.9;%34.8;
River.flag(r) = 3;
filen='RivDis';
River.file(r,1:length(filen))=filen;
River.area(r)=30;
River.pfact(r)=1;
River.flow_mean(r)= avgflow   ;
River.trans(r,1:14)=season_cycle.*River.flow_mean(r);

r = 11;
Name = 'SumjinRiv.';
season_cycle=[24.25	37.5 30 53 52.75 92.5 82.5 172.5 65  26.25 20.75 15.25 24.25 15.25] ; % clim
%            season_cycle(7:8) = season_cycle(7:8) + 500 ;
avgflow=mean(season_cycle); season_cycle=season_cycle./avgflow;
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 127.56;%126.6;
River.lat(r) = 35.18;%34.8;
River.flag(r) = 3;
filen='RivDis';
River.file(r,1:length(filen))=filen;
River.area(r)=30;
River.pfact(r)=1;
River.flow_mean(r)= avgflow   ;
River.trans(r,1:14)=season_cycle.*River.flow_mean(r);


r = 12;
Name = 'NakdongRiv.';
season_cycle=[61.45	83.75	130.25	218.25	149.2	119.65	608.9	404.9	531.65	188.7	121.35	83.3 61.45 83.3] ; % clim
%            season_cycle(7:8) = season_cycle(7:8) + 500 ;
avgflow=mean(season_cycle); season_cycle=season_cycle./avgflow;
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 128.95;  %126.6;
River.lat(r) = 35.15;   %35.4; 
River.flag(r) = 3;
filen='RivDis';
River.file(r,1:length(filen))=filen;
River.area(r)=30;
River.pfact(r)=1;
River.flow_mean(r)= avgflow   ;
River.trans(r,1:14)=season_cycle.*River.flow_mean(r);

r = 13;
Name = 'TumenRiv.';
season_cycle=[30, 20,	50, 160, 200, 340, 380, 570, 410, 160, 80, 10, 30, 10] ; % clim
%            season_cycle(7:8) = season_cycle(7:8) + 500 ;
avgflow=mean(season_cycle); season_cycle=season_cycle./avgflow;
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 130.6989; 
River.lat(r) = 42.2928;   
River.flag(r) = 3;
filen='RivDis';
River.file(r,1:length(filen))=filen;
River.area(r)=30;
River.pfact(r)=1;
River.flow_mean(r)= avgflow   ;
River.trans(r,1:14)=season_cycle.*River.flow_mean(r);
save eas_rivers_discharge_YS_EJS River

% This script will create a river runoff FORCING NetCDF file.
%
% Paul Goodman's modification of useast_rivers.m (J. Wilkin)
% BJ Choi's modification AUG 2, 2004.

%add_ocean_toolboxes;

IWRITE=1;
IPLOT=0;

%-----------------------------------------------------------------------
% specify the input files
%-----------------------------------------------------------------------

% roms grid
% grdname  = 'D:\add2_ini_bry_grd\grid\roms_grid2_ADD_08_2_ep.nc'; % new 1/10
% grdname= 'D:\MEPL\project\NWP\make_init\grid\roms_grid_combine2_test34.nc';
% grdname is read from romstools_param;

% grdname = ['/data1/kimyy/Model/ROMS/roms_nwp/nwp_1_10/input/',testname,'/spinup/roms_grid_nwp_1_10_',testname,'.nc'];
grdname = ['E:\Data\Model\ROMS\nwp_1_10\input\',testname,'\roms_grid_nwp_1_10_',testname,'.nc'];

grd_file = grdname ;

N = 40; % Number of vertical levels


% the new netcdf forcing file - river data will be appended
%Fname = 'E:\Working\roms\1_4\input\ecco\roms_eas4_river.nc';

head=[outputdir,'roms_river_nwp_1_10_'];
mid=num2str(year(i));
foot='.nc';
if (year==0)
    Fname = [head,'climate',foot];
else
    Fname = [head,mid,foot];
end

% load the river flow data structure
Rname = 'eas_rivers_discharge_YS_EJS';
load(Rname)
if ~exist('River')
  error([ 'Structure ''River'' does not exist in ' Rname])
end

% get the roms grid
% we need this to find the grid indices corresponding to the lon/lat
% locations of the dat

% scoord = [7.0 2.0 250 40];
grd = roms_get_grid(Vtransform, Vstretching, grd_file, [theta_s,theta_b,hc,N]);

%-----------------------------------------------------------------------
%  Set river name, location, and runoff direction.
%-----------------------------------------------------------------------
%
%  Notice that for now the river identification number is that of the
%  river count.  This variable is needed to fool ROMS generic IO
%  interphase; it can be any real number since will never used.
%
%  River.lon  = nominal longitude of river mouth (this routine finds i,j)
%  River.lat  = nominal latitude of river mouth (this routine finds i,j)
%
%  River.dir  = 0 for flow entering cell through u-direction face
%             = 1 for flow entering cell through v-direction face
%
%  River.Xpos,Ypos = the i,j index of the u-point or v-point location on 
%             the ROMS C-grid thyat defines the face of the cell through
%             which the river enters. 
%             The river input is added to the appropriate ubar, vbar, u, 
%             v, and tracer fluxes on this face.
%
%  River.sens = 1 for river flow in the positive u or v direction
%             = -1 for river flow in the negative u or v direction, 
%             This factor multiplies the river flow rate Q (>0) (see
%             sources.h for why this is done)
%
%  River.flag can have any of the following values:
%             = 0,  All Tracer source/sink are off.
%             = 1,  Only temperature is on.
%             = 2,  Only salinity is on.
%             = 3,  Both temperature and salinity are on. 
%
%  This is documented further in ROMS sources.h

Nrivers = length(River.lon);

%-----------------------------------------------------------------------
% Find the ij coordinates of the river lat/lon points
%-----------------------------------------------------------------------

for r = 1:Nrivers
  
  River.num(r) = r;

  % The logic below checks whether the u-face or v-face is the better 
  % approximation to the nominal location of the river mouth, and assigns 
  % River.dir accordingly. 
  % The sign of diff(rhomask) is used to determine whether the flow 
  % enters from the NSEW

  % case flow enters cell through u-direction face
  glon = grd.lon_u;
  glat = grd.lat_u;
  % keep only u-faces that are coastline
  drhomask = diff(grd.mask_rho')';
  notcoast = find(drhomask==0);
  glon(notcoast) = NaN;
  glat(notcoast) = NaN;
  % find the u-face closest to the river mouth
  [J,I,du] = closest(glon,glat,River.lon(r),River.lat(r));

  % assume this is best option until we test v-face
  River.Xpos(r) = I;
  River.Ypos(r) = J-1; % ROMS u points start from index j=0
  River.dir(r)  = 0;
  River.sens(r) = drhomask(J,I);
  River.glon(r) = glon(J,I);
  River.glat(r) = glat(J,I);

  % case flow enters cell through v-direction face
  glon = grd.lon_v;
  glat = grd.lat_v;
  % keep only v-faces that are coastline
  drhomask = diff(grd.mask_rho);
  notcoast = find(drhomask==0);
  glon(notcoast) = NaN;
  glat(notcoast) = NaN;
  % find the v-face closest to the river mouth
  [J,I,dv] = closest(glon,glat,River.lon(r),River.lat(r));
  
  if dv < du 
    % overwrite because v-face result is closer to river mouth
    River.Xpos(r) = I-1; % ROMS v points start from index i=0
    River.Ypos(r) = J;   
    River.dir(r)  = 1;
    River.sens(r) = drhomask(J,I);
    River.glon(r) = glon(J,I);
    River.glat(r) = glat(J,I);
  end
 
end

if IPLOT
  clf
  % plot the original and selected lon/lat to check that the lookup was done 
  % sensibly
  pcolorjw(grd.lon_rho,grd.lat_rho,grd.mask_rho_nan)
% amerc;
  % shading faceted
% plotnenacoast(2)
  hold on
  han = plot(River.lon,River.lat,'r^');
  set(han,'markersize',10,'MarkerFaceColor',get(han,'color'))
  han = plot(River.glon,River.glat,'bd');
  set(han,'markersize',10,'MarkerFaceColor',get(han,'color'))
  han = plot([River.lon; River.glon],[River.lat; River.glat],'r-');
  set(han,'linewidth',2)
  for r=1:length(River.lon)
    han = text(River.lon(r),River.lat(r),...
	[ ' (' int2str(r) ') ' deblank(River.Name(r,:))]);
    set(han,'fontsize',12)
  end  
  hold off
end

for r=1:Nrivers  
  River.vshape(r,1:N)=1/N;
  River.trans(r,:) = River.sens(r)*River.trans(r,:);
end

  detailstr = ['The ' int2str(Nrivers) ' Rivers are : '];
for r=1:Nrivers;  
  detailstr = [detailstr int2str(r) '. ' strcat(River.Name(r,:)) ', '];
end
  detailstr = [detailstr(1:end-2) '.'];
  

%  River.time is evenly spaced through the year of 365.25 days
River.time =time;
cycle=425.25;
River.time = [15:30:405];
River.time_units = 'days';


% temp/salt get used, or not, according to the value of River.flag 
salt = 0;
%temp = 0;
Ytemp_max=28;
Ytemp_min=8;
Yamp=(Ytemp_max-Ytemp_min)/2;
Htemp_max=25;
Htemp_min=3;
Hamp=(Htemp_max-Htemp_min)/2;
factor=2*pi/12;

 for i=1:length(River.time)
 Htemp(i)=Hamp*( sin(factor*(i-5)) + 1 )+Htemp_min;
 Ytemp(i)=Yamp*( sin(factor*(i-5)) + 1 )+Ytemp_min;
 end


% for i=1:length(River.time)
%   River.temp(1,:,i)=ones([1 N]).*Ytemp(i);
%   River.temp(2,:,i)=ones([1 N]).*Htemp(i);
%   River.salt(:,:,i)=ones([Nrivers N]).*salt;
% end

% aaaaaaaaaaaaaaaaaa

Rtemp_max = [28.23 27.52  27.19 27.09 26.86 27    27.23 27.54 27.83  28.02  27.97 27.97 20.49] ;
Rtemp_min = [6.01  0.12  -1.48  -2.01 -3.51 -2.55 -1.27 0.23  1.73    2.9   2.6   2.63 1.09] ;

% % Tuman river temp(OISST) -> [1.51, 1.09, 1.68, 3.90, 7.87, 13.00, 17.54, 20.49, 18.91, 14.65, 9.37, 4.32]

for i=1:length(River.time)
    for r=1:Nrivers 
        %   River.temp(1,:,i)=ones([1 N]).*Ytemp(i);
        Htemp_max= Rtemp_max(r) ;
        Htemp_min= Rtemp_min(r) ;
        Hamp=(Htemp_max-Htemp_min)/2;
        factor=2*pi/12;
        if (i<=12)
            Htemp(i)=Hamp*( sin(factor*(i-5)) + 1 )+Htemp_min;
            River.temp(r,:,i)=ones([1 N]).*Htemp(i);
        elseif (i==13)
            Htemp(i)=Hamp*( sin(factor*(1-5)) + 1 )+Htemp_min;
            River.temp(r,:,i)=ones([1 N]).*Htemp(13);
        elseif (i==14)
            Htemp(i)=Hamp*( sin(factor*(12-5)) + 1 )+Htemp_min;
            River.temp(r,:,i)=ones([1 N]).*Htemp(14);
        end
    end
  River.salt(:,:,i)=ones([Nrivers N]).*salt;
end
 
%River.salt(41,:,:) = 28.9;

%-----------------------------------------------------------------------
%  Create empty river data FORCING NetCDF file.
%-----------------------------------------------------------------------

  disp([ 'Creating ' Fname '...'])

      create_empty_EAS_rivers_Y

%-----------------------------------------------------------------------
%  Write river data into existing FORCING NetCDF file.
%-----------------------------------------------------------------------

  disp([ 'Appending rivers data to ' Fname '...'])
%  disp('paused...');pause
%  [Vname,status]=wrt_rivers(Fname,River);

%     write_rivers
     
     nc = netcdf(Fname,'write');
     
theVarname = 'river';
nc{theVarname}(:) = River.num;

theVarname = 'river_Xposition';
nc{theVarname}(:) = River.Xpos;

theVarname = 'river_Eposition';
nc{theVarname}(:) = River.Ypos;

theVarname = 'river_direction';
nc{theVarname}(:) = River.dir;

theVarname = 'river_flag';
nc{theVarname}(:) = River.flag;

theVarname = 'river_Vshape';
nc{theVarname}(:,:) = River.vshape;

theVarname = 'river_time';
nc{theVarname}(:) = River.time;

theVarname = 'river_transport';
nc{theVarname}(:,:) = River.trans';

theVarname = 'river_temp';
nc{theVarname}(:,:,:) = permute(River.temp,[3 2 1]);

theVarname = 'river_salt';
nc{theVarname}(:,:,:) = permute(River.salt,[3 2 1]);

result = close(nc);  
 
% time=time+360;     
%time=time+365.25;     
 end

 
  
 
 