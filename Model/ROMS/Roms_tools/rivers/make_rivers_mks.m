% Make river 1992 ~ 2002
%
%  River.flag can have any of the following values:
%             = 0,  All Tracer source/sink are off.
%             = 1,  Only temperature is on.
%             = 2,  Only salinity is on.
%             = 3,  Both temperature and salinity are on. 

clear all
close all
clc

% addpath e:\working\roms\roms_tool\data_analysis\

season_cycle=[9650 9257.5 11632.5 21150 32450 41425 ...
	      47300 36925 33025 27000 19350 11215];

year=[1992:2002];

%cycle for 360
 cycle=360;
 time=[15:30:345];
% cycle for 365.25
% cycle=365.25;
% time=[1/24:1/12:23/24]*365.25;

 for i=1:length(year)
% for i=1:1
disp(' ')
disp(['Year is ',num2str(year(i))])
disp(' ')
   switch year(i)
       case 1992
           dis=[15877.0 12048.0 14160.0 39917.0 28541.0 34698.0 ...
                   41284.0 27360.0 24472.0 26016.0 19410.0 13372.0];
       case 1993
           dis=[14951.0 17125.0 22780.0 27089.0 25381.0 37336.0 ...
                   39304.0 49309.0 37186.0 26415.0 22558.0 23235.0];
       case 1994
           dis=[12478.0 14863.0 17800.0 18786.0 32471.0 26562.0 ...
                   40010.0 41306.0 33391.0 32542.0 23916.0 20490.0];
       case 1995
           dis=[16897.0 17446.0 16498.0 17797.0 35981.0 36067.0 ...
                   53968.0 47684.0 38738.0 20251.0 32647.0 12401.0];
       case 1996
           dis=[12748.0 19299.0 11732.0 27175.0 20677.0 31979.0 ...
                   53573.0 71507.0 42672.0 24483.0 25156.0 21083.0];
       case 1997
           dis=[10487.0 16134.0 19884.0 17368.0 25595.0 30261.0 ...
                   37827.0 43122.0 22475.0 25982.0 19069.0 23962.0];
       case 1998
           dis=[19697.0 22207.0 15089.0 26035.0 37194.0 40313.0 ...
                   55575.0 69106.0 42974.0 23043.0 23372.0 10601.0];
       case 1999
           dis=[13780.0 13638.0 10761.0 21968.0 43161.0 43064.0 ...
                   61477.0 41113.0 51229.0 30299.0 27379.0 17498.0];
       case 2000
           dis=[12380.0 22674.0 15064.0 18525.0 26115.0 30300.0 ...
                   39215.0 33924.0 37024.0 33622.0 29034.0 19227.0];
       case 2001
           dis=[14035.0 20892.0 14925.0 17833.0 33148.0 27650.0 ...
                   43563.0 26776.0 34081.0 19971.0 27059.0 16781.0];
       case 2002
           dis=[18154.0 15298.0 17473.0 25617.0 47192.0 45848.0 ...
                   54038.0 47352.0 42957.0 20511.0 22123.0 22616.0];
     end

avgflow=mean(season_cycle); %25031.7
season_cycle=season_cycle./avgflow;

r = 1;
Name = 'ChangJiang+Huai';
River.Name(r,1:length(Name)) = Name;
River.lon(r) = 120.0;
River.lat(r) = 31.9;
River.flag(r) = 3;
filen='RivDis';
River.file(r,1:length(filen))=filen;
River.area(r)=1794+244;
River.pfact(r)=1;
River.flow_mean(r)=mean(dis);
River.trans(r,1:12)=dis;

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
River.flow_mean(r)=1490.0;
River.trans(r,1:12)=season_cycle.*River.flow_mean(r);

save eas_rivers_discharge River

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
% grd_file ='E:\Working\roms\1_4\input\roms_eas_grid_jeju.nc';
% grd_file ='D:\matlab\roms\input\roms_eas_grid_jeju.nc';
grd_file ='D:\Roms\06-12-25(1_4)\roms_grid_4degree.nc';
N = 20; % Number of vertical levels

% the new netcdf forcing file - river data will be appended
%Fname = 'E:\Working\roms\1_4\input\ecco\roms_eas4_river.nc';
head='roms_';
mid=num2str(year(i));
foot='_river.nc';
Fname = [head,mid,foot];

% load the river flow data structure
Rname = 'eas_rivers_discharge';
load(Rname)
if ~exist('River')
  error([ 'Structure ''River'' does not exist in ' Rname])
end

% get the roms grid
% we need this to find the grid indices corresponding to the lon/lat
% locations of the dat
grd = roms_get_grid(grd_file);

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
% cycle=360.;
% River.time = [15:30:345];
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


for i=1:length(River.time)
  River.temp(1,:,i)=ones([1 N]).*Ytemp(i);
  River.temp(2,:,i)=ones([1 N]).*Htemp(i);
  River.salt(:,:,i)=ones([Nrivers N]).*salt;
end

%keyboard
%River.salt(41,:,:) = 28.9;

%-----------------------------------------------------------------------
%  Create empty river data FORCING NetCDF file.
%-----------------------------------------------------------------------

  disp([ 'Creating ' Fname '...'])

      create_empty_EAS_rivers

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
 
time=time+360;     
% time=time+365.25;     
end