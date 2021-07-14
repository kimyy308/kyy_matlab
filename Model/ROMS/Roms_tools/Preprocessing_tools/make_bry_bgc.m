%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Build a ROMS boundary file
%
%  Extrapole and interpole temperature and salinity from a
%  climatology to get boundary conditions for
%  ROMS (boundary netcdf file) .
%  Get the velocities and sea surface elevation via a 
%  geostrophic computation.
%
%  Data input format (netcdf):
%     temperature(T, Z, Y, X)
%     T : time [Months]
%     Z : Depth [m]
%     Y : Latitude [degree north]
%     X : Longitude [degree east]
%
%  Data source : IRI/LDEO climate Data Library (World Ocean Atlas 1998)
%    http://ingrid.ldgo.columbia.edu/
%    http://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NODC/.WOA98/
% 
%  Further Information:  
%  http://www.brest.ird.fr/Roms_tools/
%  
%  This file is part of ROMSTOOLS
%
%  ROMSTOOLS is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation; either version 2 of the License,
%  or (at your option) any later version.
%
%  ROMSTOOLS is distributed in the hope that it will be useful, but
%  WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
%  MA  02111-1307  USA
%
%  Copyright (c) 2005-2006 by Pierrick Penven 
%  e-mail:Pierrick.Penven@ird.fr  
%
%  Updated    1-Sep-2006 by Pierrick Penven
%  Pierrick Penven, IRD, 2005.                                    %
%  Olivier Aumont the master, IRD, 2006.                          %
%  Patricio Marchesiello, chief, IRD, 2007.                       %
%  Christophe Eugene Raoul Menkes, the slave, IRD, 2007.          %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  WARNING !!!!!! THIS ASSUMES THAT THE TIME FOR PISCES INITIAL.
%  IS THE SAME AS THE CLIM T AND S. ELSE, CHANGE THE PROGRAM
%
clear all
close all
%%%%%%%%%%%%%%%%%%%%% USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%
% Common parameters
%
romstools_param
Roa;
%
% Set times and cycles: monthly climatology for all data
%
%
%  Data climatologies file names:
%
%
% NO3_seas_data  = [woapisces_dir,'no3_seas.cdf'];
% NO3_ann_data   = [woapisces_dir,'no3_ann.cdf'];
% po4_seas_data  = [woapisces_dir,'po4_seas.cdf'];
% po4_ann_data   = [woapisces_dir,'po4_ann.cdf'];
% o2_seas_data   = [woapisces_dir,'o2_seas.cdf'];
% o2_ann_data    = [woapisces_dir,'o2_ann.cdf'];
% sio3_seas_data = [woapisces_dir,'sio3_seas.cdf'];
% sio3_ann_data  = [woapisces_dir,'sio3_ann.cdf'];
%
NO3_seas_data    = [woa_dir,'no3_seas.cdf'];
NO3_ann_data     = [woa_dir,'no3_ann.cdf'];
chlo_seas_data   = [chla_dir,'chla_seas.cdf'];
chlo_ann_data    = [chla_dir,'chla_ann.cdf'];
NO3min=1;
%
%
%%%%%%%%%%%%%%%%%%% END USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%
%
% Title
%
disp(' ')
disp([' Making the file: ',bryname])
disp([' Adding the PISCES variables'])
disp(' ')
disp([' Title: ',ROMS_title])
%
% Read in the grid
%
disp(' ')
disp(' Read in the grid...')
nc=netcdf(grdname);
lon=nc{'lon_rho'}(:);
lat=nc{'lat_rho'}(:);
Lp=length(nc('xi_rho'));
Mp=length(nc('eta_rho'));
hmax=max(max(nc{'h'}(:)));
result=close(nc);
% %
%Get the time of data
%
nc=netcdf(NO3_seas_data);
time_NO3=nc{'T'}(:);
close(nc)
time_NO3=(time_NO3-1)*30;
%
nc=netcdf(chlo_seas_data);
time_chlo=nc{'T'}(:);
close(nc)
time_chlo=(time_chlo-1)*30;
time_zoop=time_chlo;
time_phyto=time_chlo;
time_detritus=time_chlo;
cycle=360;
%
% Redefine the boundary file
%
if (makebry)
  disp(' ')
  disp(' Redefine the boundary file...')
  add_bry_bgc(bryname,obc,time_NO3,time_zoop,time_phyto,time_chlo,time_detritus,cycle,'write');
end
%
% Redefine the boundary file in Z-coordinates
%
if (makeZbry)
  disp(' ')
  disp(' Redefine the boundary Z-file...')
%
% get Z
%
%pause
  nc=netcdf(NO3_ann_data);
  Z=nc{'Z'}(:);
  kmax=max(find(Z<hmax))-1;
  Z=Z(1:kmax);
  close(nc)
  add_bry_bgc_Z(Zbryname,obc,Z,time_NO3,time_zoop,time_phyto,time_chlo,time_detritus,cycle,'write');
  disp(' ')
  disp(' Horizontal extrapolations')
%
% Loop on the lateral boundaries 
%
  for obcndx=1:4
    if obc(obcndx)==1
      if obcndx==1
        disp(' Processing southern boundary...')
		disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
	suffix='_south';
      elseif obcndx==2
        disp(' Processing eastern boundary...')
		disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
	suffix='_east';
      elseif obcndx==3
        disp(' Processing northern boundary...')
		disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
	suffix='_north';
      elseif obcndx==4
        disp(' Processing western boundary...')
		disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
	suffix='_west';
	  end
      %
	  disp('============================================')
	  disp('  Nitrate...')
	  cff=1;
      bry_interp_bgc(Zbryname,lon,lat,NO3_seas_data,NO3_ann_data,...
               'nitrate',['NO3',suffix],obcndx,Roa);        
      %
	  disp('============================================')
	  disp('  Chlorophylle...')
	  cff=1;
	  bry_interp_bgc_chloro(bryname,grdname,clmname,lon,lat,chlo_seas_data,chlo_ann_data,...
		'Chlorophylle',['chlo',suffix],obcndx,cff,Roa);
      %
	  disp('============================================')
	  disp('  Phytoplankton...')
	  cff=0.5;       
	  bry_interp_bgc_chloro(bryname,grdname,clmname,lon,lat,chlo_seas_data,chlo_ann_data,...
               'Phytoplancton,',['PHYTO',suffix],obcndx,cff,Roa);        
      %
	  disp('============================================')
      disp('  Zooplankton...')
	  cff=0.2;     
	  bry_interp_bgc_chloro(bryname,grdname,clmname,lon,lat,chlo_seas_data,chlo_ann_data,...
               'Zooplankton',['zoop',suffix],obcndx,cff,Roa);   
           
      disp('============================================')
      disp('  detritus...')
	  cff=0.2;     
	  bry_interp_bgc_chloro(bryname,grdname,clmname,lon,lat,chlo_seas_data,chlo_ann_data,...
               'Zooplankton',['zoop',suffix],obcndx,cff,Roa);  
    end
  end
end
%
% Vertical interpolations 
%
if (makebry)
  disp(' ')
  disp(' Vertical interpolations')
%
% Loop on the lateral boundaries 
%
  for obcndx=1:4
    if obc(obcndx)==1
      if obcndx==1
        disp(' Processing southern boundary...')
	suffix='_south';
      elseif obcndx==2
        disp(' Processing eastern boundary...')
	suffix='_east';
      elseif obcndx==3
        disp(' Processing northern boundary...')
	suffix='_north';
      elseif obcndx==4
        disp(' Processing western boundary...')
	suffix='_west';
      end
      disp(' ')
      disp('  Nitrate...')
      vinterp_bry_bgc(bryname,grdname,Zbryname,'NO3_time',['NO3',suffix],obcndx);
%       disp(' ')
%       disp('  Chlorophylle...')
%       vinterp_bry(bryname,grdname,Zbryname,['chlo',suffix],obcndx);
%       disp(' ')
%       disp('  Phytoplankton...')
%       vinterp_bry(bryname,grdname,Zbryname,['PHYTO',suffix],obcndx);
%       disp(' ')
%       disp('  Zooplankton ...')
%       vinterp_bry(bryname,grdname,Zbryname,['ZOO',suffix],obcndx);
    end
  end
end
%
% Make a few plots
%
if makeplot==1
disp(' ')
disp(' Make a few plots...')
test_bry(bryname,grdname,'NO3',1,obc)
figure
test_bry(bryname,grdname,'chlo',1,obc)
figure
test_bry(bryname,grdname,'phyt',1,obc)
figure
test_bry(bryname,grdname,'zoop',1,obc)
figure
test_bry(bryname,grdname,'detritus',1,obc)

end
%
% End
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
