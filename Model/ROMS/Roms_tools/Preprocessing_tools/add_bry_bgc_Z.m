function add_bry_pisces_Z(zbryname,obc,Z,time_NO3,time_zoop,time_phyto,time_chlo,time_detritus,cycle,clobber);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                 %
%   function add_bry_pisces_Z(zbryname,obc,...                    %
%                             Z,time,cycle,clobber);              %
%                                                                 %
%   This function create the header of a Netcdf climatology       %
%   file.                                                         %
%                                                                 %
%   Input:                                                        %
%                                                                 %
%   zbryname     Netcdf climatology file name (character string). %
%   obc          open boundaries flag (1=open , [S E N W]).       %
%   Z            Depth of vertical levels.(Vector)                %
%   time         time.(vector)                                    %
%   cycle        Length (days) for cycling the climatology.(Real) %
%   clobber      Switch to allow or not writing over an existing  %
%                file.(character string)                          %
%                                                                 %
%  Pierrick Penven, IRD, 2005.                                    %
%  Olivier Aumont the master, IRD, 2006.                          %
%  Patricio Marchesiello, chief, IRD, 2007.                       %
%  Christophe Eugene Raoul Menkes, the slave, IRD, 2007.          %
%                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp([' Adding PISCES data into file : ',zbryname])
disp(' ')
%
%  Create the boundary file
%
nc = netcdf(zbryname,clobber);
result = redef(nc);
%
%  Create dimensions
%
nc('NO3_time') = length(time_NO3);
nc('chlo_time') = length(time_chlo);
nc('zoop_time') = length(time_zoop);
nc('phyt_time') = length(time_phyto);
nc('detritus_time')   = length(time_detritus);
nc('one') = 1;
%
% BGC
nc{'ZNO3'} = ncdouble('ZNO3') ;
nc{'ZNO3'}.long_name = ncchar('Depth');
nc{'ZNO3'}.long_name = 'Depth';
nc{'ZNO3'}.units = ncchar('m');
nc{'ZNO3'}.units = 'm';
%
%  Create variables and attributes
%
nc{'NO3_time'} = ncdouble('NO3_time') ;
nc{'NO3_time'}.long_name = ncchar('time for NO3 climatology')
nc{'NO3_time'}.long_name = 'time for NO3 climatology';
nc{'NO3_time'}.units = ncchar('day');
nc{'NO3_time'}.units = 'day';
nc{'NO3_time'}.cycle_length = cycle;%
%
nc{'chlo_time'} = ncdouble('chlo_time') ;
nc{'chlo_time'}.long_name = ncchar('time for chlo climatology');
nc{'chlo_time'}.long_name = 'time for chlo climatology';
nc{'chlo_time'}.units = ncchar('day');
nc{'chlo_time'}.units = 'day';
nc{'chlo_time'}.cycle_length = cycle;%
%
nc{'phyt_time'} = ncdouble('phyt_time') ;
nc{'phyt_time'}.long_name = ncchar('time for phyt climatology');
nc{'phyt_time'}.long_name = 'time for phyt climatology';
nc{'phyt_time'}.units = ncchar('day');
nc{'phyt_time'}.units = 'day';
nc{'phyt_time'}.cycle_length = cycle;%
%
nc{'zoop_time'} = ncdouble('zoop_time') ;
nc{'zoop_time'}.long_name = ncchar('time for zoop climatology');
nc{'zoop_time'}.long_name = 'time for zoop climatology';
nc{'zoop_time'}.units = ncchar('day');
nc{'zoop_time'}.units = 'day';
nc{'zoop_time'}.cycle_length = cycle;
%
nc{'detritus_time'} = ncdouble('detritus_time') ;
nc{'detritus_time'}.long_name = ncchar('time for detritus climatology');
nc{'detritus_time'}.long_name = 'time for detritus climatology';
nc{'detritus_time'}.units = ncchar('day');
nc{'detritus_time'}.units = 'day';
nc{'detritus_time'}.cycle_length = cycle;%
%
if obc(1)==1
%
%   Southern boundary
%
  nc{'NO3_south'} = ncdouble('NO3_time','Z','xi_rho') ;
  nc{'NO3_south'}.long_name = ncchar('southern boundary NO3');
  nc{'NO3_south'}.long_name = 'southern boundary NO3';
  nc{'NO3_south'}.units = ncchar('mMol N m-3');
  nc{'NO3_south'}.units = 'mMol N m-3';
%
  nc{'chlo_south'} = ncdouble('chlo_time','Z','xi_rho') ;
  nc{'chlo_south'}.long_name = ncchar('southern boundary chlo');
  nc{'chlo_south'}.long_name = 'southern boundary chlo';
  nc{'chlo_south'}.units = ncchar('mMol N m-3');
  nc{'chlo_south'}.units = 'mMol N m-3';
%
  nc{'zoop_south'} = ncdouble('zoop_time','Z','xi_rho') ;
  nc{'zoop_south'}.long_name = ncchar('southern boundary zoop');
  nc{'zoop_south'}.long_name = 'southern boundary zoop';
  nc{'zoop_south'}.units = ncchar('mMol N m-3');
  nc{'zoop_south'}.units = 'mMol N m-3';
%
  nc{'phyt_south'} = ncdouble('phyt_time','Z','xi_rho') ;
  nc{'phyt_south'}.long_name = ncchar('southern boundary phyt');
  nc{'phyt_south'}.long_name = 'southern boundary phyt';
  nc{'phyt_south'}.units = ncchar('mMol N m-3');
  nc{'phyt_south'}.units = 'mMol N m-3';
  %
  nc{'detritus_south'} = ncdouble('detritus_time','Z','xi_rho') ;
  nc{'detritus_south'}.long_name = ncchar('southern boundary detritus');
  nc{'detritus_south'}.long_name = 'southern boundary detritus';
  nc{'detritus_south'}.units = ncchar('mMol N m-3');
  nc{'detritus_south'}.units = 'mMol N m-3';
  
%
end
%
if obc(2)==1
%
%   Eastern boundary
%
  nc{'NO3_east'} = ncdouble('NO3_time','Z','eta_rho') ;
  nc{'NO3_east'}.long_name = ncchar('eastern boundary NO3');
  nc{'NO3_east'}.long_name = 'eastern boundary NO3';
  nc{'NO3_east'}.units = ncchar('mMol N m-3');
  nc{'NO3_east'}.units = 'mMol N m-3';
%
  nc{'chlo_east'} = ncdouble('chlo_time','Z','eta_rho') ;
  nc{'chlo_east'}.long_name = ncchar('eastern boundary chlo');
  nc{'chlo_east'}.long_name = 'eastern boundary chlo';
  nc{'chlo_east'}.units = ncchar('mMol N m-3');
  nc{'chlo_east'}.units = 'mMol N m-3';
%
  nc{'zoop_east'} = ncdouble('zoop_time','Z','eta_rho') ;
  nc{'zoop_east'}.long_name = ncchar('eastern boundary zoop');
  nc{'zoop_east'}.long_name = 'eastern boundary zoop';
  nc{'zoop_east'}.units = ncchar('mMol N m-3');
  nc{'zoop_east'}.units = 'mMol N m-3';
%
  nc{'phyt_east'} = ncdouble('phyt_time','Z','eta_rho') ;
  nc{'phyt_east'}.long_name = ncchar('eastern boundary phyt');
  nc{'phyt_east'}.long_name = 'eastern boundary phyt';
  nc{'phyt_east'}.units = ncchar('mMol N m-3');
  nc{'phyt_east'}.units = 'mMol N m-3';
%
%
  nc{'detritus_east'} = ncdouble('detritus_time','Z','xi_rho') ;
  nc{'detritus_east'}.long_name = ncchar('east boundary detritus');
  nc{'detritus_east'}.long_name = 'east boundary detritus';
  nc{'detritus_east'}.units = ncchar('mMol N m-3');
  nc{'detritus_east'}.units = 'mMol N m-3';
  
%
end
%
if obc(3)==1
%
%   Northern boundary
%
  nc{'NO3_north'} = ncdouble('NO3_time','Z','xi_rho') ;
  nc{'NO3_north'}.long_name = ncchar('northern boundary NO3');
  nc{'NO3_north'}.long_name = 'northern boundary NO3';
  nc{'NO3_north'}.units = ncchar('mMol N m-3');
  nc{'NO3_north'}.units = 'mMol N m-3';
%
  nc{'chlo_north'} = ncdouble('chlo_time','Z','xi_rho') ;
  nc{'chlo_north'}.long_name = ncchar('northern boundary chlo');
  nc{'chlo_north'}.long_name = 'northern boundary chlo';
  nc{'chlo_north'}.units = ncchar('mMol N m-3');
  nc{'chlo_north'}.units = 'mMol N m-3';
%
  nc{'zoop_north'} = ncdouble('zoop_time','Z','xi_rho') ;
  nc{'zoop_north'}.long_name = ncchar('northern boundary zoop');
  nc{'zoop_north'}.long_name = 'northern boundary zoop';
  nc{'zoop_north'}.units = ncchar('mMol N m-3');
  nc{'zoop_north'}.units = 'mMol N m-3';
%
  nc{'phyt_north'} = ncdouble('phyt_time','Z','xi_rho') ;
  nc{'phyt_north'}.long_name = ncchar('northern boundary phyt');
  nc{'phyt_north'}.long_name = 'northern boundary phyt';
  nc{'phyt_north'}.units = ncchar('mMol N m-3');
  nc{'phyt_north'}.units = 'mMol N m-3';
%
%
  nc{'detritus_north'} = ncdouble('detritus_time','Z','xi_rho') ;
  nc{'detritus_north'}.long_name = ncchar('north boundary detritus');
  nc{'detritus_north'}.long_name = 'north boundary detritus';
  nc{'detritus_north'}.units = ncchar('mMol N m-3');
  nc{'detritus_north'}.units = 'mMol N m-3';
  
%
end
%
if obc(4)==1
%
%   Western boundary
%
  nc{'NO3_west'} = ncdouble('NO3_time','Z','eta_rho') ;
  nc{'NO3_west'}.long_name = ncchar('western boundary NO3');
  nc{'NO3_west'}.long_name = 'western boundary NO3';
  nc{'NO3_west'}.units = ncchar('mMol N m-3');
  nc{'NO3_west'}.units = 'mMol N m-3';
%
  nc{'chlo_west'} = ncdouble('chlo_time','Z','eta_rho') ;
  nc{'chlo_west'}.long_name = ncchar('western boundary chlo');
  nc{'chlo_west'}.long_name = 'western boundary chlo';
  nc{'chlo_west'}.units = ncchar('mMol N m-3');
  nc{'chlo_west'}.units = 'mMol N m-3';
%
  nc{'zoop_west'} = ncdouble('zoop_time','Z','eta_rho') ;
  nc{'zoop_west'}.long_name = ncchar('western boundary zoop');
  nc{'zoop_west'}.long_name = 'western boundary zoop';
  nc{'zoop_west'}.units = ncchar('mMol N m-3');
  nc{'zoop_west'}.units = 'mMol N m-3';
%
  nc{'phyt_west'} = ncdouble('phyt_time','Z','eta_rho') ;
  nc{'phyt_west'}.long_name = ncchar('western boundary phyt');
  nc{'phyt_west'}.long_name = 'western boundary phyt';
  nc{'phyt_west'}.units = ncchar('mMol N m-3');
  nc{'phyt_west'}.units = 'mMol N m-3';
%
%
  nc{'detritus_west'} = ncdouble('detritus_time','Z','xi_rho') ;
  nc{'detritus_west'}.long_name = ncchar('west boundary detritus');
  nc{'detritus_west'}.long_name = 'west boundary detritus';
  nc{'detritus_west'}.units = ncchar('mMol N m-3');
  nc{'detritus_west'}.units = 'mMol N m-3';
  
%
end
%
% Leave define mode
%
result = endef(nc);
%
% Write variables
%
nc{'detritus_time'}(:) = time_detritus;
nc{'NO3_time'}(:) = time_NO3;
nc{'chlo_time'}(:) = time_chlo;
nc{'zoop_time'}(:) = time_zoop;
nc{'phyt_time'}(:) = time_phyto;
if obc(1)==1
  nc{'NO3_south'}(:)  =  0;
  nc{'chlo_south'}(:)  =  0;
  nc{'zoop_south'}(:) =  0;
  nc{'phyt_south'}(:)   =  0;
  nc{'detritus_south'}(:)   =  0;
end
if obc(2)==1 
  nc{'NO3_east'}(:)  =  0;
  nc{'chlo_east'}(:)  =  0;
  nc{'zoop_east'}(:) =  0;
  nc{'phyt_east'}(:)   =  0;
  nc{'detritus_east'}(:)   =  0;
end 
if obc(3)==1 
  nc{'NO3_north'}(:)  =  0;
  nc{'chlo_north'}(:)  =  0;
  nc{'zoop_north'}(:) =  0;
  nc{'phyt_north'}(:)   =  0;
  nc{'detritus_north'}(:)   =  0;
end 
if obc(4)==1 
  nc{'NO3_west'}(:)  =  0;
  nc{'chlo_west'}(:)  =  0;
  nc{'zoop_west'}(:) =  0;
  nc{'phyt_west'}(:)   =  0;
  nc{'detritus_west'}(:)   =  0;
end 
close(nc)
return


