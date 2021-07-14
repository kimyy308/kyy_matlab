function nc_make_tides(fname,Ntides,start_tide_mjd,components)
%
% % % This File is made by Kim. Y. Y. at 2017-07-13
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nc=netcdf(fname,'write');
nc=netcdf(fname,'clobber');
% redef(nc);
%
%  Add dimension
%
nc('tide_period')=Ntides;
%
%  Add variables and attributes
%
nc{'tide_period'} = ncdouble('tide_period');
nc{'tide_period'}.long_name = ncchar('Tide angular period');
nc{'tide_period'}.long_name = 'Tide angular period';
nc{'tide_period'}.units = ncchar('Hours');
nc{'tide_period'}.units = 'Hours';

nc{'tide_Ephase'} = ncdouble('tide_period', 'eta_rho', 'xi_rho');
nc{'tide_Ephase'}.long_name = ncchar('Tidal elevation phase angle');
nc{'tide_Ephase'}.long_name = 'Tidal elevation phase angle';
nc{'tide_Ephase'}.units = ncchar('Degrees');
nc{'tide_Ephase'}.units = 'Degrees';

nc{'tide_Eamp'} = ncdouble('tide_period', 'eta_rho', 'xi_rho');
nc{'tide_Eamp'}.long_name = ncchar('Tidal elevation amplitude');
nc{'tide_Eamp'}.long_name = 'Tidal elevation amplitude';
nc{'tide_Eamp'}.units = ncchar('Meter');
nc{'tide_Eamp'}.units = 'Meter';

nc{'tide_Cmin'} = ncdouble('tide_period', 'eta_rho', 'xi_rho');
nc{'tide_Cmin'}.long_name = ncchar('Tidal current ellipse semi-minor axis');
nc{'tide_Cmin'}.long_name = 'Tidal current ellipse semi-minor axis';
nc{'tide_Cmin'}.units = ncchar('Meter second-1');
nc{'tide_Cmin'}.units = 'Meter second-1';

nc{'tide_Cmax'} = ncdouble('tide_period', 'eta_rho', 'xi_rho');
nc{'tide_Cmax'}.long_name = ncchar('Tidal current, ellipse semi-major axis');
nc{'tide_Cmax'}.long_name = 'Tidal current, ellipse semi-major axis';
nc{'tide_Cmax'}.units = ncchar('Meter second-1');
nc{'tide_Cmax'}.units = 'Meter second-1';

nc{'tide_Cangle'} = ncdouble('tide_period', 'eta_rho', 'xi_rho');
nc{'tide_Cangle'}.long_name = ncchar('Tidal current inclination angle');
nc{'tide_Cangle'}.long_name = 'Tidal current inclination angle';
nc{'tide_Cangle'}.units = ncchar('Degrees between semi-major axis and East');
nc{'tide_Cangle'}.units = 'Degrees between semi-major axis and East';

nc{'tide_Cphase'} = ncdouble('tide_period', 'eta_rho', 'xi_rho');
nc{'tide_Cphase'}.long_name = ncchar('Tidal current phase angle');
nc{'tide_Cphase'}.long_name = 'Tidal current phase angle';
nc{'tide_Cphase'}.units = ncchar('Degrees');
nc{'tide_Cphase'}.units = 'Degrees';

nc.date = ncchar(date);
nc.date = date;
nc.start_tide_mjd=start_tide_mjd;
nc.components = ncchar(components);
nc.components = components;

close(nc)
