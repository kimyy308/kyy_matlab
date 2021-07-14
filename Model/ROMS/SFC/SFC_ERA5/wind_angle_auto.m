clear all; clc

disp(' ')
disp('Rotating wind vector')
disp(' ')

ufile = ['asan_frc_KMA_2019Uwind_1hourly.nc'];
vfile = ['asan_frc_KMA_2019Vwind_1hourly.nc'];

g = grd('asan');
cosa = cos(g.angle);
sina = sin(g.angle);

unc = netcdf(ufile, 'w');
u = unc{'Uwind'}(:);

vnc = netcdf(vfile, 'w');
v = vnc{'Vwind'}(:);

cosa_mat = repmat(cosa, [1,1,length(u)]);
cosa_mat = permute(cosa_mat, [3 1 2]);

sina_mat = repmat(sina, [1, 1, length(v)]);
sina_mat = permute(sina_mat, [3 1 2]);

urot = u.*cosa_mat + v.*sina_mat;
vrot = v.*cosa_mat - u.*sina_mat;

unc{'Uwind'}(:) = urot;
vnc{'Vwind'}(:) = vrot;

close(unc);
close(vnc);
