function grd = gg(location)

% if nargin == 0
%   location = 'eas'; % default
% end
%    scoord = [5 0.4 50 20];

switch location

    case 'NWP_1_20'
        vert_param
        grd_file ='/scratch/snu02/roms_nwp/output/smooth13_vtvs/roms_grid_combine2_smooth13.nc' ;
        scoord = [7.0 2.0 250 40]; % theta_s theta_b hc N
        
        disp(' ')
        disp([ 'Loading ROMS grd for application: ' location])
        disp([ 'using grid file ' grd_file])
        disp(' ')
        grd = roms_get_grid(grd_file,scoord);

        
end