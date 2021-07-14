function grd = gg(location, tyear)

% if nargin == 0
%   location = 'eas'; % default
% end
%    scoord = [5 0.4 50 20];

switch location
    case 'N_model'
%         grd_file = '../02_grid_depth/smoothing/grid_sumjin_estuary_v09.nc' ;
%         Vtransform = 2;
%         Vstretching = 4;
%         scoord = [1 1 1 20]; % theta_s theta_b hc N
%         
%         disp(' ')
%         disp([ 'Loading ROMS grd for application: ' location])
%         disp([ 'using grid file ' grd_file])
%         disp(' ')
%         grd = roms_get_grid(Vtransform, Vstretching, grd_file,scoord);
        grd_file = 'E:\Data\Model\ROMS\ysecs\input\roms_grd.nc' ;
        Vtransform = 1;
        Vstretching = 1;
        scoord = [5 0.4 4 20]; % theta_s theta_b hc N
        
        disp(' ')
        disp([ 'Loading ROMS grd for application: ' location])
        disp([ 'using grid file ' grd_file])
        disp(' ')
        grd = roms_get_grid(Vtransform, Vstretching, grd_file,scoord);

    case 'L_model'
%       grd_file = 'L_model/roms_grd.nc';
%       Vtransform = 2;
%       Vstretching = 4;
%       scoord = [5 0.4 4 20]; % theta_s theta_b hc N
% %--- ÁöÈÆÀÇ NWP ¸ðµ¨ -------------------
% %         grd_file = 'F:\ROMS\Sumjin\08_bounday_ts\make_bndy\L_model\NWP2010_byJJH\monthly\roms_grid_NWP.nc' ;
% %         scoord = [10 0 250 40]; % theta_s theta_b hc N
% %--------------------------------------
%         disp(' ')
%         disp([ 'Loading ROMS grd for application: ' location])
%         disp([ 'using grid file ' grd_file])
%         disp(' ')
%         grd = roms_get_grid(Vtransform, Vstretching, grd_file,scoord);
      
        grd_file = ['E:\Data\Model\ROMS\nwp_1_10\test06\DA\', num2str(tyear,'%04i'),'\test06_monthly_',num2str(tyear,'%04i'), '_01.nc'];
        Vtransform = 2;
        Vstretching = 4;
        scoord = [10 1 250 40]; % theta_s theta_b hc N

        disp(' ')
        disp([ 'Loading ROMS grd for application: ' location])
        disp([ 'using grid file ' grd_file])
        disp(' ')
        grd = roms_get_grid(Vtransform, Vstretching, grd_file,scoord);
        
end