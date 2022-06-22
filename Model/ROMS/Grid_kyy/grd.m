function grd = gg(location, testname)
% %  Updated 21-May-2018   by Yong-Yub Kim

% if nargin == 0
%   location = 'eas'; % default
% end
%    scoord = [5 0.4 50 20];
if nargin == 1
    switch location  
    case 'test06'
%         vert_param
        grd_file ='/data1/kimyy/Model/ROMS/roms_nwp/nwp_1_20/grids/roms_grid_nwp_1_20_test06.nc' ;
        scoord = [5.0 0.4 5.0 40] % theta_s theta_b hc N
        Vtransform = 1;
        Vstretching = 1;
        disp(' ')
        disp([ 'Loading ROMS grd for application: ' location])
        disp([ 'using grid file ' grd_file])
        disp(' ')
        grd = roms_get_grid(Vtransform,Vstretching,grd_file,scoord);
        
    case 'NWP_1_20'
%         vert_param
%         grd_file ='/data1/kimyy/Model/ROMS/roms_nwp/nwp_1_20/grids/roms_grid_nwp_1_20_test37.nc' ;
%         grd_file ='E:\Data\Model\ROMS\nwp_1_20\input\test48\roms_grid_nwp_1_20_test48.nc' ;
%             grd_file ='D:\Data\Model\ROMS\nwp_1_20\input\test53\roms_grid_nwp_1_20_test53.nc' ;
            grd_file ='D:\Data\Model\ROMS\nwp_1_20\input\test2117\roms_grid_nwp_1_20_test2117.nc' ;

        scoord = [10.0 1.0 250.0 40]; % theta_s theta_b hc N
        Vtransform = 2;
        Vstretching = 4;
        disp(' ')
        disp([ 'Loading ROMS grd for application: ' location])
        disp([ 'using grid file ' grd_file])
        disp(' ')
        grd = roms_get_grid(Vtransform,Vstretching,grd_file,scoord);
   
    case 'ES_1_40'
        grd_file ='/data1/kimyy/Model/ROMS/roms_nwp/nwp_1_20/grids/etopo1_Eastsea_40.nc' ;
        scoord = [10.0 1.0 250.0 40] % theta_s theta_b hc N
        Vtransform = 2;
        Vstretching = 4;
        disp(' ')
        disp([ 'Loading ROMS grd for application: ' location])
        disp([ 'using grid file ' grd_file])
        disp(' ')
        grd = roms_get_grid(Vtransform,Vstretching,grd_file,scoord);
    
    case 'seo_NWP_reanalysis'    
        grd_file ='/data2/kimyy/Reanalysis/nwp_1_10_seo/ens_mean_monthly/roms_grid_10km_new.nc' ;
        scoord = [5.0 0.4 5.0 20] % theta_s theta_b hc N
        Vtransform = 1;
        Vstretching = 1;
        disp(' ')
        disp([ 'Loading ROMS grd for application: ' location])
        disp([ 'using grid file ' grd_file])
        disp(' ')
        grd = roms_get_grid(Vtransform,Vstretching,grd_file,scoord);
     case 'ES_convection'    
        grd_file ='E:\Data\Model\ROMS\ES_conv\input\roms_grid_ES_conv_test01.nc' ;
        scoord = [10.0 0.5 1100.0 208] % theta_s theta_b hc N
        Vtransform = 2;
        Vstretching = 4;
        disp(' ')
        disp([ 'Loading ROMS grd for application: ' location])
        disp([ 'using grid file ' grd_file])
        disp(' ')
        grd = roms_get_grid(Vtransform,Vstretching,grd_file,scoord);
    end
end

if nargin == 2
    switch location
    case 'NWP_1_10'
%         vert_param
%         grd_file ='/data1/kimyy/Model/ROMS/roms_nwp/nwp_1_20/grids/roms_grid_nwp_1_20_test38.nc' ;
        grd_file =['/scratch/snu01/kimyy/roms_nwp/nwp_1_10/input/',testname,'/spinup/2017/roms_grid_nwp_1_20_',testname,'.nc'] ;
        scoord = [10.0 1.0 250.0 40] % theta_s theta_b hc N
        Vtransform = 2;
        Vstretching = 4;
        disp(' ')
        disp([ 'Loading ROMS grd for application: ' location])
        disp([ 'using grid file ' grd_file])
        disp(' ')
        grd = roms_get_grid(Vtransform,Vstretching,grd_file,scoord);
    case 'nwp_1_10_linux'
%         vert_param
        grd_file =['/data1/kimyy/Model/ROMS/nwp_1_10/input/',testname,'/spinup/roms_grid_nwp_1_10_',testname,'.nc'] ;
        scoord = [10.0 1.0 250.0 40] % theta_s theta_b hc N
        Vtransform = 2;
        Vstretching = 4;
        disp(' ')
        disp([ 'Loading ROMS grd for application: ' location])
        disp([ 'using grid file ' grd_file])
        disp(' ')
        grd = roms_get_grid(Vtransform,Vstretching,grd_file,scoord); 
        
    case 'nwp_1_10_window'
%         vert_param
        grd_file =['D:\Data\Model\ROMS\nwp_1_10\input\',testname,'\roms_grid_nwp_1_10_',testname,'.nc'] ;
        scoord = [10.0 1.0 250.0 40] % theta_s theta_b hc N
        Vtransform = 2;
        Vstretching = 4;
        disp(' ')
        disp([ 'Loading ROMS grd for application: ' location])
        disp([ 'using grid file ' grd_file])
        disp(' ')
        grd = roms_get_grid(Vtransform,Vstretching,grd_file,scoord); 

    case 'NWP_1_20'
%         vert_param
%         grd_file ='/data1/kimyy/Model/ROMS/roms_nwp/nwp_1_20/grids/roms_grid_nwp_1_20_test38.nc' ;
%         grd_file ='E:\Data\Model\ROMS\nwp_1_20\input\test38\roms_grid_nwp_1_20_test38.nc' ;
%         grd_file ='E:\Data\Model\ROMS\nwp_1_20\input\test48\roms_grid_nwp_1_20_test48.nc' ;
        grd_file = ['D:\Data\Model\ROMS\nwp_1_20\input\',testname,'\roms_grid_nwp_1_20_',testname,'.nc'] ;
        scoord = [10.0 1.0 250.0 40] % theta_s theta_b hc N
        Vtransform = 2;
        Vstretching = 4;
        disp(' ')
        disp([ 'Loading ROMS grd for application: ' location])
        disp([ 'using grid file ' grd_file])
        disp(' ')
        grd = roms_get_grid(Vtransform,Vstretching,grd_file,scoord);
    case 'NWP_1_20_linux'
%         vert_param
        grd_file =['/data1/kimyy/Model/ROMS/roms_nwp/nwp_1_20/input/',testname,'/spinup/roms_grid_nwp_1_20_',testname,'.nc'] ;
        scoord = [10.0 1.0 250.0 40] % theta_s theta_b hc N
        Vtransform = 2;
        Vstretching = 4;
        disp(' ')
        disp([ 'Loading ROMS grd for application: ' location])
        disp([ 'using grid file ' grd_file])
        disp(' ')
        grd = roms_get_grid(Vtransform,Vstretching,grd_file,scoord); 
    case 'NWP_1_40_linux'
%         vert_param
        grd_file =['/data1/kimyy/Model/ROMS/roms_nwp/nwp_1_40/input/',testname,'/spinup/roms_grid_nwp_1_40_',testname,'.nc'] ;
        scoord = [5.0 0.4 5.0 40] % theta_s theta_b hc N
        Vtransform = 1;
        Vstretching = 1;
        disp(' ')
        disp([ 'Loading ROMS grd for application: ' location])
        disp([ 'using grid file ' grd_file])
        disp(' ')
        grd = roms_get_grid(Vtransform,Vstretching,grd_file,scoord);
    
        
    case 'NWP_1_50_linux'
%         vert_param
        grd_file =['/data1/kimyy/Model/ROMS/roms_nwp/nwp_1_50/input/',testname,'/spinup/roms_grid_nwp_1_50_',testname,'.nc'] ;
        scoord = [5.0 0.4 5.0 40] % theta_s theta_b hc N
        Vtransform = 1;
        Vstretching = 1;
        disp(' ')
        disp([ 'Loading ROMS grd for application: ' location])
        disp([ 'using grid file ' grd_file])
        disp(' ')
        grd = roms_get_grid(Vtransform,Vstretching,grd_file,scoord);
        
        
    case 'NWP_1_100_linux'
%         vert_param
        grd_file =['/data1/kimyy/Model/ROMS/roms_nwp/nwp_1_100/input/',testname,'/spinup/roms_grid_nwp_1_100_',testname,'.nc'] ;
        scoord = [5.0 0.4 5.0 40] % theta_s theta_b hc N
        Vtransform = 1;
        Vstretching = 1;
        disp(' ')
        disp([ 'Loading ROMS grd for application: ' location])
        disp([ 'using grid file ' grd_file])
        disp(' ')
        grd = roms_get_grid(Vtransform,Vstretching,grd_file,scoord);    
        
    case 'NWP_1_20_DGIST'
        fid=fopen('/scratch/snu01/kimyy/roms_nwp/nwp_1_20/output/modelinfo');
        modelinfo=textscan(fid,'%s');
        fclose(fid);
        grd_file =['/scratch/snu01/kimyy/roms_nwp/nwp_1_20/input/',testname,'/spinup/roms_grid_nwp_1_20_',testname,'.nc'] ; %% same grid is used at both spinup and run case
        scoord(1) = str2num(modelinfo{1,1}{9,1})
        scoord(2) = str2num(modelinfo{1,1}{10,1})
        scoord(3) = str2num(modelinfo{1,1}{11,1})
        scoord(4) = str2num(modelinfo{1,1}{6,1})
        Vtransform = str2num(modelinfo{1,1}{7,1})
        Vstretching = str2num(modelinfo{1,1}{8,1})
        disp(' ')
        disp([ 'Loading ROMS grd for application: ' location])
        disp([ 'using grid file ' grd_file])
        disp(' ')
        grd = roms_get_grid(Vtransform,Vstretching,grd_file,scoord);
    case 'ES_1_40'
        grd_file ='/data1/kimyy/Model/ROMS/roms_nwp/nwp_1_20/grids/etopo1_Eastsea_40.nc' ;
        scoord = [7.0 2.0 250.0 40] % theta_s theta_b hc N
        Vtransform = 2;
        Vstretching = 4;
        disp(' ')
        disp([ 'Loading ROMS grd for application: ' location])
        disp([ 'using grid file ' grd_file])
        disp(' ')
        grd = roms_get_grid(Vtransform,Vstretching,grd_file,scoord);
    
    case 'seo_NWP_reanalysis'    
        grd_file ='/data2/kimyy/Reanalysis/nwp_1_10_seo/ens_mean_monthly/roms_grid_10km_new.nc' ;
        scoord = [5.0 0.4 5.0 20] % theta_s theta_b hc N
        Vtransform = 1;
        Vstretching = 1;
        disp(' ')
        disp([ 'Loading ROMS grd for application: ' location])
        disp([ 'using grid file ' grd_file])
        disp(' ')
        grd = roms_get_grid(Vtransform,Vstretching,grd_file,scoord);
    end
end
end