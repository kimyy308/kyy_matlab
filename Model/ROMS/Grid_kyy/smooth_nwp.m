clear all
close all
dropboxpath='C:\Users\KYY\Dropbox';
addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));

romstools_param
warning off

% smoothname=['D:\MEPL\project\NWP\make_init\grid\roms_grid_combine2_smooth10.nc'];
rawname=['D:\MEPL\project\NWP\make_init\grid\roms_grid_combine2_edited_nwp_kyy.nc'];


disp(' ')
disp([' Start smoothing the grid: ',smoothname])
disp(' ')
disp([' Title: ',ROMS_title])
disp(' ')
disp([' Resolution: 1/',num2str(1/dl),' deg'])


%
%  Smooth the topography
%
% read not smoothed grid
nc2=netcdf('D:\MEPL\project\NWP\make_init\grid\roms_grid_combine2_not_smooth.nc','write');
h_raw=nc2{'h'}(:);
lon_rho=nc2{'lon_rho'}(:);
lat_rho=nc2{'lat_rho'}(:);
close(nc2);
nc3=netcdf(rawname);
h_0_8=nc3{'h'}(:);
close(nc3);

% % nc4=netcdf('D:\MEPL\project\NWP\make_init\grid\merged_etopo.nc');
% % h_merged=nc4{'depth_merged'}(:);
% % lon_merged=nc4{'lon'}(:);
% % lat_merged=nc4{'lat'}(:);
% % h_merged2 = griddata(lon_rho,lat_rho,h_merged,lon_merged,lat_merged);
% % close(nc4);

% load('D:\MEPL\project\NWP\KorBathy30s\etopo_merged.mat');
% depth_merged(find(depth_merged(:,:)<-10000))=depth_merged(find(depth_merged(:,:)<-10000))/1000.;
% lon_merged=xx;
% lat_merged=yy;
% h_merged =
% griddata(lon_merged,lat_merged,depth_merged,lon_rho(1,:),lat_rho(:,1));

nc4=netcdf('D:\MEPL\project\NWP\make_init\grid\merged_interp_roms.nc');
h_merged2=nc4{'H_MERGED4'}(:);
% lon_merged=nc4{'lon'}(:);
% lat_merged=nc4{'lat'}(:);
% h_merged2 = griddata(lon_rho,lat_rho,h_merged,lon_merged,lat_merged);
close(nc4);


nc=netcdf(smoothname,'write');
% h=nc{'h'}(:);
h=h_raw;
maskr=nc{'mask_rho'}(:);

% % for smooth 1
% h=smoothgrid(h,maskr,hmin,hmax_coast,hmax,...
%              rtarget,n_filter_deep_topo,n_filter_final);

% % for smooth 2
% % % h(y,x)
% h1=h(1:200,1:200);
% h2=h(701:900,701:900);
% maskr1=maskr(1:200,1:200);
% maskr2=maskr(701:900,701:900);
% for i=1:5
% h1=smoothgrid(h1,maskr1,hmin,hmax_coast,hmax,...
%              rtarget,n_filter_deep_topo,n_filter_final);
% h2=smoothgrid(h2,maskr2,hmin,hmax_coast,hmax,...
%              rtarget,n_filter_deep_topo,n_filter_final);
% end
% h(1:200,1:200)=h1;
% h(701:900,701:900)=h2;

% % for smooth3         
% h=smoothgrid(h,maskr,hmin,hmax_coast,hmax,...
%              rtarget,n_filter_deep_topo,n_filter_final);
% h(650:920,650:980)=h_0_95(650:920,650:980);

% % for smooth4         
% h=smoothgrid(h,maskr,hmin,hmax_coast,hmax,...
%              rtarget,n_filter_deep_topo,n_filter_final);
% h(:,:)=h_0_95(:,:);

% % % % for smooth5     
% % rtarget=0.5;
% % h=smoothgrid(h,maskr,hmin,hmax_coast,hmax,...
% %              rtarget,n_filter_deep_topo,n_filter_final);
% % for i=1:10
% %     rtarget=0.2
% %     h(1:170,1:250)=smoothgrid(h(1:170,1:250),maskr(1:170,1:250),hmin,hmax_coast,hmax,...
% %                  rtarget,n_filter_deep_topo,n_filter_final);
% %     h(701:920,701:980)=smoothgrid(h(701:920,701:980),maskr(701:920,701:980),hmin,hmax_coast,hmax,...
% %                  rtarget,n_filter_deep_topo,n_filter_final);
% %      h(701:800,701:850)=smoothgrid(h(701:800,701:850),maskr(701:800,701:850),hmin,hmax_coast,hmax,...
% %      rtarget,n_filter_deep_topo,n_filter_final);
% % end

% % for smooth6     

% % h=h_merged2;
% % rtarget=0.5;
% % h=smoothgrid(h,maskr,hmin,hmax_coast,hmax,...
% %              rtarget,n_filter_deep_topo,n_filter_final);
% % for i=1:10
% %     rtarget=0.2
% %     h(1:170,1:250)=smoothgrid(h(1:170,1:250),maskr(1:170,1:250),hmin,hmax_coast,hmax,...
% %                  rtarget,n_filter_deep_topo,n_filter_final);
% %     h(701:920,701:980)=smoothgrid(h(701:920,701:980),maskr(701:920,701:980),hmin,hmax_coast,hmax,...
% %                  rtarget,n_filter_deep_topo,n_filter_final);
% %      h(701:800,701:850)=smoothgrid(h(701:800,701:850),maskr(701:800,701:850),hmin,hmax_coast,hmax,...
% %      rtarget,n_filter_deep_topo,n_filter_final);
% % end

% % % % for smooth7     
% % h=h_merged2;
% % % % korea strait
% % h(380:510,260:330)=h(380:510,260:330)+20;
% % h(420:480,330:380)=h(420:480,330:380)+20;
% % % % tsugaru strait
% % maskr(604,508)=1;  maskr(604,509)=1;
% % maskr(611,506)=1;  maskr(611,507)=1; maskr(611,508)=1;
% % maskr(612,520)=1;
% % maskr(611,519)=1;  maskr(611,520)=1; maskr(611,521)=1;
% % maskr(617,520)=1;  maskr(617,521)=1;
% % nc{'mask_rho'}(:) = maskr
% % h(590:620,500:530)=h(590:620,500:530)+20;
% % % % soya strait
% % h(715:735,533:546)=h(715:735,533:546)+20;
% % % % whole domain smoothing
% % rtarget=0.8;
% % h=smoothgrid(h,maskr,hmin,hmax_coast,hmax,...
% %              rtarget,n_filter_deep_topo,n_filter_final);
% % % % oyashio
% % for i=1:5
% %     rtarget=0.2
% %     h(701:920,701:980)=smoothgrid(h(701:920,701:980),maskr(701:920,701:980),hmin,hmax_coast,hmax,...
% %                  rtarget,n_filter_deep_topo,n_filter_final);
% %      h(701:800,701:850)=smoothgrid(h(701:800,701:850),maskr(701:800,701:850),hmin,hmax_coast,hmax,...
% %      rtarget,n_filter_deep_topo,n_filter_final);
% % end
% % for i=1:3
% %     rtarget=0.2
% %     h(1:170,1:250)=smoothgrid(h(1:170,1:250),maskr(1:170,1:250),hmin,hmax_coast,hmax,...
% %                  rtarget,n_filter_deep_topo,n_filter_final);
% % end
% % % % korea strait smoothing
% % rtarget=0.2
% % h(450:510,280:400)=smoothgrid(h(450:510,280:400),maskr(450:510,280:400),hmin,hmax_coast,hmax,...
% %                rtarget,n_filter_deep_topo,n_filter_final);
% % h(420:510,280:400)=smoothgrid(h(420:510,280:400),maskr(420:510,280:400),hmin,hmax_coast,hmax,...
% %                rtarget,n_filter_deep_topo,n_filter_final);
% % h(360:510,260:330)=smoothgrid(h(360:510,260:330),maskr(360:510,260:330),hmin,hmax_coast,hmax,...
% %                rtarget,n_filter_deep_topo,n_filter_final);



% % % % for smooth8
% % h=h_merged2;
% % % % korea strait
% % h(380:530,260:330)=h(380:530,260:330)+20;
% % h(420:480,330:380)=h(420:480,330:380)+20;
% % % % tsugaru strait
% % maskr(604,508)=1;  maskr(604,509)=1;
% % maskr(611,506)=1;  maskr(611,507)=1; maskr(611,508)=1;
% % maskr(612,520)=1;
% % maskr(611,519)=1;  maskr(611,520)=1; maskr(611,521)=1;
% % maskr(617,520)=1;  maskr(617,521)=1;
% % nc{'mask_rho'}(:) = maskr
% % h(590:620,500:530)=h(590:620,500:530)+20;
% % % % soya strait
% % h(715:735,533:546)=h(715:735,533:546)+20;
% % % % whole domain smoothing
% % rtarget=0.8;
% % h=smoothgrid(h,maskr,hmin,hmax_coast,hmax,...
% %              rtarget,n_filter_deep_topo,n_filter_final);
% % % % oyashio
% % for i=1:5
% %     rtarget=0.2
% %     h(701:920,701:980)=smoothgrid(h(701:920,701:980),maskr(701:920,701:980),hmin,hmax_coast,hmax,...
% %                  rtarget,n_filter_deep_topo,n_filter_final);
% %      h(701:800,701:850)=smoothgrid(h(701:800,701:850),maskr(701:800,701:850),hmin,hmax_coast,hmax,...
% %      rtarget,n_filter_deep_topo,n_filter_final);
% % end
% % for i=1:3
% %     rtarget=0.2
% %     h(1:170,1:250)=smoothgrid(h(1:170,1:250),maskr(1:170,1:250),hmin,hmax_coast,hmax,...
% %                  rtarget,n_filter_deep_topo,n_filter_final);
% % end
% % % % korea strait smoothing
% % rtarget=0.2
% % h(450:530,280:400)=smoothgrid(h(450:530,280:400),maskr(450:530,280:400),hmin,hmax_coast,hmax,...
% %                rtarget,n_filter_deep_topo,n_filter_final);
% % h(420:530,280:400)=smoothgrid(h(420:530,280:400),maskr(420:530,280:400),hmin,hmax_coast,hmax,...
% %                rtarget,n_filter_deep_topo,n_filter_final);
% % h(360:530,260:330)=smoothgrid(h(360:530,260:330),maskr(360:530,260:330),hmin,hmax_coast,hmax,...
% %                rtarget,n_filter_deep_topo,n_filter_final);
           
           

% % % % for smooth13
% % h=h_merged2;
% % 
% % % % tsugaru strait
% % maskr(604,508)=1;  maskr(604,509)=1;
% % maskr(611,506)=1;  maskr(611,507)=1; maskr(611,508)=1;
% % maskr(612,520)=1;
% % maskr(611,519)=1;  maskr(611,520)=1; maskr(611,521)=1;
% % maskr(617,520)=1;  maskr(617,521)=1;
% % nc{'mask_rho'}(:) = maskr
% % 
% % % % whole domain smoothing
% % rtarget=0.8;
% % h=smoothgrid(h,maskr,hmin,hmax_coast,hmax,...
% %              rtarget,n_filter_deep_topo,n_filter_final);
% % % % oyashio
% % for i=1:1
% %     rtarget=0.2
% %     h(701:920,701:980)=smoothgrid(h(701:920,701:980),maskr(701:920,701:980),hmin,hmax_coast,hmax,...
% %                  rtarget,n_filter_deep_topo,n_filter_final);
% %      h(701:800,701:850)=smoothgrid(h(701:800,701:850),maskr(701:800,701:850),hmin,hmax_coast,hmax,...
% %      rtarget,n_filter_deep_topo,n_filter_final);
% % end
% % for i=1:5
% %     rtarget=0.2
% %     h(740:770,710:730)=smoothgrid(h(740:770,710:730),maskr(740:770,710:730),hmin,hmax_coast,hmax,...
% %                  rtarget,n_filter_deep_topo,n_filter_final);
% %      h(770:810,750:780)=smoothgrid(h(770:810,750:780),maskr(770:810,750:780),hmin,hmax_coast,hmax,...
% %      rtarget,n_filter_deep_topo,n_filter_final);
% %      h(810:860,770:810)=smoothgrid(h(810:860,770:810),maskr(810:860,770:810),hmin,hmax_coast,hmax,...
% %      rtarget,n_filter_deep_topo,n_filter_final);
% %       h(890:910,840:860)=smoothgrid(h(890:910,840:860),maskr(890:910,840:860),hmin,hmax_coast,hmax,...
% %      rtarget,n_filter_deep_topo,n_filter_final);
% % end
% % for i=1:1
% %     rtarget=0.2
% %     h(1:170,1:250)=smoothgrid(h(1:170,1:250),maskr(1:170,1:250),hmin,hmax_coast,hmax,...
% %                  rtarget,n_filter_deep_topo,n_filter_final);
% % end
% % 
% % % % korea strait
% % h(380:510,260:330)=h(380:510,260:330)+20;
% % % % tsugaru strait
% % h(590:620,500:530)=h(590:620,500:530)+20;
% % % % soya strait
% % h(715:735,533:546)=h(715:735,533:546)+20;
% % % % kuroshio path
% % h(430:500,490:560)=h_raw(430:500,490:560);


% % % % for smooth13
% % 
% % % % tsugaru strait
% % maskr(604,508)=1;  maskr(604,509)=1;
% % maskr(611,506)=1;  maskr(611,507)=1; maskr(611,508)=1;
% % maskr(612,520)=1;
% % maskr(611,519)=1;  maskr(611,520)=1; maskr(611,521)=1;
% % maskr(617,520)=1;  maskr(617,521)=1;
% % nc{'mask_rho'}(:) = maskr
% % 
% % % % whole domain smoothing
% % rtarget=0.8;
% % h=smoothgrid(h,maskr,hmin,hmax_coast,hmax,...
% %              rtarget,n_filter_deep_topo,n_filter_final);
% % % % oyashio
% % for i=1:1
% %     rtarget=0.2
% %     h(701:920,701:980)=smoothgrid(h(701:920,701:980),maskr(701:920,701:980),hmin,hmax_coast,hmax,...
% %                  rtarget,n_filter_deep_topo,n_filter_final);
% %      h(701:800,701:850)=smoothgrid(h(701:800,701:850),maskr(701:800,701:850),hmin,hmax_coast,hmax,...
% %      rtarget,n_filter_deep_topo,n_filter_final);
% % end
% % for i=1:5
% %     rtarget=0.2
% %     h(740:770,710:730)=smoothgrid(h(740:770,710:730),maskr(740:770,710:730),hmin,hmax_coast,hmax,...
% %                  rtarget,n_filter_deep_topo,n_filter_final);
% %      h(770:810,750:780)=smoothgrid(h(770:810,750:780),maskr(770:810,750:780),hmin,hmax_coast,hmax,...
% %      rtarget,n_filter_deep_topo,n_filter_final);
% %      h(810:860,770:810)=smoothgrid(h(810:860,770:810),maskr(810:860,770:810),hmin,hmax_coast,hmax,...
% %      rtarget,n_filter_deep_topo,n_filter_final);
% %       h(890:910,840:860)=smoothgrid(h(890:910,840:860),maskr(890:910,840:860),hmin,hmax_coast,hmax,...
% %      rtarget,n_filter_deep_topo,n_filter_final);
% % end
% % for i=1:1
% %     rtarget=0.2
% %     h(1:170,1:250)=smoothgrid(h(1:170,1:250),maskr(1:170,1:250),hmin,hmax_coast,hmax,...
% %                  rtarget,n_filter_deep_topo,n_filter_final);
% % end
% % 
% % % % korea strait
% % h(380:510,260:330)=h(380:510,260:330)+20;
% % % % tsugaru strait
% % h(590:620,500:530)=h(590:620,500:530)+20;
% % % % soya strait
% % h(715:735,533:546)=h(715:735,533:546)+20;
% % % % kuroshio path
% % h(430:500,490:560)=h_raw(430:500,490:560);



% % % % for test17
% % 
% % % % tsugaru strait
% % maskr(604,508)=1;  maskr(604,509)=1;
% % maskr(611,506)=1;  maskr(611,507)=1; maskr(611,508)=1;
% % maskr(612,520)=1;
% % maskr(611,519)=1;  maskr(611,520)=1; maskr(611,521)=1;
% % maskr(617,520)=1;  maskr(617,521)=1;
% % nc{'mask_rho'}(:) = maskr
% % 
% % % % korea strait
% % h(380:510,260:330)=h(380:510,260:330)+20;
% % % % tsugaru strait
% % h(590:620,500:530)=h(590:620,500:530)+20;
% % % % soya strait
% % h(715:735,533:546)=h(715:735,533:546)+20;
% % 
% % 
% % % % korea strait smoothing
% % rtarget=0.2
% % h(450:530,280:400)=smoothgrid(h(450:530,280:400),maskr(450:530,280:400),hmin,hmax_coast,hmax,...
% %                rtarget,n_filter_deep_topo,n_filter_final);
% % h(420:530,280:400)=smoothgrid(h(420:530,280:400),maskr(420:530,280:400),hmin,hmax_coast,hmax,...
% %                rtarget,n_filter_deep_topo,n_filter_final);
% % h(360:530,260:330)=smoothgrid(h(360:530,260:330),maskr(360:530,260:330),hmin,hmax_coast,hmax,...
% %                rtarget,n_filter_deep_topo,n_filter_final);
% % 
% % % % oyashio
% % for i=1:1
% %     rtarget=0.2
% %     h(701:920,701:980)=smoothgrid(h(701:920,701:980),maskr(701:920,701:980),hmin,hmax_coast,hmax,...
% %                  rtarget,n_filter_deep_topo,n_filter_final);
% %      h(701:800,701:850)=smoothgrid(h(701:800,701:850),maskr(701:800,701:850),hmin,hmax_coast,hmax,...
% %      rtarget,n_filter_deep_topo,n_filter_final);
% % end
% % for i=1:5
% %     rtarget=0.2
% %     h(740:770,710:730)=smoothgrid(h(740:770,710:730),maskr(740:770,710:730),hmin,hmax_coast,hmax,...
% %                  rtarget,n_filter_deep_topo,n_filter_final);
% %      h(770:810,750:780)=smoothgrid(h(770:810,750:780),maskr(770:810,750:780),hmin,hmax_coast,hmax,...
% %      rtarget,n_filter_deep_topo,n_filter_final);
% %      h(810:860,770:810)=smoothgrid(h(810:860,770:810),maskr(810:860,770:810),hmin,hmax_coast,hmax,...
% %      rtarget,n_filter_deep_topo,n_filter_final);
% %       h(890:910,840:860)=smoothgrid(h(890:910,840:860),maskr(890:910,840:860),hmin,hmax_coast,hmax,...
% %      rtarget,n_filter_deep_topo,n_filter_final);
% % end
% % for i=1:1
% %     rtarget=0.2
% %     h(1:170,1:250)=smoothgrid(h(1:170,1:250),maskr(1:170,1:250),hmin,hmax_coast,hmax,...
% %                  rtarget,n_filter_deep_topo,n_filter_final);
% % end
% % 
% % % % whole domain smoothing
% % rtarget=0.8;
% % h=smoothgrid(h,maskr,hmin,hmax_coast,hmax,...
% %              rtarget,n_filter_deep_topo,n_filter_final);
% % 
% % % % kuroshio path
% % h(430:500,490:560)=h_raw(430:500,490:560);




% % % for test18
% 
% % % tsugaru strait
% maskr(604,508)=1;  maskr(604,509)=1;
% maskr(611,506)=1;  maskr(611,507)=1; maskr(611,508)=1;
% maskr(612,520)=1;
% maskr(611,519)=1;  maskr(611,520)=1; maskr(611,521)=1;
% maskr(617,520)=1;  maskr(617,521)=1;
% nc{'mask_rho'}(:) = maskr
% 
% % % korea strait
% h(380:510,260:330)=h(380:510,260:330)+25;
% % % tsugaru strait
% h(590:620,500:530)=h(590:620,500:530)+25;
% % % soya strait
% h(715:735,533:546)=h(715:735,533:546)+25;
% 
% 
% % % korea strait smoothing
% rtarget=0.2
% h(450:530,280:400)=smoothgrid(h(450:530,280:400),maskr(450:530,280:400),hmin,hmax_coast,hmax,...
%                rtarget,n_filter_deep_topo,n_filter_final);
% h(420:530,280:400)=smoothgrid(h(420:530,280:400),maskr(420:530,280:400),hmin,hmax_coast,hmax,...
%                rtarget,n_filter_deep_topo,n_filter_final);
% h(360:530,260:330)=smoothgrid(h(360:530,260:330),maskr(360:530,260:330),hmin,hmax_coast,hmax,...
%                rtarget,n_filter_deep_topo,n_filter_final);
% 
% % % oyashio
% for i=1:1
%     rtarget=0.2
%     h(701:920,701:980)=smoothgrid(h(701:920,701:980),maskr(701:920,701:980),hmin,hmax_coast,hmax,...
%                  rtarget,n_filter_deep_topo,n_filter_final);
%      h(701:800,701:850)=smoothgrid(h(701:800,701:850),maskr(701:800,701:850),hmin,hmax_coast,hmax,...
%      rtarget,n_filter_deep_topo,n_filter_final);
% end
% for i=1:5
%     rtarget=0.2
%     h(740:770,710:730)=smoothgrid(h(740:770,710:730),maskr(740:770,710:730),hmin,hmax_coast,hmax,...
%                  rtarget,n_filter_deep_topo,n_filter_final);
%      h(770:810,750:780)=smoothgrid(h(770:810,750:780),maskr(770:810,750:780),hmin,hmax_coast,hmax,...
%      rtarget,n_filter_deep_topo,n_filter_final);
%      h(810:860,770:810)=smoothgrid(h(810:860,770:810),maskr(810:860,770:810),hmin,hmax_coast,hmax,...
%      rtarget,n_filter_deep_topo,n_filter_final);
%       h(890:910,840:860)=smoothgrid(h(890:910,840:860),maskr(890:910,840:860),hmin,hmax_coast,hmax,...
%      rtarget,n_filter_deep_topo,n_filter_final);
% end
% for i=1:1
%     rtarget=0.2
%     h(1:170,1:250)=smoothgrid(h(1:170,1:250),maskr(1:170,1:250),hmin,hmax_coast,hmax,...
%                  rtarget,n_filter_deep_topo,n_filter_final);
% end
% 
% % % whole domain smoothing
% rtarget=0.5;
% h=smoothgrid(h,maskr,hmin,hmax_coast,hmax,...
%              rtarget,n_filter_deep_topo,n_filter_final);
% h=smoothgrid(h,maskr,hmin,hmax_coast,hmax,...
%              rtarget,n_filter_deep_topo,n_filter_final);
% h=smoothgrid(h,maskr,hmin,hmax_coast,hmax,...
%              rtarget,n_filter_deep_topo,n_filter_final);
% % % kuroshio path
% h(430:500,490:560)=h_raw(430:500,490:560);



% % % % for test19
% % 
% % % % tsugaru strait
% % maskr(604,508)=1;  maskr(604,509)=1;
% % maskr(611,506)=1;  maskr(611,507)=1; maskr(611,508)=1;
% % maskr(612,520)=1;
% % maskr(611,519)=1;  maskr(611,520)=1; maskr(611,521)=1;
% % maskr(617,520)=1;  maskr(617,521)=1;
% % nc{'mask_rho'}(:) = maskr
% % 
% % % % korea strait
% % h(380:510,260:330)=h(380:510,260:330)+20;
% % % % tsugaru strait
% % h(590:620,500:530)=h(590:620,500:530)+0; % there are too much transport in the tsugaru strait and negative transport in the soya strait
% % % % soya strait
% % h(715:735,533:546)=h(715:735,533:546)+25;
% % 
% % 
% % % % korea strait smoothing
% % rtarget=0.2
% % h(450:530,280:400)=smoothgrid(h(450:530,280:400),maskr(450:530,280:400),hmin,hmax_coast,hmax,...
% %                rtarget,n_filter_deep_topo,n_filter_final);
% % h(420:530,280:400)=smoothgrid(h(420:530,280:400),maskr(420:530,280:400),hmin,hmax_coast,hmax,...
% %                rtarget,n_filter_deep_topo,n_filter_final);
% % h(360:530,260:330)=smoothgrid(h(360:530,260:330),maskr(360:530,260:330),hmin,hmax_coast,hmax,...
% %                rtarget,n_filter_deep_topo,n_filter_final);
% % 
% % % % oyashio
% % for i=1:1
% %     rtarget=0.2
% %     h(701:920,701:980)=smoothgrid(h(701:920,701:980),maskr(701:920,701:980),hmin,hmax_coast,hmax,...
% %                  rtarget,n_filter_deep_topo,n_filter_final);
% %      h(701:800,701:850)=smoothgrid(h(701:800,701:850),maskr(701:800,701:850),hmin,hmax_coast,hmax,...
% %      rtarget,n_filter_deep_topo,n_filter_final);
% % end
% % for i=1:5
% %     rtarget=0.2
% %     h(740:770,710:730)=smoothgrid(h(740:770,710:730),maskr(740:770,710:730),hmin,hmax_coast,hmax,...
% %                  rtarget,n_filter_deep_topo,n_filter_final);
% %      h(770:810,750:780)=smoothgrid(h(770:810,750:780),maskr(770:810,750:780),hmin,hmax_coast,hmax,...
% %      rtarget,n_filter_deep_topo,n_filter_final);
% %      h(810:860,770:810)=smoothgrid(h(810:860,770:810),maskr(810:860,770:810),hmin,hmax_coast,hmax,...
% %      rtarget,n_filter_deep_topo,n_filter_final);
% %       h(890:910,840:860)=smoothgrid(h(890:910,840:860),maskr(890:910,840:860),hmin,hmax_coast,hmax,...
% %      rtarget,n_filter_deep_topo,n_filter_final);
% % end
% % for i=1:1
% %     rtarget=0.2
% %     h(1:170,1:250)=smoothgrid(h(1:170,1:250),maskr(1:170,1:250),hmin,hmax_coast,hmax,...
% %                  rtarget,n_filter_deep_topo,n_filter_final);
% % end
% % 
% % % % whole domain smoothing
% % rtarget=0.5;
% % h=smoothgrid(h,maskr,hmin,hmax_coast,hmax,...
% %              rtarget,n_filter_deep_topo,n_filter_final);
% % % % kuroshio path
% % h(430:500,490:560)=h_raw(430:500,490:560);





% % % for test23
% 
% % % tsugaru strait
% % maskr(604,508)=1;  maskr(604,509)=1;
% % maskr(611,506)=1;  maskr(611,507)=1; maskr(611,508)=1;
% % maskr(612,520)=1;
% % maskr(611,519)=1;  maskr(611,520)=1; maskr(611,521)=1;
% % maskr(617,520)=1;  maskr(617,521)=1;
% % nc{'mask_rho'}(:) = maskr
% 
% % % korea strait
% h(380:510,260:330)=h(380:510,260:330)+10;
% % % tsugaru strait
% h(590:620,500:530)=h(590:620,500:530)+0; % there are too much transport in the tsugaru strait and negative transport in the soya strait
% % % soya strait
% h(715:735,533:546)=h(715:735,533:546)+25;
% 
% 
% % % korea strait smoothing
% rtarget=0.2
% h(450:530,280:400)=smoothgrid(h(450:530,280:400),maskr(450:530,280:400),hmin,hmax_coast,hmax,...
%                rtarget,n_filter_deep_topo,n_filter_final);
% h(420:530,280:400)=smoothgrid(h(420:530,280:400),maskr(420:530,280:400),hmin,hmax_coast,hmax,...
%                rtarget,n_filter_deep_topo,n_filter_final);
% h(360:530,260:330)=smoothgrid(h(360:530,260:330),maskr(360:530,260:330),hmin,hmax_coast,hmax,...
%                rtarget,n_filter_deep_topo,n_filter_final);
% 
% % % oyashio
% for i=1:1
%     rtarget=0.2
%     h(701:920,701:980)=smoothgrid(h(701:920,701:980),maskr(701:920,701:980),hmin,hmax_coast,hmax,...
%                  rtarget,n_filter_deep_topo,n_filter_final);
%      h(701:800,701:850)=smoothgrid(h(701:800,701:850),maskr(701:800,701:850),hmin,hmax_coast,hmax,...
%      rtarget,n_filter_deep_topo,n_filter_final);
% end
% for i=1:5
%     rtarget=0.2
%     h(740:770,710:730)=smoothgrid(h(740:770,710:730),maskr(740:770,710:730),hmin,hmax_coast,hmax,...
%                  rtarget,n_filter_deep_topo,n_filter_final);
%      h(770:810,750:780)=smoothgrid(h(770:810,750:780),maskr(770:810,750:780),hmin,hmax_coast,hmax,...
%      rtarget,n_filter_deep_topo,n_filter_final);
%      h(810:860,770:810)=smoothgrid(h(810:860,770:810),maskr(810:860,770:810),hmin,hmax_coast,hmax,...
%      rtarget,n_filter_deep_topo,n_filter_final);
%       h(890:910,840:860)=smoothgrid(h(890:910,840:860),maskr(890:910,840:860),hmin,hmax_coast,hmax,...
%      rtarget,n_filter_deep_topo,n_filter_final);
% end
% for i=1:1
%     rtarget=0.2
%     h(1:170,1:250)=smoothgrid(h(1:170,1:250),maskr(1:170,1:250),hmin,hmax_coast,hmax,...
%                  rtarget,n_filter_deep_topo,n_filter_final);
% end
% 
% % % whole domain smoothing
% rtarget=0.5;
% h=smoothgrid(h,maskr,hmin,hmax_coast,hmax,...
%              rtarget,n_filter_deep_topo,n_filter_final);
% % % kuroshio path
% h(430:500,490:560)=h_raw(430:500,490:560);




% % % for test27
% 
% % % tsugaru strait
% % maskr(604,508)=1;  maskr(604,509)=1;
% % maskr(611,506)=1;  maskr(611,507)=1; maskr(611,508)=1;
% % maskr(612,520)=1;
% % maskr(611,519)=1;  maskr(611,520)=1; maskr(611,521)=1;
% % maskr(617,520)=1;  maskr(617,521)=1;
% % nc{'mask_rho'}(:) = maskr
% 
% % % korea strait
% h(380:510,260:330)=h(380:510,260:330)+5;
% % % tsugaru strait
% h(590:620,500:530)=h(590:620,500:530)+0; % there are too much transport in the tsugaru strait and negative transport in the soya strait
% % % soya strait
% h(715:735,533:546)=h(715:735,533:546)+25;
% 
% 
% % % korea strait smoothing
% rtarget=0.2
% h(450:530,280:400)=smoothgrid(h(450:530,280:400),maskr(450:530,280:400),hmin,hmax_coast,hmax,...
%                rtarget,n_filter_deep_topo,n_filter_final);
% h(420:530,280:400)=smoothgrid(h(420:530,280:400),maskr(420:530,280:400),hmin,hmax_coast,hmax,...
%                rtarget,n_filter_deep_topo,n_filter_final);
% h(360:530,260:330)=smoothgrid(h(360:530,260:330),maskr(360:530,260:330),hmin,hmax_coast,hmax,...
%                rtarget,n_filter_deep_topo,n_filter_final);
% 
% % % onshore region smoothing
% h(160:400,140:340)=smoothgrid(h(160:400,140:340),maskr(160:400,140:340),hmin,hmax_coast,hmax,...
%                rtarget,n_filter_deep_topo,n_filter_final);
%            
% % % oyashio
% for i=1:1
%     rtarget=0.2
%     h(701:920,701:980)=smoothgrid(h(701:920,701:980),maskr(701:920,701:980),hmin,hmax_coast,hmax,...
%                  rtarget,n_filter_deep_topo,n_filter_final);
%      h(701:800,701:850)=smoothgrid(h(701:800,701:850),maskr(701:800,701:850),hmin,hmax_coast,hmax,...
%      rtarget,n_filter_deep_topo,n_filter_final);
% end
% for i=1:5
%     rtarget=0.2
%     h(740:770,710:730)=smoothgrid(h(740:770,710:730),maskr(740:770,710:730),hmin,hmax_coast,hmax,...
%                  rtarget,n_filter_deep_topo,n_filter_final);
%      h(770:810,750:780)=smoothgrid(h(770:810,750:780),maskr(770:810,750:780),hmin,hmax_coast,hmax,...
%      rtarget,n_filter_deep_topo,n_filter_final);
%      h(810:860,770:810)=smoothgrid(h(810:860,770:810),maskr(810:860,770:810),hmin,hmax_coast,hmax,...
%      rtarget,n_filter_deep_topo,n_filter_final);
%       h(890:910,840:860)=smoothgrid(h(890:910,840:860),maskr(890:910,840:860),hmin,hmax_coast,hmax,...
%      rtarget,n_filter_deep_topo,n_filter_final);
% end
% for i=1:1
%     rtarget=0.2
%     h(1:170,1:250)=smoothgrid(h(1:170,1:250),maskr(1:170,1:250),hmin,hmax_coast,hmax,...
%                  rtarget,n_filter_deep_topo,n_filter_final);
% end
% 
% % % whole domain smoothing
% rtarget=0.5;
% h=smoothgrid(h,maskr,hmin,hmax_coast,hmax,...
%              rtarget,n_filter_deep_topo,n_filter_final);
% 
%          
%          
% % % kuroshio path 139.4 ~ 142.9E (490:560) , 34.36 ~ 37.20N (430:500)
% h(430:500,490:560)=h_raw(430:500,490:560);




% % % for test30
% 
% 
% % % tsugaru strait
% % maskr(604,508)=1;  maskr(604,509)=1;
% % maskr(611,506)=1;  maskr(611,507)=1; maskr(611,508)=1;
% % maskr(612,520)=1;
% % maskr(611,519)=1;  maskr(611,520)=1; maskr(611,521)=1;
% % maskr(617,520)=1;  maskr(617,521)=1;
% % nc{'mask_rho'}(:) = maskr
% 
% % % korea strait
% h(380:510,260:330)=h(380:510,260:330)+5;
% % % tsugaru strait
% h(590:620,500:530)=h(590:620,500:530)+0; % there are too much transport in the tsugaru strait and negative transport in the soya strait
% % % soya strait
% h(715:735,533:546)=h(715:735,533:546)+25;
% 
% 
% % % korea strait smoothing
% rtarget=0.2
% h(450:530,280:400)=smoothgrid(h(450:530,280:400),maskr(450:530,280:400),hmin,hmax_coast,hmax,...
%                rtarget,n_filter_deep_topo,n_filter_final);
% h(420:530,280:400)=smoothgrid(h(420:530,280:400),maskr(420:530,280:400),hmin,hmax_coast,hmax,...
%                rtarget,n_filter_deep_topo,n_filter_final);
% h(360:530,260:330)=smoothgrid(h(360:530,260:330),maskr(360:530,260:330),hmin,hmax_coast,hmax,...
%                rtarget,n_filter_deep_topo,n_filter_final);
% 
% % % onshore region smoothing
% for i=1:2
% h(160:400,140:340)=smoothgrid(h(160:400,140:340),maskr(160:400,140:340),hmin,hmax_coast,hmax,...
%                rtarget,n_filter_deep_topo,n_filter_final);
% end   
% % % oyashio
% for i=1:1
%     rtarget=0.2
%     h(701:920,701:980)=smoothgrid(h(701:920,701:980),maskr(701:920,701:980),hmin,hmax_coast,hmax,...
%                  rtarget,n_filter_deep_topo,n_filter_final);
%      h(701:800,701:850)=smoothgrid(h(701:800,701:850),maskr(701:800,701:850),hmin,hmax_coast,hmax,...
%      rtarget,n_filter_deep_topo,n_filter_final);
% end
% for i=1:5
%     rtarget=0.2
%     h(740:770,710:730)=smoothgrid(h(740:770,710:730),maskr(740:770,710:730),hmin,hmax_coast,hmax,...
%                  rtarget,n_filter_deep_topo,n_filter_final);
%      h(770:810,750:780)=smoothgrid(h(770:810,750:780),maskr(770:810,750:780),hmin,hmax_coast,hmax,...
%      rtarget,n_filter_deep_topo,n_filter_final);
%      h(810:860,770:810)=smoothgrid(h(810:860,770:810),maskr(810:860,770:810),hmin,hmax_coast,hmax,...
%      rtarget,n_filter_deep_topo,n_filter_final);
%       h(890:910,840:860)=smoothgrid(h(890:910,840:860),maskr(890:910,840:860),hmin,hmax_coast,hmax,...
%      rtarget,n_filter_deep_topo,n_filter_final);
% end
% for i=1:1
%     rtarget=0.2
%     h(1:170,1:250)=smoothgrid(h(1:170,1:250),maskr(1:170,1:250),hmin,hmax_coast,hmax,...
%                  rtarget,n_filter_deep_topo,n_filter_final);
% end
% 
% 
% % % kuroshio separation region smoothing
% for i=1:5
% %     rtarget = 0.07
%         rtarget = 0.12 - 0.01*i
% %     h(430:500,500:700)=smoothgrid(h(430:500,500:700),maskr(430:500,500:700),hmin,hmax_coast,hmax,...
% %             rtarget,n_filter_deep_topo,n_filter_final);
% %     h(430:460,530:570)=smoothgrid(h(430:460,530:570),maskr(430:460,530:570),hmin,hmax_coast,hmax,...
% %             rtarget,n_filter_deep_topo,n_filter_final);
% %     h(441:460,500:570)=smoothgrid(h(441:460,500:570),maskr(441:460,500:570),hmin,hmax_coast,hmax,...
% %     rtarget,n_filter_deep_topo,n_filter_final);   %%fine  (iter=10)
%     h(443+i*3:460,500:570)=smoothgrid(h(443+i*3:460,500:570),maskr(443+i*3:460,500:570),hmin,hmax_coast,hmax,...
%     rtarget,n_filter_deep_topo,n_filter_final);  
% end
% 
% % % whole domain smoothing
% rtarget=0.5;
% h=smoothgrid(h,maskr,hmin,hmax_coast,hmax,...
%              rtarget,n_filter_deep_topo,n_filter_final);
% 
%          
%          
% % % % kuroshio path 139.4 ~ 142.9E (490:560) , 34.36 ~ 37.20N (430:500)
% % h(430:500,490:560)=h_raw(430:500,490:560);


% % % for test34
% 
% % % tsugaru strait
% % maskr(604,508)=1;  maskr(604,509)=1;
% % maskr(611,506)=1;  maskr(611,507)=1; maskr(611,508)=1;
% % maskr(612,520)=1;
% % maskr(611,519)=1;  maskr(611,520)=1; maskr(611,521)=1;
% % maskr(617,520)=1;  maskr(617,521)=1;
% % nc{'mask_rho'}(:) = maskr
% 
% % % korea strait
% h(380:510,260:330)=h(380:510,260:330)+10;
% % % tsugaru strait
% h(590:620,500:530)=h(590:620,500:530)+0; % there are too much transport in the tsugaru strait and negative transport in the soya strait
% % % soya strait
% h(715:735,533:546)=h(715:735,533:546)+25;
% 
% 
% % % korea strait smoothing
% rtarget=0.2
% h(450:530,280:400)=smoothgrid(h(450:530,280:400),maskr(450:530,280:400),hmin,hmax_coast,hmax,...
%                rtarget,n_filter_deep_topo,n_filter_final);
% h(420:530,280:400)=smoothgrid(h(420:530,280:400),maskr(420:530,280:400),hmin,hmax_coast,hmax,...
%                rtarget,n_filter_deep_topo,n_filter_final);
% h(360:530,260:330)=smoothgrid(h(360:530,260:330),maskr(360:530,260:330),hmin,hmax_coast,hmax,...
%                rtarget,n_filter_deep_topo,n_filter_final);
% 
% % % oyashio
% for i=1:1
%     rtarget=0.2
%     h(701:920,701:980)=smoothgrid(h(701:920,701:980),maskr(701:920,701:980),hmin,hmax_coast,hmax,...
%                  rtarget,n_filter_deep_topo,n_filter_final);
%      h(701:800,701:850)=smoothgrid(h(701:800,701:850),maskr(701:800,701:850),hmin,hmax_coast,hmax,...
%      rtarget,n_filter_deep_topo,n_filter_final);
% end
% for i=1:5
%     rtarget=0.2
%     h(740:770,710:730)=smoothgrid(h(740:770,710:730),maskr(740:770,710:730),hmin,hmax_coast,hmax,...
%                  rtarget,n_filter_deep_topo,n_filter_final);
%      h(770:810,750:780)=smoothgrid(h(770:810,750:780),maskr(770:810,750:780),hmin,hmax_coast,hmax,...
%      rtarget,n_filter_deep_topo,n_filter_final);
%      h(810:860,770:810)=smoothgrid(h(810:860,770:810),maskr(810:860,770:810),hmin,hmax_coast,hmax,...
%      rtarget,n_filter_deep_topo,n_filter_final);
%       h(890:910,840:860)=smoothgrid(h(890:910,840:860),maskr(890:910,840:860),hmin,hmax_coast,hmax,...
%      rtarget,n_filter_deep_topo,n_filter_final);
% end
% 
% % izu-ogasawara ridge smoothing
% rtarget=0.3
%   h(340:440,460:540)=smoothgrid(h(340:440,460:540),maskr(340:440,460:540),hmin,hmax_coast,hmax,...
%  rtarget,n_filter_deep_topo,n_filter_final);
% % Luzon strait smoothing
% for i=1:1
%     rtarget=0.4
%     h(1:170,1:250)=smoothgrid(h(1:170,1:250),maskr(1:170,1:250),hmin,hmax_coast,hmax,...
%                  rtarget,n_filter_deep_topo,n_filter_final);
% end
% 
% % % whole domain smoothing
% rtarget=0.5;
% h=smoothgrid(h,maskr,hmin,hmax_coast,hmax,...
%              rtarget,n_filter_deep_topo,n_filter_final);
% 
% % % kuroshio path
% h(430:550,490:560)=h_raw(430:550,490:560);
% 
% %
% %  Write it down
% %
% disp(' ')
% disp(' Write it down...')
% nc{'h'}(:)=h;
% close(nc);




% % % for test35
% 
% % % korea strait
% h(380:510,260:330)=h(380:510,260:330)+10;
% % % tsugaru strait
% h(590:620,500:530)=h(590:620,500:530)+0; % there are too much transport in the tsugaru strait and negative transport in the soya strait
% % % soya strait
% h(715:735,533:546)=h(715:735,533:546)+25;
% 
% 
% % % korea strait smoothing
% rtarget=0.2
% h(450:530,280:400)=smoothgrid(h(450:530,280:400),maskr(450:530,280:400),hmin,hmax_coast,hmax,...
%                rtarget,n_filter_deep_topo,n_filter_final);
% h(420:530,280:400)=smoothgrid(h(420:530,280:400),maskr(420:530,280:400),hmin,hmax_coast,hmax,...
%                rtarget,n_filter_deep_topo,n_filter_final);
% h(360:530,260:330)=smoothgrid(h(360:530,260:330),maskr(360:530,260:330),hmin,hmax_coast,hmax,...
%                rtarget,n_filter_deep_topo,n_filter_final);
% 
% % % oyashio
% for i=1:1
%     rtarget=0.2
%     h(701:920,701:980)=smoothgrid(h(701:920,701:980),maskr(701:920,701:980),hmin,hmax_coast,hmax,...
%                  rtarget,n_filter_deep_topo,n_filter_final);
%      h(701:800,701:850)=smoothgrid(h(701:800,701:850),maskr(701:800,701:850),hmin,hmax_coast,hmax,...
%      rtarget,n_filter_deep_topo,n_filter_final);
% end
% for i=1:5
%     rtarget=0.2
%     h(740:770,710:730)=smoothgrid(h(740:770,710:730),maskr(740:770,710:730),hmin,hmax_coast,hmax,...
%                  rtarget,n_filter_deep_topo,n_filter_final);
%      h(770:810,750:780)=smoothgrid(h(770:810,750:780),maskr(770:810,750:780),hmin,hmax_coast,hmax,...
%      rtarget,n_filter_deep_topo,n_filter_final);
%      h(810:860,770:810)=smoothgrid(h(810:860,770:810),maskr(810:860,770:810),hmin,hmax_coast,hmax,...
%      rtarget,n_filter_deep_topo,n_filter_final);
%       h(890:910,840:860)=smoothgrid(h(890:910,840:860),maskr(890:910,840:860),hmin,hmax_coast,hmax,...
%      rtarget,n_filter_deep_topo,n_filter_final);
% end
% 
% % izu-ogasawara ridge smoothing (south japan)
% rtarget=0.1
%   h(340:460,370:540)=smoothgrid(h(340:460,370:540),maskr(340:460,370:540),hmin,hmax_coast,hmax,...
%  rtarget,n_filter_deep_topo,n_filter_final);
% 
% % Luzon strait smoothing
% for i=1:1
%     rtarget=0.4
%     h(1:170,1:250)=smoothgrid(h(1:170,1:250),maskr(1:170,1:250),hmin,hmax_coast,hmax,...
%                  rtarget,n_filter_deep_topo,n_filter_final);
% end
% 
% for i=1:1
%     rtarget=0.2
%     h(60:170,100:150)=smoothgrid(h(60:170,100:150),maskr(60:170,100:150),hmin,hmax_coast,hmax,...
%                  rtarget,n_filter_deep_topo,n_filter_final);
% end
% 
% 
% 
% % % whole domain smoothing
% rtarget=0.5;
% h=smoothgrid(h,maskr,hmin,hmax_coast,hmax,...
%              rtarget,n_filter_deep_topo,n_filter_final);
% 
% % % kuroshio path
% h(460:550,490:560)=h_raw(460:550,490:560);
% 
% %
% %  Write it down
% %
% disp(' ')
% disp(' Write it down...')
% nc{'h'}(:)=h;
% close(nc);

% % % for test36
% 
% % % korea strait
% h(380:510,260:330)=h(380:510,260:330)+10;
% % % tsugaru strait
% h(590:620,500:530)=h(590:620,500:530)+0; % there are too much transport in the tsugaru strait and negative transport in the soya strait
% % % soya strait
% h(715:735,533:546)=h(715:735,533:546)+25;
% 
% 
% % % korea strait smoothing
% rtarget=0.2
% h(450:530,280:400)=smoothgrid(h(450:530,280:400),maskr(450:530,280:400),hmin,hmax_coast,hmax,...
%                rtarget,n_filter_deep_topo,n_filter_final);
% h(420:530,280:400)=smoothgrid(h(420:530,280:400),maskr(420:530,280:400),hmin,hmax_coast,hmax,...
%                rtarget,n_filter_deep_topo,n_filter_final);
% h(360:530,260:330)=smoothgrid(h(360:530,260:330),maskr(360:530,260:330),hmin,hmax_coast,hmax,...
%                rtarget,n_filter_deep_topo,n_filter_final);
% 
% % % oyashio
% for i=1:1
%     rtarget=0.2
%     h(701:920,701:980)=smoothgrid(h(701:920,701:980),maskr(701:920,701:980),hmin,hmax_coast,hmax,...
%                  rtarget,n_filter_deep_topo,n_filter_final);
%      h(701:800,701:850)=smoothgrid(h(701:800,701:850),maskr(701:800,701:850),hmin,hmax_coast,hmax,...
%      rtarget,n_filter_deep_topo,n_filter_final);
% end
% for i=1:5
%     rtarget=0.2
%     h(740:770,710:730)=smoothgrid(h(740:770,710:730),maskr(740:770,710:730),hmin,hmax_coast,hmax,...
%                  rtarget,n_filter_deep_topo,n_filter_final);
%      h(770:810,750:780)=smoothgrid(h(770:810,750:780),maskr(770:810,750:780),hmin,hmax_coast,hmax,...
%      rtarget,n_filter_deep_topo,n_filter_final);
%      h(810:860,770:810)=smoothgrid(h(810:860,770:810),maskr(810:860,770:810),hmin,hmax_coast,hmax,...
%      rtarget,n_filter_deep_topo,n_filter_final);
%       h(890:910,840:860)=smoothgrid(h(890:910,840:860),maskr(890:910,840:860),hmin,hmax_coast,hmax,...
%      rtarget,n_filter_deep_topo,n_filter_final);
% end
% 
% % izu-ogasawara ridge smoothing (south japan)
% rtarget=0.7
% n_filter_deep_topo=30;
% %   h(340:460,370:540)=smoothgrid(h(340:460,370:540),maskr(340:460,370:540),hmin,hmax_coast,hmax,...
% %  rtarget,n_filter_deep_topo,n_filter_final);
% for i=470:520
%     for j = 430:440
%          if (h(j,i)<2000) 
%              h(j,i)= 2000;
%          end
%     end
% end
% for i=470:520
%     for j = 400:430
%          if (h(j,i)<2000 + (430-j)/30*500) 
%              h(j,i)= 2000 + (430-j)/30*500;
%          end
%     end
% end
% for i=470:520
%     for j = 330:400
%          if (h(j,i)<2500 + (400-j)/70*500) 
%              h(j,i)= 2500 + (400-j)/70*500;
%          end
%     end
% end
% for i=470:520
%     for j = 250:330
%          if (h(j,i)<2000 + (j-250)/80*1000) 
%              h(j,i)= 2000 + (j-250)/80*1000;
%          end
%     end
% end
%   h(340:440,460:540)=smoothgrid(h(340:440,460:540),maskr(340:440,460:540),hmin,hmax_coast,hmax,...
%  rtarget,n_filter_deep_topo,n_filter_final);
% n_filter_deep_topo=0;
% 
% 
% % Luzon strait smoothing
% n_filter_deep_topo=50;
% for i=1:1
%     rtarget=0.5
%     h(70:130,110:140)=smoothgrid(h(70:130,110:140),maskr(70:130,110:140),hmin,hmax_coast,hmax,...
%                  rtarget,n_filter_deep_topo,n_filter_final);
% end
% for i=1:1
%     rtarget=0.5
%     h(130:170,90:150)=smoothgrid(h(130:170,90:150),maskr(130:170,90:150),hmin,hmax_coast,hmax,...
%                  rtarget,n_filter_deep_topo,n_filter_final);
% end
% % n_filter_deep_topo=1000;
% % for i=1:1
% %     rtarget=0.5
% %     h(80:120,110:140)=smoothgrid(h(80:120,110:140),maskr(80:120,110:140),hmin,hmax_coast,hmax,...
% %                  rtarget,n_filter_deep_topo,n_filter_final);
% % end
% for i=110:140
%     for j = 85:150
%         if (h(j,i)<2000)
%             h(j,i) = 2000;
%         end
%     end
% end
% 
% n_filter_deep_topo=50;
% for i=1:1
%     rtarget=0.5
%     h(60:170,110:150)=smoothgrid(h(60:170,110:150),maskr(60:170,110:150),hmin,hmax_coast,hmax,...
%                  rtarget,n_filter_deep_topo,n_filter_final);
% end
% n_filter_deep_topo=0;
% for i=1:1
%     rtarget=0.5
%     h(1:170,1:250)=smoothgrid(h(1:170,1:250),maskr(1:170,1:250),hmin,hmax_coast,hmax,...
%                  rtarget,n_filter_deep_topo,n_filter_final);
% end
% 
% % % whole domain smoothing
% rtarget=0.6;
% h=smoothgrid(h,maskr,hmin,hmax_coast,hmax,...
%              rtarget,n_filter_deep_topo,n_filter_final);
% 
% % % kuroshio path
% h(460:550,490:560)=h_raw(460:550,490:560);
% 
% %
% %  Write it down
% %
% disp(' ')
% disp(' Write it down...')
% nc{'h'}(:)=h;
% close(nc);


% % for test 37

% % korea strait
h(380:510,260:330)=h(380:510,260:330)+10;
% % tsugaru strait
h(590:620,500:530)=h(590:620,500:530)+0; % there are too much transport in the tsugaru strait and negative transport in the soya strait
% % soya strait
h(715:735,533:546)=h(715:735,533:546)+25;


% % korea strait smoothing
rtarget=0.2
h(450:530,280:400)=smoothgrid(h(450:530,280:400),maskr(450:530,280:400),hmin,hmax_coast,hmax,...
               rtarget,n_filter_deep_topo,n_filter_final);
h(420:530,280:400)=smoothgrid(h(420:530,280:400),maskr(420:530,280:400),hmin,hmax_coast,hmax,...
               rtarget,n_filter_deep_topo,n_filter_final);
h(360:530,260:330)=smoothgrid(h(360:530,260:330),maskr(360:530,260:330),hmin,hmax_coast,hmax,...
               rtarget,n_filter_deep_topo,n_filter_final);

% % oyashio
n_filter_final = 5;
n_filter_deep_topo=30;
for i=1:1
    rtarget=0.2
    h(701:920,701:980)=smoothgrid(h(701:920,701:980),maskr(701:920,701:980),hmin,hmax_coast,hmax,...
                 rtarget,n_filter_deep_topo,n_filter_final);
     h(701:800,701:850)=smoothgrid(h(701:800,701:850),maskr(701:800,701:850),hmin,hmax_coast,hmax,...
     rtarget,n_filter_deep_topo,n_filter_final);
end
for i=1:1
    rtarget=0.2
    h(740:770,710:730)=smoothgrid(h(740:770,710:730),maskr(740:770,710:730),hmin,hmax_coast,hmax,...
                 rtarget,n_filter_deep_topo,n_filter_final);
     h(770:810,750:780)=smoothgrid(h(770:810,750:780),maskr(770:810,750:780),hmin,hmax_coast,hmax,...
     rtarget,n_filter_deep_topo,n_filter_final);
     h(810:860,770:810)=smoothgrid(h(810:860,770:810),maskr(810:860,770:810),hmin,hmax_coast,hmax,...
     rtarget,n_filter_deep_topo,n_filter_final);
      h(890:910,840:860)=smoothgrid(h(890:910,840:860),maskr(890:910,840:860),hmin,hmax_coast,hmax,...
     rtarget,n_filter_deep_topo,n_filter_final);
end
n_filter_deep_topo=30;
n_filter_final = 5;


% izu-ogasawara ridge smoothing (south japan)
rtarget=0.7
n_filter_deep_topo=30;
n_filter_final=3;
%   h(340:460,370:540)=smoothgrid(h(340:460,370:540),maskr(340:460,370:540),hmin,hmax_coast,hmax,...
%  rtarget,n_filter_deep_topo,n_filter_final);
for i=470:500
    for j = 400:415
         if (h(j,i)>300) 
             h(j,i)= h(j,i)+200.0;
         end
    end
end
for i=470:500
    for j = 415:427
%          if (h(j,i)<2000 + (430-j)/30*500) 
%              h(j,i)= 2000 + (430-j)/30*500;
%          if (h(j,i)>300) 
%              h(j,i)= h(j,i)+2000.0;
             h(j,i)=2000.0;
%          end
%          end
    end
end
% for i=470:520
%     for j = 330:400
%          if (h(j,i)<2500 + (400-j)/70*500) 
%              h(j,i)= 2500 + (400-j)/70*500;
%          end
%     end
% end
% for i=470:520
%     for j = 250:330
%          if (h(j,i)<2000 + (j-250)/80*1000) 
%              h(j,i)= 2000 + (j-250)/80*1000;
%          end
%     end
% end
  h(340:440,460:540)=smoothgrid(h(340:440,460:540),maskr(340:440,460:540),hmin,hmax_coast,hmax,...
 rtarget,n_filter_deep_topo,n_filter_final);
n_filter_deep_topo=0;
n_filter_final=3;



% Luzon strait smoothing
n_filter_final = 10;
for i=100:150
    for j = 75:150
        if (h(j,i)<2000)
            h(j,i) = h(j,i) + 50;
        end
    end
end

n_filter_deep_topo=100;
for i=1:5
    rtarget=0.8
    h(70:125,70:140)=smoothgrid(h(70:125,70:140),maskr(70:125,70:140),hmin,hmax_coast,hmax,...
                 rtarget,n_filter_deep_topo,n_filter_final);
end
for i=1:1
    rtarget=0.8
    h(130:170,70:150)=smoothgrid(h(130:170,70:150),maskr(130:170,70:150),hmin,hmax_coast,hmax,...
                 rtarget,n_filter_deep_topo,n_filter_final);
end

% n_filter_deep_topo=1000;
% for i=1:1
%     rtarget=0.5
%     h(80:120,110:140)=smoothgrid(h(80:120,110:140),maskr(80:120,110:140),hmin,hmax_coast,hmax,...
%                  rtarget,n_filter_deep_topo,n_filter_final);
% end



n_filter_deep_topo=50;
for i=1:5
    rtarget=0.7
    h(60:170,70:150)=smoothgrid(h(60:170,70:150),maskr(60:170,70:150),hmin,hmax_coast,hmax,...
                 rtarget,n_filter_deep_topo,n_filter_final);
end
n_filter_deep_topo=0;
for i=1:1
    rtarget=0.5
    h(1:170,1:250)=smoothgrid(h(1:170,1:250),maskr(1:170,1:250),hmin,hmax_coast,hmax,...
                 rtarget,n_filter_deep_topo,n_filter_final);
end

n_filter_final = 5;


% % whole domain smoothing
rtarget=0.6;
h=smoothgrid(h,maskr,hmin,hmax_coast,hmax,...
             rtarget,n_filter_deep_topo,n_filter_final);

% % kuroshio path
h(460:550,490:560)=h_raw(460:550,490:560);
h(250:350,470:520)=h_raw(250:350,470:520);

%
%  Write it down
%
disp(' ')
disp(' Write it down...')
nc{'h'}(:)=h;
close(nc);









% 
%  plot
% 
% pcolor(h)

% figure;

% pcolor(h(700:920,700:980))
% pcolor(h(700:800,700:850))


lonlat = [132 143 27 35];  %%kuro, japan south
m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
m_grid('fontsize',20, 'box', 'fancy', 'tickdir', 'in');   %% for nwp = 25, for es = 20
hold on;
m_pcolor(lon_rho,lat_rho,h);
shading interp;
m_gshhs_i('color','k')
m_gshhs_i('patch',[.8 .8 .8]);  
% set colorbar 
ccc= colorbar;
load D:\MEPL\project\SSH\중간보고\smooth13_vtvs\kyy_plot_subroutine\jet_mod
colormap(jet_mod);
caxis([0 5000]);
hold off;

% figure;
% lonlat = [117 124 16 23]; %% luzon
% m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
% m_grid('fontsize',20, 'box', 'fancy', 'tickdir', 'in');   %% for nwp = 25, for es = 20
% hold on;
% m_pcolor(lon_rho,lat_rho,h);
% shading interp;
% % set colorbar 
% ccc= colorbar;
% load D:\MEPL\project\SSH\중간보고\smooth13_vtvs\kyy_plot_subroutine\jet_mod
% colormap(jet_mod);
% caxis([0 5000]);
% hold off;


% figure;
% lonlat = [127 144 33 52]; %% East Sea
% m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
% m_grid('fontsize',20, 'box', 'fancy', 'tickdir', 'in');   %% for nwp = 25, for es = 20
% hold on;
% m_pcolor(lon_rho,lat_rho,h);
% shading interp;
% % set colorbar 
% ccc= colorbar;
% load D:\MEPL\project\SSH\중간보고\smooth13_vtvs\kyy_plot_subroutine\jet_mod
% colormap(jet_mod);
% caxis([0 5000]);
% hold off;


% figure;
% lonlat = [115 164 15 52]; %% nwp
% m_proj('mercator','lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
% m_grid('fontsize',20, 'box', 'fancy', 'tickdir', 'in');   %% for nwp = 25, for es = 20
% hold on;
% m_pcolor(lon_rho,lat_rho,h);
% shading interp;
% % set colorbar 
% ccc= colorbar;
% load D:\MEPL\project\SSH\중간보고\smooth13_vtvs\kyy_plot_subroutine\jet_mod
% colormap(jet_mod);
% caxis([0 5000]);
% hold off;



figure; 
pcolor(450:550,350:450,h(350:450,450:550))  %% izu ridge
shading interp;
colormap jet;
colorbar;
% 
% figure;
% pcolor(h(1:170,1:250))  %% luzon
% shading interp
% colormap jet;
% colorbar;
% 
% figure;
% pcolor(h);
% shading interp;
% colormap jet;
% colorbar;