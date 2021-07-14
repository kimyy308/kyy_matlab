%% example_animation.m
%-- Create an example animation
%--     and then compress it to an mp4 file with h264 codec.
clearvars; close all;
x = (0:0.1:10);
omega = 0.25;

imgSize = [0 0 4 3];
imgRes  = 200;
scrDPI  = get(groot,'ScreenPixelsPerInch');
vidFPS  = 12;
vidPath = '/home/kimyy/';
vidName = fullfile(vidPath,'animation_file');

%% open video object
if ~exist(vidPath,'dir'),  mkdir(vidPath);    end
vidObj = VideoWriter(vidName,'Uncompressed AVI');
% vidObj = VideoWriter(vidName,'Motion JPEG AVI');	% compressed format
% vidObj.Quality = vidQty;                          % set video quality if compressed format is used
vidObj.FrameRate = vidFPS;
open(vidObj);

%% draw and write frames
figure('PaperPosition',imgSize,'Position',imgSize*scrDPI);
for t = 1:48
	plot(x,sin(x-omega*t))
	axis([0 10 -1 1])
    imgDat = print(gcf,'-RGBImage',['-r',num2str(imgRes)]);  % save image to variable 'image'
%     print(gcf,'foo','-dpng',['-r',num2str(imgRes)]);  % in case, save image to 'foo.png'
    writeVideo(vidObj,imgDat);
end
%% close video object
close(vidObj);

%% compress video file
%-- CAUTION: If you want to compress your video, image height (pixel) must be divisible by 2.
shared_lib_path = '/usr/local/lib64/lib';
[status,cmdout]= system(['LD_LIBRARY_PATH=',shared_lib_path,' ffmpeg -y -i ',vidName,'.avi -c:v h264 -an -pix_fmt yuv420p ',vidName,'.mp4']);

if status ~= 0, display(cmdout);
% else,           delete([out_name,'.avi']);   % delete uncompressed file if compression was successful
end

