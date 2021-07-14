function avi=png2avi( name, directory )
% avi=pn2avi( name, directory) creates avi object from all png files
%in the directory, and saves it with name 'name'
% AUTOSTOP

% get all png files
%if exist( 'directory' ) == 1
%	cd( directory );
%end

pngfiles = dir([directory,'/*png']);
numfiles = length(pngfiles);

% create avi
avi = avifile( name,'FPS',1.5);

for i = 1:numfiles
	f = pngfiles(i);
	% load frame
	disp( sprintf('%d/%d (%s)', i, numfiles, f.name ) );
	frame = imread( [directory,'/',f.name] );
	% add frame
	avi = addframe( avi, frame );
end

% finish
avi = close(avi);


