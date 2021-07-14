%----This is the code to make cseof.c script file for same grid but
%differnet variables

% %---------STEP1. set the variables-----------------------------------%
% name_layer = char('uv','air','hgt','omg'); %variable name to n layer
% nn = char('slp','sst2','ssh2'); %variable name to one layer
% 
% for i = 1:4 %number of name layer
%     n1 = n(i,:);
%     for j = 1 : 10 %layer number
%         name = strcat(n1,'_',num2str(j));
%         nn = char(nn,name);
%     end
% end

load('STEP2_1_output.mat')

%----------STEP2. make script file-----------------------------------%
name = cell(28);
name{1} = '#!/bin/csh';     %use csh shell
% for i = 1:length(nn);
name{2} = char('');
name{3} = char( ['gfortran -fno-backtrace -o cseof ', workdir, '/programs/eigen/cseof.f']); %compile the eigen file
name{4} = char( '' );
name{5} = char( 'cat >! cseof.com <<ENDc'  );
name{6} = char('0');  %% job number (only 0)
name{7} = char( strcat('./../cseofs/pct_', 'NorESM1-M_', varname, '.dat')  );     %load variable file
name{8} = char( '(5e13.5)' );
name{9} = char( [maxmodenum, ' 1']  ); %max.mode number(same to eigen.c)
name{10} = char( num2str(tlen)); %time number
name{11} = char(  '1'    ); %% time index(=1) or space index(=2) first
LVnumber = length(time)/length(inputyear);
name{12} = char(  num2str(LVnumber) );   %want period(LV number)
spectralpoints = floor(LVnumber/2);
name{13} = char(  num2str(spectralpoints)  );   %want period/2 (number of spectral points)
name{14} = char( '1              ' );  %% interval subdivisions for integrations
name{15} = char( '0' );  %% cycle for detrending
name{16} = char( num2str(tlen));    %time number (size of cov matrix)
name{17} = char( '99.99' );   %want explain %
cseofmodenum=10;
name{18} = char( num2str(cseofmodenum) );    %want max.mode number
name{19} = char( '1.');   %% eof scale factor (Do not change)
name{20} = char( '2');  %% 1: rc ts    2 : cov
name{21} = char( 'ENDc' );
name{22} = char( [workdir, '/script/cseof < ', workdir, '/script/cseof.com'] );
name{23} = char( strcat('mv -f inform.d  ../cseofs/cinf_', 'NorESM1-M_', varname, '.d'));
name{24} = char(  strcat('mv -f Bloch.d  ../cseofs/blo_', 'NorESM1-M_', varname,'.d') );
name{25} = char(  strcat('mv -f pcts.d  ../cseofs/cpct_', 'NorESM1-M_', varname,'.d') );
name{26} = char(  strcat('mv -f eigen.d  ../cseofs/ceig_', 'NorESM1-M_', varname,'.d') );
name{27} = char( 'rm -f cseof cseof.com' );
name{28} = char( '' );
% end

fid = fopen([workdir, '/script/cseof.c'], 'w+')
for nline=1:size(name,1)
    fprintf(fid, '%s\n', name{nline});
end
fclose(fid);
                      

save('STEP2_2_output.mat')


