%----This is the code to make eigen.c script file for same grid but
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

load('STEP2_2_output.mat')

%----------STEP2. make script file-----------------------------------%
%------case1. one layer-----------------%
name = cell(23);
name{1} = '#!/bin/csh';     %use csh shell

% for i = 2:2;
name{2} = char('');
name{3} = char('gfortran -fno-backtrace -o recast ../programs/util/recastx.f');
name{4} = char( '' );
name{5} = char( 'cat >! recast.com <<ENDc'  );
name{6} = char( strcat('../cseofs/eof_', 'NorESM1-M_', varname, '.dat  ')  ); %predictor eof file
name{7} = char( 'DIR' );
name{8} = char( 'nofile'  );
name{9} = char( num2str(cseofmodenum));                     %max. eof mode number
name{10} = char(  num2str(slen) );   %space number
name{11} = char(  '1.      ' );                              
name{12} = char(  [num2str(slen), ' 1']  );     %space structure
% name{13} = char( strcat('../regress/blo_reg_', 'NorESM1-M_', varname,'.d') ); %regressed blo file
name{13} = char( strcat('../cseofs/blo_', 'NorESM1-M_', varname,'.d') ); %blo file
name{14} = char( 'SEQ' );
name{15} = char( '0');
sizeLV = LVnumber * cseofmodenum;
name{16} = char( num2str(sizeLV) );  %LV number * cseof mode number (365*10)
name{17} = char( '1' );
name{18} = char( strcat('../cseofs/LV_', 'NorESM1-M_', varname,'.dat'));     %LV file name
name{19} = char( 'DIR');
name{20} = char( 'ENDc' );
name{21} = char( [workdir, '/script/recast < ', workdir, '/script/recast.com'] );
name{22} = char( 'rm -f recast recast.com' );
name{23} = char( '' );
% end

fid = fopen([workdir, '/script/recast.c'], 'w+')
for nline=1:size(name,1)
    fprintf(fid, '%s\n', name{nline});
end
fclose(fid);
               

% %----case 2. multi layer---------------%                      
% for i = 4:13;
%     name{} = char('');
%     name{} = char('gfortran -o recast ../programs/util/recastx.f');
%     name{} = char( '' );
%     name{} = char( 'cat >! recast.com <<ENDc'  );
%     name{} = char( strcat('../cseofs/eof_',nn(i,:),'.dat')  );
%     name{} = char( 'DIR' );
%     name{} = char( 'nofile'  );
%     name{} = char( '20');
%     name{} = char(  '3250 '    );
%     name{} = char(  '1.      ' );
%     name{} = char(  '65 50   '  );
%     name{} = char( strcat('../regress/blo_reg_',nn(i,:),'.d') );
%     name{} = char( 'SEQ' );
%     name{} = char( '0');
%     name{} = char( '36' );
%     name{} = char( '1' );
%     name{} = char( strcat('../cseofs/LV_',nn(i,:),'.dat'));
%     name{} = char( 'DIR');
%     name{} = char( 'ENDc' );
%     name{} = char( './recast < recast.com' );
%     name{} = char( 'rm recast recast.com' );
%     name{} = char( '' );
%     name = char(name,name{});
% end
           

save('STEP2_5_output.mat')
