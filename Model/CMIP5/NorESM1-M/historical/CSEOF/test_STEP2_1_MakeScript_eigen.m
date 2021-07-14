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

name= cell(28);
%----------STEP2. make script file-----------------------------------%
name{1} = '#!/bin/csh';     %use csh shell
% for i = 1:length(nn);
    name{2} = char('');
    name{3} = char(['gfortran -fno-backtrace -o eigen ', workdir, '/programs/eigen/eigen.f']); %compile the eigen file
    name{4} = char('');
    name{5} = char('cat >! eof.com <<ENDc'  );
%     name1 = char(name1, strcat('../data/',nn(i,:),'.data')  );         %load variable file
    name{6} = char('./../data/test_NorESM1-M_tas1.data' );  %load variable file
    name{7} = char(['((', num2str(size(data,2)), 'e16.7))'] );  %ascii file read, (time_num)e16.7
    name{8} = char([num2str(size(data,1)), ' 1']); %space number
    name{9} = char(num2str(size(data,2)));  %time number
    name{10} = char('2            '    );  %space*time matrix or time*space matrix (1 or 2)
    name{11} = char('0            ' );                            
    name{12} = char('1          '  );
    name{13} = char('0              ' ); %% area adjustment (0: No)
    name{14} = char('15. 37.' );  %% starting latitude and increment for area adjustment
    name{15} = char('99.');  %want to modes to explain 95% (percent variance)
    name{16} = char('1.' );  %% EOF scaling factor
    maxmodenum = '5';
    name{17} = char(maxmodenum);  %% number of EOFs to be printed                              
    name{18} = char('0');  %% PC normalization
%     name{} = char(name{}, strcat('../cseofs/eof_',nn(i,:),'.dat')   );  %save eof LV file
    name{19} = char(strcat('../cseofs/eof_', 'NorESM1-M_', varname, '.dat')   );  %save eof LV file    
    name{20} = char('DIR');                                       %DIR : binary
%     name{} = char(name{}, strcat('../cseofs/pct_',nn(i,:),'.dat')   );  %save eof pct file
    name{21} = char(strcat('../cseofs/pct_', 'NorESM1-M_', varname, '.dat')   );  %save eof pct file
    name{22} = char('(5e13.5)' );     %% ASCII format "in parenthesis"                               %save option: 6e13.5 acsii
    name{23} = strtrim(char('ENDc'));
    name{24} = char(['./eigen < ./eof.com'] );
%     name{} = char(name{}, [workdir, '/script/eigen < ', workdir, '/script/eof.com'] );
    name{25} = char(strcat('mv -f inform.d',' ../cseofs/inf_', 'NorESM1-M_', varname,'.d')); %save information file
    name{26} = char(strcat('mv -f avg.d',' ../cseofs/avg_', 'NorESM1-M_', varname,'.d') );  %save average file
    name{27} = char('rm -f eigen eof.com' );
    name{28} = char('' );
%     name = char(name,strtrim(name1));
% end

fid = fopen([workdir, '/script/eigen.c'], 'w+')
for nline=1:length(name)
    fprintf(fid, '%s\n', name{nline});
end
fclose(fid);