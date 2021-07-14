%----This is the code to make eigen.c script file for same grid but
%differnet variables

%---------STEP1. set the variables-----------------------------------%
name_layer = char('uv','air','hgt','omg'); %variable name to n layer
nn = char('slp','sst2','ssh2'); %variable name to one layer

for i = 1:4 %number of name layer
    n1 = n(i,:);
    for j = 1 : 10 %layer number
        name = strcat(n1,'_',num2str(j));
        nn = char(nn,name);
    end
end

%----------STEP2. make script file-----------------------------------%
name = '#!/bin/csh';     %use csh shell
for i = 1:length(nn);
    name1 = char('');
    name1 = char(name1,'gfortran -o combin ../programs/util/combin_new.f');
    name1 = char(name1, '' );
    name1 = char(name1, 'cat >! combin.com << ENDc '  );
    name1 = char(name1, strcat('../cseofs/blo_',nn(i,:),'.d')  );     %predictor blo file
    name1 = char(name1, 'SEQ' );
     name1 = char(name1, strcat('../cseofs/cinf_',nn(i,:),'.d') );    %predictor inf file
    name1 = char(name1, '(22x,e15.7,////,10(40x,e16.7,/))'  );        %format, 10 should be changed to regress mode number
    name1 = char(name1, '10');                                        %10 should be changed to regress mode number
    name1 = char(name1,  '12            '    );                       %LV number
    name1 = char(name1,  '20 1          '  );                         %max.EOF mode number of predictor :20
    name1 = char(name1, '3              ' );                          %regressed mode number:3
    name1 = char(name1, '0' );                  
    name1 = char(name1, strcat('../regress/blo_reg_',nn(i,:),'.d'));  %predictor,regressed blo file
    name1 = char(name1, 'SEQ' );
    name1 = char(name1, '2' );
    name1 = char(name1, strcat('../regress/reg_',nn(i,:),'.d'));      %target-predictor regress file name
    name1 = char(name1, '(///,10(7X,E13.5,/),/)');                    %10 is regress mode number
    name1 = char(name1, '../cseofs/cinf_sst2.d');                     %target cinf file name
    name1 = char(name1, '(22x,e15.7,////,10(41x,e15.7,/))');          %10 is target mode number
    name1 = char(name1, 'ENDc' );
    name1 = char(name1, './combin < combin.com  ' );
    name1 = char(name1, '.rm combin combin.com  ' );
    name = char(name,name1);
end


                      




