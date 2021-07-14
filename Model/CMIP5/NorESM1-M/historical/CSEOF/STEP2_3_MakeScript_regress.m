%----This is the code to make regress.c script file for same grid but
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
    %Target 1 mode to predictor 10 modes
    name1 = char('');
    name1 = char(name1,'gfortran -o regress ../programs/util/regress_new_mean.f');
    name1 = char(name1, '' );
    name1 = char(name1, 'cat >! regress.com << END'  );
    name1 = char(name1,'768 10    ');                     %Want regress PC time series (time,predicted mode)
    name1 = char(name1, '../cseofs/cpct_sst2.d' );        %target PC file name
    name1 = char(name1, '(6e13.5)' );                     %format
    name1 = char(name1, '0          '  );                 %skip amount        
    name1 = char(name1, '1 768     ');                    %data interval for target time series
    name1 = char(name1,  strcat('../cseofs/cpct_',nn(i,:),'.d' )    ); %predictor time series
    name1 = char(name1,  '(6e13.5)            ' );                     %format
    name1 = char(name1,  '0                  '  );                     %skip amout
    name1 = char(name1, '1 768                   ' );                  %data interval for predictor time series
    name1 = char(name1, '1 768    ' );                                 %regression interval
    name1 = char(name1, 'regress.d');
    name1 = char(name1, '1         ' );
    name1 = char(name1, '2        ' );                                %confidence interval(0:no,2:90%,3:95%,4:99%)
    name1 = char(name1, '2          ');
    name1 = char(name1, '../cseofs/cinf_sst2.d');                    %target cinf file
    name1 = char(name1, '(41x,e15.7)' );                             %format
    name1 = char(name1, strcat('../cseofs/cinf_',nn(i,:),'.d' ) );   %predictor cinf file
    name1 = char(name1, '(41x,e15.7) ' );                            %format
    name1 = char(name1, 'END' );
    name1 = char(name1, './regress < regress.com' );
    name1 = char(name1, 'mv regress.d regress.s');
    name1 = char(name1, '' ); 
    
    %Target 2 mode 
    name1 = char(name1,'gfortran -o regress ../programs/util/regress_new_mean.f');
    name1 = char(name1, '' );
    name1 = char(name1, 'cat >! regress.com << END'  );
    name1 = char(name1,'768 10   ');
    name1 = char(name1, '../cseofs/cpct_sst2.d' );
    name1 = char(name1, '(6e13.5)' );
    name1 = char(name1, '128          '  );       %skip amount sholud be changed : (regress interval)/6 
    name1 = char(name1, '1 768     ');
    name1 = char(name1,  strcat('../cseofs/cpct_',nn(i,:),'.d' )    );
    name1 = char(name1,  '(6e13.5)            ' );
    name1 = char(name1,  '0                 '  );
    name1 = char(name1, '1 768                ' );
    name1 = char(name1, '1 768     ' );
    name1 = char(name1, 'regress.d');
    name1 = char(name1, '1         ' );
    name1 = char(name1, '2          ' );
    name1 = char(name1, '2          ');
    name1 = char(name1, '../cseofs/cinf_sst2.d');
    name1 = char(name1, '(/,41x,e15.7)' );       %read 2mode inf, so should input /
    name1 = char(name1, strcat('../cseofs/cinf_',nn(i,:),'.d' ) );
    name1 = char(name1, '(41x,e15.7) ' );
    name1 = char(name1, 'END' );
    name1 = char(name1, './regress < regress.com' );
    name1 = char(name1, 'cat regress.d >> regress.s');
    name1 = char(name1, '' );

    %Target 3mode
    name1 = char(name1,'gfortran -o regress ../programs/util/regress_new_mean.f');
    name1 = char(name1, '' );
    name1 = char(name1, 'cat >! regress.com << END'  );
    name1 = char(name1,'768 10    ');
    name1 = char(name1, '../cseofs/cpct_sst2.d' );
    name1 = char(name1, '(6e13.5)' );
    name1 = char(name1, '256    '       ); %skip amount sholud be changed : (regress interval)/6 *2
    name1 = char(name1, '1 768     ');
    name1 = char(name1,  strcat('../cseofs/cpct_',nn(i,:),'.d' )    );
    name1 = char(name1,  '(6e13.5)            ' );
    name1 = char(name1,  '0                '  );
    name1 = char(name1, '1 768                  ' );
    name1 = char(name1, '1 768     ' );
    name1 = char(name1, 'regress.d');
    name1 = char(name1, '1         ' );
    name1 = char(name1, '2         ' );
    name1 = char(name1, '2          ');
    name1 = char(name1, '../cseofs/cinf_sst2.d');
    name1 = char(name1, '(//,41x,e15.7)' ); %read 3mode inf, so should input //
    name1 = char(name1, strcat('../cseofs/cinf_',nn(i,:),'.d' ));
    name1 = char(name1, '(41x,e15.7) ' );
    name1 = char(name1, 'END' );
    name1 = char(name1, './regress < regress.com' );
    name1 = char(name1, 'cat regress.d >> regress.s');
    %save script
    name1 = char(name1, strcat('mv regress.s ../regress/reg_',nn(i,:),'.d') );
    name1 = char(name1, 'rm regress regress.com regress.d' );
    name1 = char(name1, 'rm Y_err.d Y_est.d' );
    
    name = char(name,name1);
end

















                      




