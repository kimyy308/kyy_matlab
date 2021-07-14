function line=ext_line_auto(line, upperline, gridd)
    tempdata=line;
    ismask=isnan(tempdata);
    if (sum(ismask-1)==0)
        line=upperline;
    elseif (sum(ismask)~=0)
        save('abc.mat')
        tempdata(ismask) = ...
            interp1(gridd(~ismask),tempdata(~ismask), ...
            gridd(ismask), 'nearest', 'extrap');
        line=tempdata;
    else 
        line=line;    
    end
return