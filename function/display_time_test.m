

lap_time_j=tic;
jlen=10;
for j=1:jlen
elapsed=toc(lap_time_j);
nchar(1)= fprintf(' %.0f sec. (%.0f sec. left)', elapsed, elapsed*(jlen-j+1)/j);
pause(2)
fprintf(repmat('\b',1,sum(nchar)))  % remove printed time
end


% fprintf(repmat('\b',1,1))