% x=-50:50;
x=0:100;
N=length(x);
se = std(x)/sqrt(N);
conf=0.95;
alpha = 1 - conf;
pLo = alpha/2;
pUp = 1-alpha/2;
crit=tinv([pLo pUp], N-1);
xbar= mean(x);
ci=xbar + crit*se

pd = fitdist(x', 'Normal');

% paramci(pd)

% [h, p, ci, stats]=ttest(x, mean(x), 'Alpha', 0.01)