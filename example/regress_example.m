


%% How to regress Y on X using matlab?

Y=[81  55 80  24 78 52  88 45 50 69 66 45 24 43 38 72  41  48 52  52  66 89];

X=[124 49 181 4  22 152 75 54 43 41 17 22 16 10 63 170 125 15 222 171 97 254];

b=regress(Y.', X.');

close all;
plot(X);
hold on
plot(Y);
plot(X*b);
% plot(Y*b);
legend


% example 2
x = 1:10;
y = 1:10;
X = [x; ones(1,length(x))];
b = regress(y.',X.');