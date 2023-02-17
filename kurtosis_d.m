function S = kurtosis_d(x,xdata)

q = xdata(:,1);
Delta = xdata(1,2);
delta = xdata(2,2);

D = x(1)* 1000;
K = x(2);

S = exp(- q.^2 *(Delta - delta/3) * D + 1/6 * q.^4 * (Delta - delta/3)^2 * D^2 * K); 



