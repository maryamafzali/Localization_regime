function S = localization_code(x,q)

r = x(1);
u = x(2);
v = x(3);
w = x(4);
p = x(5);

S = exp( -(q*r).^(2/3) .* ( (q*u).^(4/3* exp(-q.^6 * p^6)) - (q*v).^(10/3 * exp(-q.^6 *p^6)) + (q*w).^(16/3 * exp(-q.^6 *p^6)) ) );

