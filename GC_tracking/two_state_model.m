function two_state_model(t, v);

%setting paramters
t = 0:0.001:1;
t0 = 0.5;
k = 0.5;
nu_0 = 0.1;

s = sinh(k*(t-t0));
l = nu_0/k;

%calc of p(t=0)
s0 = s(1);
x = 0:s0/1000:s0;
y = exp(2.*l.*x)*(1-x/sqrt(1+x.^2));
p0 = l/(exp(2*l*s0)-1)*trapz(x,y);

%calc of probability p(t)
p = zeros(length(t), 1);
p(1) = p0;
for i=2:length(t)
    s_t = s(i); 
    x = 0:s_t/(i*100):s_t;
    y = exp(-2.*l.*(s_t-x))*(1-x/sqrt(1+x.^2));
    p(i) = exp(-2*l*s_t)*p0 + l*trapz(x,y);
end

