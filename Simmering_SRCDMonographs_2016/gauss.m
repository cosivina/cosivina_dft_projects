function y = gauss(n,n0,G0,sigma,k)

x= 0:n;

y = G0 * exp(-0.5 * (x-n0).^2/sigma.^2) + k;

