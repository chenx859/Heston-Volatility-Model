function ret = HestonP(kappa,theta,sigma,v0,r,T,s0,K,N,u)
trunc=ceil(N^(1/3));
delta=trunc/N;
a=0;
P=0;
for i=1:N;
    P=P+HestonPIntegrand(a,kappa,theta,sigma,v0,r,T,s0,K,u)*delta;
    P(isnan(P))= 0;
    a=a+delta;
end
ret = 0.5 + 1/pi*P;
end

function ret = HestonPIntegrand(phi,kappa,theta,sigma,v0,r,T,s0,K,u)
ret = real(exp(-i*phi*log(K))*Hestf(phi,kappa,theta,sigma,v0,r,T,s0,u)/(i*phi));
end

function f = Hestf(phi,kappa,theta,sigma,v0,r,T,s0,u)
x = log(s0);
d = sqrt(kappa^2-sigma^2*(2*u*phi*i-phi^2));
g = (kappa + d)/(kappa - d);
C = r*phi*i*T + kappa*theta/sigma^2*((kappa + d)*T -2*log((1-g*exp(d*T))/(1-g)));
D = (kappa + d)/sigma^2*((1-exp(d*T))/(1-g*exp(d*T)));
f = exp(C + D*v0 + i*phi*x);
end

