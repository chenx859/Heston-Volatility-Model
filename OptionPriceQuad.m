s0=100;
u1=0.5;
u2=-0.5;
kappa=2;
r=0;
sigma=0.8;
theta=0.15;
v0=0.12;
K=[70 90 100 110 130];
T=[1/4 1/2 1 5];
N=[10 50 100 1000];
C=zeros(4,20);


for n=1:4;
    for k=1:5;
        for t=1:4;
            c=s0*HestonP(kappa,theta,sigma,v0,r,T(t),s0,K(k),N(n),u1)-...
                K(k)*exp(-r*T(t))*HestonP(kappa,theta,sigma,v0,r,T(t),s0,K(k),N(n),u2);
            C(n,t+4*(k-1))=c;
        end
    end
end


CN=zeros(1,20);
for k=1:5;
    for t=1:4;
        c=s0*HestonP(kappa,theta,sigma,v0,r,T(t),s0,K(k),100000,u1)-...
            K(k)*exp(-r*T(t))*HestonP(kappa,theta,sigma,v0,r,T(t),s0,K(k),100000,u2);
        CN(1,t+4*(k-1))=c;
    end
end

CNN=[CN;CN;CN;CN];
MSE=mean((C-CNN).^2,2);
plot(N,MSE),xlabel('N'),ylabel('MSE');
hold on
plot(N,N.^(-1/3));
hold off