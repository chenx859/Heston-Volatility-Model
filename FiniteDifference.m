% Solutions to Heston PDE using explicit scheme.
function [Un,U]=hestonexeuro(kappa,theta,sigma,rho,V,J,r,T,dt,S,I,K)
nt=ceil(T/dt);
% dt = T/nt;
% the lower bounds of s and v are both 030 
sbound = S*2;
vbound = V;
ds = S/I; % step length of s
dv = V/J; % step length of v
svalue = 0:ds:sbound;
vvalue = 0:dv:vbound;
ns = size(svalue,2)-1;
nv = size(vvalue,2)-1;
iline = 0:1:ns;
jline = 0:1:nv;
temp = max(0, svalue-K); % as usual, the payoff of European call option.
uinitial = temp'*ones(1,nv+1); % to construct the initial condition of u.
U = uinitial; % u in dimension of (ns+1)*(nv+1).
for n=1:nt
 
 % the interior elements. Cross term part.
 Upp = U(3:ns+1,3:nv+1);
 Umm = U(1:ns-1,1:nv-1);
 Ump = U(1:ns-1,3:nv+1);
 Upm = U(3:ns+1,1:nv-1);
 Ut = Upp+Umm-Ump-Upm;
 G0 = rho*sigma*iline'*jline/4;
 G0t = G0(2:ns,2:nv);
 U1t = G0t.*Ut*dt;
 U1 = [zeros(1,nv+1);zeros(ns-1,1) U1t zeros(ns-1,1);zeros(1,nv+1)];
 
 % the elements in S direction. 
 B1 = iline(2:ns).*iline(2:ns)/2; B1 = [B1';0;0];
 B2 = -iline(2:ns).*iline(2:ns); B2 = [0;B2';0];
 B3 = iline(2:ns).*iline(2:ns)/2; B3 = [0;0;B3'];
 B = spdiags([B1 B2 B3],[-1 0 1],ns+1,ns+1);
 C1 = -r*iline(2:ns)/2; C1 = [C1';0;0];
 C2 = -r*ones(ns+1,1)/2;
 C3 = r*iline(2:ns)/2; C3 = [0;0;C3'];
 C = spdiags([C1 C2 C3],[-1 0 1],ns+1,ns+1);
 for j=2:nv
 A = U(:,j);
 G1 = vvalue(j)*B+C;
 U2temp = dt*G1*A; % the interior elements along the j-th sub-vector.
 U2temp(ns+1) = 0; % where s=S;
 U2temp(1) = 0; % where s=0;
 if j==2
 U2 = U2temp; 
 else U2 = [U2 U2temp];
 end
 end
 
 U2 = [zeros(ns+1,1) U2 zeros(ns+1,1)];
 
% elements in V direction
 D1 = sigma*sigma*vvalue(2:nv)/(2*dv*dv)-kappa*(theta-vvalue(2:nv))/(2*dv); 
 D1 = [D1';0;0];
 D2 = -sigma*sigma*vvalue(2:nv)/(dv*dv)-r/2;
 D2 = [0;D2';0];
 D3 = sigma*sigma*vvalue(2:nv)/(2*dv*dv)+kappa*(theta-vvalue(2:nv))/(2*dv); 
 D3 = [0;0;D3'];
 D = spdiags([D1 D2 D3],[-1 0 1],nv+1,nv+1);
 for i = 2:ns
 A = U(i,:)';
 U3temp = dt*D*A;
 U3temp(1) = 0;
 U3temp(nv+1) = 0;
 if i==2
 U3 = U3temp';
 else U3 = [U3; U3temp'];
 end
 end
 
 U3 = [zeros(1,nv+1);U3;zeros(1,nv+1)];
 
 Uvtank = U(:,1); % the first column of the previous U
 U = U + U1 + U2 + U3 ;
 U(1,:) = zeros(1,nv+1);
 U(:,nv+1) = svalue';
 % elements where v=0
 for i=2:ns
 U(i,1) = 2*U(i,2)-U(i,3);
% U(i,1)=(Uvtank(i)-r*i*dt*U(i-1,1)+kappa*theta*dt*U(i,2)/dv)/ ...
% (1-r*i*dt+kappa*theta*dt/dv+r*dt);
 end
 %U(ns+1,2:nv) = (2*ds+4*U(ns,2:nv)-U(ns-1,2:nv))/3;
 U(ns+1,2:nv) = ds+U(ns,2:nv);
 
end

 Un = U(1:I+1,1:J+1);

 U = Un(100/ds+1,ceil(0.12/dv))*(ceil(0.12/dv)*dv-0.12)/dv+Un(100/ds+1,ceil(0.12/dv)+1)*(0.12-(ceil(0.12/dv)-1)*dv)/dv;


end