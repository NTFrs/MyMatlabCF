function [Price]=Asian_BS_FD_Call(S0,r,sigma,T,M,N,Mode)



if (nargin<7)
    Mode=0;
end



Smin=S0*exp( (r-sigma^2/2)*T-sigma*sqrt(T)*6);
Smax=S0*exp( (r-sigma^2/2)*T+sigma*sqrt(T)*6);
Rmax=T*Smax/Smin;
R=linspace(0,Rmax,N+1)'; dR=R(2)-R(1);

t=linspace(0,T,M+1)'; dt=t(2)-t(1);

nodes=R(2:end-1);

under=0.5*sigma^2*nodes.^2/dR^2 -(1 -r*nodes)/(2*dR);
BCa=under(1);
mid=-1/dt-sigma^2*nodes.^2/dR^2;
upper=0.5*sigma^2*nodes.^2/dR^2 +(1 -r*nodes)/(2*dR);
BCend=upper(end);

tridiag=[ under, mid, upper ];
A=spdiags(tridiag,1:-1:-1,N-1,N-1)';

AA=sparse(N,N);
AA(2:end,2:end)=A;
AA(1,1)=-1/dt-1/dR; AA(1,2)=1/dR; AA(2,1)=BCa;

V=max(1-R(1:end-1)/T,0);

for j=1:M

    V=AA\(-V/dt);
    
end

V=[V;0];
Price=S0*V(1);
end