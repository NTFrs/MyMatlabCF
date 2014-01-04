function Price=EuroBS_FD(S0,K,r,d,sigma,T,M,N,FlagCP,theta)
%EUROBS_FD Computes price of an european option using DF
%
%Computes price of a vanilla option  in B&S framework using a FD method. By
%default it uses Crank-Nicolson algorithm for time discretization, thus
%being unconditionnaly stable.
%
%Output:
%Price at time 0 and S value S0
%
%Inputs:
%S0: Initial stock value
%K: Strike
%r: risk free constant interest rate
%d: dividend rate
%sigma: volatility
%T: maturity
%M: Number of time intervals
%N: Number of "space" intervals
%FlagCP: 1 for Call, -1 (or anything else) for Put
%theta: (OPTIONAL) theta coefficient:
%        theta=0 Implicit Euler
%        theta=1 Explicit Euler


if (nargin<10)
    theta=0.5;
end

ctheta=1-theta;
%Needed Parameters
q=2*r/sigma^2; qd=2*(r-d)/sigma^2;

%Creating domain
Smin=S0*exp( (r-sigma^2/2)*T-6*sigma*sqrt(T));
Smax=S0*exp( (r-sigma^2/2)*T+6*sigma*sqrt(T));
x=linspace(log(Smin/K),log(Smax/K),N+1)'; dx=x(2)-x(1);
tau=linspace(0,T*sigma^2/2,M+1); dtau=tau(2)-tau(1);

nodi=x(2:end-1);

%Stability Parameter
lambda=dtau/dx^2;
if theta > 0.5
if lambda>0.5
    disp('May have stability issues')
end
end

tridiag=[ones(N-1,1)*(-ctheta)*lambda,ones(N-1,1)*(1+2*ctheta*lambda),...
    ones(N-1,1)*(-ctheta)*lambda];
A=spdiags(tridiag,-1:1,N-1,N-1);
tridiag=[ones(N-1,1)*theta*lambda,ones(N-1,1)*(1-2*theta*lambda),...
    ones(N-1,1)*theta*lambda];
B=spdiags(tridiag,-1:1,N-1,N-1);
clear tridiag;

if FlagCP==1
    r1=@(j) 0;
    r2=@(j) exp(0.5*(qd+1)*x(end)+0.25*(qd+1)^2*tau(j));
    Y=max(exp(nodi/2*(qd+1))-exp(nodi/2*(qd-1)),0);
else
    r1=@(j) exp(0.5*(qd-1)*x(1)+0.25*(qd-1)^2*tau(j));
    r2=@(j) 0;
    Y=max(exp(nodi/2*(qd-1))-exp(nodi/2*(qd+1)),0);
end

BC=sparse(N-1,1);

for j=2:M+1
    BC(1)=ctheta*lambda*r1(j)+theta*lambda*r1(j-1);
    BC(end)=ctheta*lambda*r2(j)+theta*lambda*r2(j-1);
    Y=A\(B*Y+BC);
end
Y=[r1(M+1);Y;r2(M+1)];

V=K*exp(-0.5*(qd-1)*x-(0.25*(qd-1)^2+q)*tau(end)).*Y;
S=K*exp(x);
plot(S,V);
title('Prezzo dell''opzione al variare di S');
Price=interp1(S,V,S0,'spline');


end