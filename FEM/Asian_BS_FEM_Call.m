function [price]=Asian_BS_FEM_Call(S0,r,sigma,T,M,N)
%ASIAN_BS_FEM_CALL does nothing for now, doesn't work

Smin=S0*exp( (r-sigma^2/2)*T-sigma*sqrt(T)*6);
Smax=S0*exp( (r-sigma^2/2)*T+sigma*sqrt(T)*6);
Rmax=T*Smax/Smin;
R=linspace(0,Rmax,N+1)'; dR=R(2)-R(1);
t=linspace(0,T,M+1)'; dt=t(2)-t(1);

%assembling B
outer=dR/6*ones(N,1);
diag=[dR/3 ; 2*dR/3*ones(N-1,1)];
B=spdiags([outer diag outer],-1:1,N,N);
clear outer diag

%assembling A
aux=(3*dR^2); lost =aux;
upper = R(1:end-1).^3/aux-R(2:end).^3/aux;
diag=R(3:end).^3/aux-R(1:end-2).^3/aux; diag=[R(2).^3/aux; diag ];
under=R(1:end-2).^3/aux-R(2:end-1).^3/aux; under=[lost; under];
tridiag=[under, diag, upper];
A=spdiags(tridiag,1:-1:-1,N,N)';
clear aux lost upper diag under tridiag

%assembling C
C=spdiags([ 0.5*ones(N,1) zeros(N,1) -0.5*ones(N,1)],-1:1,N,N); C(1,1)=-0.5;

%assembling D
diag= R(1:end-2)/2-2*dR/3-R(3:end)/2; diag=[-dR/6 ; diag];
upper=-R(1:end-1)-dR/3;
under=R(1:end-1)/2-dR/3;
D=spdiags([under, diag, upper],1:-1:-1,N,N)';
clear diag upper under

Mat=-(sigma^2+r)*D-0.5*sigma^2*A+C;
H=B'/dt-Mat';
V=max(1-R(1:end-1)/T,0);



for k=M:-1:1
%   disp('INHERE')
    V=H\(B'*V/dt);
end

price=S0*V(1);


end