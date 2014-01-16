function [Price]=FFT_Call(S0,Strike,r,T,M,charExp)

%risk neutral function definition
neutralExp=@(u) charExp(u)-charExp(-1i)*1i*u;
charfunc=@(u) exp(T*neutralExp(u));

N=2^M;
A=600; %think this

%fourier space discretization
eta=A/N; 
v=[0:eta:(N-1)*A/N]; v(1)=1e-22;

%integrand preparation
Z_k=exp(1i*v*r*T).*(charfunc(v-1i)-1)./(1i*v.*(1i*v+1));

%weights
w=ones(1,N); w(1)=0.5;w(end)=0.5;

%transform
trasf=w.*eta.*Z_k.*exp(1i*pi*(0:N-1));
z_k=real(fft(trasf))/pi;

%real space discretization
lambda=2*pi/(N*eta);
k=-lambda*N/2+lambda*(0:N-1);

C=S0*(z_k+max(1-exp(k-r*T),0));
K=S0*exp(k);

%Output formatting
index=find(K>0.1*S0 & K<3*S0);
C=C(index);K=K(index);

Price=interp1(K,C,Strike,'spline');





end