function [Prezzo]=albero_Europea(S0,K,r,T,sigma,N,flag)

tic

dt=T/N;
u=exp(sigma*sqrt(dt)); d=1/u;
q=(exp(r*dt)-d)/(u-d);

Albero=zeros(N+1,N+1);
Albero(:,1)=max(flag*(S0*u.^[N:-2:-N]-K),0);

for j=2:N+1
    Albero(1:end-j+1,j)=exp(-r*dt)*(q*Albero(1:end-j+1,j-1)+...
        (1-q)*(Albero(2:end-j+2,j-1)));
end

Prezzo=Albero(1,end);

toc
end
    