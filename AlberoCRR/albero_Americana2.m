function [Prezzo]=albero_Americana2(S0,K,r,T,sigma,N,flag)

tic

dt=T/N;
u=exp(sigma*sqrt(dt)); d=1/u;
q=(exp(r*dt)-d)/(u-d);


Prezzo=max(flag*(S0*u.^[N:-2:-N]-K),0);

for j=2:N+1
    Prezzo=max(exp(-r*dt)*(q*Prezzo(1:end-1)+...
        (1-q)*(Prezzo(2:end))),flag*(S0*u.^[N+1-j:-2:-N-1+j]-K) );
end

toc
end
    