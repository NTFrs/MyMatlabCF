function Prezzo=Exact_Euro_BS(St,K,r,t,T,sigma,flag)

%EXACT_EURO_BS Compute exact price of a Call (1) or a Put (-1) in B&S

%Computes the exact value in Black And Scholes model for a Call if flag =1,
%for a put if flag =-1
%
%St price at time t
%K Strike
%r risk free constant interest rate
%t evaluation time
%T maturity
%sigma volatility
%flag 1 for Call -1 for Put



d1=(log(St/K)+(r+0.5*sigma^2)*(T-t))/(sigma*sqrt(T-t));
d2=d1-sigma*sqrt(T-t);

if flag==1
    Prezzo=St*normcdf(d1)-K*exp(-r*(T-t))*normcdf(d2);
elseif flag==-1
    Prezzo=normcdf(-d2)*K*exp(-r*(T-t))-normcdf(-d1)*St;
else
    disp('Flag not recognized, use 1 for Call and -1 for Put');
end

end