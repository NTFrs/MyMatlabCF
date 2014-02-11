function [Price,CI]=EuroBS_MC(S0,K,r,T,sigma,N,flag,alpha,Mode)
%EUROBS_MC Compute call or put price in BS using a Montecarlo Methods
%
%INPUTS:
%S0 spot price
%K Strike
%r risk free constant interest rate
%T maturity
%sigma volatility
%flag 1 for Call -1 for Put
%alpha Confidence level (CI given at 100*(1-alpha)%), default 0.05
%Mode Evaluation mode (Optional)
%   Mode=0 for classical MC (default)
%   Mode=1 for Antitetic variables
%   Mode=2 for Control Variables
%
%OUTPUTS:
%Price price at spot price
%CI confidence interval at 95%

if nargin<8
    alpha=0.05;
    Mode=0;
end
if nargin<9
    Mode=0;
end

if Mode==0


    disp('Normal');
    ST=S0*exp( (r-sigma^2/2)*T+sigma*sqrt(T)*randn(N,1) );
    
    DiscPayoff=max(flag*(ST-K),0)*exp(-r*T);
    
    [Price,Dummy,CI]=normfit(DiscPayoff,alpha);
    
elseif Mode==1
    disp('anti');
    temp=randn(N,1);
    
    ST=S0*exp( (r-sigma^2/2)*T+sigma*sqrt(T)*temp );
    ST_AV=S0*exp( (r-sigma^2/2)*T+sigma*sqrt(T)*(-temp) );
    
    DiscPayoff=0.5*max(flag*(ST-K),0)*exp(-r*T)+...
        0.5*max(flag*(ST_AV-K),0)*exp(-r*T);
       
    [Price,Dummy,CI]=normfit(DiscPayoff,alpha);
elseif Mode==2
    disp('Control');
    Ef=S0*exp(r*T); %media variabile di controllo
    smallN=N/100; %parametro simulazioni
    
    f=S0*exp( (r-sigma^2/2)*T+sigma*sqrt(T)*randn(smallN,1) ); %var control
    g=max(flag*(f-K),0)*exp(-r*T);%function to evaluate
    
    VC=cov(f,g);
    
    par=-VC(1,2)/VC(1,1);
    
    f=S0*exp( (r-sigma^2/2)*T+sigma*sqrt(T)*randn(N,1));
    g=max(flag*(f-K),0)*exp(-r*T);
    
    [Price,Dummy,CI]=normfit(g+par*(f-Ef),alpha);
    
else
    disp('Mode not reconigzed, consult help for more information');
end


end