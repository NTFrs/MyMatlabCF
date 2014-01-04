function [Price,CI]=AsianBS_MC(S0,K,r,sigma,T,Nsim,FlagCP,FlagFloat,varargin)
%ASIANBS_MC Computes price of an asian option in BS framework with MC
%
%Computes price and confidence interval for an asiatic call or put using a
%montecarlo method. User can specify if floating or fixed strike which type
%of mean to use and what Montecarlo Method.
%
%OUTPUTS:
%Price Option price
%CI confidence interval for option price at level alpha (default 0.05)
%
%INPUTS:
%S0 spot price
%K strike (dummy if floating strike)
%r constant risk free interest rate
%sigma volatility
%T maturity
%Nsim Number of simulations
%FlagCP 1 for call, -1 for put
%FlagFloat 1 for fixed strike -1 for floating strike
%
%Optional INPUTS:
%..,'ConfidenceLevel',value,.... to specify confidence level
%...,'Mean',string,... to specify if 'Ari' or 'Geo' mean
%OUTPUT: Price and confidence interval at 95%
%
%EXAMPLE:
%   AsianBS_MC(100,105,0.03,0.4,1,1e6,1,-1,'Mode','Anti','Mean','Geo')
%computes the price of a floating strike call with geometrical mean



indexM=find(strcmp('Mode',varargin));
indexA=find(strcmp('ConfidenceLevel',varargin));
indexT=find(strcmp('Mean',varargin));

%%%%%%%%%%%%%%%%Mode determination%%%%%%%%%%%%%%%%%%%%%%

Mode=0;
if(isempty(indexM)==0)
    if(strcmp(varargin{indexM+1},'Anti'));
        Mode=1;
    elseif(strcmp(varargin{indexM+1},'Control'));
        Mode=2;
    else
        disp('Mode not recognized, using default');
    end
end

%%%%%%%%%%%%%%Confidence Level Determination%%%%%%%%%%%%%%%%%%

alpha=0.05;
if(isempty(indexA)==0)
    alpha=varargin{indexA+1};
end

%%%%%%%%%%%%%%Mean type determination%%%%%%%%%%%%%%%%%
Geo=0;
ComputeMean=@(X) mean(X(:,2:end),2);
if(isempty(indexT)==0)
    if strcmp(varargin{indexT+1},'Geo')
        ComputeMean=@(X) geomean(X(:,2:end),2);
        Geo=1;
    else
        disp('Type of mean not reognized, using default (Ari)');
    end
end

%%%%%%%%%%%%%%%%%%Different Modes%%%%%%%%%%%%%%%%%
Nsteps=252*T;

if Mode==0
    
    S=underlyngBS(S0,r,sigma,T,Nsim,Nsteps);
    
    if(FlagFloat==1)
        
        DiscPayoff=max(FlagCP*(ComputeMean(S)-K),0)*exp(-r*T);
        
    elseif(FlagFloat==-1)
        
        DiscPayoff=max(FlagCP*(S(:,end)-ComputeMean(S)),0)*exp(-r*T);
        
    else
        disp('FlagFloat not recognized');
        return
    end
    
    
elseif Mode==1 %%%%%%%%%%%%Anti%%%%%%%%%%%%%%%
    
    [S,S_AV]=underlyngBS(S0,r,sigma,T,Nsim,Nsteps);
    
    if(FlagFloat==1)
        
        DiscPayoff_S=max(FlagCP*(ComputeMean(S)-K),0)*exp(-r*T);
        DiscPayoff_AV=max(FlagCP*(ComputeMean(S_AV)-K),0)*exp(-r*T);
        
    elseif(FlagFloat==-1)
        
        DiscPayoff_S=max(FlagCP*(S(:,end)-ComputeMean(S)),0)*exp(-r*T);
        DiscPayoff_AV=max(FlagCP*(S_AV(:,end)-ComputeMean(S_AV)),0)*exp(-r*T);
        
    else
        disp('FlagFloat not recognized');
        return
    end
    
    DiscPayoff=0.5*DiscPayoff_S+0.5*DiscPayoff_AV;
    
else %%%%%%%%%%control%%%%%%%%%%%%%
    
    smallN=Nsims/100;
    if FlagFloat==1
        
        %mean(S)(T) is control variable
        if strcmp(varargin{indexT+1},'Geo')
            Ef=1/T*exp(log(S0)*(r-0.5*sigma^2)*T^2/2+sigma^2*T^3/3);%%<---see note
        else
            Ef=S0*mean(exp(r*T/Nsteps*(1:Nsteps)));
        end
        
        S=underlyngBS(S0,r,sigma,T,smallN,Nsteps);
        f=ComputeMean(S)-K;
        g=max(f,0)*exp(-r*T);
        VC=cov(f,g);
        par=-VC(1,2)/VC(1,1);
        
        S=underlyngBS(S0,r,sigma,T,Nsim,Nsteps);
        f=ComputeMean(S)-K;
        g=max(f,0)*exp(-r*T);
        
        DiscPayoff=g-par*(f-EF);
        
    elseif FlagFloat==-1
        
        %S(T)-mean(S)(T) is control variable
        if Geo
            Ef=S0*exp(r*T)-1/T*...
                exp(log(S0)*(r-0.5*sigma^2)*T^2/2+sigma^2*T^3/3); %%%<---see note
        else
            Ef=S0*exp(r*T)-S0*mean(exp(r*T/Nsteps*(1:Nsteps)));
        end
        
        S=underlyngBS(S0,r,sigma,T,smallN,Nsteps);
        f=S(:,end)-ComputeMean(S);
        g=max(f,0)*exp(-r*T);
        VC=cov(f,g);
        par=-VC(1,2)/VC(1,1);
        
        S=underlyngBS(S0,r,sigma,T,Nsim,Nsteps);
        f=S(:,end)-ComputeMean(S);
        g=max(f,0)*exp(-r*T);
        
        DiscPayoff=g-par*(f-EF);
        
        
    else
        disp('FlagFloat not recognized');
        return
    end
    
end


[Price,Dummy,IC]=normfit(DiscPayoff,alpha);


end



%Contrl vaiate choice for geometrical mean: since in a time continuous
%setting the geometrical mean is given by GM(T)=1/T*exp(I(T)), with
%I(T)=int_0^T ln(S(t))dt, we have that I(T) is a normal variable with mean
%ln(S0)*T+(r-0.5*sigma^2)*T^2/2 and variance sigma^2*T^3/3 we obtain that the
%expected value of GM(T) is 1/T*exp(ln(S0)*T+(r-0.5*sigma^2)*T^2/2+sigma^2*T^3/3)