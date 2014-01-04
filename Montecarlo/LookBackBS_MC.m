function [Price,CI]=LookBackBS_MC(S0,K,r,sigma,T,Nsim,FlagCP,FlagFixed,varargin)
%LOOKBACKBS_MC Computes the price and c.i. for a LookBack using MC
%
%Computes Lookback Options price and confidence interval, using
%different Montecarlo methods
%
%INPUTS:
%
%S0 spot price
%K Strike
%r risk free constant interest rate
%sigma volatility
%T maturity
%Nsim Number of simulations
%FLagCP 1 for Call and -1 for Put
%FlagFixed 1 for LookBack options and -1 for forward LookBack options
%
%Other arguments:
%..,'ConfidenceLevel',value,.... to specify confidence level
%...,'Mode',string,...) to specify Method
%   'Anti' for anthitethic variables
%
%
%OUTPUT: Price and confidence interval at 95%
%
%EXAMPLE:
%   LookBackBS_MC(100,105,0.03,0.4,1,1e6,-1,1,'Mode','Anti',)
%Computes the price of an Lookback Put using Antithetic variables method


indexM=find(strcmp('Mode',varargin));
indexA=find(strcmp('ConfidenceLevel',varargin));

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


Nsteps=252*T;

%%%%%%%%%%Different modes%%%%%%%%%%%%%%

if Mode==0 %%%%%%%%Normal Montecarlo%%%%%%%%%%%%%%%
    
    %Asset simulation
    S=underlyngBS(S0,r,sigma,T,Nsim,Nsteps);
    
    %LookBack or Forward LookBack
    if(FlagFixed==1) %Out opton
        
        
        %Put or call
        if FlagCP==1
            DiscPayoff=max(S(:,end)-min(S,[],2),0)*exp(-r*T);
        else
            DiscPayoff=max(max(S,[],2)-S(:,end),0)*exp(-r*T);
        end
        
    elseif(FlagFixed==-1)
        
        %Put or call
        if FlagCP==1
            DiscPayoff=max(max(S,[],2)-K,0)*exp(-r*T);
        else
            DiscPayoff=max(K-min(S,[],2),0)*exp(-r*T);
        end
        
        
    else
        disp('Wrong Flag for type of option');
    end
    
    
    
elseif Mode==1 %%%%%%%%%%%Antithetic Variables%%%%%%%%%%%%%
    
    %Asset simulatio
    [S,S_AV]=underlyngBS(S0,r,sigma,T,Nsim,Nsteps);
    
    %LookBack or Forward LookBack
    if(FlagFixed==1)
        
        
        %Put or call
        if FlagCP==1
            DiscPayoff=max(S(:,end)-min(S,[],2),0)*exp(-r*T);
            DiscPayoff_AV=max(S_AV(:,end)-min(S_AV,[],2),0)*exp(-r*T);
        else
            DiscPayoff=max(max(S,[],2)-S(:,end),0)*exp(-r*T);
            DiscPayoff_AV=max(max(S_AV,[],2)-S_AV(:,end),0)*exp(-r*T);
        end
        
    elseif(FlagType==-1)
        
        %Put or call
        if FlagCP==1
            DiscPayoff=max(max(S,[],2)-K,0)*exp(-r*T);
            DiscPayoff_AV=max(max(S_AV,[],2)-K,0)*exp(-r*T);
        else
            DiscPayoff=max(K-min(S,[],2),0)*exp(-r*T);
            DiscPayoff_AV=max(K-min(S_AV,[],2),0)*exp(-r*T);
        end
        
    else
        disp('Wrong Flag for type of option');
        return
    end
    
    DiscPayoff=0.5*DiscPayoff_S+0.5*DiscPayoff_AV;
    
    
    
% else %%%%%%%%%%%%%%%%%Control Variables%%%%%%%%%%%%%%%%
%     
%     %we use S(T)-min(S)(T)(or max) as control variable
%     
%     %small Asset simulation
%     
%     smallN=Nsim/100;
%     S=underlyngBS(S0,r,sigma,T,smallN,Nsteps);
%     
%     
%     %LookBack or Forward Lookback
%     if(FlagType==1) %lookback
%         
% %         if FlagCP==1
% %             
% %         EF=S0*exp(r*T)
% %         
% %         else
% %             
% %         end
%         
%     elseif(FlagType==-1) %In option
% 
%         
% %         if FlagCP==1
% %             
% %         EF=S0*exp(r*T)
% %         
% %         else
% %             
% %         end
%         
%         
%     else
%         disp('Wrong Flag for type of option');
%         return;
%     end
%     
%     %DiscPayoff=g+par*(f-Ef);
    
    
end

[Price,Dummy,CI]=normfit(DiscPayoff,alpha);

end
