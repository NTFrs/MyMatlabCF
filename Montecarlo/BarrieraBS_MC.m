function [Price,CI]=BarrieraBS_MC(S0,K,r,sigma,T,Nsim,FlagCP,FlagType,varargin)
%BARRIERABS_MC Computes the price and c.i. for a Barrier using MC
%
%Computes Barrier Options price and confidence interval, using
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
%FlagType 1 for Out options and -1 for In options
%
%Other arguments:
%..,'ConfidenceLevel',value,.... to specify confidence level
%...,'Low',Number,... to specify lower Barrier
%...,'Upp',Number,... to specify Upper Barrier
%...,'Mode',string,...) to specify Method
%   'Anti' for anthitethic variables
%   'Control' for control variable
%
%OUTPUT: Price and confidence interval at 95%
%
%EXAMPLE:
%   BarrieraBS_MC(100,105,0.03,0.4,1,1e6,-1,1,'Mode','Anti','Low',91)
%Computes the price of an Down&Out Put with barrier 91 using Antithetic
% variables method


indexU=find(strcmp('Upp',varargin));
indexL=find(strcmp('Low',varargin));
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

%%%%%%%%%%%%%Barriers determinarion%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(isempty(indexU)==0)
    U=varargin{indexU+1};
end

if(isempty(indexL)==0)
    L=varargin{indexL+1};
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
    
    %Distinction between OUT or IN
    if(FlagType==1) %Out opton
        
        
        %Barrier use
        if(exist('L','var'))
            tempL=(min(S,[],2)>L);
        else
            tempL=ones(Nsim,1);
        end
        
        if(exist('U','var'))
            tempU=(max(S,[],2)<U);
        else
            tempU=ones(Nsim,1);
        end
        
        DiscPayoff=max(FlagCP*(S(:,end)-K),0)*exp(-r*T).*tempL.*tempU;
        
    elseif(FlagType==-1) %In option
        
        %Barrier use
        if(exist('L','var'))
            tempL=(min(S,[],2)<=L);
        else
            tempL=ones(Nsim,1);
        end
        
        if(exist('U','var'))
            tempU=(max(S,[],2)>=U);
        else
            tempU=ones(Nsim,1);
        end
        
        DiscPayoff=max(FlagCP*(S(:,end)-K),0)*exp(-r*T).*tempL.*tempU;
        
    else
        disp('Wrong Flag for type of option');
    end
    
    
    
elseif Mode==1 %%%%%%%%%%%Antithetic Variables%%%%%%%%%%%%%
    
    %Asset simulatio
    [S,S_AV]=underlyngBS(S0,r,sigma,T,Nsim,Nsteps);
    
    %Distinction between OUT or IN
    if(FlagType==1) %Out opton
        
        %Barrier use
        if(exist('L','var'))
            tempL=(min(S,[],2)>L);
            tempL_AV=(min(S_AV,[],2)>L);
        else
            tempL=ones(Nsim,1);
            tempL_AV=ones(Nsim,1);
        end
        
        if(exist('U','var'))
            tempU=(max(S,[],2)<U);
            tempU_AV=(max(S_AV,[],2)<U);
        else
            tempU=ones(Nsim,1);
            tempU_AV=ones(Nsim,1);
        end
        
        DiscPayoff_S=max(FlagCP*(S(:,end)-K),0)*exp(-r*T).*tempL.*tempU;
        DiscPayoff_AV=max(FlagCP*(S_AV(:,end)-K),0)*exp(-r*T).*tempL_AV.*tempU_AV;
        
        DiscPayoff=0.5*DiscPayoff_S+0.5*DiscPayoff_AV;
        
    elseif(FlagType==-1) %In option
        
        %Barrier use
        if(exist('L','var'))
            tempL=(min(S,[],2)<=L);
            tempL_AV=(min(S_AV,[],2)<=L);
        else
            tempL=ones(Nsim,1);
            tempL_AV=ones(Nsim,1);
        end
        
        if(exist('U','var'))
            tempU=(max(S,[],2)>=U);
            tempU_AV=(max(S_AV,[],2)>=U);
        else
            tempU=ones(Nsim,1);
            tempU_AV=ones(Nsim,1);
        end
        
        DiscPayoff_S=max(FlagCP*(S(:,end)-K),0)*exp(-r*T).*tempL.*tempU;
        DiscPayoff_AV=max(FlagCP*(S_AV(:,end)-K),0)*exp(-r*T).*tempL_AV.*tempU_AV;
        
        DiscPayoff=0.5*DiscPayoff_S+0.5*DiscPayoff_AV;
        
    else
        disp('Wrong Flag for type of option');
    end
    
    
else %%%%%%%%%%%%%%%%%Control Variables%%%%%%%%%%%%%%%%
    
    %small Asset simulation
    Ef=S0*exp(r*T);
    smallN=Nsim/100;
    S=underlyngBS(S0,r,sigma,T,smallN,Nsteps);
    f=S(:,end);
    
    %Distinction between OUT or IN
    if(FlagType==1) %Out opton
        
        %Barrier use
        if(exist('L','var'))
            tempL=(min(S,[],2)>L);
        else
            tempL=ones(smallN,1);
        end
        
        if(exist('U','var'))
            tempU=(max(S,[],2)<U);
        else
            tempU=ones(smallN,1);
        end
        
        %parameter determination
        g=max(FlagCP*(S(:,end)-K),0)*exp(-r*T).*tempL.*tempU;
        
        VC=cov(f,g);
        par=-VC(1,2)/VC(1,1);

        %real simulation
        S=underlyngBS(S0,r,sigma,T,Nsim,Nsteps);
    
        %Barrier use
        if(exist('L','var'))
            tempL=(min(S,[],2)>L);
        else
            tempL=ones(Nsim,1);
        end
        
        if(exist('U','var'))
            tempU=(max(S,[],2)<U);
        else
            tempU=ones(Nsim,1);
        end
        
        f=S(:,end);
        
        g=max(FlagCP*(S(:,end)-K),0)*exp(-r*T).*tempL.*tempU;
        
        DiscPayoff=g+par*(f-Ef);
        
        
    elseif(FlagType==-1) %In option
        
        %Barrier use
        if(exist('L','var'))
            tempL=(min(S,[],2)<=L);
        else
            tempL=ones(smallN,1);
        end
        
        if(exist('U','var'))
            tempU=(max(S,[],2)>=U);
        else
            tempU=ones(smallN,1);
        end
        
        g=max(FlagCP*(S(:,end)-K),0)*exp(-r*T).*tempL.*tempU;
        
        %parameter determination
        VC=cov(f,g);
        par=-VC(1,2)/VC(1,1);

        %real simulation
        S=underlyngBS(S0,r,sigma,T,Nsim,Nsteps);
        
        %Barrier use
        if(exist('L','var'))
            tempL=(min(S,[],2)<=L);
        else
            tempL=ones(Nsim,1);
        end
        
        if(exist('U','var'))
            tempU=(max(S,[],2)>=U);
        else
            tempU=ones(Nsim,1);
        end
        
        f=S(:,end);
        
        g=max(FlagCP*(S(:,end)-K),0)*exp(-r*T).*tempL.*tempU;
        
        DiscPayoff=g+par*(f-Ef);
        
    else
        disp('Wrong Flag for type of option');
    end
    
    
end

[Price,Dummy,CI]=normfit(DiscPayoff,alpha);

end
