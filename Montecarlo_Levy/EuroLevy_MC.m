function [Price,CI]=EuroLevy_MC(S0,K,r,T,FlagCP,N,model,varargin)



if strcmpi(model,'kou')
    
    gamma=varargin{1}; sigma=varargin{2}; lambda=varargin{3};
    lambda_u=varargin{4}; lambda_d=varargin{5}; p_up=varargin{6};
    X=MCL_Kou(gamma,sigma,lambda,lambda_u,lambda_d,p_up,T,1,N)
    char=@(u) charExp_Kou(u,gamma,sigma,lambda,lambda_u,lambda_d,p_up);
    
elseif strcmpi(model,'merton')
    
    gamma=varargin{1}; sigma=varargin{2}; lambda=varargin{3};
    mu=varargin{4}; delta=varargin{5};
    X=MCL_Merton(gamma,sigma,lambda,mu,delta,T,1,N);
    char=@(u) charExp_Merton(u,gamma,sigma,lambda,mu,delta);
    
elseif strcmpi(model,'vg')
    
    theta=varargin{1}; sigma=varargin{2}; kappa=varargin{3};
    X=MCL_VG(theta,sigma,kappa,T,1,N);
    char=@(u) charExp_VG(u,theta,sigma,kappa);
    
elseif strcmpi(model,'nig')
    
    theta=varargin{1}; sigma=varargin{2}; kappa=varargin{3};
    X=MCL_NIG(theta,sigma,kappa,T,1,N);
    char=@(u) charExp_NIG(u,theta,sigma,kappa);
    
end
char(-1i)
ST=getExpRiskNeutral(S0,r,T,X,char);

DiscPayoff=max(FlagCP*(ST(:,end)-K),0)*exp(-r*T);

[Price,Dummy,CI]=normfit(DiscPayoff);

end