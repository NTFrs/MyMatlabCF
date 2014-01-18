function [Price,ExerciseTime]=AmericanBS_MC(S0,K,r,d,sigma,T,FlagCP,Nsims,Nsteps)

Paths=underlyngBS_div(S0,r,d,sigma,T,Nsims,Nsteps);
Paths=Paths(:,2:end);
ExerciseTime=Nsteps*ones(Nsims,1);
dt=T/Nsteps;
CashFlows=max(FlagCP*(Paths(:,end)-K),0);

for j=Nsteps-1:-1:1
    InMoney=find(FlagCP*(Paths(:,j)-K)>0);
    S=Paths(InMoney,j);
    
    %regression part: we find CV with a regression
    RegrMat=[ones(length(S),1), S, S.^2];
    Ydata=CashFlows(InMoney).*exp(-r*dt*(ExerciseTime(InMoney)-j));
    alpha=RegrMat\Ydata;
    
    %Intrinsec and continuation values
    IV=FlagCP*(S-K);
    CV=RegrMat*alpha;
    
    Index=find(IV>CV);
    ExercisePaths=InMoney(Index);
    
    %update CashFlows
    CashFlows(ExercisePaths)=IV(Index);
    
    %update exercise time
    ExerciseTime(ExercisePaths)=j;
    
end

Price=max(FlagCP*(S0-K),mean(CashFlows.*exp(-r*dt*ExerciseTime)));

end