% Monte Carlo Simulation, Estimator for Windstream's values in the RDOF (Working paper)
% Ignacio Nunez, Economics Department, UT-Austin, ijnunez@utexas.edu

%% Create Sample

clear all;
close all;

N=1000; %sample size
x=rand(N,1); %common covariates
y=rand(N,1); %competitor covariates
d=rand(N,1); %Windstream covariates

intercept=ones(N,1);
%competitors parameters
mu_1=0; 
mu_1x=0.5;
mu_1y=2;
st_1=0.5;

%Windstream parameters
mu_2=0;
mu_2x=0.5;
mu_2d=2;
st_2=0.5;
kappa=0.5;

%Set number of simulations
simulations=1;
Results_coef=zeros(simulations,9);

for s=1:simulations
    %Create sample
    xi=randn(N,1)*st_1;
    epsilon=xi*kappa+randn(N,1)*st_2;
    X1 = intercept*mu_1+x*mu_1x+y*mu_1y+xi; %Competitors
    X2 = intercept*mu_2+x*mu_2x+d*mu_2d+epsilon; %Windstream
    X2_censored = zeros(N,1);  
    dummy_censored = zeros(N,1);    
    for i=1:N
        if X1(i)<=X2(i)
        X2_censored(i)=X2(i);  %uncensored, observe WT value
        else
        X2_censored(i)=X1(i);  %censored, observe competitors' value
        dummy_censored(i)=1;
        end
    end


    global data data2
    data=[x,y,d,X2_censored,dummy_censored];
    count_censored=sum(data(:,5));
    data2=zeros(count_censored,5);
    count=0;
    for i=1:N
    if data(i,5)==1    
    count=count+1;   
    data2(count,:)=data(i,:);
    end
    end
    
    global simdraws
    sim=10000;
    pd = makedist('Normal','mu',0,'sigma',1);
    simdraws = random(pd,sim,1);

    %Estimate parameters
    param=[0,0,0,0,0,0,0,1,1];
    real=[mu_2,mu_2x,mu_2d,mu_1,mu_1x,mu_1y,kappa,st_2,st_1];
    %param=[mu_1,mu_1x,mu_1y,mu_2,mu_2x,mu_2d,k,sigma2,sigma_xi];
    options = optimset('Display','iter','MaxFunEvals',1e3,'MaxIter',1e3);
    solution=fminsearch('Obj_Heckman_full',param,options);
    solution=fminunc('Obj_Heckman_full',solution,options);
    solution=fminsearch('Obj_Heckman_full',solution,options);
    %solution=fminsearch('Obj_Heckman_full',solution,options);
    %solution=fminsearch('Obj_Heckman_full',solution,options);
    %solution=fminsearch('Obj_Heckman_full',solution,options);

    Obj_Heckman_full(param)
    Obj_Heckman_full(real)
    Obj_Heckman_full(solution)

    Results_coef(s,:)=solution;    
    
end




