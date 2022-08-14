% Monte Carlo - Tobit
% Ignacio Nunez, ijnunez@utexas.edu
% Constant   Constant-Urban	Area	Area-Urban	Connections	Connections-Urban

%% Create Sample

clear all;
close all;

N=1000;
x=rand(N,1); %common covariates
y=rand(N,1); %competitor covariates
d=rand(N,1); %Windstream covariates

intercept=ones(N,1);
%competitors parameters
mu_1=1; 
mu_1x=0.5;
mu_1y=0.5;
st_1=1;

%Windstream parameters
mu_2=0.5;
mu_2x=1;
mu_2d=0.3;
st_2=1;
kappa=-0.5;

simulations=50;
Results_coef=zeros(simulations,15);

for s=1:simulations
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
param=[1,1,1,1,1,1,1];
h=0.5547;

count_censored=sum(data(:,5));
data2=zeros(count_censored,5);
count=0;
for i=1:N
if data(i,5)==1    
count=count+1;   
data2(count,:)=data(i,:);
end
end

Results_coef(s,:)=Obj_MonteCarlo(param);    

% current=param;
% Obj_current=Obj_MonteCarlo_tobit(current);
% for j=1:grid
% for k=1:grid
%     test_mu=mu_min+(mu_max-mu_min)/grid*j;
%     test_sigma=sigma_min+(sigma_max-sigma_min)/grid*j;
%     test=[test_mu,test_sigma];
%     aux=Obj_MonteCarlo_tobit(test);
%     if aux<Obj_current
%     Obj_current=aux;
%     current=test;
%     current
%     end
% end
% end
% Results_coef_global(s,:)=current;    
end

Obj_MonteCarlo_tobit([1,2])

fig=figure;
subplot(1,2,1)
histogram(Results_coef(:,1))
xlabel('\mu')
title ('\mu=1, \sigma=2, S=100, N=1000')
subplot(1,2,2)
histogram(Results_coef(:,2))
xlabel('\sigma')
set(fig,'PaperOrientation','landscape');
print('Monte_Carlo_Starlink','-dpdf')
close all;

Results=[mean(Results_coef(:,1)) std(Results_coef(:,1));mean(Results_coef(:,2)) std(Results_coef(:,2))];



