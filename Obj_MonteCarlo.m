function ret=Obj_MonteCarlo(param0)

global data beta_wt data2

%% Estimate probit model
param_num=3;
beta_wt=transpose([0.28,-0.28,0.28,-1.11]);
param=beta_wt;

%Probit_Heckman(beta_wt)
options = optimset('Display','iter','MaxFunEvals',4e3,'MaxIter',4e3);
%options = optimset('MaxFunEvals',4e3,'MaxIter',4e3);
solution=fminunc('Probit_Heckman',param,options);
%solution=fminsearch('Probit_Heckman',solution,options);
beta_wt=solution;
param=beta_wt;
x=data(:,1);
y=data(:,2);
z=data(:,3);


clear param solution
alpha_wt=transpose([0.5,1,-0.5,1.8,1]); %alpha, ratio, ratio2 ,sigma_uncond,
%Obj_Heckman_correlated_particle(alpha_wt)
param=alpha_wt;

%alpha_wt=transpose([1,0.1,-1,3,0.5]); %alpha, ratio, ratio2 ,sigma_uncond,
%Obj_Heckman_correlated_particle(alpha_wt)

% Obj_Heckman_correlated(alpha_wt)
% Obj_Heckman_correlated_particle(alpha_wt)
% fun = @(x)Obj_Heckman_correlated_particle(x);
% lb = [-5,-5,-5,-5,-5];
% ub = [5,5,5,5,5];
% %options = optimoptions('particleswarm','SwarmSize',500,'Display','iter','PlotFcn','pswplotbestf','FunctionTolerance',1e-7,'HybridFcn',@fmincon);
% options = optimoptions('particleswarm','SwarmSize',1000,'FunctionTolerance',1e-6,'HybridFcn',@fmincon,'Display','iter','PlotFcn','pswplotbestf');
% nvars = 5;
% solution_par = particleswarm(fun,nvars,lb,ub,options);
% like1=Obj_Heckman_correlated_particle(solution_par);
% solution_par2 = particleswarm(fun,nvars,lb,ub,options);
% like2=Obj_Heckman_correlated_particle(solution_par);
% if like1>like2
% param=transpose(solution_par); 
% else
% param=transpose(solution_par2);
% end

options = optimset('Display','iter','MaxFunEvals',4e3,'MaxIter',4e3);
%options = optimset('MaxFunEvals',4e3,'MaxIter',4e3);
%solution=fminunc('Obj_Heckman_correlated',alpha_wt,options);
solution=fminsearch('Obj_Heckman_correlated_particle',param,options);
solution=fminunc('Obj_Heckman_correlated_particle',solution,options);

param=solution;
k=min(param(3),0.99);
h=max(abs(param(4)),0.001);
alpha_wt=zeros(3,1);
alpha_wt(1)=param(1);
alpha_wt(2)=param(2);
alpha_wt(3)=-beta_wt(4)*h;
delta=zeros(3,1);
delta(1)=beta_wt(1)*h+alpha_wt(1);
delta(2)=beta_wt(2)*h+alpha_wt(2);
delta(3)=beta_wt(3)*h;
sigma2 = min(max(abs(param(5)),0.001),h-0.0001);
sigma_xi=((h^2-sigma2^2)^(0.5))/abs(1-k);
sigma1=max(((k^2)*(sigma_xi^2)+sigma2^2)^0.5,0.001);
ret=[beta_wt(1),beta_wt(2),beta_wt(3),beta_wt(4),alpha_wt(1),alpha_wt(2),alpha_wt(3),delta(1),delta(2),delta(3),h,k,sigma2,sigma_xi,sigma1];




param=[alpha_wt(1),alpha_wt(2),alpha_wt(3),delta(1),delta(2),delta(3),k,sigma2,sigma_xi];
options = optimset('Display','iter','MaxFunEvals',4e3,'MaxIter',4e3);
%options = optimset('MaxFunEvals',4e3,'MaxIter',4e3);
%solution=fminunc('Obj_Heckman_correlated',alpha_wt,options);
Obj_Heckman_full(param)
solution=fminsearch('Obj_Heckman_full',param,options);
solution=fminunc('Obj_Heckman_full',solution,options);





% 
% 
% 
% current_sol=param;
% current_sign=0;
% current_loglike=10000;
% 
% sign_k=0;
% solution=fminunc('Obj_Heckman_correlated',param,options);
% solution=fminsearch('Obj_Heckman_correlated',solution,options);
% solution=fminunc('Obj_Heckman_correlated',solution,options);
% if logll<current_loglike
% current_loglike=logll;
% current_sol=solution;
% end
% 
% sign_k=1;
% solution=fminunc('Obj_Heckman_correlated',param,options);
% solution=fminsearch('Obj_Heckman_correlated',solution,options);
% solution=fminsearch('Obj_Heckman_correlated',solution,options);
% solution=fminunc('Obj_Heckman_correlated',solution,options);
% if logll<current_loglike
% current_loglike=logll;
% current_sol=solution;
% current_sign=1;
% end
% 
% sign_k=-1;
% solution=fminunc('Obj_Heckman_correlated',param,options);
% solution=fminsearch('Obj_Heckman_correlated',solution,options);
% solution=fminunc('Obj_Heckman_correlated',solution,options);
% if logll<current_loglike
% current_loglike=logll;
% current_sol=solution;
% current_sign=2;
% end
% 
% param=solution;
% alpha=param;
% 
% 
% clear param solution
% param_num=5;
% delta=transpose([1,0.5,0.5,0.66,1.2]);
% param=delta;
% 
% solution=fminunc('Obj_Heckman_correlated_w',param,options);
% solution=fminsearch('Obj_Heckman_correlated_w',solution,options);
% 
% param=solution;
% delta=param;

%ret=[beta(1),beta(2),beta(3),alpha(1),alpha(2),alpha(3),alpha(4),delta(1),delta(2),delta(3),delta(4),delta(5)];
disp('done_solver')
end