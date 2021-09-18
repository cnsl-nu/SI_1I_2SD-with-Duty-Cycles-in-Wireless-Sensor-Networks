clear; clc; close all;
addpath(genpath('./helpers/'));
addpath(genpath('./GEMF/'));
rng(19)
global RunTime
% # Initial Setup
r0 = 3;
monte_rounds = 100;

dim = [30,30];
N = dim(1)*dim(2);
RunTime = 50;

% generate network
Net1 = NetGen_GeoUniform(N,r0,dim,1);
Net2 = NetGen_GeoRandom(N,r0,dim);
NetUni = NetCmbn({Net1, Net1});
NetRnd = NetCmbn({Net2, Net2});

% ## Parameters and initial conditions
% % R1<1; R2<1
% alpha = [0.04, 0.05]; % infect rate 
% mu = [0.04, 0.04]; % sleep s 
% gamma = [0.3, 0.3]; % rec rate 
% lambda = [0.05, 0.2]; % sleep I1 
% kappa = [0.04, 0.1]; % sleep I2    % [I2_s -> I2_a, I2_a -> I2_s]

% % R1>1; R2<1
% alpha = [0.25, 0.05]; % infect rate 
% mu = [0.04, 0.04]; % sleep s
% gamma = [0.25, 0.5]; % rec rate 
% lambda = [0.05, 0.12]; % sleep I1 
% kappa = [0.04, 0.15]; % sleep I2    % [I2_s -> I2_a, I2_a -> I2_s]

% % R1<1; R2>1
% alpha = [0.04, 0.15]; % infect rate 
% mu = [0.04, 0.04]; % sleep s 
% gamma = [0.25, 0.07]; % rec rate 
% lambda = [0.06, 0.12]; % sleep I1 
% kappa = [0.05, 0.15]; % sleep I2    % [I2_s -> I2_a, I2_a -> I2_s]

% R1>1  R2>1
alpha = [0.25, 0.35];   % infect rate
mu = [0.04, 0.04];     % sleep s
gamma = [0.04, 0.35];  % rec rate
lambda = [0.06, 0.07];  % sleep I1
kappa = [0.05, 0.06];  % sleep I2    % [I2_s -> I2_a, I2_a -> I2_s]
% initial infection
I1_a_initial_uni=10;
I2_a_initial_uni=10; 
I1_a_initial_rnd=10;
I2_a_initial_rnd=10;
% struct that hold network parameters
Para = Para_active_sleep_SI1I2S(alpha, mu, gamma, lambda, kappa); 
M = Para{1}; StopCond={'RunTime', RunTime};

x0_uni = Initial_Cond_Gen(N,'Population',[3,5],[I1_a_initial_uni,I2_a_initial_uni]);
x0_rnd = Initial_Cond_Gen(N,'Population',[3,5],[I1_a_initial_rnd,I2_a_initial_rnd]);
x0 = {x0_uni, x0_rnd};

init_infection_uni = [I1_a_initial_uni, I2_a_initial_uni];
init_infection_rnd = [I1_a_initial_rnd, I2_a_initial_rnd];
init_infection = {init_infection_uni, init_infection_rnd};
paramet = {alpha, mu, gamma, lambda, kappa};
%
[R0_uni,R1_uni,R2_uni] = calc_R0(NetUni, alpha, mu, gamma, lambda, kappa, N);
[R0_rnd,R1_rnd,R2_rnd] = calc_R0(NetRnd, alpha, mu, gamma, lambda, kappa, N);

[R0_uni, R1_uni, R2_uni; R0_rnd ,R1_rnd, R2_rnd]
%% Visualize network 
renderNetwork(Net1,21, 'Uniform', x0_uni);
renderNetwork(Net2,22, 'Random', x0_rnd);
%% Homogen
[tt, sol_values] = ode_homogenious(alpha, mu, gamma, lambda, kappa,...
    RunTime,N,I1_a_initial_uni,I2_a_initial_uni,r0);

%% Monte Carlo
rng(42)
[t, uni_, rnd_] = MonteCarlo(monte_rounds,N,Para,NetUni,NetRnd,x0,StopCond,init_infection);
% plot_MonteCarlo(paramet,uni_,rnd_,N,t,NetUni,NetRnd, false)

%% ODE
[t, Xuni, Xrnd] = ode(N,Para,NetUni,NetRnd,x0,StopCond);
% plot_ode(paramet,Xuni,Xrnd,N,t,NetUni,NetRnd,false)

%%
tit = {'Sa','Ss','Infected_a 1 vs Time','Infected_s 1 vs Time',...
    'Infected_a 2 vs Time','Infected_s 2 vs Time'};
i = 1;
variable1 = RunTime/.05+1;
shape = [4,6,variable1];
export = zeros(shape);
for z=[3,5,4,6]  % [3,5,4,6,1,2]
    export(i,1,:) = t;
    export(i,2,:) = uni_(z,:)./N;
    export(i,3,:) = rnd_(z,:)./N;
    export(i,4,:) = Xuni(z,:)./N;
    export(i,5,:) = Xrnd(z,:)./N;
    export(i,6,:) = sol_values(z,:)./N;
    
    figure(i)
    plot(t,uni_(z,:)./N, t,rnd_(z,:)./N); 
    hold on
    plot(t,Xuni(z,:)./N,'-.b','linewidth',1);
    plot(t,Xrnd(z,:)./N,'--r','linewidth',1);
    plot(t,sol_values(z,:)./N,'--k','linewidth',1)
    title(tit{z})
    legend('I1_a uni','I1_a rnd','ODE uni','ODE rnd','homogen','Location','northwest');
    xlim([0, RunTime])
%     ylim([0,1])
    grid on
    hold off
    i = i + 1;
end;
%%
export_i1a=reshape(export(1,:,:),6,variable1)';
export_i2a=reshape(export(2,:,:),6,variable1)';
export_i1s=reshape(export(3,:,:),6,variable1)';
export_i2s=reshape(export(4,:,:),6,variable1)';

plot(t,export_i1a(:,[2,3,4,5,6]))
% dlmwrite('1export_i1a.txt',export_i1a(1:15:end,:),'delimiter','\t')
% dlmwrite('1export_i2a.txt',export_i2a(1:15:end,:),'delimiter','\t')
% dlmwrite('1export_i1s.txt',export_i1s(1:15:end,:),'delimiter','\t')
% dlmwrite('1export_i2s.txt',export_i2s(1:15:end,:),'delimiter','\t')
