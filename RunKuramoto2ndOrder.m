%A software package for optimizing synchronization of coupled oscillators with high-order networks
%(c) 2021 Ying Tang, Dinghua Shi, Linyuan Lv
%All rights reserved. 
%This MATLAB code package optimizes network topology for synchronization of coupled oscillators 
%with high-order interactions. The current focus is the system with Kuramoto-type coupling function 
%for identical oscillators, the second-order interactions (triangle). The optimization is realized by 
%minimizing the eigenratios or the spread of eigenvalues for the generalized Laplacian matrices. For the undirected network, 
%we rewire the triangle interactions and use simulated annealing to optimize the network synchronizability. 
%For the directed network, we selectively remove directional triangle interactions to optimize synchronizability, 
%and investigate asymmetry for the optimized directed network.

%A detailed description on the scripts is in README file. 
%Contact: Ying Tang, jamestang23@gmail.com

%%
clear;
global omega  kappa n Tensor
delta = 100;  % The range of time for figure slide change
n = 10;% Number of oscillators.
FreqMean=0; % Mean of omega distribution
beta =1;% Gaussian distribution's variance (0 will make the initial uniformly distributed state as a steady state), or Half-width of omega interval: use 0 for now %.660*kappa;
%omega = FreqMean+beta*tan(pi*(rand(n,1)-1/2));%draw random number from a Cauchy Distribution
omega = randn(n,1);%draw random number from a Cauchy Distribution
load('Rewire_1_InitialNet_2.mat')
%load('Triangles.mat')
Summary.KcScan=[0.3:-0.005:0.1]*n^2;%Works good for delta 100 and n 100; [5000:-1000:1000];%[10000:-2:0];%10;%[0.01:0.5:1];
step = 0.1;  % Step size.
s = max(step,.01);
SteadyTime=round(delta/step*0.2);
SteadyTime2=round(delta/step*0.5);

for Model=2:3%; % 1 is all-to-all, 2 is initial network, 3 is optimized network

%%
if Model==1
    Tensor=1;
    Name='All';
elseif  Model==2
    Tensor=Summary.AdjTensorInitial{n};
    Name='Initial';
else
   Tensor=Summary.AdjTensorOptimal{n};
    Name='Optimized';
end
filename=[pwd,'\',Name,'_T',num2str(delta),'n',num2str(n),'MeanF',num2str(FreqMean),'VarF',num2str(beta)];



tspan = 0:s:delta;
theta = zeros(n,1);%(1:n)'/n*2*pi; 
tic;
for i=1:length(Summary.KcScan)
    kappa=Summary.KcScan(i);
    [t,theta] = ode45(@ode,tspan,theta);%,options);
    for j = 1:length(t)
    z = 1/n*sum(exp(1i*theta(j,:))); %psi = angle(z);%z = abs(z);
    H.psi(j) = angle(z);
    H.order(j) = abs(z);
    z2 = 1/n*sum(exp(2i*theta(j,:))); %psi = angle(z);%z = abs(z);
    H.psi2(j) = angle(z2);
    H.order2(j) = abs(z2);
    end
    
    Summary.OrderPara{Model}(i)=mean(H.order(end-SteadyTime:end));
    Summary.OrderPara_T2{Model}(i)=mean(H.order(end-SteadyTime2:end));
    Summary.OrderPara2{Model}(i)=mean(H.order2(end-SteadyTime:end));
    Summary.OrderPara2_T2{Model}(i)=mean(H.order2(end-SteadyTime2:end));
    
    figure
 plot(H.order)
 ylim([0 1.1]);
% pause(1);
 if i==1
 figurename2=[filename,'_loss.jpg'];
saveas(gcf,figurename2);
 end
 close all;
 toc;
 
 theta=theta(end,:)';
 

end



%%
 figure('Position',[100 00 800 600]);
plot(Summary.KcScan,Summary.OrderPara{Model},'bo','markerFaceColor','b');hold on;
plot(Summary.KcScan,Summary.OrderPara_T2{Model},'ro','markerFaceColor','r');
plot(Summary.KcScan,Summary.OrderPara2{Model},'bx','markerFaceColor','b');
plot(Summary.KcScan,Summary.OrderPara2_T2{Model},'rx','markerFaceColor','r');
xlabel('Coupling constant')
xlabel('Order parameter');
ylim([0 1]);
xlim([min(Summary.KcScan), max(Summary.KcScan)]);
set(gca,'FontSize',12);
figurename2=[filename,'.jpg'];
saveas(gcf,figurename2);
     
fileoutput=[filename,'.mat'];
save(fileoutput,'Summary');
 
close all;
end
%%
function theta_dot =  ode(~,theta)%,n,kappa,beta,FreqMean,Tensor)%(~,theta)
global omega kappa n Tensor
       theta_dot = omega  - kappa/n^2*sum(sum(Tensor.*sin(repmat(2*theta-theta',[1 1 length(theta)])...
            -permute(repmat(theta,[1 length(theta) length(theta)]),[3,2,1])),3),2);
       
end