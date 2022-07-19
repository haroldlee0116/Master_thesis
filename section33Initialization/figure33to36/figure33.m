%% Figure 3-3: Algebriac and random initializaiton (one trial)
close all; clear variables; clc; 
% download the 'recording1' and 'recording2' data set from the following
% link
% https://drive.google.com/drive/folders/1us0RHSgVMUb1_0ROEOWwL3ohBlylcLWB?usp=sharing
load recording1
load recording2
subplot(2,3,1)
plot(0:recording1.Algorithm.iterations,recording1.Algorithm.fval,'r-*')
hold on
plot(0:recording2.Algorithm.iterations,recording2.Algorithm.fval,'b-*')
set(gca, 'YScale', 'log')
title('Objective function value','FontSize', 12)
xlabel('Iteration','FontSize', 12)
ylabel('f(z_{i})','FontSize', 12)
legend('Algebraic Initialization','Random Initialization')
hold off

set(gca,'FontSize',18) % Creates an axes and sets its FontSize to 18
subplot(2,3,2)
plot(0:recording1.Algorithm.iterations,recording1.Algorithm.fval*frob(y),'r-*')
hold on
plot(0:recording2.Algorithm.iterations,recording2.Algorithm.fval*frob(y),'b-*')
set(gca, 'YScale', 'log')
title('Training error','FontSize', 12)
xlabel('Iteration','FontSize', 12)
ylabel('r_{t}(i)','FontSize', 12)
set(gca,'FontSize',18) % Creates an axes and sets its FontSize to 18
hold off

subplot(2,3,3)
plot(1:recording1.Algorithm.iterations,recording1.Algorithm.relfval,'r-*')
hold on
plot(1:recording2.Algorithm.iterations,recording2.Algorithm.relfval,'b-*')
set(gca, 'YScale', 'log')
title('Relative objective function values','FontSize', 12)
xlabel('Iteration','FontSize', 12)
ylabel('f(z_{i-1})-f(z_{i})','FontSize', 12)
set(gca,'FontSize',18) % Creates an axes and sets its FontSize to 18
hold off

subplot(2,3,4)
plot(1:recording1.Algorithm.iterations,recording1.Algorithm.relstep,'r-*')
hold on
plot(1:recording2.Algorithm.iterations,recording2.Algorithm.relstep,'b-*')
set(gca, 'YScale', 'log')
title('Newton step size','FontSize', 12)
xlabel('Iteration','FontSize', 12)
ylabel('p_{i}','FontSize', 12)
set(gca,'FontSize',18) % Creates an axes and sets its FontSize to 18
hold off


subplot(2,3,5)
plot(recording1.Algorithm.rho,'r-*')
hold on
plot(recording2.Algorithm.rho,'b-*')
set(gca, 'YScale', 'log')
title('Trust region','FontSize', 12)
xlabel('Iteration','FontSize', 12)
ylabel('\delta_{i}','FontSize', 12)
set(gca,'FontSize',18) % Creates an axes and sets its FontSize to 18

hold off


subplot(2,3,6)
for i=1:recording1.Algorithm.iterations
    T_est= cpdgen([repmat(recording1.Algorithm.z{1,i}{1,1}(1),1,d),recording1.Algorithm.z{1,i}{1,1}(2)]);
    recording1.Algorithm.cof(i)= norm(symcheck(reshape(T_est,I*ones(1,d))));
end
for i=1:recording2.Algorithm.iterations
    T_est= cpdgen([repmat(recording2.Algorithm.z{1,i}{1,1}(1),1,d),recording2.Algorithm.z{1,i}{1,1}(2)]);
    recording2.Algorithm.cof(i)= norm(symcheck(reshape(T_est,I*ones(1,d))));
end
plot(1:recording1.Algorithm.iterations,recording1.Algorithm.cof,'r-*')
hold on
plot(1:recording2.Algorithm.iterations,recording2.Algorithm.cof,'b-*')
set(gca, 'YScale', 'log')
title('Symmetric coefficient','FontSize', 12)
xlabel('Iteration','FontSize', 12)
ylabel('s_{i}','FontSize', 12)
set(gca,'FontSize',18) % Creates an axes and sets its FontSize to 18
hold off