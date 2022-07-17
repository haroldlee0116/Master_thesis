%% Ploting recording2 from random initialization:
% download the 'recording2' data set from the following link:
% https://drive.google.com/file/d/1JuksjwmnNBnkdfPvhlnRjlUbBW_HvFYC/view?usp=sharing
close all; clear variables; clc; 
load recording2
subplot(2,3,1)
plot(0:recording2.Algorithm.iterations,recording2.Algorithm.fval,'-*')
set(gca, 'YScale', 'log')
title('Objective function value','FontSize', 12)
xlabel('Iteration','FontSize', 12)
ylabel('f(z_{i})','FontSize', 12)
xlim([0,recording2.Algorithm.iterations])
set(gca,'FontSize',18) % Creates an axes and sets its FontSize to 18
subplot(2,3,2)
plot(0:recording2.Algorithm.iterations,recording2.Algorithm.fval/frob(y),'-*')
set(gca, 'YScale', 'log')
title('Training error','FontSize', 12)
xlabel('Iteration','FontSize', 12)
ylabel('r_{t}(i)','FontSize', 12)
xlim([0,recording2.Algorithm.iterations])
set(gca,'FontSize',18) % Creates an axes and sets its FontSize to 18
subplot(2,3,3)
plot(1:recording2.Algorithm.iterations,recording2.Algorithm.relfval,'-*')
set(gca, 'YScale', 'log')
title('Relative objective function values','FontSize', 12)
xlabel('Iteration','FontSize', 12)
ylabel('f(z_{i-1})-f(z_{i})','FontSize', 12)
xlim([1,recording2.Algorithm.iterations])
set(gca,'FontSize',18) % Creates an axes and sets its FontSize to 18
subplot(2,3,4)
plot(1:recording2.Algorithm.iterations,recording2.Algorithm.relstep,'-*')
set(gca, 'YScale', 'log')
title('Newton step size','FontSize', 12)
xlabel('Iteration','FontSize', 12)
ylabel('p_{i}','FontSize', 12)
set(gca,'FontSize',18) % Creates an axes and sets its FontSize to 18
xlim([1,recording2.Algorithm.iterations])
subplot(2,3,5)
plot(recording2.Algorithm.rho,'-*')
set(gca, 'YScale', 'log')
title('Trust region','FontSize', 12)
xlabel('Iteration','FontSize', 12)
ylabel('\delta_{i}','FontSize', 12)
set(gca,'FontSize',18) % Creates an axes and sets its FontSize to 18
xlim([1,recording2.Algorithm.iterations])
subplot(2,3,6)
for i=1:recording2.Algorithm.iterations
    T_est= cpdgen([repmat(recording2.Algorithm.z{1,i}{1,1}(1),1,d),recording2.Algorithm.z{1,i}{1,1}(2)]);
    recording2.Algorithm.cof(i)= norm(symcheck(reshape(T_est,I*ones(1,d))));
end
plot(1:recording2.Algorithm.iterations,recording2.Algorithm.cof,'-*')
set(gca, 'YScale', 'log')
title('Symmetric coefficient','FontSize', 12)
xlabel('Iteration','FontSize', 12)
ylabel('s_{i}','FontSize', 12)
set(gca,'FontSize',18) % Creates an axes and sets its FontSize to 18
xlim([1,recording2.Algorithm.iterations])