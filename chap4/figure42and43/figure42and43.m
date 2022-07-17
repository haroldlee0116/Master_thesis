% download the 'recording' data set from:
% https://drive.google.com/file/d/17kT2B3BW2XB8qL6ju9lv0YJsqzU4NoFz/view?usp=sharing
close all
clear all
clc
load recording
%%  Figure 4_2
figure(1)
for d=1:3
subplot(3,2,2*d-1)
plot(recording.y_est,'r')
hold on
grid on
plot(recording.y_hat_est(:,d),'b')
hold off
grid off
xlabel('Sample')
ylabel('Output value')
vaf_est=(1-var(recording.y_est-recording.y_hat_est(:,d))/var(recording.y_est))*100;
legend({'$$\mathbf{y}^{7}_{t}$$','$$\hat{\mathbf{y}}^{7}_{t}$$'},'Interpreter','latex')
title(['Estimation performance: VAF_{t}= ' num2str(vaf_est) '% (d= ' num2str(d) ')'])
xlim([1,1200]);

subplot(3,2,2*d)
plot(recording.y_val,'black')
hold on
grid on
plot(recording.y_hat_val(:,d),'m')
hold off
grid off
xlabel('Sample')
ylabel('Output value')
vaf_val= (1-var(recording.y_val-recording.y_hat_val(:,d))/var(recording.y_val))*100;
legend({'$$\mathbf{y}^{7}_{v}$$','$$\hat{\mathbf{y}}^{7}_{v}$$'},'Interpreter','latex')
title(['Validation performance: VAF_{v}= ' num2str(vaf_val) '% (d= ' num2str(d) ')'])
xlim([1,200]);
end
%% Figure 4_3
figure(2)
subplot(3,2,1)
autocorr(recording.y_val-recording.y_hat_val(:,1));
title('Auto-correlation for d=1')
ylabel('E[r(n)r(n-\tau)]')
xlabel('Lag \tau')
ylim([-1,1])
set(gca,'FontSize',20)

subplot(3,2,3)
autocorr(recording.y_val-recording.y_hat_val(:,2));
title('Auto-correlation for d=2')
ylabel('E[r(n)r(n-\tau)]')
xlabel('Lag \tau')
ylim([-1,1])
set(gca,'FontSize',20)

subplot(3,2,5)
autocorr(recording.y_val-recording.y_hat_val(:,3));
title('Auto-correlation for d=3')
ylabel('E[r(n)r(n-\tau)]')
xlabel('Lag \tau')
ylim([-1,1])
set(gca,'FontSize',20)

subplot(3,2,2)
crosscorr(recording.y_val-recording.y_hat_val(:,1),u_est(26:end,1))
xlabel('Lag \tau')
ylabel('E[r(n)u(n-\tau)]')
title('Cross-correlation for d=1')
set(gca,'FontSize',20)
subplot(3,2,4)
crosscorr(recording.y_val-recording.y_hat_val(:,2),u_est(26:end,1))
xlabel('Lag \tau')
ylabel('E[r(n)u(n-\tau)]')
title('Cross-correlation for d=2')
set(gca,'FontSize',20)
subplot(3,2,6)
crosscorr(recording.y_val-recording.y_hat_val(:,3),u_est(26:end,1))
xlabel('Lag \tau')
ylabel('E[r(n)u(n-\tau)]')
title('Cross-correlation for d=3')
set(gca,'FontSize',20)