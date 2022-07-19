%% Figure 4_4

% 'recording2' data set comes from the 187th row of 'exp5application'.
% It records the validation perofrmance of the developed framework on EEG 
% dataset for differnt M.

% download the  'recording2' data set from:
% https://drive.google.com/file/d/1R0vymGQqNizshtd5coL08WRqDV_mm3-F/view?usp=sharing

close all
clear all
clc
load recording2
figure(4)
for i=1:4
    subplot(4,1,i)
    plot(recording2.y_val);
    hold on
    plot(recording2.y_hat_val(:,i));
    vaf_val= (1-var(recording2.y_val-recording2.y_hat_val(:,i)))/var(recording2.y_val)*100;
    hold off
    xlabel('Sample')
    ylabel('Output')
    title(['M = ' num2str(i*10), ', VAF =' num2str(vaf_val) '%'])
    set(gca,'FontSize',20)
end
%% Figure 4_5
figure(6)
subplot(4,2,1)
autocorr(recording2.y_val-recording2.y_hat_val(:,1));
title('Auto-correlation for M=10')
ylabel('E[r(n)r(n-\tau)]')
xlabel('Lag \tau')
ylim([-1,1])
set(gca,'FontSize',20)

subplot(4,2,3)
autocorr(recording2.y_val-recording2.y_hat_val(:,2));
title('Auto-correlation for M=20')
ylabel('E[r(n)r(n-\tau)]')
xlabel('Lag \tau')
ylim([-1,1])
set(gca,'FontSize',20)

subplot(4,2,5)
autocorr(recording2.y_val-recording2.y_hat_val(:,3));
title('Auto-correlation for M=30')
ylabel('E[r(n)r(n-\tau)]')
xlabel('Lag \tau')
ylim([-1,1])
set(gca,'FontSize',20)

subplot(4,2,7)
autocorr(recording2.y_val-recording2.y_hat_val(:,4));
title('Auto-correlation for M=40')
ylabel('E[r(n)r(n-\tau)]')
xlabel('Lag \tau')
ylim([-1,1])
set(gca,'FontSize',20)


subplot(4,2,2)
crosscorr(recording2.y_val-recording2.y_hat_val(:,4),u_est(26:end,1))
xlabel('Lag \tau')
ylabel('E[r(n)u(n-\tau)]')
title('Cross-correlation for M=10')
set(gca,'FontSize',20)

subplot(4,2,4)
crosscorr(recording2.y_val-recording2.y_hat_val(:,2),u_est(26:end,1))
xlabel('Lag \tau')
ylabel('E[r(n)u(n-\tau)]')
title('Cross-correlation for M=20')
set(gca,'FontSize',20)

subplot(4,2,6)
crosscorr(recording2.y_val-recording2.y_hat_val(:,3),u_est(26:end,1))
xlabel('Lag \tau')
ylabel('E[r(n)u(n-\tau)]')
title('Cross-correlation for M=30')
set(gca,'FontSize',20)

subplot(4,2,8)
crosscorr(recording2.y_val-recording2.y_hat_val(:,1),u_est(26:end,1))
xlabel('Lag \tau')
ylabel('E[r(n)u(n-\tau)]')
title('Cross-correlation for M=40')
set(gca,'FontSize',20)