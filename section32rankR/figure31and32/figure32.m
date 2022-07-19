%% Figure 3_2: Rank value performance boxplots
close all; clear variables; clc;

% 'data' is collected from the 87th to 102th row of 'exp2RankR'. It
% contains the identification performance for differentt combination of R
% and I.

% download the 'data' file for ploting from the following link:
% https://drive.google.com/file/d/1Sqy2k21jWze_VnHFMAJ7jBQqsNY0VeYB/view?usp=sharing

load data
figure(2)
hold on
time = {data.exp1_R.relerr_optimization(:,1:8:end),data.exp1_R.relerr_optimization(:,2:8:end),data.exp1_R.relerr_optimization(:,3:8:end),data.exp1_R.relerr_optimization(:,4:8:end),data.exp1_R.relerr_optimization(:,5:8:end),data.exp1_R.relerr_optimization(:,6:8:end),data.exp1_R.relerr_optimization(:,7:8:end),data.exp1_R.relerr_optimization(:,8:8:end)};         
h = boxplotGroup(time, 'PrimaryLabels', {'R=1','R=2','R=3','R=4','R=5*','R=6','R=7','R=8'}, ...
  'SecondaryLabels',{'I=5' 'I=10' 'I=15','I=20'}, ....
'InterGroupSpace',3, 'MedianStyle','target', 'OutlierSize',5,...
'Symbol','s','Widths',.5,'Colors',lines(8),...
'ExtremeMode','compress','jitter',0.9,'Whisker',1.2,'GroupLines',true);
ylabel('Training error');
title('Influence of R for different I on training error')
h.axis.FontSize = 15;
h.axis2.FontSize = 15;
set(gca,'yscale','log')
hold off

figure(3)
hold on
time = {data.exp1_R.t_total(:,1:8:end),data.exp1_R.t_total(:,2:8:end),data.exp1_R.t_total(:,3:8:end),data.exp1_R.t_total(:,4:8:end),data.exp1_R.t_total(:,5:8:end),data.exp1_R.t_total(:,6:8:end),data.exp1_R.t_total(:,7:8:end),data.exp1_R.t_total(:,8:8:end)};         
t = boxplotGroup(time, 'PrimaryLabels', {'R=1','R=2','R=3','R=4','R=5*','R=6','R=7','R=8'}, ...
  'SecondaryLabels',{'I=5' 'I=10' 'I=15','I=20'}, ....
'InterGroupSpace',3, 'MedianStyle','target', 'OutlierSize',5,...
'Symbol','s','Widths',.5,'Colors',lines(8),...
'ExtremeMode','compress','jitter',0.9,'Whisker',1.2,'GroupLines',true);
ylabel('Training run time (s)');
title('Influence of R for different I on training run time (s)')
t.axis.FontSize = 15;
t.axis2.FontSize = 15;
set(gca,'yscale','log')
hold off

figure(4)
hold on
time = {data.exp1_R.relerr_validation(:,1:8:end),data.exp1_R.relerr_validation(:,2:8:end),data.exp1_R.relerr_validation(:,3:8:end),data.exp1_R.relerr_validation(:,4:8:end),data.exp1_R.relerr_validation(:,5:8:end),data.exp1_R.relerr_validation(:,6:8:end),data.exp1_R.relerr_validation(:,7:8:end),data.exp1_R.relerr_validation(:,8:8:end)};         
t = boxplotGroup(time, 'PrimaryLabels', {'R=1','R=2','R=3','R=4','R=5*','R=6','R=7','R=8'}, ...
  'SecondaryLabels',{'I=5' 'I=10' 'I=15','I=20'}, ....
'InterGroupSpace',3, 'MedianStyle','target', 'OutlierSize',5,...
'Symbol','s','Widths',.5,'Colors',lines(8),...
'ExtremeMode','compress','jitter',0.9,'Whisker',1.2,'GroupLines',true);
ylabel('Validation error');
title('Influence of R for different I on validation error')
t.axis.FontSize = 15;
t.axis2.FontSize = 15;
set(gca,'yscale','log')
hold off

figure(5)
hold on
time = {data.exp1_R.s(:,1:8:end),data.exp1_R.s(:,2:8:end),data.exp1_R.s(:,3:8:end),data.exp1_R.s(:,4:8:end),data.exp1_R.s(:,5:8:end),data.exp1_R.s(:,6:8:end),data.exp1_R.s(:,7:8:end),data.exp1_R.s(:,8:8:end)};         
t = boxplotGroup(time, 'PrimaryLabels', {'R=1','R=2','R=3','R=4','R=5*','R=6','R=7','R=8'}, ...
  'SecondaryLabels',{'I=5' 'I=10' 'I=15','I=20'}, ....
'InterGroupSpace',3, 'MedianStyle','target', 'OutlierSize',5,...
'Symbol','s','Widths',.5,'Colors',lines(8),...
'ExtremeMode','compress','jitter',0.9,'Whisker',1.2,'GroupLines',true);
ylabel('Symmetric coefficient');
title('Influence of R for different I on symmetric coefficient')
t.axis.FontSize = 15;
t.axis2.FontSize = 15;
set(gca,'yscale','log')
hold off
ylim([e-22,e-15])