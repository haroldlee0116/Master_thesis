%% Figure 3_8: Different methods comparision (with budget)
% download the 'recording' and 'recording2' data set from:
% https://drive.google.com/drive/folders/1RoLjXPPN0uIxnvRUNz2AuH8DKFIT1SDV?usp=sharing

clear all
close all
clc
load recording
figure(2)
subplot(2,2,1)
time = {recording.method1.train_error',recording.method2.train_error',recording.method3.train_error',recording.method5.train_error'};         
h = boxplotGroup(time, 'PrimaryLabels', {'A','B','C','D'}, ...
  'SecondaryLabels',{'SNR=+Inf' 'SNR=40' 'SNR=20','SNR=-20'}, ....
'InterGroupSpace',2, 'MedianStyle','target', 'OutlierSize',5,...
'Symbol','s','Widths',.5,'Colors',lines(8),...
'ExtremeMode','compress','jitter',0.9,'Whisker',1.2,'GroupLines',true);
ylabel('Traning error on y');
h.axis.FontSize = 18;
h.axis2.FontSize = 18;
grid on
set(gca,'yscale','log')

subplot(2,2,2)
time = {recording.method1.val_error',recording.method2.val_error',recording.method3.val_error',recording.method5.val_error'};         
h = boxplotGroup(time, 'PrimaryLabels', {'A','B','C','D'}, ...
  'SecondaryLabels',{'SNR=+Inf' 'SNR=40' 'SNR=20','SNR=-20'}, ....
'InterGroupSpace',2, 'MedianStyle','target', 'OutlierSize',5,...
'Symbol','s','Widths',.5,'Colors',lines(8),...
'ExtremeMode','compress','jitter',0.9,'Whisker',1.2,'GroupLines',true);
ylabel('Validation error on y');
h.axis.FontSize = 18;
h.axis2.FontSize = 18;
grid on
set(gca,'yscale','log')

subplot(2,2,3)
time = {recording.method1.time',recording.method2.time',recording.method3.time',recording.method5.time'};         
h = boxplotGroup(time, 'PrimaryLabels', {'A','B','C','D'}, ...
  'SecondaryLabels',{'SNR=+Inf' 'SNR=40' 'SNR=20','SNR=-20'}, ....
'InterGroupSpace',2, 'MedianStyle','target', 'OutlierSize',5,...
'Symbol','s','Widths',.5,'Colors',lines(8),...
'ExtremeMode','compress','jitter',0.9,'Whisker',1.2,'GroupLines',true);
ylabel('Running time (s)');
h.axis.FontSize = 18;
h.axis2.FontSize = 18;
grid on

subplot(2,2,4)
time = {recording.method1.symmetric',recording.method2.symmetric',recording.method3.symmetric',recording.method5.symmetric'};         
h = boxplotGroup(time, 'PrimaryLabels', {'A','B','C','D'}, ...
  'SecondaryLabels',{'SNR=+Inf' 'SNR=40' 'SNR=20','SNR=-20'}, ....
'InterGroupSpace',2, 'MedianStyle','target', 'OutlierSize',5,...
'Symbol','s','Widths',.5,'Colors',lines(8),...
'ExtremeMode','compress','jitter',0.9,'Whisker',1.2,'GroupLines',true);
ylabel('Symmetric coefficient');
h.axis.FontSize = 18;
h.axis2.FontSize = 18;
grid on
set(gca,'yscale','log')
%% Figure 3_9: Different methods comparision (without budget)
clear all 
close all 
clc
load recording2
figure(2)
subplot(2,2,1)
time = {recording.method1.train_error',recording.method2.train_error',recording.method3.train_error',recording.method5.train_error'};         
h = boxplotGroup(time, 'PrimaryLabels', {'A','B','C','D'}, ...
  'SecondaryLabels',{'SNR=+Inf' 'SNR=40' 'SNR=20','SNR=-20'}, ....
'InterGroupSpace',2, 'MedianStyle','target', 'OutlierSize',5,...
'Symbol','s','Widths',.5,'Colors',lines(8),...
'ExtremeMode','compress','jitter',0.9,'Whisker',1.2,'GroupLines',true);
ylabel('Traning error on y');
h.axis.FontSize = 18;
h.axis2.FontSize = 18;
grid on
set(gca,'yscale','log')

subplot(2,2,2)
time = {recording.method1.val_error',recording.method2.val_error',recording.method3.val_error',recording.method5.val_error'};         
h = boxplotGroup(time, 'PrimaryLabels', {'A','B','C','D'}, ...
  'SecondaryLabels',{'SNR=+Inf' 'SNR=40' 'SNR=20','SNR=-20'}, ....
'InterGroupSpace',2, 'MedianStyle','target', 'OutlierSize',5,...
'Symbol','s','Widths',.5,'Colors',lines(8),...
'ExtremeMode','compress','jitter',0.9,'Whisker',1.2,'GroupLines',true);
ylabel('Validation error on y');
h.axis.FontSize = 18;
h.axis2.FontSize = 18;
grid on
set(gca,'yscale','log')

subplot(2,2,3)
time = {recording.method1.time',recording.method2.time',recording.method3.time',recording.method5.time'};         
h = boxplotGroup(time, 'PrimaryLabels', {'A','B','C','D'}, ...
  'SecondaryLabels',{'SNR=+Inf' 'SNR=40' 'SNR=20','SNR=-20'}, ....
'InterGroupSpace',2, 'MedianStyle','target', 'OutlierSize',5,...
'Symbol','s','Widths',.5,'Colors',lines(8),...
'ExtremeMode','compress','jitter',0.9,'Whisker',1.2,'GroupLines',true);
ylabel('Running time (s)');
h.axis.FontSize = 18;
h.axis2.FontSize = 18;
grid on

subplot(2,2,4)
time = {recording.method1.symmetric',recording.method2.symmetric',recording.method3.symmetric',recording.method5.symmetric'};         
h = boxplotGroup(time, 'PrimaryLabels', {'A','B','C','D'}, ...
  'SecondaryLabels',{'SNR=+Inf' 'SNR=40' 'SNR=20','SNR=-20'}, ....
'InterGroupSpace',2, 'MedianStyle','target', 'OutlierSize',5,...
'Symbol','s','Widths',.5,'Colors',lines(8),...
'ExtremeMode','compress','jitter',0.9,'Whisker',1.2,'GroupLines',true);
ylabel('Symmetric coefficient');
h.axis.FontSize = 18;
h.axis2.FontSize = 18;
grid on
set(gca,'yscale','log')