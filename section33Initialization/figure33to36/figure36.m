%% Figure 3_6: PCR Visualization

% 'pcr' is the result from the 66th-71th row of 'exp3PCR'. It contains the
% PCR value for different tensor size (d and I).

% down load the 'pcr' data set from the following link:
% https://drive.google.com/file/d/1nqt-irx457cMIGTcYd2S1cFfnt4FwGK6/view?usp=sharing
close all; clear variables; clc; 
load pcr
figure(4)
data = {[pcr(1:30),pcr(91:120),pcr(271:300)],[pcr(31:60),pcr(121:150), pcr(301:330)],[pcr(61:90),pcr(151:180),pcr(331:360)]};    
h = boxplotGroup(data,...
  'PrimaryLabels', {'I=9','I=15','I=25'}, ...
  'SecondaryLabels',{'d=3' 'd=4' 'd=5'},....
'InterGroupSpace',2, 'MedianStyle','target', 'OutlierSize',5,...
'Symbol','s','Widths',.5,'Colors','krg',...
'ExtremeMode','compress','jitter',0.9,'Whisker',1.2,'grouplines',true);
ylabel('PCR (%)');
h.axis.FontSize = 18;
h.axis2.FontSize = 18;