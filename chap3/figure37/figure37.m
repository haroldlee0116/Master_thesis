%% Figure 3_7: Visualization of storage complexity
clear all 
close all
clc

d=3;
p=1;
l=1;
M=5;
I=p*M+1;
N=1000;

for r=1:50
    s_1(r)=N*l+(r^2)*(N*l*I+2*I+1);
    s_2(r)=N*l+N*l*r*r*I*I + r*r*I*I + 2*r*r*I;
    s_3(r)=N*l+N*r*I + I*r*l+2*(d-1)*r*r*I+N*r+r*r+r*r*I;
    s_4(r)=r*(I+1)*(3+N)+N*(I+2)+r*r*(I+1)^2;
end

figure(3)
semilogy(s_1,'LineWidth',2)
hold on
semilogy(s_2,'LineWidth',2)
semilogy(s_3,'LineWidth',2)
semilogy(s_4,'LineWidth',2)
s_6=ones(50,1)*(10^5);
semilogy(s_6,'-.r','LineWidth',3)
hold off
grid on
grid minor
legend('ALS','MALS','MPP','GN (Proposed)','Threshold')
xlim([1,50])
xlabel('Rank R')
ylabel('Storage Complexity (Byte)')
title('Maximum rank value within the storage complexity threshold')