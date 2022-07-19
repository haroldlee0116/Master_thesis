%% Figure 2-4: Volterra tensor three formats
figure(2)
for x = 2:50
y1(x) = (5).^x;
y2(x) = (10).^x;
y3(x) = (20).^x;
y11(x) = nchoosek(5-1+x,5-1);
y21(x) = nchoosek(10-1+x,10-1);
y31(x) = nchoosek(20-1+x,20-1);
y12(x) = (5+1)*10;
y22(x) = (10+1)*10;
y32(x) = (20+1)*10;
y13(x) = (5+1)*50;
y23(x) = (10+1)*50;
y33(x) = (20+1)*50;
y14(x) = (5+1)*100;
y24(x) = (10+1)*100;
y34(x) = (20+1)*100;
end
subplot(1,3,1)
semilogy(y1,'r:o')
hold on
grid on
semilogy(y11,'b:+')
semilogy(y12,'g-p')
semilogy(y13,'m-p')
semilogy(y14,'k-p')
hold off
xlim([2,50])
title('I=5')
legend({'$\mathcal{V}$','$\mathcal{V}_{sym}$','CPD($\mathcal{V}_{sym}$) (R=10)', 'CPD($\mathcal{V}_{sym}$) (R=50)', 'CPD($\mathcal{V}_{sym}$) (R=100)'},'Interpreter','latex')
ylabel('Storage complexity')
xlabel('Order d')
set(gca,'FontSize',18)

subplot(1,3,2)
semilogy(y2,'r:o')
hold on
grid on
semilogy(y21,'b:+')
semilogy(y22,'g-p')
semilogy(y23,'m-p')
semilogy(y24,'k-p')
hold off
xlim([2,50])
title('I=10')
ylabel('Storage complexity')
xlabel('Order d')
set(gca,'FontSize',18)

subplot(1,3,3)
semilogy(y3,'r:o')
hold on
grid on
semilogy(y31,'b:+')
semilogy(y32,'g-p')
semilogy(y33,'m-p')
semilogy(y34,'k-p')
hold off
xlim([2,50])
title('I=20')
%legend({'$\mathcal{V}$','$\mathcal{V}_{sym}$','CPD($\mathcal{V}_{sym}$) (R=5)', 'CPD($\mathcal{V}_{sym}$) (R=20)', 'CPD($\mathcal{V}_{sym}$) (R=40)'},'Interpreter','latex')
ylabel('Storage complexity')
xlabel('Order d')
set(gca,'FontSize',18)