%% Figure 1-2: Voltera kernel coefficient numbers
figure(1)
for x = 1:100
y1(x) = 5^x;
y2(x) = 10^x;
y3(x) = 20^x;
y4(x) = 40^x;
y11(x) = nchoosek(5-1+x,x);
y21(x) = nchoosek(10-1+x,x);
y31(x) = nchoosek(20-1+x,x);
y41(x) = nchoosek(40-1+x,x);
end
subplot(2,1,1)
loglog(y1,'r:>')
hold on
grid on
loglog(y2,'b:>')
loglog(y3,'g:>')
loglog(y4,'m:>')
ylabel('Number of coefficients for h_{d}')
xlabel('Order d')
legend('M=5','M=10','M=20','M=40')

subplot(2,1,2)
loglog(y11,'r-*')
hold on
grid on
loglog(y21,'b-*')
loglog(y31,'g-*')
loglog(y41,'m-*')
hold off
legend('M=5','M=10','M=20','M=40')
ylabel('Number of coefficients for h_{d}^{sym} ')
xlabel('Order d')