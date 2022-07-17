%% Figure 4_6: Different R
% download the 'recording4' data set from:
% https://drive.google.com/file/d/1d0g9KVznKy8YLx2ind3OhhbTNIaNtcvl/view?usp=sharing
load recording4
figure(7)
R=[10,15,20,30,40,70];
for i=1:6
    subplot(3,2,i)
    plot(recording4.y_val);
    hold on
    plot(recording4.y_hat_val(:,i));
    vaf_val= (1-var(recording4.y_val-recording4.y_hat_val(:,i)))/var(recording4.y_val)*100;
    hold off
    legend({'$$\mathbf{y}^{7}_{v}$$','$$\hat{\mathbf{y}}^{7}_{v}$$'},'Interpreter','latex')
    xlabel('Sample')
    ylabel('Output')
    title(['R = ' num2str(R(i)), ', VAF =' num2str(vaf_val) '%'])
    set(gca,'FontSize',20)
end