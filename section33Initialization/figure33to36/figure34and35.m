%%  Complexity compensation rate visualization
% figure 3-4
close all; clear variables; clc;  
N = 2000;
a=1;
for I=1:30
        b=1;
        for R=1:30
            d=5;
            Cg=5*(N*(R^2)*(I^2)+(R^3)*(I^3));
            Ca=0.5*N*I*(nchoosek(I+d-1,I-1)^2)+(7/6)*(nchoosek(I+d-1,I-1)^3)+N*nchoosek(I+d-1,I-1);
            alphair5c(a,b) = 100*(Ca-Cg)/Cg;
            Sa=nchoosek(I+d-1,I-1)*(I^d);
            Sg=R*(I+1)*(I+4+N)+N;
            alphair5s(a,b) = 100*(Sa-Sg)/Sg;
            d=10;
            Cg=5*(N*(R^2)*(I^2)+(R^3)*(I^3));
            Ca=0.5*N*I*(nchoosek(I+d-1,I-1)^2)+(7/6)*(nchoosek(I+d-1,I-1)^3)+N*nchoosek(I+d-1,I-1);
            alphair10c(a,b) = 100*(Ca-Cg)/Cg;
            Sa=nchoosek(I+d-1,I-1)*(I^d);
            Sg=R*(I+1)*(I+4+N)+N;
            alphair10s(a,b) = 100*(Sa-Sg)/Sg;
            d=20;
            Cg=5*(N*(R^2)*(I^2)+(R^3)*(I^3));
            Ca=0.5*N*I*(nchoosek(I+d-1,I-1)^2)+(7/6)*(nchoosek(I+d-1,I-1)^3)+N*nchoosek(I+d-1,I-1);
            alphair20c(a,b) = 100*(Ca-Cg)/Cg;
            Sa=nchoosek(I+d-1,I-1)*(I^d);
            Sg=R*(I+1)*(I+4+N)+N;
            alphair20s(a,b) = 100*(Sa-Sg)/Sg;
            b=b+1;
        end
        a=a+1;
end
figure(1)
subplot(1,2,1)
surf(alphair5c,'FaceColor','r', 'FaceAlpha',0.9, 'EdgeColor','g')
hold on
surf(alphair10c,'FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor','r');
surf(alphair20c,'FaceColor','b', 'FaceAlpha',0.2, 'EdgeColor','b');
shading interp;
xlabel('I')
ylabel('R')
zlabel('$$CRR_{c}$$','interpreter','latex')
set(gca,'zscale','log')
set(gca,'FontSize',20) 
hold off
subplot(1,2,2)
surf(alphair5s,'FaceColor','r', 'FaceAlpha',0.9, 'EdgeColor','g')
hold on
surf(alphair10s,'FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor','r');
surf(alphair20s,'FaceColor','b', 'FaceAlpha',0.2, 'EdgeColor','b');
shading interp;
xlabel('I')
ylabel('R')
zlabel('$$CRR_{s}$$','interpreter','latex')
set(gca,'zscale','log')
set(gca,'FontSize',20) 
hold off

% figure 3-5
a=1;
for d=1:30
        b=1;
        for R=1:30
            I=5;
            Cg=5*(N*(R^2)*(I^2)+(R^3)*(I^3));
            Ca=0.5*N*I*(nchoosek(I+d-1,I-1)^2)+(7/6)*(nchoosek(I+d-1,I-1)^3)+N*nchoosek(I+d-1,I-1);
            alphadr5c(a,b) = 100*(Ca-Cg)/Cg;
            Sa=nchoosek(I+d-1,I-1)*(I^d);
            Sg=R*(I+1)*(I+4+N)+N;
            alphadr5s(a,b) = 100*(Sa-Sg)/Sg;
            I=10;
            Cg=5*(N*(R^2)*(I^2)+(R^3)*(I^3));
            Ca=0.5*N*I*(nchoosek(I+d-1,I-1)^2)+(7/6)*(nchoosek(I+d-1,I-1)^3)+N*nchoosek(I+d-1,I-1);
            alphadr10c(a,b) = 100*(Ca-Cg)/Cg;
            Sa=nchoosek(I+d-1,I-1)*(I^d);
            Sg=R*(I+1)*(I+4+N)+N;
            alphadr10s(a,b) = 100*(Sa-Sg)/Sg;
            I=20;
            Cg=5*(N*(R^2)*(I^2)+(R^3)*(I^3));
            Ca=0.5*N*I*(nchoosek(I+d-1,I-1)^2)+(7/6)*(nchoosek(I+d-1,I-1)^3)+N*nchoosek(I+d-1,I-1);
            alphadr20c(a,b) = 100*(Ca-Cg)/Cg;
            Sa=nchoosek(I+d-1,I-1)*(I^d);
            Sg=R*(I+1)*(I+4+N)+N;
            alphadr20s(a,b) = 100*(Sa-Sg)/Sg;
            b=b+1;
        end
        a=a+1;
end
figure(2)
subplot(1,2,1)
surf(alphadr5c,'FaceColor','r', 'FaceAlpha',0.9, 'EdgeColor','g')
hold on
surf(alphadr10c,'FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor','r');
surf(alphadr20c,'FaceColor','b', 'FaceAlpha',0.2, 'EdgeColor','b');
shading interp;
xlabel('d')
ylabel('R')
zlabel('$$CRR_{c}$$','interpreter','latex')
set(gca,'zscale','log')
set(gca,'FontSize',20) 
hold off
subplot(1,2,2)
surf(alphadr5s,'FaceColor','r', 'FaceAlpha',0.9, 'EdgeColor','g')
hold on
surf(alphadr10s,'FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor','r');
surf(alphadr20s,'FaceColor','b', 'FaceAlpha',0.2, 'EdgeColor','b');
shading interp;
xlabel('d')
ylabel('R')
zlabel('$$CRR_{s}$$','interpreter','latex')
set(gca,'zscale','log')
set(gca,'FontSize',20) 
hold off