clc;clear;

% res = load('voronoi_diff_network.txt');
% N_many = [32,42,48,64,128,256];
res = load('delaunay_diff_network.txt');
N_many = [32,42,48,64];
pav_f = [];
Delta = [];

for n = N_many
    p_first_time = res(res(:,1)==n,2);
    pav = mean(p_first_time);
    pav_f = [pav_f pav];
    del2 = mean(p_first_time.^2)-mean(p_first_time)^2;
    Delta = [Delta sqrt(del2)];
end

L = sqrt(N_many);

% fit from \Delta ~ L^{-1/v} to find v
lfv = fit(log(L)',log(Delta)','poly1');
pv = coeffvalues(lfv);
nu = -1.0/pv(1);
confidence = confint(lfv,0.95);
stddnu = diff(confidence(:,1))/2*1/pv(1)^2*1;
yv_fit = L.^(-1/nu)*Delta(1)/L(1)^(-1/nu);

figure;
loglog(L,Delta,'^','MarkerSize',10,'LineWidth',2);
hold on;
loglog(L,yv_fit,'LineWidth',2,'LineStyle',':');
xlabel('L','FontSize',14);
ylabel('\Delta','FontSize',14);
title(sprintf("\\nu = %.4f \\pm %.4f",nu,stddnu),'FontSize',14);

% from the exponent v we got to find p_c
x = L.^(-1/nu);
[Ppc,S] = polyfit(x,pav_f,1); 
xfit = linspace(0,max(x));
[ypc_fit,delta_pc] = polyval(Ppc,xfit,S);

figure;
%plot(x,pav_f,'o','MarkerSize',10,'LineWidth',2);
errorbar(x,pav_f,Delta,'o','MarkerSize',10,'LineWidth',2);
hold on;
plot(xfit,ypc_fit,'LineWidth',2,'LineStyle',':');
xlabel('L^{-1/\nu}','FontSize',14);
ylabel('p_{av}','FontSize',14);
title(sprintf("p_c = %.4f \\pm %.4f",ypc_fit(1),2*delta_pc(1)),'FontSize',14);

