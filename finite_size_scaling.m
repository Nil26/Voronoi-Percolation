clc;clear;

res = load('voronoi_diff_network.txt');
N_many = [32,42,48,64,128,256];
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
Pv = polyfit(log(L),log(Delta),1);
v = -1/Pv(1);
yv_fit = L.^(-1/v)*Delta(1)/L(1)^(-1/v);

figure;
loglog(L,Delta,'^');
hold on;
loglog(L,yv_fit);
xlabel('L');
ylabel('\Delta');

% from the exponent v we got to find p_c
x = L.^(-1/v);
Ppc = polyfit(x,pav_f,1);
xfit = linspace(0,max(x));
ypc_fit = xfit*Ppc(1)+Ppc(2);

figure;
plot(x,pav_f,'o');
hold on;
plot(xfit,ypc_fit);
xlabel('L^{-1/\nu}');
ylabel('p_{av}');
