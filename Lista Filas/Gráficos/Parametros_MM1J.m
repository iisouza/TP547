
clc, clear all, close all
p = 0.8; y = 200;
P0 = (1-p)/(1-p^(J+2));

buf = 50;
Eq = zeros(1, buf+1); Etq = zeros(1, buf+1); Pb = zeros(1, buf+1);

for J = 0:buf
		Pb(J+1) = p^(J+1)*P0;
		Eq(J+1) = p/(1-p) - ((J+2)*p^(J+2))/(1-p^(J+2));
		Etq(J+1) = Eq(J+1)/(y*(1-Pb(J+1)))*1000;
end, J = [1 5 10 15];

figure, plot(J, Pb(J+1), 'ro', 'linewidth', 1.25, 'HandleVisibility', 'off'), hold on
plot(0:buf, Pb, 'b', 'linewidth', 1.5), ylim([-0.05 0.45]), grid on
xlabel('\fontsize{14}{Número de servidores (J)}'), ylabel('\fontsize{14}{Probabilidade de bloqueio (P_b)}')
title('\fontsize{14}{Probabilidade de bloqueio do sistema em função do número de servidores}')
print -dpng Pb.png

figure, plot(J, Eq(J+1), 'ro', 'linewidth', 1.25, 'HandleVisibility', 'off'), hold on
plot(0:buf, Eq, 'k', 'linewidth', 1.5), ylim([0.5 4.2]), grid on
xlabel('\fontsize{14}{Número de servidores (J)}'), ylabel('\fontsize{14}{Número médio de elementos no sistema E[q]}')
title('\fontsize{14}{Número médio de elementos no sistema em função do número de servidores}')
print -dpng Eq.png

figure, plot(J, Etq(J+1), 'ro', 'linewidth', 1.25, 'HandleVisibility', 'off'), hold on
plot(0:buf, Etq, 'm', 'linewidth', 1.5), ylim([3 22]), grid on
xlabel('\fontsize{14}{Número de servidores (J)}'), ylabel('\fontsize{14}{Tempo médio de sistema E[t_q] [ms]}')
title('\fontsize{14}{Tempo médio de sistema em função do número de servidores}')
print -dpng Etq.png, fprintf('\n')

