
% Simulação da Pout de sistema Cognitivo com Improper Gaussian Signaling
clc, clear all, close all, fprintf('\n')
warning('off', 'all');
pkg load signal

% Modelo do Sistema:
N0 = 1; 			% Variância do ruído.
N = 1e6; 		% Amostras Monte Carlo.
alpha = 4; 	% Coeficiente de Perda de Percurso.

% Distâncias entre os nós para a Configuração I e
% médias dos canais entre links fonte secundária e usuário primário:
xA = 1; yA = 0; xB = 1.5; yB = 0; xS = 1; yS = 1; xD = 0.5; yD = 0.5;
lambdaAD = sqrt((xA - xD)^2 + (yA - yD)^2)^(-alpha);
lambdaAB = sqrt((xA - xB)^2 + (yA - yB)^2)^(-alpha);
lambdaSD = sqrt((xS - xD)^2 + (yS - yD)^2)^(-alpha);
lambdaSB = sqrt((xS - xB)^2 + (yS - yB)^2)^(-alpha);

% Sorteio de canais usando Rayleigh fading:
gad = (max(0, sqrt(lambdaAD).*randraw('nakagami', 1, [1 N]))).^2;
hab = (max(0, sqrt(lambdaAB).*randraw('nakagami', 1, [1 N]))).^2;
hsd = (max(0, sqrt(lambdaSD).*randraw('nakagami', 1, [1 N]))).^2;
gsb = (max(0, sqrt(lambdaSB).*randraw('nakagami', 1, [1 N]))).^2;

% SOP na rede secundária e Outage na primária em função de C_x
Ps = db2pow(10); 			% Potência da fonte primária
Ra = 1; Rs = 1; Pa_max = [5 20];
rho_values = 0:0.05:1; style = {'b--', 'b:'};

Cx = 0;
A = -mean(hsd).*Ps + N0.*(4.^Rs - 1);
B = (Cx.^2 - 1).*mean(gad).*(4.^Rs - 1);
C = (mean(hsd).^2).*(Ps.^2).*(4.^Rs)+(Cx.^2).*(4.^Rs - 1).*(-(N0 + mean(hsd).*Ps).^2 + (4.^Rs.*N0.^2));
D = ((Cx.^2 - 1).^2).*(mean(gad).^2).*(4.^Rs - 1).^2;
Pa_dag = (A./B) + sqrt(C./D); % Equação 10

for p = 1:numel(Pa_max)
		SOP = {};
		for rho = 0:0.05:1
				lambdaAE = sqrt((rho - xA)^2 + (rho - yA)^2)^(-alpha);
				hae = (max(0, sqrt(lambdaAE).*randraw('nakagami', 1, [1 N]))).^2;

				lambdaSE = sqrt((rho - xS)^2 + (rho - yS)^2)^(-alpha);
				gse = (max(0, sqrt(lambdaSE).*randraw('nakagami', 1, [1 N]))).^2;

				Pa = min(Pa_dag, db2pow(Pa_max(p))); % Equação 12

				% Informação mútua secundários e SOP na rede secundária
				PGS_ab = log2(1 + ((Pa.*hab)./(Ps.*mean(gsb)+N0)));
				IGS_ab = 0.5.*log2(1 - ((Pa.*hab.*Cx)./(Pa.*(hab) + Ps.*mean(gsb)+N0)).^2);
				PGS_ae = log2(1 + ((Pa.*mean(hae))./(Ps.*mean(gse)+N0)));
				IGS_ae =0.5.*log2((1 - ((Pa.*mean(hae).*Cx)./(Pa.*mean(hae) + Ps.*mean(gse)+N0)).^2));
				SOP{end + 1} = sum((PGS_ab + IGS_ab - PGS_ae - IGS_ae) < Ra)./N; % Equação 13
		end

		semilogy(rho_values, cell2mat(SOP), style(p), 'linewidth', 1.5), hold on
end

style = {'r-', 'r-.'}; Cx = 1;
B = (Cx.^2 - 1).*mean(gad).*(4.^Rs - 1);
C = (mean(hsd).^2).*(Ps.^2).*(4.^Rs)+(Cx.^2).*(4.^Rs - 1).*(-(N0 + mean(hsd).*Ps).^2 + (4.^Rs.*N0.^2));
D = ((Cx.^2 - 1).^2).*(mean(gad).^2).*(4.^Rs - 1).^2;
Pa_dag = (A./B) + sqrt(C./D); % Equação 10

for p = 1:numel(Pa_max)
		SOP = {};
		for rho = 0:0.05:1
				lambdaAE = sqrt((rho - xA)^2 + (rho - yA)^2)^(-alpha);
				hae = (max(0, sqrt(lambdaAE).*randraw('nakagami', 1, [1 N]))).^2;

				lambdaSE = sqrt((rho - xS)^2 + (rho - yS)^2)^(-alpha);
				gse = (max(0, sqrt(lambdaSE).*randraw('nakagami', 1, [1 N]))).^2;

				Pa = min(Pa_dag, db2pow(Pa_max(p))); % Equação 12

				% Informação mútua secundários e SOP na rede secundária
				PGS_ab = log2(1 + ((Pa.*hab)./(Ps.*mean(gsb)+N0)));
				IGS_ab = 0.5.*log2(1 - ((Pa.*hab.*Cx)./(Pa.*(hab) + Ps.*mean(gsb)+N0)).^2);
				PGS_ae = log2(1 + ((Pa.*mean(hae))./(Ps.*mean(gse)+N0)));
				IGS_ae =0.5.*log2((1 - ((Pa.*mean(hae).*Cx)./(Pa.*mean(hae) + Ps.*mean(gse)+N0)).^2));
				SOP{end + 1} = sum((PGS_ab + IGS_ab - PGS_ae - IGS_ae) < Ra)./N; % Equação 13
		end

		semilogy(rho_values, cell2mat(SOP), style(p), 'linewidth', 1.5), hold on
		legend('SOP PGS P_{a_{max}} = 5 dB', 'SOP PGS P_{a_{max}} = 20 dB', ['SOP' ...
				' IGS P_{a_{max}} = 5 dB'], 'SOP IGS P_{a_{max}} = 20 dB', 'Location', 'southwest')
		xlabel('Eve coordinates on the Cartesian plan - \rho'), xticks(0:0.1:1), ylabel('Secrecy Outage Probability - Os')
end, grid on, axis([0 1 10^-3 1])

title({'Physical Layer Security in Cognitive Radio', 'Networks Using Improper Gaussian Signaling'})
print -dpng EveCoordinates.png

