% Génération des échantillons
% ENPC 2017
clear all
close all
clc
N=1000;

% Génération des échantillons
symboles=sign(randn(1,N))+j*sign(randn(1,N));

% Gnération du coefficient H du canal
H=randn+j*randn;

% Génération du bruit
bruit=randn(1,N)+j*randn(1,N);

% Règlage du rapport signal sur bruit (SNR : signal to noise ratio)
% On peut éviter de lire ce paragraphe en première lecture du code
Ps=abs(H)^2*(symboles*symboles'/length(symboles)); % Puissance signal
Pb=(bruit*bruit'/length(bruit)); % Puissance bruit
SNR=10; % rapport puissance signal sur puissance bruit cible (en linéaire)
alpha=sqrt((Ps/Pb)/SNR); % coefficient multiplicatif à appliquer sur le bruit

% Signal reçu
r=H*symboles+alpha*bruit;

figure(1)
plot(r,'x');
grid on
xlabel('real');
ylabel('imaginary');

% Algorithme de Viterbi & Viterbi
H_estime=((-1/4)*sum(r.^4)/length(r)).^(1/4);
% Levée de l'ambiguité potentielle de exp(j*pi/2)
hypotheses=[H_estime;j*H_estime;-H_estime;-j*H_estime];
ecart=H-hypotheses;
[valeur,indice]=min(abs(ecart));
H_estime=hypotheses(indice);

disp('-----------------------');
H
H_estime
disp('-----------------------');
 




