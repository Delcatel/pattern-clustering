% G�n�ration des �chantillons
% ENPC 2017
clear all
close all
clc
N=1000;

% G�n�ration des �chantillons
symboles=sign(randn(1,N))+j*sign(randn(1,N));

% Gn�ration du coefficient H du canal
H=randn+j*randn;

% G�n�ration du bruit
bruit=randn(1,N)+j*randn(1,N);

% R�glage du rapport signal sur bruit (SNR : signal to noise ratio)
% On peut �viter de lire ce paragraphe en premi�re lecture du code
Ps=abs(H)^2*(symboles*symboles'/length(symboles)); % Puissance signal
Pb=(bruit*bruit'/length(bruit)); % Puissance bruit
SNR=10; % rapport puissance signal sur puissance bruit cible (en lin�aire)
alpha=sqrt((Ps/Pb)/SNR); % coefficient multiplicatif � appliquer sur le bruit

% Signal re�u
r=H*symboles+alpha*bruit;

figure(1)
plot(r,'x');
grid on
xlabel('real');
ylabel('imaginary');

% Algorithme de Viterbi & Viterbi
H_estime=((-1/4)*sum(r.^4)/length(r)).^(1/4);
% Lev�e de l'ambiguit� potentielle de exp(j*pi/2)
hypotheses=[H_estime;j*H_estime;-H_estime;-j*H_estime];
ecart=H-hypotheses;
[valeur,indice]=min(abs(ecart));
H_estime=hypotheses(indice);

disp('-----------------------');
H
H_estime
disp('-----------------------');
 




