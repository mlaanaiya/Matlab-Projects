%% 1)Representation temporelle

%% Données de l'exercice 
    
A = 1;
f_0 = 1100;
Fe_1 = 10000;
Te_1 = 1/Fe_1;
Fe_3 = 1000;
Te_3 = 1/Fe_3;
N = 90;

%% Question1

% Générer un 90 echantillon d'un cosinus d'amplitude A = 1 V et de fréquence f_0
% = 1100 Hz et de fréquence d'echantillonnage Fe = 10000 Hz.

% implémentation des echantillons.
T_1 = 0:Te_1:(N-1)*Te_1;
X_1 = A*cos(2*pi*f_0*T_1);

%% Question2

% Tracer le cosinus généré dans la Question1 avec une échelle temporelle en
% secondes.

subplot(3,2,1);

plot(T_1,X_1,'b');

%Recadrer le graphe
axis([-0.001 0.01 -1.1 1.1]);

% Nommer les axes
title('Representation temporelle du signal X_1');
xlabel('Temps en secondes');
ylabel('Amplitude en Volts');

%% Question3

% Générer un 90 echantillon d'un cosinus d'amplitude A = 1 V et de fréquence f_0
% = 1100 Hz et de fréquence d'echantillonnage Fe3 = 1000 Hz.

T_3 = 0:Te_3:(N-1)*Te_3;
X_3 = A*cos(2*pi*f_0*T_3);


%% Question 4

% Tracer le cosinus généré dans la Question3 avec une échelle temporelle en
% secondes

subplot(3,2,2);
plot(T_3,X_3,'r')

%Recadrer le graphe
axis([-0.01 0.1 -1.1 1.1]);

% Nommer les axes
title('Representation temporelle du signal X_3');
xlabel('Temps en secondes');
ylabel('Amplitude en Volts');

% On ne peut pas retrouver à partir du graphe la frequence f_0 parce que la
% condition de Shannon-Nyquist n'est pas respectée  "Fe > 2*f_0" ici est fausse 

%% 2)Representation frequentielle

%% Question1

% En effectuant un échantillonage temporel, on périodise le signal selon la
% fréquence Fe d'où le calcul de la transformée entre 0 et Fe.

%% Question2

% Estimer la transformée de Fourier des signaux X_1 et X_3

TF_1 = fft(X_1);
TF_3 = fft(X_3);

% Transformer les complexes en réels
aTF_1 = abs(TF_1);
aTF_3 = abs(TF_3);

% Tracer en échelle log les transformées de Fourier du signal X_1 
subplot(3,2,3);
semilogy(aTF_1,'b')

% Nommer les axes
title('Representation frequentielle du signal X_1');
xlabel('Frequence en Hz');
ylabel('Gain en dB');

% Tracer en échelle log les transformées de Fourier du signal X_3
subplot(3,2,4);
semilogy(aTF_3,'r')

% Nommer les axes
title('Representation frequentielle du signal X_3');
xlabel('Frequence en Hz');
ylabel('Gain en dB');

%% Question3

% Estimer la transformée du signal X_1 en utilisant la méthode du
% Zero Padding

TF_1_zp = fft(X_1,2^9);

% Transformer les complexes en réels
aTF_1_zp = abs(TF_1_zp);

% Tracer en échelle log les transformées de Fourier du signal X_1 en zero
% padding

subplot(3,2,5);
semilogy(aTF_1_zp,'b')

% Nommer les axes
title('Representation frequentielle du signal X_1 en Zero Padding');
xlabel('Frequence en Hz');
ylabel('Gain en dB');
