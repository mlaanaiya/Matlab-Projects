%% 1)Generation du signal à filtrer

% Données du sujet :
f_1 = 1000;
f_2 = 3000;
A = 1;
Fe = 10000;
N = 100;
Te = 1/Fe;
f_c = 2000;

%% Question1

% Générer N = 100 échantillons d'une somme de deux cosinus d'amplitude 1 V, et
% de frequences f_1 = 1000 Hz, et f_2 = 3000 Hz, échantillonés à Fe =
% 10000 Hz

% Générons le tableau des indices temporels T_12, et les signaux X_1, X_2

T_12 = 0:Te:(N-1)*Te;

X_1 = A*cos(2*pi*f_1*T_12);
X_2 = A*cos(2*pi*f_2*T_12);

% Ainsi le signal somme X_1 et X_2, X_12 est :

X_12 = X_1 + X_2;

%% Question2

% Tracer le signal obtenu en échelle temporelle
% Tracer les deux autres signaux au dessus en commentaires pour la suite

subplot(3,2,1);

%plot(T_12,X_1,'b')
%hold on
%plot(T_12,X_2,'r')
%hold on

plot(T_12,X_12,'b')

axis([-0.0005 0.0105 -2.2 2.2]);

%Nommer les axes
title('Représentation temporelle du signal X_1_2');
xlabel('Temps en secondes');
ylabel('Amplitude en Volts');

%% Question3

% Tracer une representation frequentielle du signal obtenu
TF_12 = fft(X_12);
aTF_12 = abs(TF_12);

subplot(3,2,2);

semilogy(aTF_12,'r')

%Nommer les axes
title('Représentation fréquentielle du signal X_1_2');
xlabel('Fréquence en Hz');
ylabel('Gain en dB');


%% 2)Synthèse du filtre passe-bas

%%  Question1

% La réponse impulsionnelle du filtre passe-bas qui ne laisse passer que le
% cosinus de fréquence f_1 = 1000 Hz est h = int Porte_f_c(f) * exp(j2pikf) * df de -f_0/2 à +f_0/2
% avec une fréquence de coupure f_c = 2000 Hz

I_5 = -5:5; % Intervalle discret
I_30 = -30:30;
h_5 = 2*f_c/Fe*sinc(2*I_5*(f_c/Fe)); % Réponse impulsionnelle ordre 11
h_30 = 2*f_c/Fe*sinc(2*I_30*(f_c/Fe)); % Réponse impulsionnelle ordre 61

%% Question2

% Tracer la réponse aux filtres d'ordre 11 et 61
subplot(3,2,3);

plot(I_5,h_5,'g')

%Nommer les axes
title("Représentation de la réponse impulsionnelle en ordre 11");

subplot(3,2,4);

plot(I_30,h_30,'g')

%Nommer les axes
title("Représentation de la réponse impulsionnelle en ordre 61");


%% Question3

% Estimez la réponse de fréquences des filtres d'ordre 11 et 61
Y_12_11 = filter(h_5,1,X_12);
Y_12_61 = filter(h_30,1,X_12);

TF_11 = fft(Y_12_11);

TF_61 = fft(Y_12_61);

aTF_11 = abs(TF_11);

aTF_61 = abs(TF_61);

subplot(3,2,5);

semilogy(aTF_11)

title("Représentation de la réponse filtre synthétisé d'ordre 11");

subplot(3,2,6);

semilogy(aTF_61)
title("Représentation de la réponse filtre synthétisé d'ordre 61");


