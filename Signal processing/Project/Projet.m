%load DonneesBinome1.mat;
%load DonneesBinome2.mat;
%load DonneesBinome3.mat;
%load DonneesBinome4.mat;
%load DonneesBinome5.mat;
load DonneesBinome6.mat;
%% NRZ %%
% Donnees = randi([0,1],1,100);
Nombre_bits = length(bits);
F1 = 980;%1080-100
F0 = 1180;%1080+100
Fe = 48*10^3;
Te = 1/Fe;
Ts = 1/300;
Ns = floor(Ts*Fe);
Nombre_Tot = Ns * Nombre_bits;
t = 0 : Te : (Nombre_Tot-1)*Te;
Taille_echantillon = ones(1,Ns);
F_0 = 6*10^3;
F_1 = 2*10^3;
VA_0 = rand*2*pi;
VA_1 = rand*2*pi;

%% X %% 
NRZ = kron(bits , Taille_echantillon); 
t = 0:Te:(length(NRZ)-1)*Te;
frequence = 0:Fe/(length(NRZ)-1):Fe;
% plot(t,NRZ);
% axis([0 0.5 -1 2]);
% title("NRZ");
% xlabel("Tps(s)");
% ylabel("NRZ(t)");
Freq_nrz = 2^nextpow2(length(NRZ));
Tf_nrz = fft(NRZ, Freq_nrz);
Tf_nrz = fftshift(Tf_nrz);
Snrz = 1/Nombre_Tot * abs(Tf_nrz).^2;
x = (1-NRZ).*cos(2*pi*F_0*t + VA_0) + NRZ.*cos(2*pi*F_1*t + VA_1);
x_fsk = (1-NRZ).*cos(2*pi*F0*t + VA_0) + NRZ.*cos(2*pi*F1*t + VA_1);
%figure,semilogy(linspace(-Fe/2,Fe/2,Freq_nrz),Snrz) 

%%Implantation du bruit
SNR = 50;
Px = mean(abs(x).^2);
Pb = Px/(10^(SNR/10));
somme = sqrt(Pb);
bruit = (somme)*randn(1,length(x));
x_sans_bruit = x; %signale sans bruit
x_bruit = x + bruit; %signal avec bruit

%figure,plot(t,x),
%axis([0 1.5*Nombre_Tot*Te -2 2]);
Freq_X = 2^nextpow2(length(x_bruit));
Tf_X = fft(NRZ, Freq_X);
Tf_X = fftshift(Tf_X);
Sx = 1/Nombre_Tot * abs(Tf_X).^2;
%figure,semilogy(linspace(-Fe/2,Fe/2,Freq_X),Sx)

%***************************************************************************%

%% Filtrage
%%3.3.1 Filtre passe bas 

ordre=100;
coeff= -ordre:1:ordre; % Intervalle discret
fc=4000;
h_1 = (2*fc/Fe)*sinc(2*(fc/Fe)*coeff);

%***************************************************************************%

%%3.3.2 Filtre passe haut

Feh = 10000;
Teh = 1/Feh; 
h_2 = -h_1;
h_2(floor(length(coeff)/2)+1) = 1 + h_2(floor(length(coeff)/2)+1);

%***************************************************************************%
%%Retard

retard = zeros(1,(length(h_1)-1)/2);
signal_retard_bruit = [retard x_bruit];
signal_retard_sans_bruit = [retard x_sans_bruit];
   
%%3.3.3 Filtrage
%Signal bruité
passe_bas = h_1;
passe_haut = h_2;
y_bas = filter(passe_bas,1,signal_retard_bruit);
y_haut = filter(passe_haut,1,signal_retard_bruit);
y_bas = y_bas(1:Ns*length(bits));
y_haut = y_haut(1:Ns*length(bits));

%Signal non bruité
y_bas_1 = filter(passe_bas,1,signal_retard_sans_bruit);
y_haut_1 = filter(passe_haut,1,signal_retard_sans_bruit);
y_bas_1 = y_bas_1(1:Ns*length(bits));
y_haut_1 = y_haut_1(1:Ns*length(bits));



%***************************************************************************%

% %%Tracés à réaliser 
% plot(t,x);
% axis([0 0.01 -2 2]);
% title("Le signal");
% figure;
% plot(passe_bas);
% title("réponse impulsionnelle du passe-bas");
% figure;
% freq_bas = fft(passe_bas, length(Sx));
% plot(abs(freq_bas).*2);
% title("réponse frenquentielle du passe-bas");
% figure;
% plot(Sx,'g');
% hold on
% plot(abs(freq_bas).*2,'r');
% hold off
% figure;
% Syb = (abs(fft(y_bas)).^2)/length(y_bas);
% plot(y_bas);
% title("Signal filtré par passe-bas");
% hold on
% plot(Syb);
% title("Densité spectrale du signal filtré en basses fréquences");
% hold off
% figure;
% 
% plot(passe_haut);
% title("réponse impulsionnelle du passe-haut");
% figure;
% freq_haut = fft(passe_haut, length(Sx));
% plot(abs(freq_haut).*2);
% title("réponse frenquentielle du passe-haut");
% figure;
% plot(Sx,'g');
% title ("densité spectrale du signal");
% hold on
% plot(abs(freq_haut).*2,'r');
% title ("densité spectrale du signal");
% hold off
% figure;
% Syh = (abs(fft(y_haut)).^2)/length(y_haut);
% plot(y_haut);
% title("signal filtré en hautes fréquences");
% hold on
% plot(Syh);
% title("densité spectrale du signal filtré par passe-haut");
% hold off

%***************************************************************************%
%%Détection d´énergie

%Taux d´erreur : voir après la reconstitution de l´image

%Avec bruit

%Pour le passe-bas
X_bas = reshape(y_bas,Ns,length(bits)); % Ns*length(bits) = length(bits)
E_bas = zeros(length(X_bas),1);
E_bas = sum(X_bas.^2);
K_bas = (max(E_bas)-min(E_bas))/2; %Seuil
%Reconstitution en bas
E_bas(E_bas<=K_bas) = 0;
E_bas(E_bas>K_bas) = 1;


%Pour le passe-haut
X_haut = reshape(y_haut,Ns,length(bits));
E_haut = zeros(length(X_haut),1);
E_haut = sum(X_haut.^2);
K_haut = (max(E_haut)-min(E_haut))/2;%Seuil
%Reconstitution en haut
E_haut(E_haut<=K_haut) = 1;
E_haut(E_haut>K_haut) = 0;

%Reconstitution image en bas
pcode reconstitution_image;
reconstitution_image(E_bas);
which reconstitution_image;
title("On passe aux images :D");
title("image reconstituée en passe-bas avec le signal bruité");

%Reconstitution image en haut
pcode reconstitution_image;
reconstitution_image(E_haut);
which reconstitution_image;
title("image reconstituée en passe-haut avec le signal bruité");


%Sans bruit

%Pour le passe_bas
X_bas_1 = reshape(y_bas_1,Ns,length(bits)); 
E_bas_1 = zeros(length(X_bas_1),1);
E_bas_1 = sum(X_bas_1.^2);
K_bas_1 = (max(E_bas_1)-min(E_bas_1))/2;%Seuil
%Reconstitution en bas
E_bas_1(E_bas_1<=K_bas_1) = 0;
E_bas_1(E_bas_1>K_bas_1) = 1;


%Pour le passe_haut
X_haut_1 = reshape(y_haut_1,Ns,length(bits));
E_haut_1 = zeros(length(X_haut_1),1);
E_haut_1 = sum(X_haut_1.^2);
K_haut_1 = (max(E_haut_1)-min(E_haut_1))/2;%Seuil
%Reconstitution en haut
E_haut_1(E_haut_1<=K_haut_1) = 1;
E_haut_1(E_haut_1>K_haut_1) = 0;

%Reconstitution image en bas
pcode reconstitution_image;
reconstitution_image(E_bas_1);
which reconstitution_image;
title("image reconstituée en passe-bas avec le signal non bruité");

%Reconstitution image en haut
pcode reconstitution_image;
reconstitution_image(E_haut_1);
which reconstitution_image;
title("image reconstituée en passe-haut avec le signal non bruité");

%%Calcul du taux d´erreur

%Signal avec bruit
tauxderreur = sum(abs(E_bas - bits))/length(bits);
%Signal sans bruit
taux_d_erreur_1 = sum(abs(E_bas_1 - bits))/length(bits);

%***************************************************************************%
                               %FSK
%***************************************************************************%
%% Démodulateur FSK
x_reshape = reshape(x_fsk,Ns,length(bits));
t_reshape = reshape(t,Ns,length(bits));
x0 = x_reshape.*cos(2*pi*F0*t_reshape + VA_0); %Le produit du signal avec le cos en F0 = 1180
x1 = x_reshape.*cos(2*pi*F1*t_reshape + VA_1); %Le produit du signal avec le cos en F1 = 980
H = trapz(Ts,x1) - trapz(Ts,x0);%difference des deux intégrale
bits_reconstitue = H > 0;
pcode reconstitution_image;
reconstitution_image(bits_reconstitue);
which reconstitution_image;
title("image reconstitué en FSK");

%***************************************************************************%

%%4.2/Erreur de synchronisation
%1/
x0_phase = x_reshape.*cos(2*pi*F0*t_reshape + VA_0 + 2*pi*rand); %Déphasage aléatoire ajouté
x1_phase = x_reshape.*cos(2*pi*F1*t_reshape + VA_1 + 2*pi*rand); %Déphasage aléatoire ajouté
H_1 = trapz(Ts,x1_phase) - trapz(Ts,x0_phase);%difference des deux intégrale
bits_reconstitue_1 = H_1 < 0;
pcode reconstitution_image;
reconstitution_image(bits_reconstitue_1);
which reconstitution_image;
title("image avec déphasage");
%Remarque : Image affichée non nette

%2/
x0_h = x_reshape.*cos(2*pi*F0*t_reshape);
x0_sin = x_reshape.*sin(2*pi*F0*t_reshape);
H0_h = (trapz(Ts,x0_h)).^2 + (trapz(Ts,x0_sin)).^2;
x1_h = x_reshape.*cos(2*pi*F1*t_reshape);
x1_sin = x_reshape.*sin(2*pi*F1*t_reshape);
H1_h = (trapz(Ts,x1_h)).^2 + (trapz(Ts,x1_sin)).^2;
H_tot = H1_h - H0_h;
bits_restitue_tot = H_tot > 0;
pcode reconstitution_image;
reconstitution_image(bits_restitue_tot);
which reconstitution_image;
title("image reconstitué en FSK même avec le déphasage");

% %
% Y_1 = filter(h_1,1,x);
% TF_1 = fft(Y_1);
% ShTD =fftshift(TF_1);
% aTF_1 = abs(ShTD);
% 
% 
% figure, semilogy(linspace(-Fe/2,Fe/2,length(h_1)),abs(fftshift(fft(h_1)))),
% title("Réponse en fréquence du filtre")
% figure,plot(h_1)
% title("Réponse Impulsionnelle");
% 
% figure, semilogy(linspace(-Fe/2,Fe/2,length(Y_1)),aTF_1),
% hold on, semilogy(linspace(-Fe/2,Fe/2,length(x)),abs(fftshift(fft(x)))),
% title("Densité spectrale de puissance de X, et réponse filtre implanté");
% 
% 
% figure, plot(Y_1),
% title("Signal en sortie du filtre");
% figure, semilogy(linspace(-Fe/2,Fe/2,length(Y_1)),(abs(fftshift(fft(Y_1))).^2)/length(Y_1)),
% title("Densité spectrale en puissance du signal en sortie");
% 
% figure,semilogy(length(aTF_2),aTF_2);
% %4.1 Plan qui prend beaucoup de temps

% Démodulateur FSK
% demodulation = [];
% t1 = 0:Ts/Ns:Ts;
% Fc_0 = 1180;%Fc+100 = 1080 + 100
% Fc_1 = 980; %Fc-100 = 1080 -100 
% cos0 = cos(2*pi*Fc_0*t1 + VA_0); % pour un binaire = 0
% cos1 = cos(2*pi*Fc_1*t1 + VA_1); % pour un binaire = 1
% mc0 = cos0*reshape(x,length(x)/length(cos0),length(cos0));
% mc1 = cos1*reshape(x,length(x)/length(cos0),length(cos0));
% integ0 = trapz(Ts,mc0);
% integ1 = trapz(Ts,mc1);
% u = integ1 - integ0;
% for i = 1:length(u)
%     if u(i) > 0
%         a = 1;
%     else
%         a = 0;
%     end
%     demodulation = [demodulation a];
% end
% x_reshape = reshape(x,Ns,length(bits));
% % 
% demodulation = [];
% x_reshape = reshape(x,Ns,length(bits));
% for n = 1:length(bits)
%     mc0 = cos0.*x(n); %produit
%     mc1 = cos1.*x(n); %produit
%     integ0 = trapz(Ts,mc0); %Integrale
%     integ1 = trapz(Ts,mc1); %Integrale
%     u = integ0 - integ1;
%     if (u > 0)
%         a=0;
%     elseif (u < 0)
%         a=1;
%     end
%     demodulation = [demodulation a];
% end
% % Démodulation en numérique
% numerique = [];
% for n = 1:length(demodulation)
%     if (demodulation(n) == 1)
%         digit = ones(1,Ns);
%     else (demodulation(n) == 0);
%         digit = zeros(1,Ns);
%     end
%     numerique = [numerique digit];
% end
% Temps_Periode = Ts/Ns:Ts/Ns:Ts*length(bits);
% subplot(3,1,2);
% plot(Temps_Periode,demodulation);
% axis([0 1 -0.5 1.5]);
% xlabel("Temps");
% ylabel("Amplitude");
% title("Demodulation FSK");
%%Retour vers NRZ
% demodulation = transpose(bits_restitue).*ones(1,Ns);
% Temps_periode = Ts/Ns:Ts/Ns:Ts*length(bits);
% demod = reshape(demodulation,1,length(NRZ));
%xlabel("Temps");
%ylabel("Amplitude");
%title("Demodulation FSK");
% subplot(3,1,2);
% plot(Temps_periode,demod);
% axis([0 1 -0.5 1.5]);