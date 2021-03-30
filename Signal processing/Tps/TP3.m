%% Filtrage

f_1 = 6000;
A=1; 
Fe = 48000;
Te = 1/Fe;
Ns = 300;
Ts = Ns*Te;
f_c = 6000;

Var_0 = rand*2*pi;
Var_1 = rand*2*pi;
Donnee = randi([0,100],1);
Fe = 48*10^3;
ordre = 101;
fc = 4000;
Te = 1/Fe;
Debit = 300;
Ns = 1/(Debit*Te);
Nombre_Tot = Ns*length(Donnee);
t = 0:Te:(Nombre_Tot-1)*Te;
Taille_echantillon = ones(1,Ns);
NRZ = kron(Donnee , Taille_echantillon);
T_1 = 0:Te:(Nombre_Tot-1)*Te;
coef = -(ordre-1)/2:1:(ordre-1)/2;
X_1 = (1 - NRZ).*cos(2*pi*f_1*T_1 + Var_0 ) + NRZ.* cos(2*pi*f_1*T_1 + Var_1 );

I = -Nombre_Tot*Te:Nombre_Tot*Te; % Intervalle discret

h_1 = 2*f_c/Fe * sinc(2*coef*(f_c/Fe)); 

Y_1 = filter(h_1,1,X_1);

TF_1 = fft(Y_1);

aTF_1 = abs(TF_1);

semilogy(aTF_1)


title("Représentation de la réponse filtre");


%%filtre passe-haut
h_2 = -h_1;
h_2(floor(length(coef)/2)+1) = 1 + h_2(floor(length(coef)/2));
y_haut = filter(h_2,1,X_1);
TF_2 = fft(y_haut);
aTF_2 = abs(TF_2);
semilogy(aTF_2);

title("filtre_haut");
%DSP



