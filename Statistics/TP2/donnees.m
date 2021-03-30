clear;
close all;

taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);

% Fenetre d'affichage :
figure('Name','Points situes au voisinage d''une droite', ...
	'Position',[0.4*L,0,0.6*L,0.8*H]);
axis equal;
hold on;
set(gca,'FontSize',20);
hx = xlabel('$x$','FontSize',30);
set(hx,'Interpreter','Latex');
hy = ylabel('$y$','FontSize',30);
set(hy,'Interpreter','Latex');

% Bornes d'affichage des donnees centrees en (0,0) :
taille = 20;
bornes = [-taille taille -taille taille];

% Parametres de la droite :
theta_D = 2*pi*(rand-0.5);
cos_theta_D = cos(theta_D);
sin_theta_D = sin(theta_D);
marge = 5;
rho_D = (taille-marge)*rand;

% Affichage de la droite :
droite_horizontale = abs(cos_theta_D)<abs(sin_theta_D);
pas = 0.01;
if droite_horizontale
	x_D = -taille:pas:taille;
	y_D = (rho_D-cos_theta_D*x_D)/sin_theta_D;
else
	y_D = -taille:pas:taille;
	x_D = (rho_D-sin_theta_D*y_D)/cos_theta_D;
end
plot(x_D,y_D,'r-','LineWidth',3);

% Donnees non bruitees :
n = 50;
if droite_horizontale
	x_donnees = taille*(2*rand(1,n)-1);
	y_donnees = (rho_D-cos_theta_D*x_donnees)/sin_theta_D;
else
	y_donnees = taille*(2*rand(1,n)-1);
	x_donnees = (rho_D-sin_theta_D*y_donnees)/cos_theta_D;
end

% Donnees bruitees :
sigma = 2;
x_donnees_bruitees = x_donnees+sigma*randn(1,n);
y_donnees_bruitees = y_donnees+sigma*randn(1,n);
indices_visibles = find(x_donnees_bruitees>-taille & ...
			x_donnees_bruitees<taille & ...
			y_donnees_bruitees>-taille & ...
			y_donnees_bruitees<taille);
x_donnees_bruitees = x_donnees_bruitees(indices_visibles);
y_donnees_bruitees = y_donnees_bruitees(indices_visibles);
n = length(indices_visibles);

% Affichage des donnees bruitees :
plot(x_donnees_bruitees,y_donnees_bruitees,'k+','MarkerSize',10,'LineWidth',2);
axis(bornes);
lg = legend(' Droite', ...
	' Donnees bruitees', ...
	'Location','Best');
grid on;
