donnees;

n_tests = 500;

% Estimation de la position du centre :
[C_estime,R_moyen] = estimation_1(x_donnees_bruitees,y_donnees_bruitees,n_tests);

% Affichage du cercle estime :
n_points_cercle = 100;
theta_cercle = 2*pi/n_points_cercle:2*pi/n_points_cercle:2*pi;
x_cercle_estime = C_estime(1)+R_moyen*cos(theta_cercle);
y_cercle_estime = C_estime(2)+R_moyen*sin(theta_cercle);
plot(x_cercle_estime([1:end 1]),y_cercle_estime([1:end 1]),'b','LineWidth',3);
lg = legend(' Cercle initial', ...
		' Donnees bruitees', ...
		' Cercle estime', ...
		'Location','Best');
