donnees;

n_tests = 500;

% Estimation du rayon et de la position du centre :
[C_estime,R_estime] = estimation_2(x_donnees_bruitees,y_donnees_bruitees,n_tests);

% Affichage du cercle estime :
n_points_cercle = 100;
theta_cercle = 2*pi/n_points_cercle:2*pi/n_points_cercle:2*pi;
x_cercle_estime = C_estime(1)+R_estime*cos(theta_cercle);
y_cercle_estime = C_estime(2)+R_estime*sin(theta_cercle);
plot(x_cercle_estime([1:end 1]),y_cercle_estime([1:end 1]),'b','LineWidth',3);
lg = legend(' Cercle initial', ...
		' Donnees bruitees', ...
		' Cercle estime', ...
		'Location','Best');
