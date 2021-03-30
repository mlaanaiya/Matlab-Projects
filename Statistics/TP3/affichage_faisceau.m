function affichage_faisceau(rho,theta,limites_affichages,couleur)

% Limites de l'affichage du faisceau :
x_min = limites_affichages(1);
x_max = limites_affichages(2);
y_min = limites_affichages(3);
y_max = limites_affichages(4);

% Trace des droites en boucle :
for k = 1:length(rho)

	% Parametres de la droite D_k :
	rho_k = rho(k);
	theta_k = theta(k);
	sin_theta_k = sin(theta_k);
	cos_theta_k = cos(theta_k);

	% Coordonnees des extremites du segment a afficher :
	x_extremites = zeros(1,2);
	y_extremites = zeros(1,2);

	% Cas ou D_k coupe le bord gauche :
	i = 1;
	x = x_min;
	y = (rho_k-x*cos_theta_k)/sin_theta_k;
	if y>=y_min & y<=y_max
		x_extremites(i) = x;
		y_extremites(i) = y;
		i = i+1;
	end

	% Cas ou D_k coupe le bord droit :
	x = x_max;
	y = (rho_k-x*cos_theta_k)/sin_theta_k;
	if y>=y_min & y<=y_max
		x_extremites(i) = x;
		y_extremites(i) = y;
		i = i+1;
	end

	% Cas ou D_k coupe le bord superieur :
	y = y_min;
	x = (rho_k-y*sin_theta_k)/cos_theta_k;
	if x>=x_min & x<=x_max
		x_extremites(i) = x;
		y_extremites(i) = y;
		i = i+1;
	end

	% Cas ou D_k coupe le bord inferieur :
	y = y_max;
	x = (rho_k-y*sin_theta_k)/cos_theta_k;
	if x>=x_min & x<=x_max
		x_extremites(i) = x;
		y_extremites(i) = y;
	end

	% Affichage de la droite D_k :
	plot(x_extremites,y_extremites,couleur,'LineWidth',2);
end
