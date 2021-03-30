donnees;

n_tests = 100;

% Estimation de la droite de regression par le maximum de vraisemblance :
[a_DYX_1,b_DYX_1] = estimation_1(x_donnees_bruitees,y_donnees_bruitees,n_tests);

% Affichage de la droite de regression estimee par le maximum de vraisemblance :
if abs(a_DYX_1)<1
	x_DYX_1 = x_D;
	y_DYX_1 = a_DYX_1*x_DYX_1+b_DYX_1;
else
	y_DYX_1 = y_D;
	x_DYX_1 = (y_DYX_1-b_DYX_1)/a_DYX_1;
end
plot(x_DYX_1,y_DYX_1,'b','LineWidth',3);
axis(bornes);
lg = legend('~Droite', ...
	'~Donnees bruitees', ...
	'~$D_{YX}$ (maximum de vraisemblance)', ...
	'Location','Best');
set(lg,'Interpreter','Latex');

% Calcul et affichage de l'ecart angulaire :
theta_DYX_1 = atan2(b_DYX_1,-a_DYX_1*b_DYX_1);
EA_DYX_1 = abs(theta_DYX_1-theta_D);
fprintf('D_YX (maximum de vraisemblance) : ecart angulaire = %.2f degres\n',EA_DYX_1/pi*180);
