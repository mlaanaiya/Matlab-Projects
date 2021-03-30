donnees;

n_tests = 100;

% Estimation de la droite de regression par le maximum de vraisemblance :
[theta_Dorth_1,rho_Dorth_1] = estimation_3(x_donnees_bruitees,y_donnees_bruitees,n_tests);

% Affichage de la droite de regression estimee par le maximum de vraisemblance :
cos_theta_Dorth_1 = cos(theta_Dorth_1);
sin_theta_Dorth_1 = sin(theta_Dorth_1);
if abs(cos_theta_Dorth_1)<abs(sin_theta_Dorth_1)
	x_Dorth_1 = x_D;
	y_Dorth_1 = (rho_Dorth_1-x_Dorth_1*cos_theta_Dorth_1)/sin_theta_Dorth_1;
else
	y_Dorth_1 = y_D;
	x_Dorth_1 = (rho_Dorth_1-y_Dorth_1*sin_theta_Dorth_1)/cos_theta_Dorth_1;
end
plot(x_Dorth_1,y_Dorth_1,'b','LineWidth',3);
axis(bornes);
lg = legend('~Droite', ...
	'~Donnees bruitees', ...
	'~$D_\perp$ (maximum de vraisemblance)', ...
	'Location','Best');
set(lg,'Interpreter','Latex');

% Calcul et affichage de l'ecart angulaire :
EA_Dorth_1 = abs(theta_Dorth_1-theta_D);
fprintf('D_perp (maximum de vraisemblance) : ecart angulaire = %.2f degres\n',EA_Dorth_1/pi*180);
