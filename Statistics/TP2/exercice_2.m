exercice_1;

% Estimation de la droite de regression par resolution du systeme AX = B :
[a_DYX_2,b_DYX_2] = estimation_2(x_donnees_bruitees,y_donnees_bruitees);

% Affichage de la droite de regression estimee par resolution du systeme AX = B :
if abs(a_DYX_2)<1
	x_DYX_2 = x_D;
	y_DYX_2 = a_DYX_2*x_DYX_2+b_DYX_2;
else
	y_DYX_2 = y_D;
	x_DYX_2 = (y_DYX_2-b_DYX_2)/a_DYX_2;
end
plot(x_DYX_2,y_DYX_2,'g','LineWidth',3);
lg = legend('~Droite', ...
	'~Donnees bruitees', ...
	'~$D_{YX}$ (maximum de vraisemblance)', ...
	'~$D_{YX}$ (moindres carres)', ...
	'Location','Best');
set(lg,'Interpreter','Latex');

% Calcul et affichage de l'ecart angulaire :
theta_DYX_2 = atan2(b_DYX_2,-a_DYX_2*b_DYX_2);
EA_DYX_2 = abs(theta_DYX_2-theta_D);
fprintf('D_YX (moindres carres) : ecart angulaire = %.2f degres\n',EA_DYX_2/pi*180);
