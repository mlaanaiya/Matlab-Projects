function [theta_Dorth_2,rho_Dorth_2] = estimation_4(x_donnees_bruitees, y_donnees_bruitees)
x_G = mean(x_donnees_bruitees);
y_G = mean(y_donnees_bruitees);
x_centre = x_donnees_bruitees - x_G;
y_centre = y_donnees_bruitees - y_G;
C = [x_centre; y_centre];
C = C';
%calcul des vecteurs propres
[V, D] = eig(C'*C);
VP = diag(D);
[vp_min, I_min] = min(VP);
Y = V(:, I_min);%vecteur propre associé à vp_min
theta_Dorth_2 = atan2( Y(2), Y(1) ); %arctan(cos(theta)/sin(theta))
rho_Dorth_2 = x_G*cos(theta_Dorth_2) + y_G*sin(theta_Dorth_2);
end