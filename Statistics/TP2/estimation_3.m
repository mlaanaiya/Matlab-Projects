function [theta_Dorth_1,rho_Dorth_1] = estimation_3(x_donnees_bruitees,y_donnees_bruitees,n_tests)
 theta_tests = pi*rand(n_tests,1);
 x_G = mean(x_donnees_bruitees);
 y_G = mean(y_donnees_bruitees);
 x_centre=x_donnees_bruitees-x_G;
 y_centre=y_donnees_bruitees-y_G;
 sigma=sum((x_centre.*cos(theta_tests)+y_centre.*sin(theta_tests)).^2,2);
 [min_sig,I_min]=min(sigma);
 theta_Dorth_1 = theta_tests(I_min);
 rho_Dorth_1=x_G*cos(theta_Dorth_1)+y_G*sin(theta_Dorth_1);
end
