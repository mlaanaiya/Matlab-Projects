function [a_DYX_1,b_DYX_1] = estimation_1(x_donnees_bruitees,y_donnees_bruitees,n_tests)
ab_tests = -pi/2+pi*rand(n_tests,1);
n=size(x_donnees_bruitees,2);
m=size(ab_tests,1);    
x_G=mean(x_donnees_bruitees);
y_G=mean(y_donnees_bruitees);
x_centre=x_donnees_bruitees-x_G;
y_centre=y_donnees_bruitees-y_G;
Y=repmat(y_centre,[m,1]);
sigma=sum((Y-tan(ab_tests).*x_centre).^2,2);
[min_sig,I_min]=min(sigma);
ab_min = ab_tests(I_min);
a_DYX_1= tan(ab_min);
b_DYX_1 = y_G - a_DYX_1*x_G;
end

