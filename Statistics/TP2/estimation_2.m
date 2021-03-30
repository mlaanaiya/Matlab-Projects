function [a_DYX_2,b_DYX_2] = estimation_2(x_donnees_bruitees,y_donnees_bruitees);
n=size(x_donnees_bruitees,2);
A=[x_donnees_bruitees',ones(n,1)];
B = y_donnees_bruitees';
X=pinv(A)*B;%A\B 
a_DYX_2=X(1);
b_DYX_2=X(2);
end

