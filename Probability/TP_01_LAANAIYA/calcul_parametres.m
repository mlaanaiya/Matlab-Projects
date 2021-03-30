function [r,a,b] = calcul_parametres(X,Y)
sigmaX2=1/length(X) *sum(X.^2)-mean(X)^2;
sigmaY2=1/length(Y) *sum(Y.^2)-mean(Y)^2;
sigmaXY=1/length(X) *sum(X.*Y)-mean(X)*mean(Y);
r=sigmaXY/(sqrt(sigmaX2)*sqrt(sigmaY2));
a=sigmaXY/sigmaX2;
b=mean(Y)-a*mean(X);
end