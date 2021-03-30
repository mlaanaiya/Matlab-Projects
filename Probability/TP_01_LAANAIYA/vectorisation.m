function [X,Y] = vectorisation(I);
A=I(:,1:end-1);
B=I(:,2:end);
X=A(:);
Y=B(:);
end
