function[C_x,C_y,M]=matrice_inertie(E_x,E_y,G_norme_E)
n = size(G_norme_E,1);
P = ones(1,n) * G_norme_E;
C_x = 0;
C_y = 0;
M = zeros(2,2);
for i = 1:n
    C_x = C_x + G_norme_E(i,1)*E_x(i,1);
    C_y = C_y + G_norme_E(i,1)*E_y(i,1);
end
C_x = C_x/sum(P);
C_y = C_y/sum(P);
for i = 1:n
    M(1,1) = M(1,1) + G_norme_E(i,1)*(E_x(i,1)-C_x).^2;
    M(1,2) = M(1,2) + G_norme_E(i,1)*(E_x(i,1)-C_x)*(E_y(i,1)-C_y);
    M(2,1) = M(2,1) + G_norme_E(i,1)*(E_x(i,1)-C_x)*(E_y(i,1)-C_y);
    M(2,2) = M(2,2) + G_norme_E(i,1)*(E_y(i,1)-C_y).^2;
end
M=M./sum(P);
end