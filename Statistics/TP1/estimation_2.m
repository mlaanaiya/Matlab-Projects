function [C_estime , R_moyen] = estimation_2(x_donnees_bruitees,y_donnees_bruitees,n_tests)
    barycentre_x = mean(x_donnees_bruitees);
    G1 = repmat(barycentre_x,1,length(x_donnees_bruitees));
    barycentre_y = mean(y_donnees_bruitees); 
    G2 = repmat(barycentre_y,1,length(y_donnees_bruitees));
    Gx = x_donnees_bruitees - G1;
    Gy = y_donnees_bruitees - G2;
    R_moyen = mean(sqrt(Gx.^2 + Gy.^2));
    B = [barycentre_x barycentre_y];
    G = repmat(B,n_tests,1);
    C = G + rand(n_tests,2)*R_moyen - repmat(R_moyen*0.5,n_tests,2);
    repmat(C(:,1),1,size(x_donnees_bruitees,2));
    repmat(x_donnees_bruitees,n_tests,1);
    repmat(C(:,2),1,size(y_donnees_bruitees,2));
    repmat(y_donnees_bruitees,n_tests,1);
    R = repmat(0.5,n_tests,1).*R_moyen + rand(n_tests,1)*R_moyen;
    repmat(R,1,size(x_donnees_bruitees,2));
    L = sqrt((C(:,1) - x_donnees_bruitees).^2 + (C(:,2) - y_donnees_bruitees).^2);
    L = (L - R).^2;
    Somme = sum(transpose(L));
    [m, b] = min(Somme);
    C_estime = C(b, :);
    R_estime = R(b);
end 