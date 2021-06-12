function [TEB_bruit, mapping, X_ech] = Calcul_TEB_Egaliseur(precision, eb_N0_db, bits, Ns, hc, h, M, heg)

% Tampon pour calculer le TEB avec précision
tmp = zeros(1,length(eb_N0_db));
for i = 1 : precision
    Donnees = bits;
    for k = 1 : length(eb_N0_db)
        % Mapping
        mapping = 2*Donnees - 1;
        % Positionnement des zéros entre les ak
        Kronr = kron(mapping, [1 zeros(1, Ns-1)]);
        % Filtrage
        xe = filter(h, 1, Kronr);
        % Filtre Hc
        ye = filter(hc, 1, xe);
        % Bruit
        Px = mean(abs(ye).^2);
        sigma2 = Px*Ns/(log2(M)*10^(eb_N0_db(k)/10));
        bruit = sqrt(sigma2)*randn(1, length(ye));
        % Filtrage de récéption
        X = filter(h, 1, ye+bruit);
        % Échantillonage
        X_ech = X(Ns:Ns:end);
        % Égalisation
        X_ech = filter(heg, 1, X_ech);
        % Décision
        decision = sign(X_ech);
        % Démapping
        demapping = (decision + 1)/2;
        % Calcul du TEB
        test = Donnees - demapping;
        TEB_bruit(k) = length(find(test~=0))/length(Donnees);
    end
    tmp = tmp + TEB_bruit;
end
% TEB calculé avec precision
TEB_bruit = tmp/precision;