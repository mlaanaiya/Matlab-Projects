function [TEB, mapping, X, X_ech] = BPSK(bits, Ns, hc, h, canal)

% Mapping
mapping = 2*bits - 1;
% Positionnement des zéros entre les ak
Kronr = kron(mapping, [1 zeros(1, Ns-1)]);
% Filtrage
xe = filter(h, 1, Kronr);
% Filtre Hc
if(canal == 1) % S'il y a le canal
    ye = filter(hc, 1, xe);
else % S'il n'y a pas le canal
    ye = xe;
end
% Filtrage de récéption
X = filter(h, 1, ye);
% Échantillonage
X_ech = X(Ns:Ns:end);
% Décision
decision = sign(X_ech);
% Demapping
demapping = (decision + 1)/2;
% Calcul du TEB
test = bits - demapping ;
TEB = length(find(test~=0))/length(bits);