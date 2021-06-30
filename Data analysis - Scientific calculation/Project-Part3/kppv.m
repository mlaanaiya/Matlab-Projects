%--------------------------------------------------------------------------
% ENSEEIHT - 1SN - Analyse de donnees
% TP4 - Reconnaissance de chiffres manuscrits par k plus proches voisins
% fonction kppv.m
%
% Données :
% DataA      : les données d'apprentissage (connues) nos 4 visages et  4
% posture
% visages avec postures
% LabelA     : les labels des données d'apprentissage 1-4 et 1-4
%
% DataT      : les données de test (on veut trouver leur label) notre image
% aléatoire
% Nt_test    : nombre de données tests qu'on veut labelliser 1
%
% K          : le K de l'algorithme des k-plus-proches-voisins 1
% ListeClass : les classes possibles (== les labels possibles) 1-32 et 1-6
%
% Résultat :
% Partition : pour les Nt_test données de test, le label calculé
%
%--------------------------------------------------------------------------
function [Partition] = kppv(DataA, labelA, DataT, Nt_test, K, ListeClass)

[Na,~] = size(DataA);

% Initialisation du vecteur d'étiquetage des images tests

disp(['Classification des images test dans ' num2str(length(ListeClass)) ' classes'])
disp(['par la methode des ' num2str(K) ' plus proches voisins:'])

% Boucle sur les vecteurs test de l'ensemble de l'évaluation
for i = 1:Nt_test
    
    disp(['image test n°' num2str(i)])

    % Calcul des distances entre les vecteurs de test 
    % et les vecteurs d'apprentissage (voisins)
    % À COMPLÉTER
    distances = sqrt(sum((repmat(DataT,Na,1)-DataA).^2,2));
    % On ne garde que les indices des K + proches voisins
    % À COMPLÉTER
    [~,index_sorted] = sort(distances);
    index_K = index_sorted(1:K);

    individu_trouve = labelA(index_K);

    % Assignation de l'étiquette correspondant à la classe trouvée au point 
    % correspondant à la i-ème image test dans le vecteur "Partition" 
    % À COMPLÉTER
    Partition = individu_trouve;    
end

