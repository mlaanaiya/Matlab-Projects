clear;
close all;

load eigenfaces;

% Tirage aleatoire d'une image de test :
personne = randi(nb_personnes);
posture = randi(nb_postures);
% si on veut tester/mettre au point, on fixe l'individu
individu = 17;
posture = 5;

ficF = strcat('./Data/', liste_personnes{individu}, liste_postures{posture}, '-300x400.gif')
img1 = imread(ficF);
image_test = double(transpose(img(:)));
 
% Pourcentage d'information 
per = 0.95;

% Nombre N de composantes principales a prendre en compte 
% [dans un second temps, N peut etre calcule pour atteindre le pourcentage
% d'information avec N valeurs propres (contraste)] :
N = find(C > per,1,'first'); %% Extraction de l'indice du contraste supérieur au pourcentage requis

CP = X*W; % Toutes 
CP = CP(:,1:N);

%...
% Détermination de l'individu le plus proche
%% A la recherche de l'individu
DataA  = CP;
ficF = strcat('./Data/', liste_personnes{individu}, liste_postures{posture}, '-300x400.gif');
img = imread(ficF);
DataT = double(transpose(img(:))); 
DataT = DataT*W; % Calcul de la composante principale totale
DataT = DataT(1:N); % N composantes principales
labelA = repmat(liste_personnes_base,4,1);
labelA = labelA(:);
Nt_test = 1;
K = 1;
ListeClass = repmat(liste_personnes,6,1);
[Partition] = kppv(DataA, labelA, DataT, Nt_test, K, ListeClass);
personne_proche = Partition;
if(strcmp(personne_proche{1}(1),'f'))
    indice_personne_proche = str2double(personne_proche{1}(2:end));
else
    indice_personne_proche = sum(str2double(personne_proche{1}(2:end))+16);
end
%% A la recherche de sa posture
Y = [];
for k = liste_postures_base
    ficF = strcat('./Data/', liste_personnes{indice_personne_proche}, liste_postures{k}, '-300x400.gif');
    img = imread(ficF);
    Y = [Y ; double(transpose(img(:)))];
end
% Remplissage de la matrice Y :
DataA = Y*W; % Composante principales
DataA = DataA(:,1:N); % N premières composantes principales
labelA = liste_postures_base;
Nt_test = 1;
K = 1;
ListeClass = liste_postures;
[Partition] = kppv(DataA, labelA, DataT, Nt_test, K, ListeClass);
posture_proche = Partition;
indice_posture_proche = posture_proche;

figure('Name','Image tiree aleatoirement','Position',[0.2*L,0.2*H,0.8*L,0.5*H]);

subplot(1, 2, 1);
% Affichage de l'image de test :
colormap gray;
imagesc(img1);
title({['posture ' num2str(posture) ' de ', liste_personnes{individu}]}, 'FontSize', 20);
axis image;


ficF = strcat('./Data/', liste_personnes{indice_personne_proche}, liste_postures{indice_posture_proche}, '-300x400.gif')
img = imread(ficF);
        
subplot(1, 2, 2);
imagesc(img);
title({['posture ' num2str(posture_proche) ' de ', liste_personnes{indice_personne_proche}]}, 'FontSize', 20);
axis image;
