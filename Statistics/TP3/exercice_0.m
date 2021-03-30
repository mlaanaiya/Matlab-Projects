clear;
close all;
taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);

% Parametres :
alpha = pi/180*12;			% Seuil sur l'orientation du gradient (en radians)
cos_alpha = cos(alpha);
seuil_norme = 2/sin(alpha);		% Seuil sur la norme du gradient (decoule d'une etude theorique)

% Chargement des donnees :
load Donnees/parquet;
%load Donnees/bateau;
[n_lignes,n_colonnes] = size(I);
limites_affichages = [1 n_colonnes 1 n_lignes];

figure('Name','Detection des alignements','Position',[0.2*L,0,0.8*L,H]);

% Affichage de l'image :
subplot(2,2,1);
imagesc(I);
axis equal;
axis off;
colormap gray;
title('Image originale','FontSize',30);
axis(limites_affichages);

% Fenetre d'affichage des pixels de contour :
subplot(2,2,3);
imagesc(I);
axis equal;
axis off;
colormap gray;
hold on;

% Gradient du niveau de gris (x vers la droite, y vers le bas) :
I = double(I);
[G_x,G_y] = gradient(I);
G_norme = sqrt(G_x.^2+G_y.^2);

% Selection des pixels de contour :
contour = G_norme>seuil_norme;

% Pas de pixel de contour sur le bord de l'image => Traitement simplifie !
contour([1 n_lignes],:) = 0;
contour(:,[1 n_colonnes]) = 0;

% Tri des pixels de contour :
indices_contour = find(contour);
[~,indices] = sort(G_norme(indices_contour),'descend');
indices_contour = indices_contour(indices);

% Affichage d'une petite fleche sur les pixels de contour :
[i,j] = ind2sub(size(I),indices_contour);
quiver(j,i,G_x(indices_contour),G_y(indices_contour),'r');	% Attention : x = j et y = i
axis equal;
title('Pixels de contour','FontSize',30);
axis(limites_affichages);
drawnow;

% Affichage des ensembles E (la fonction label2rgb donne a chaque ensemble E une couleur differente) :
subplot(2,2,2);
imagesc(I_resultat);
axis(limites_affichages);
axis equal;
axis off;
title('Ensembles candidats','FontSize',30);

% Affichage de l'esquisse :
subplot(2,2,4);
imagesc(120*ones(size(I)),[0,255]);
axis equal;
axis off;
title('Alignements','FontSize',30);
hold on;
axis(limites_affichages);

% Boucle sur les ensembles E :
for k = 1:size(extremites,3)
	% Affichage du segment :
	line(extremites(1,:,k),extremites(2,:,k),'Color','w','LineWidth',2);
end

save exercice_0;
