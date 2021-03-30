clear
close all

liste_postures = {'v1e1','v3e1','v1e2','v3e2','v1e3','v3e3'};

% taille des images
nb_lignes = 400;
nb_colonnes = 300;

nb_personnes = 4; % 2 hommes - 2 femmes

nb_postures = 4;
postures_seance = [1 2 3 4];

X=[];
liste_base=[];

% Dimensions du masque (question 3)
%ligne_min = 200;
%ligne_max = 350;
%colonne_min = 60;
%colonne_max = 290;
%X_masque = [];

taille_ecran = get(0,'ScreenSize');
L = taille_ecran(3);
H = taille_ecran(4);
figure('Name','Individus','Position',[0,0,0.67*L,0.67*H]);
colormap(gray(256));

% Lecture des fichiers Images, Constuction de X
% et Affichage des images sous forme de planche-contact 
% (une personne par ligne, une posture par colonne) :

% Femmes
for j = 1:nb_personnes/2
    no_posture = 0;
	for k = postures_seance
        no_posture = no_posture + 1;
        
        ficF=strcat('./Data/','f',num2str(j),liste_postures{k},'-300x400.gif');
        liste_base=[liste_base;ficF];
        img=imread(ficF);
        % Remplissage de la matrice X :
        X= [X; double(transpose(img(:)))];
        
        % Dégradation de l'image (question 3)
        %img(ligne_min:ligne_max,colonne_min:colonne_max)=0;
        % Remplissage de la matrice X_masque :
        %X_masque= [X_masque; double(transpose(img(:)))];
        
        % Affichage
		subplot(nb_personnes,nb_postures,2*(j-1)*nb_postures+no_posture);
		imagesc(img);
		hold on;
		axis image;
		title(['Personne ' num2str(j) ', posture ' num2str(k)]);
        
	end
end

% Hommes
for j = 1:nb_personnes/2
    no_posture = 0;
	for k = postures_seance
        no_posture = no_posture + 1;
        
        ficM=strcat('./Data/','m',num2str(j),liste_postures{k},'-300x400.gif');
        liste_base=[liste_base;ficM];
        img=imread(ficM);
        % Remplissage de la matrice X :
        X= [X; double(transpose(img(:)))];  
        
        % Dégradation de l'image (question 3)
        %img(ligne_min:ligne_max,colonne_min:colonne_max)=0;
        % Remplissage de la matrice X_masque :
        %X_masque= [X_masque; double(transpose(img(:)))];
    
        % Affichage
		subplot(nb_personnes,nb_postures,(2*(j-1)+1)*nb_postures+no_posture);
		imagesc(img);
		hold on;
		axis image;
		title(['Personne ' num2str(j+16) ', posture ' num2str(k)]);    
        
	end
end

%%%% Partie à compléter %%%%%%%

% nombre d'individus
n = size(X,1);

% Calcul de l'individu moyen :
individu_moyen = mean(X);

% Centrage de la matrice X :
X_centre = X - individu_moyen;

% Calcul de la matrice de covariance (impossible de calculer ainsi à cause de sa taille) :
% Sigma = transpose(X_centre)*X_centre/n

% Calcul de la matrice résultant du calcul commuté (voir annexe du sujet) :
%
Sigma2 = X_centre*X_centre'/n;

% Calcul des vecteurs/valeurs propres de la matrice Sigma2
[V2,D2] = eig(Sigma2);

% Les vecteurs propres de Sigma (les eigenfaces) se déduisent de ceux de Sigma2 :
V = X_centre'*V2;

% Tri par ordre décroissant des valeurs propres de Sigma2 :
D2_diag = diag(D2);
[D2_diag_trie,indices2_tries] = sort(D2_diag,'descend');

% Tri des eigenfaces dans le même ordre (on enlève la dernière eigenface, qui appartient au noyau de Sigma) :
V_tries = V(:,indices2_tries)';
%V_tries = [individu_moyen;V_tries(1:(n-1),:)];

% Normalisation des eigenfaces :
V_tries_normalisee = V_tries'./(sqrt(sum(V_tries'.^2)));
W = V_tries_normalisee(:,1:n-1)';
% Affichage de l'individu moyen et des eigenfaces sous la forme de "pseudo-images" 
% (leurs coordonnees sont interprétees comme des niveaux de gris) :

figure('Name','Individu moyen et eigenfaces','Position',[0,0,0.67*L,0.67*H]);
colormap(gray(256)); 
img = reshape(individu_moyen,nb_lignes,nb_colonnes);
subplot(nb_personnes,nb_postures,1)
imagesc(img); 
hold on; 
axis image; 
title('Individu moyen');
for k = 1:n-1
	img = reshape(W(k,:),nb_lignes,nb_colonnes);
	subplot(nb_personnes,nb_postures,k+1);
	imagesc(img); 
	hold on; 
	axis image; 
	title(['Eigenface ',num2str(k)]);
end

% Sauvegarde des données pour l'exercice 2
save eigenfaces
