clear;
close all;
load eigenfaces;
figure('Name',['Plan des trois premieres composantes principales'],'Position',[0,0,0.67*L,0.67*H]);

% Composantes principales des donnees d'apprentissage
C = X_centre*W;			

% Affichage des deux premieres composantes principales des donnees d'apprentissage :
hold on;
for i = 1:4
	lignes = (i-1)*nb_postures_base+1:i*nb_postures_base;
	plot3(C(lignes,1),C(lignes,2),C(lignes,3),['*' 0.5*(rand(1,3)+1)],'MarkerSize',10,'LineWidth',2);
end
rotate3d on
grid on
proportion = 0.5;
ecart_x = max(C(:,1))-min(C(:,1));
min_x = min(C(:,1))-proportion*ecart_x;
max_x = max(C(:,1))+proportion*ecart_x;
ecart_y = max(C(:,2))-min(C(:,2));
min_y = min(C(:,2))-proportion*ecart_y;
max_y = max(C(:,2))+proportion*ecart_y;
ecart_z = max(C(:,3))-min(C(:,3));
min_z = min(C(:,3))-proportion*ecart_z;
max_z = max(C(:,3))+proportion*ecart_z;
axis([min_x max_x min_y max_y min_z max_z]);
axis equal;
set(gca,'FontSize',20);
xlabel('$C_1$','FontSize',30,'Interpreter','Latex');
ylabel('$C_2$','FontSize',30,'Interpreter','Latex');
zlabel('$C_3$','FontSize',30,'Interpreter','Latex');
