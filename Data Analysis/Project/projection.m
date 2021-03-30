clear
close all

load eigenfaces;

h = figure('Position',[0,0,0.67*L,0.67*H]);
figure('Name','RMSE en fonction du nombre de composantes principales','Position',[0.67*L,0,0.33*L,0.3*L]);

% Calcul de la RMSE entre images originales et images reconstruites :
RMSE_max = 0;

% Composantes principales des données d'apprentissage
C = (X_centre*W')';

for q = 0:n-1
    Q = C(1:q,:);		% q premières composantes principales
    E = W(1:q,:);		% q premières eigenfaces
    X_reconstruit = X_centre*(E)'*E+ones(n,1)*individu_moyen;
    figure(1);
    set(h,'Name',['Utilisation des ' num2str(q) ' premieres composantes principales']);
    colormap gray;
    hold off;
    for k = 1:n
        subplot(nb_personnes, nb_postures,k);
        img = reshape(X_reconstruit(k,:),nb_lignes,nb_colonnes);
        imagesc(img);
        hold on;
        axis image;
        axis off;
    end
    
    figure(2);
    hold on;
    Y_predicted = X_reconstruit(:);
    Y_actual = X(:);
    diff = (Y_predicted-Y_actual).^2;
    RMSE = sqrt(sum(diff)/length(Y_actual));
    RMSE_max = max(RMSE,RMSE_max);

    plot(q,RMSE,'r+','MarkerSize',8,'LineWidth',2);
    axis([0 n-1 0 1.1*RMSE_max]);
    set(gca,'FontSize',20);
    hx = xlabel('$q$','FontSize',30);
    set(hx,'Interpreter','Latex');
    ylabel('RMSE','FontSize',30);
    
    pause(0.01);
end

% save projection;
