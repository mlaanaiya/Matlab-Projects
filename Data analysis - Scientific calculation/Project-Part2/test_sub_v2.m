%
% PARAMÈTRES
%%%%%%%%%%%%

% taille de la matrice symétrique
n = 500;

% type de la matrice (voir matgen_csad)
% imat == 1 valeurs propres D(i) = i
% imat == 2 valeurs propres D(i) = random(1/cond, 1) avec leur logarithmes
%                                  uniformément répartie, cond = 1e10
% imat == 3 valeurs propres D(i) = cond**(-(i-1)/(n-1)) avec cond = 1e5
% imat == 4 valeurs propres D(i) = 1 - ((i-1)/(n-1))*(1 - 1/cond) avec cond = 1e2
imat = 2;

% tolérance
eps = 1e-8;
% nombre d'itérations max pour atteindre la convergence
maxit = 10000;

% on génère la matrice (1) ou on lit dans un fichier (0)
genere = 1;

% méthode de calcul
v = 0; % subspace iteration v0

% nombre de valeurs propres cherchées (v0)
m = 20;

% puissance p
p = [1,5,10,20,40,100];

% percentage
percentage = 0.4;


fprintf('\n******* calcul avec subspace iteration v2 ******\n');

% Génération d'une matrice rectangulaire aléatoire symétrique définie
    % positive A de taille (n x n)
    % A matrice
    % D ses valeurs propres
    fprintf('\n******* création de la matrice ******\n');
    % appel à eig de matlab : calcul de toutes les valeurs propres
    t_v =  cputime;
    [A, D, ~] = matgen_csad(imat,n);
    t_v = cputime-t_v;
    fprintf('\nTemps de création de la matrice = %0.3e\n',t_v)
    save(['A_' num2str(n) '_' num2str(imat)], 'A', 'D', 'imat', 'n');
    
    % appel à la version 2 de la subspace iteration method
        
        t_v2 =  cputime;
        % W1 valeurs propres
        % V1 vecteurs propres
        for i=1:length(p)
          fprintf('\nTest pour p = %d\n',p(i))
          [ W2, V2, n_ev2, it2, itv2, flag2 ] = subspace_iter_v2( A, m, percentage, p(i), eps, maxit );
          t_v2 = cputime-t_v2;
          
          if(flag2 == 0)
              
              [q2, qv2] = verification_qualite(A, D, W2, V2, n_ev2);
              
              fprintf('\nTemps subspace iteration v2 = %0.3e\n',t_v2);
              fprintf('Nombre d''itérations : %d\n', it2);
              fprintf('Nombre de valeurs propres pour attendre le pourcentage = %d\n', n_ev2);

          else
              if(flag2 == 1)
                  fprintf('subspace iteration v2 : pourcentage %0.3e non atteint avec %d valeurs propres\n', percentage, m);
              else
                  fprintf('subspace iteration v2 : convergence non atteinte: %d\n', it2)
              end
              
              W = 0;
              V = 0;
          end
          
          flag = flag2;
        end