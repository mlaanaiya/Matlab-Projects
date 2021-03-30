function [rho_F_1,theta_F_1] = RANSAC_2(rho,theta,parametres)
     n = length(rho);
     res_inf = Inf;
     for i = 1:parametres(3)
        a = randperm(n,2);
        [rho_F , theta_F] = estimation_F(rho(a) , theta(a));
        res = abs(rho - rho_F*cos(theta - theta_F));
        Bool = (res <=parametres(1));
        droites_conformes = sum(Bool)/length(Bool);
        if  droites_conformes >= parametres(2) 
             res = sum(1/length(Bool)*abs(rho(Bool) - rho_F*cos(theta(Bool) - theta_F)));
              if res <= res_inf
                  res_inf = res;
                  rho_F_1 = rho_F;
                  theta_F_1 = theta_F;
              end
        end 
     end 
     



end

