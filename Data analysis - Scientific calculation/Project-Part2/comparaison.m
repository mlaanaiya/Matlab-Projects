
[A,D,info] = matgen_csad(1,500);
subplot(2,2,1);
plot(D);
title("Distritbution spectrale d'une matrice de type 1");

[A,D,info] = matgen_csad(2,500);
subplot(2,2,2);
plot(D);
title("Distritbution spectrale d'une matrice de type 2");

[A,D,info] = matgen_csad(3,500);
subplot(2,2,3);
plot(D);
title("Distritbution spectrale d'une matrice de type 3");

[A,D,info] = matgen_csad(4,500);
subplot(2,2,4);
plot(D);
title("Distritbution spectrale d'une matrice de type 4");

