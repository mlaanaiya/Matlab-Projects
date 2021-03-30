function frequences = calcul_frequences(texte,alphabet);
frequences = zeros(length(alphabet));
for i=1:length(alphabet)
    frequences(i)=(1/length(texte))*length(find(texte==alphabet(i)));
end
end
    
