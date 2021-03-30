function [selection_frequences,selection_alphabet] = selection(frequences,alphabet);
ind = find(frequences>0);
selection_frequences = zeros(1,length(ind));
selection_alphabet=zeros(1,length(ind));
for i=1:length(ind)
    selection_frequences(i)=frequences(ind(i));
    selection_alphabet(i)=alphabet(ind(i));   
end 
end 
    
    
        
    