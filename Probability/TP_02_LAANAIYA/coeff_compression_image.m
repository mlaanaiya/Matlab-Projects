function coeff = coeff_compression_image(histogramme,dico)
ap_compression=0;
av_compression=0;
for i=1:length(histogramme) 
    ap_compression = ap_compression + length(dico{i,2})*(histogramme(i)/sum(histogramme));
    av_compression=av_compression+histogramme(i)*8/sum(histogramme);
end
coeff = av_compression/ap_compression;
end

