function binaryChi = binaryDesign(chi,chiMat,threshold,part,M)
    switch part
        case 'abs'
            material = find(abs(chiMat) >= threshold*chi);
        case 'real'
            material = find(real(chiMat) >= threshold*chi);
    end
    binaryChi = transpose(repelem(1e-7,M));
    binaryChi(material) = chi;
end