function dens = get_IntraInterHemDensity_fromAgreement(CC, hem)

dens = zeros(2);

for i=[1 2]
    idx = logical(hem==i);
    for j=[1 2]
        jdx = logical(hem==j);
        tmp_cc = CC(idx,jdx);
        if i==j
            tmp_cc = tmp_cc-eye(size(tmp_cc,1));
        end
        dens(i,j) = mean(tmp_cc(:), 1, 'omitnan');
    end
end