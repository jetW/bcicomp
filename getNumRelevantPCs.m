%returns the ceiled average principle components that encomposs 95%
%variability of the data.
function cp = getNumRelevantPCs(NC)
    cp = 0;
    for i  = 1:size(NC, 3)
        E = NC(:,:,i);  
        [COEFF, SCORE, LATENT] = pca(E');

        p = cumsum(LATENT)/sum(LATENT);
        pi = find(p >= 0.95);
        pi = pi(1);
        fprintf('%d) %d\n', i, pi);
        cp = cp+ pi;
    end
    cp = ceil(cp/i);
end
