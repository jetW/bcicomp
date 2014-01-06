%%%%%%%%%%
%% Returns the 1,2,last-1, and last rows of W, i.e 'stationary spatial filters'
%%%%%%%%%%
function W = getCSPVectors(M1, M2)

    R1 = zeros(size(M1, 1));
    R2 = zeros(size(M1, 1));
    
    for i = 1:size(M1,3)
       E = M1(:,:,i);
       [chans,T]=size(E);

       E = E-repmat(mean(E,2),[1 T]);

       tmpC = (E*E');
       R1 = R1 + tmpC./trace(tmpC);
    end
    R1 = R1/i;

    for i = 1:size(M2,3)
       E = M2(:,:,i);
       [chans,T]=size(E);

       E = E-repmat(mean(E,2),[1 T]);

       tmpC = (E*E');
       R2 = R2 + tmpC./trace(tmpC);
    end
    R2 = R2/i;

    %just to see if there were any differences....too much overlap to notice
    %anything....
    % for i = 1:size(NC,3)
    %    E = NC(:,:,i);
    %    E = E-repmat(mean(E,2),[1 T]);
    %    
    %    tmpC = (E*E');
    %    Rnc = Rnc + tmpC./trace(tmpC);
    % end
    % Rnc = Rnc/i;

    Ccomposite=R1+R2;

    [Ucomposite,Lambdacomposite] = eig(Ccomposite);
    [Lambdacomposite,ind] = sort(diag(Lambdacomposite),'descend');
    Ucomposite = Ucomposite(:,ind);
    P=sqrt(inv(diag(Lambdacomposite)))*Ucomposite';

    S{1}=P*R1*P';
    S{2}=P*R2*P';

    [B,D] = eig(S{1},S{2});
    [D,ind] = sort(diag(D)); B = B(:,ind);

    W=(B'*P);
    for i=1:length(ind), W(i,:)=W(i,:)./norm(W(i,:)); end
    
    W = W([1, 2, end-1, end], :);
end