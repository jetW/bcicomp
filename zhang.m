%%%%%%%%%%%%%%%%%%
%Zhang's BCI method, FBCSP + PCA decomposition for NC
%%%%%%%%%%%%%%%%%%

fs = EEG.srate;
fc = [8, 9.75, 11.89, 14.49, 17.67, 21.53, 26.25, 32]; %center frequencies
qf = 0.33; % q factor (how do you use a q factor with a fc at 32 with qf 0.33?)
forder = 4;

f1 = @(f0,Q) f0*(sqrt(1+1/(4*Q^2))- 1/(2*Q));
f2 = @(f0,Q) f0*(sqrt(1+1/(4*Q^2))+ 1/(2*Q));

fpi = 1;
Wn = [ f1(fc(fpi), qf), f2(fc(fpi), qf) ] / (fs/2);
[b,a] =cheby2(forder,30, Wn); %need to know r value, for now setting to 20
%s = dfilt.df2(b,a);
%fvtool(s);
%fres = 100/length(EEG.data(1,:));
%fa = 0:fres:100-fres;
%plot(fa,abs(fft(EEG.data(1,:))));

rmpath(genpath('../Libraries/eeglab12_0_2_5b'));
EEG.data = single(filtfilt(b,a,double(EEG.data')))';
addpath(genpath('../Libraries/eeglab12_0_2_5b'));

ALLEEG(1) = EEG;
eeglab redraw
%hold on;
%plot(fa,abs(fft(EEG.data(1,:))), 'r');

%% generate epochs
windowSize = 1.5; %length of window in seconds

wl = windowSize*fs;
nEvents = size(EEG.event,2);
classLabels = zeros(1, EEG.pnts);

for i = 1:nEvents
    currEvent = EEG.event(i);
    lat = currEvent.latency;
    type = currEvent.type;
    classLabels(lat:lat+4*fs) =currEvent.type;
end
fCL = 1/10; %the fraction of classlabels to use for training 
nNC = length(find(classLabels(wl:fCL*end) == 0)); %number of NC trials
nM1 = length(find(classLabels(wl:fCL*end) == 1));
nM2 = length(find(classLabels(wl:fCL*end) == -1));

NC = zeros(59,wl,nNC); M1 = zeros(59, wl,nM1) ; M2 = zeros(59, wl,nM2); %memory initilization

nci = 1;
m1i = 1;
m2i = 1;
for i = wl:length(classLabels)*fCL
    switch classLabels(i)
        case 0
            NC(:,:,nci) = EEG.data(:, (i-wl+1):i);
            nci = nci+1;
        case 1
            M1(:,:,m1i) = EEG.data(:, (i-wl+1):i);
            m1i = m1i+1;
        case -1
            M2(:,:,m2i) = EEG.data(:, (i-wl+1):i);
            m2i = m2i+1;
    end
end
% NC(:,:,1) = [];M1(:,:,1) = [];M2(:,:,1) = [];

%% generate csp vectors

R1 = zeros(size(M1, 1));
R2 = zeros(size(M1, 1));
Rnc = zeros(size(M1,1));
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

%just to see if there were any differences....too much overlap too notice
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
%% Finding the PCA substates

%cp = getNumRelevantPCs(NC) %returns 4.0866 for wl of 1.5s
cp = 5;
ncu = zeros(cp,size(NC, 3)); %power X epochs
for i = 1:size(NC,3)
    E = NC(:,:,i); 
    [chans,T]=size(E);
    
    E = E-repmat(mean(E,2),[1 T]);
    [COEFF, SCORE, LATENT] = pca(E');

    Q = COEFF(:, 1:cp);
    
    %zhang equation 6
    ncu(:,i) = diag((Q'*E)*(Q'*E)');
    disp(i);
end
k = 2; %the paper says 2 clusters are the best, we will cv to make sure later on
[IDX, C] = kmeans(ncu', 2);
NCS1 = NC(:,:, IDX==1);
NCS2 = NC(:,:, IDX==2);

%% feature selection


%%
subplot(131); imagesc(R1); colormap bone
axis square
subplot(132); imagesc(R2); colormap bone
axis square
subplot(133); imagesc(Rnc); colormap bone
axis square