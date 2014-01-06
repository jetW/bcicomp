%%%%%%%%%%%%%%%%%%
%Zhang's BCI method, FBCSP + PCA decomposition for NC
%%%%%%%%%%%%%%%%%%

fs = EEG.srate;

fc = [8, 9.75, 11.89, 14.49, 17.67, 21.53, 26.25,32]; %center frequencies of each bank
qf = 1/0.33; % q factor (how do you use a q factor with a fc at 32 with qf 0.33?)
forder = 4; %filter order

f1 = @(f0,Q) f0*(sqrt(1+1/(4*Q^2))- 1/(2*Q));
f2 = @(f0,Q) f0*(sqrt(1+1/(4*Q^2))+ 1/(2*Q));

a = zeros(length(fc),2*forder+1);
b = zeros(length(fc),2*forder+1);

% get filters coeffs for each bank
for i = 1:length(fc)
    Wn = [ f1(fc(i), qf), f2(fc(i), qf) ] / (fs/2);
    [bi,ai] =cheby2(forder,20, Wn); %need to know r value, for now setting to 20
    b(i, :) = bi;
    a(i, :) = ai;
end

%s = dfilt.df2(b,a);
%fvtool(s);
%fres = 100/length(EEG.data(1,:));
%fa = 0:fres:100-fres;
%plot(fa,abs(fft(EEG.data(1,:))));

%% generate class time series
disp('Generate class time series...');
lat = [EEG.event.latency];

nEvents = size(EEG.event,2);
classLabels = zeros(1, EEG.pnts); %class time series

for i = 1:nEvents
    currEvent = EEG.event(i);
    lat = currEvent.latency;
    type = currEvent.type;
    classLabels(lat:lat+4*fs) =currEvent.type;
end

%% generate epochs

windowSize = 1.5; %length of window in seconds

wl = windowSize*fs;

fCL = 1/10; %the fraction of classlabels to use for training 
nNC = length(find(classLabels(wl:fCL*end) == 0)); %number of NC trials
nM1 = length(find(classLabels(wl:fCL*end) == 1));
nM2 = length(find(classLabels(wl:fCL*end) == -1));

NC = zeros(59,wl,nNC); M1 = zeros(59, wl,nM1) ; M2 = zeros(59, wl,nM2); %memory initilization

nci = 1;
m1i = 1;
m2i = 1;

Wm1vm2 = zeros(4, size(M1, 1), length(fc));
Wm1vnc1 = zeros(4, size(M1, 1), length(fc));
Wm1vnc2 = zeros(4, size(M1, 1), length(fc));
Wm2vnc1 = zeros(4, size(M1, 1), length(fc));
Wm2vnc2 = zeros(4, size(M1, 1), length(fc));

for fi = 1:length(fc)
    restoredefaultpath
%     EEG2.data = filtfilt(b(fi,:), a(fi,:), double(EEG.data'))';
    addpath(genpath('../Libraries/eeglab12_0_2_5b'));
    %%Generate epochs
    fprintf('Beginning CSP of filter bank w/ fc = %0.2f...\n', fc(fi));
    for i = wl:length(classLabels)*fCL
            E = EEG.data(:, (i-wl+1):i);
            E=filtfilt(b(fi,:), a(fi,:), double(E'))';
            switch classLabels(i)
                case 0
                    NC(:,:,nci) = E;
                    nci = nci+1;
                case 1
                    M1(:,:,m1i) = E;
                    m1i = m1i+1;
                case -1
                    M2(:,:,m2i) = E;
                    m2i = m2i+1;
            end
    end

    %% Finding the PCA substates
    disp('Finding the PCA substates...');
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
    end
    
    disp('Running k-means...');
    k = 2; %the paper says 2 clusters are the best, we will cv to make sure later on
    [IDX, C] = kmeans(ncu', 2);
    NCS1 = NC(:,:, IDX==1);
    NCS2 = NC(:,:, IDX==2);

    %% generate csp vectors
    disp('Generate csp vectors...');

    Wm1vm2(:,:,fi) = getCSPVectors(M1, M2);
    Wm1vnc1(:,:,fi) = getCSPVectors(M1, NCS1);
    Wm1vnc2(:,:,fi) = getCSPVectors(M1, NCS2);
    Wm2vnc1(:,:,fi) = getCSPVectors(M2, NCS1);
    Wm2vnc2(:,:,fi) = getCSPVectors(M2, NCS2);
    
    %% grab features
    
end
%% feature selection


%%
subplot(131); imagesc(R1); colormap bone
axis square
subplot(132); imagesc(R2); colormap bone
axis square
subplot(133); imagesc(Rnc); colormap bone
axis square