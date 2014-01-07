%%%%%%%%%%%%%%%%%%
%Zhang's BCI method
%grab window, filter
%%%%%%%%%%%%%%%%%%

fs = EEG.srate;

fc = [8, 9.75, 11.89, 14.49, 17.67, 21.53, 26.25,32]; %center frequencies of each bank
qf = 1/0.33; % q factor (inverse q mentioned in paper)
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
    classLabels(lat:lat+4*fs) = currEvent.type;
end

%% generate epochs (not processed by filter bank)

windowSize = 1.5; %length of window in seconds

wl = windowSize*fs;

fCL = 1/10; %the fraction of classlabels to use for training 
nNC = length(find(classLabels(wl:fCL*end) == 0)); %number of NC trials
nM1 = length(find(classLabels(wl:fCL*end) == 1));
nM2 = length(find(classLabels(wl:fCL*end) == -1));

NC = zeros(59,wl,nNC); M1 = zeros(59, wl,nM1) ; M2 = zeros(59, wl,nM2); %memory initilization of epoch sets
tNC = zeros(nNC,1); tM1 = zeros(nM1,1) ; tM2 = zeros(nM2,1); %memory initialization of time sets

nci = 1;
m1i = 1;
m2i = 1;

for i = wl:length(classLabels)*fCL
      E = EEG.data(:, (i-wl+1):i);
        switch classLabels(i)
            case 0
                NC(:,:,nci) = E;
                tNC(nci) = i;
                nci = nci+1;
            case 1
                M1(:,:,m1i) = E;
                tM1(nci) = i;
                m1i = m1i+1;
            case -1
                M2(:,:,m2i) = E;
                tM2(nci) = i; 
                m2i = m2i+1;
        end
end
clear ALLEEG;
clear EEG;

%% Finding the PCA substates
disp('Finding the PCA substates...');


restoredefaultpath
fprintf('Number of relevant PCs...');
cp = getNumRelevantPCs(NC); %returns 4.0866 for wl of 1.5s
fprintf('%d\n', cp);
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
addpath(genpath('../Libraries/eeglab12_0_2_5b'));   

disp('Running k-means...');
k = 2; %the paper says 2 clusters are the best, we will cv to make sure later on
[IDX, C] = kmeans(ncu', 2);
NCS1 = NC(:,:, IDX==1); tNCS1 = tNC(IDX==1); nNCS1 = sum(IDX==1); 
NCS2 = NC(:,:, IDX==2); tNCS2 = tNC(IDX==2); nNCS2 = sum(IDX==2);

%% begin csp filter calculation

Wm1vm2 = zeros(4,  size(M1, 1), length(fc));
Wm1vnc1 = zeros(4, size(M1, 1), length(fc));
Wm1vnc2 = zeros(4, size(M1, 1), length(fc));
Wm2vnc1 = zeros(4, size(M1, 1), length(fc));
Wm2vnc2 = zeros(4, size(M1, 1), length(fc));

for fi = 1:length(fc)
    fprintf('Beginning CSP of filter bank w/ fc = %0.2f...\n', fc(fi));
    
    fNCS1 = zeros(59,wl,nNCS1);fNCS2 = zeros(59,wl,nNCS2); fM1 = zeros(59, wl,nM1) ; fM2 = zeros(59, wl,nM2); %memory initilization of epoch processed by filter bank
    
    restoredefaultpath;
    for i = 1:nNCS1,    fNCS1(:,:,i) = filtfilt(b(fi,:), a(fi,:), double(NCS1(:,:,i)'))'; end
    for i = 1:nNCS2,    fNCS2(:,:,i) = filtfilt(b(fi,:), a(fi,:), double(NCS2(:,:,i)'))'; end
    for i = 1:nM1,      fM1(:,:,i) = filtfilt(b(fi,:), a(fi,:), double(M1(:,:,i)'))'; end
    for i = 1:nM2,      fM2(:,:,i) = filtfilt(b(fi,:), a(fi,:), double(M2(:,:,i)'))'; end
    addpath(genpath('../Libraries/eeglab12_0_2_5b'));   
    
    %% generate csp vectors
    disp('Generate csp vectors...');

    % TODO: getCSPVectors could perform the filtering done above, this will
    % avoid memory issues.
    Wm1vm2(:,:,fi) = getCSPVectors(fM1, fM2);
    Wm1vnc1(:,:,fi) = getCSPVectors(fM1, fNCS1);
    Wm1vnc2(:,:,fi) = getCSPVectors(fM1, fNCS2);
    Wm2vnc1(:,:,fi) = getCSPVectors(fM2, fNCS1);
    Wm2vnc2(:,:,fi) = getCSPVectors(fM2, fNCS2);
    
end
%% Collect features
disp('Collecting features...');
feNCS1 = zeros(size(NCS1, 3),4*5*length(fc)); feNCS2 = zeros(size(NCS2, 3),4*5*length(fc)); feM1 = zeros(nM1, 4*5*length(fc)) ; feM2 = zeros(nM2, 4*5*length(fc)); %memory initilization, observations X features
    
% TODO: The epoch sets used in feature collection should NOT be the same as
% thosed used for calculating the csp weights (...maybe)

   restoredefaultpath
for fi = 1:length(fc)
        fprintf('Collecting features of filter bank w/ fc = %0.2f...(%d)\n', fc(fi), fi);
     
       
        r = (fi-1)*20+1:(fi-1)*20+20;
        
%         disp('ncs1 features...');
        for i = 1:nNCS1
            fNCS1 = filtfilt(b(fi,:), a(fi,:), double(NCS1(:,:,i)'))'; 
            feNCS1(i,r(1):r(4))     = diag(  Wm1vm2(:,:,fi)*fNCS1 * ( Wm1vm2(:,:,fi)*fNCS1)' )';
            feNCS1(i,r(5):r(8))     = diag( Wm1vnc1(:,:,fi)*fNCS1 * (Wm1vnc1(:,:,fi)*fNCS1)' )';
            feNCS1(i,r(9):r(12))    = diag( Wm1vnc2(:,:,fi)*fNCS1 * (Wm1vnc2(:,:,fi)*fNCS1)' )';
            feNCS1(i,r(13):r(16))   = diag( Wm2vnc1(:,:,fi)*fNCS1 * (Wm2vnc1(:,:,fi)*fNCS1)' )';
            feNCS1(i,r(17):r(20))   = diag( Wm2vnc2(:,:,fi)*fNCS1 * (Wm2vnc2(:,:,fi)*fNCS1)' )';
        end
        
%         disp('ncs2 features...');
        for i = 1:nNCS2
            fNCS2 = filtfilt(b(fi,:), a(fi,:), double(NCS2(:,:,i)'))'; 
            feNCS2(i,r(1):r(4))     = diag(  Wm1vm2(:,:,fi)*fNCS2 * ( Wm1vm2(:,:,fi)*fNCS2)' )';
            feNCS2(i,r(5):r(8))     = diag( Wm1vnc1(:,:,fi)*fNCS2 * (Wm1vnc1(:,:,fi)*fNCS2)' )';
            feNCS2(i,r(9):r(12))    = diag( Wm1vnc2(:,:,fi)*fNCS2 * (Wm1vnc2(:,:,fi)*fNCS2)' )';
            feNCS2(i,r(13):r(16))   = diag( Wm2vnc1(:,:,fi)*fNCS2 * (Wm2vnc1(:,:,fi)*fNCS2)' )';
            feNCS2(i,r(17):r(20))   = diag( Wm2vnc2(:,:,fi)*fNCS2 * (Wm2vnc2(:,:,fi)*fNCS2)' )';
        end
        
%         disp('m1 features...');
        for i = 1:nM1
            fM1 = filtfilt(b(fi,:), a(fi,:), double(M1(:,:,i)'))'; 
            feM1(i,r(1):r(4))     = diag(  Wm1vm2(:,:,fi)*fM1 * ( Wm1vm2(:,:,fi)*fM1)' )';
            feM1(i,r(5):r(8))     = diag( Wm1vnc1(:,:,fi)*fM1 * (Wm1vnc1(:,:,fi)*fM1)' )';
            feM1(i,r(9):r(12))    = diag( Wm1vnc2(:,:,fi)*fM1 * (Wm1vnc2(:,:,fi)*fM1)' )';
            feM1(i,r(13):r(16))   = diag( Wm2vnc1(:,:,fi)*fM1 * (Wm2vnc1(:,:,fi)*fM1)' )';
            feM1(i,r(17):r(20))   = diag( Wm2vnc2(:,:,fi)*fM1 * (Wm2vnc2(:,:,fi)*fM1)' )';
        end
        
%         disp('m2 features');
        for i = 1:nM2
            fM2 = filtfilt(b(fi,:), a(fi,:), double(M2(:,:,i)'))'; 
            feM2(i,r(1):r(4))     = diag(  Wm1vm2(:,:,fi)*fM2 * ( Wm1vm2(:,:,fi)*fM2)' )';
            feM2(i,r(5):r(8))     = diag( Wm1vnc1(:,:,fi)*fM2 * (Wm1vnc1(:,:,fi)*fM2)' )';
            feM2(i,r(9):r(12))    = diag( Wm1vnc2(:,:,fi)*fM2 * (Wm1vnc2(:,:,fi)*fM2)' )';
            feM2(i,r(13):r(16))   = diag( Wm2vnc1(:,:,fi)*fM2 * (Wm2vnc1(:,:,fi)*fM2)' )';
            feM2(i,r(17):r(20))   = diag( Wm2vnc2(:,:,fi)*fM2 * (Wm2vnc2(:,:,fi)*fM2)' )';
        end
end
   
addpath(genpath('../Libraries/eeglab12_0_2_5b'));
%% feature selection (using mutual information)
disp('Calculating MI....');
%MI(A, C) = MI(C,A) = H(A) - H(A|C); 
A = [feNCS1; feNCS2; feM1; feM2];
p = size(A,2);
n = size(A,1);

S = (4/((p+2)*n))^(1/p+4)*(diag(var(A))); % bowman and azzalini optimal h
Sinv = inv(S);

m = mean(A)';
%zhang equation 12
k=p;
phi = @(t) sqrt((2*pi)^k*det(S))*exp(-0.5 * t'*Sinv*t); 

HA = 0; 
for i = 1:n
    HAi = 0;
    for j= 1:n, HAi = HAi + phi(A(i, :)' - A(j, :)');disp(HAi);  end
    HA = HA + (-1/n)*log(HAi);
end
HA = -1/n*HA;

HA_C =  0;
ns = size(feNCS1,1);
for i = 1:ns
    HA_Ci = 0;
    for j= 1:ns, HA_Ci = HA_Ci + phi(feNCS1(i, :)' - feNCS1(j, :)');  end
    HA_Co = HA_Co + (-1/n)*log(HAi);
end

ns = size(feNCS2,1);
for i = 1:ns
    HA_Ci = 0;
    for j= 1:ns, HA_Ci = HA_Ci + phi(feNCS2(i, :)' - feNCS2(j, :)');  end
    HA_Co = HA_Co + (-1/n)*log(HAi);
end
Pc = (size(feNCS1,1)+size(feNCS2,1)) / size(A,1);
HA_C = HA_C + -1/n*HA_Co*Pc;

ns = size(feM1,1);
for i = 1:ns
    HA_Ci = 0;
    for j= 1:ns, HA_Ci = HA_Ci + phi(feM1(i, :)' - feM1(j, :)');  end
    HA_Co = HA_Co + (-1/n)*log(HAi);
end
Pc = size(feM1,1)/size(A,1);
HA_C = HA_C + -1/n*HA_Co*Pc;

ns = size(feM2,1);
for i = 1:ns
    HA_Ci = 0;
    for j= 1:ns, HA_Ci = HA_Ci + phi(feM2(i, :)' - feM2(j, :)');  end
    HA_Co = HA_Co + (-1/n)*log(HAi);
end
Pc = size(feM2,1)/size(A,1);
HA_C = HA_C + -1/n*HA_Co*Pc;


MI = HA-HA_C;
