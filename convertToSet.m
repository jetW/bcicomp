%%%%%%%%%%%%%%%%%%
% Convert eeg data from bcicomp to set files for easy extraction
%%%%%%%%%%%%%%%%%%
addpath(genpath('../Libraries/eeglab12_0_2_5b'));
f = ls('dataset1/*calib*.mat');
eeglab; clc;

for fid = 1:size(f,1) 
   
   fn = cat(2, 'dataset1/',strtrim(f(fid,:)));
   data= load(fn);
    
   %in order to store eeg information in EEGLAB through scripts, we must
   %first create a struct with several properties. Some of these properties
   %are known, but some are not, but they still need to be initialized.
   NE.setname = strtrim(f(fid,:));
   NE.filename = '';
   NE.filepath = '';
   NE.pnts = size(data.cnt, 1);
   NE.nbchan = size(data.cnt,2);
   NE.srate = data.nfo.fs;
   %NE.ref = '';
   NE.data = 0.1*double(data.cnt');
   NE.icawinv = [];
   NE.icasphere= [];
   NE.icaweights = [];
   NE.icaact = [];  
   NE.trials = 0;
   NE.comments =[];
   NE.xmin = 0;
   NE.xmax = NE.pnts*NE.srate;
   
   %channel information processing
   for i = 1:length(data.nfo.clab)
        NE.chanlocs(i).labels = data.nfo.clab{i};
        NE.chanlocs(i).X = data.nfo.xpos(i);
        NE.chanlocs(i).Y= data.nfo.ypos(i);
   end
    
   %event information processing
   for i = 1:length(data.mrk.y);
       NE.event(i).type = data.mrk.y(i);
       NE.event(i).latency = data.mrk.pos(i);
   end
   
   eeg_checkset(NE);
   classes=data.nfo.classes;
   
   EEG = pop_rmbase(EEG, [], []);
   
   [ALLEEG, EEG, ~] = eeg_store(ALLEEG, NE); 
   EEG = pop_saveset( EEG, 'filename',EEG.setname(1:end-4),'filepath','C:\\Users\\pt2091\\Documents\\MATLAB\\bcicomp\\sets\\');
end