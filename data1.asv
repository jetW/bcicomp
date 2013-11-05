addpath('dataset1');
addpath(genpath('../Libraries/eeglab12_0_2_5b'));
f = ls('dataset1/*calib*.mat');

% eeglab
for fid = 1%:size(f,1)
   
   fn = cat(2, 'dataset1/',strtrim(f(fid,:)));
   data= load(fn);
    
   NE.setname = num2str(fid);
   NE.filename = '';
   NE.filepath = '';
   NE.pnts = size(data.cnt, 1);
   NE.nbchans = size(data.cnt,2);
   NE.srate = data.nfo.fs;
   %NE.ref = '';
   NE.data = data.cnt';
   
   %channel information processing
   for i = 1:length(data.nfo.clab)
        NE.chanlocs(i).labels = data.nfo.clab{i};
        NE.chanlocs(i).X = data.nfo.xpos(i);
        NE.chanlocs(i).Y= data.nfo.ypos(i);
   end
    
   %event information processing
   for i = 1:length(data.mrk.y);
       
   end
   
   eeg_checkset(NE);
   disp('~~~we good~~~');
   pause
end

