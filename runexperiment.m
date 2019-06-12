%% Semantics words MEG 
% Tijl Grootswagers

%% set up psychtoolbox if needed (should only be needed for testing)
if isempty(which('PsychtoolboxVersion.m'))
    addpath(genpath('~/PHD/MatlabToolboxes/Psychtoolbox/'))
end
close all; %close all figure windows
clearvars; %clear all variables
p=struct();

%% START define parameters

%testing parameters: CHECK THESE!
p.isrealexperiment = 1; %should only be 0 for testing; skips overwrite checks and doesn't ask for subject nr
p.fullscreen = 1; %should only be 0 for testing; if 0 it does not create a fullscreen window
p.isMEGexperiment = 1; %should only be 0 for testing; does not send actual triggers
p.synctest = 1; %should only be 0 for testing; skips synchronisation tests

% stimulus parameters
p.stimuli ={'ANKLE','BLUE','BOOT','CAPE','CLAM',...
            'CREW','CYAN','EYE','FERRY','FLY',...
            'FOOT','FRUIT','GNAT','GRAY','HOT',...
            'KAYAK','KHAKI','KNEE','LAUGH','LYNX',...
            'NEWT','PEACH','PEAR','PIE','SEW',...
            'SHOE','SIGHT','SLEIGH','SUIT','TAXI',...
            'THIGH','THINK','TIE','WHITE','YACHT'};
p.nstimuli = length(p.stimuli);
p.targetstimuli = {'BED','CHAIR','DESK','COT','COUCH',...
                   'CRIB','FUTON','SOFA','STOOL','TABLE'};

% repetition numbers
p.nblocks = 10; %number of blocks we want
p.stimrepeatsperblock = 4; %how many repeats of each stimulus per block
p.targetsperblock = [8 10]; %how many targets per block

%timing parameters
p.blockstartduration= 1; %time before first stimulus in a block (secs);
p.prestimblank = 0.200; % prestim blank duration (secs);
p.stimulusduration = 0.200; %stimulus duration (secs);
p.responseduration = [1 1.2]; %variable response duration from stimulus onset (secs);
p.halfrefresh = .5/60; %to subtract from fliptimes

%guesstimates
p.trialsperblock = p.stimrepeatsperblock * p.nstimuli + mean(p.targetsperblock);
p.averagetrialduration_seconds = p.prestimblank + mean(p.responseduration);
p.averageblockduration_minutes = p.trialsperblock * p.averagetrialduration_seconds/60;
p.estimatedexperimentduration_minutes = p.nblocks * p.averageblockduration_minutes;

%trigger parameters
p.triggerstimulus=178; %stimulus is ON trigger

%display parameters
p.fixationsize=20; %diameter of fixation cross (pixels)
p.backgroundcolor=[127 127 127]; %backround (gray)
p.fontsize=50; %font size for all the things
p.windowsize=[0 0 800 600]; %windowsize if not running fullscreen

%photodiode
p.photoflashsquaresize=100; %photodiode patch size (pixels)

%responsekey parameters, 1! and 2@ are the MEG response buttons
p.responsekeys = [KbName('1') KbName('2') KbName('3') KbName('4') KbName('1!') KbName('2@')];

%% END define parameters

%% get subject info
if ~p.isrealexperiment
    p.subjectnr = 0;
    disp('    +-----------------------------------------+')
    disp('    | WARNING: p.isrealexperiment is set to 0 |')
    disp('    | If this is a real MEG run, set to 1     |')
    disp('    +-----------------------------------------+')
else
    p.subjectnr = input('\n    Subject number: ','s');
    p.subjectnr = str2double(p.subjectnr);
    if isnan(p.subjectnr)
        error('ERROR: invalid subject number');
    end
end
p.datafilename = sprintf('sub-%02i_task-words_events.mat',p.subjectnr);
p.csvdatafilename = sprintf('sub-%02i_task-words_events.csv',p.subjectnr);

%check if we are possibly overwriting data
if p.isrealexperiment
    if exist(p.datafilename,'file')
        error(['ERROR: ' p.datafilename ' exists. Move (or delete) this file or use a different subject number']);
    end
    if exist(p.csvdatafilename,'file')
        error(['ERROR: ' p.csvdatafilename ' exists. Move (or delete) this file or use a different subject number']);
    end
end

%This should only be used to test/debug. CHECK THIS!!
Screen('Preference', 'SkipSyncTests', ~p.synctest);
if ~p.synctest
    disp('    +-------------------------------------------------+')
    disp('    | WARNING: p.synctest is set to 0                 |')
    disp('    | If this is a MEG experiment, p.synctest to 1    |')
    disp('    +-------------------------------------------------+')
end

% set up seed
p.randomseed = rng(p.subjectnr);

%% open i/o port
if p.isMEGexperiment
    %create io64 interface object
    try
        p.ioObj = io64;
        % check the port status
        status = io64(p.ioObj);
    catch e
        status = 1;
        disp(['Failed to open io64: ' e.message])
    end
else
    status = 1;
end
if status == 0
    p.address = hex2dec('DFB8'); %stim2
    %p.address = hex2dec('D020'); %stim1
    display(['Functioning parallel port opened at: ' num2str(p.address)])
else
    p.ioObj = [];
    p.address = [];
end

%make a copy of the raw experiment code and store it so we can recreate the exact script
p.experimentcode = fileread('runexperiment.m');

eventlist = table();
allstimuli = [p.stimuli p.targetstimuli];

%% create trialstructure for all blocks
for blocknr = 1:p.nblocks
    sequence = [];
    for stimrep = 1:p.stimrepeatsperblock
        stimorder = []; %find a stimulus order so that the last 5 stims do not overlap
        while isempty(stimorder) || stimrep>1 && ~isempty(intersect(sequence(end-5:end),stimorder(1:5)))
            stimorder = randsample(1:length(p.stimuli),p.nstimuli);
        end
        sequence = [sequence;stimorder'];
    end
    % add targets
    nt = randi(p.targetsperblock);
    istarget = false(length(sequence)+nt,1);
    targetlocs = round(linspace(1,length(istarget),nt+2))';
    d = floor(min(diff(targetlocs))/2)-1;
    targetlocs = targetlocs(2:end-1);
    targetlocs = targetlocs + randi([-d, d],nt,1); %jitter
    istarget(targetlocs) = true;
    fullsequence = zeros(size(istarget));
    fullsequence(~istarget) = sequence;
    fullsequence(istarget)  = randsample(length(p.stimuli)+(1:length(p.targetstimuli)),nt);
        
    % add sequence to trialstruct
    r=struct();
    r.eventnumber = size(eventlist,1) + (1:length(fullsequence))';
    r.blocknumber = blocknr*ones(size(fullsequence));
    r.presentationnumber = (1:length(fullsequence))';
    r.stimnumber = fullsequence;
    r.stim = allstimuli(fullsequence)';
    r.istarget = fullsequence > p.nstimuli;
    %concat
    eventlist = [eventlist;struct2table(r)];
end

%% setup abort function
abort = ['sca;Priority(0);ListenChar(0);ShowCursor();fprintf(''\n\nABORTING...\n'');save(''aborted.mat'');sca;'...
    'error(struct(''stack'',struct(''file'',''runexperiment'',''name'',''runexperiment'',''line'',0),'...
    '''message'',''############################## EXPERIMENT ABORTED BY USER ##############################''));'];

%% disable keyboard input
ListenChar(2);
KbName('UnifyKeyNames') %cause old version of psychtoolbox
KbCheck(); % make sure kbcheck is compiled so it is always accurate on first call

%% summarize parameters
while KbCheck()
    WaitSecs(0.01);
end
disp('    +------------------------+')
disp('    | Experiment parameters: |')
disp('    +------------------------+')
disp(p);
disp('If this is correct, press <SPACE> to start.')
[~, keycode, ~] = KbWait(); %wait for a keypress
if ~keycode(KbName('space'));eval(abort);end % check for return key
% Wait for key release

%% open window, and wait for the photodiode setup
p.screennumber=min (Screen('Screens'));
if p.screennumber>1 && p.isrealexperiment
    ListenChar(0);error('Found too many screens! Exit matlab, setup mirrored-display mode, and try again')
end
p.black=BlackIndex(p.screennumber);
p.gray=GrayIndex(p.screennumber);
p.white=WhiteIndex(p.screennumber);
if p.fullscreen
    p.windowsize=[];
elseif ismac && max(Screen('Screens',1))==2
    p.windowsize=[];
end
[p.window, p.windowrect]=Screen('OpenWindow',p.screennumber, p.backgroundcolor,p.windowsize,32,2);
Screen('TextFont',p.window, 'Arial');
Screen('TextSize',p.window,p.fontsize);
Screen('BlendFunction',p.window,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

%photoflash fuction
p.photoflashrect = [p.windowrect(3)-p.photoflashsquaresize p.windowrect(4)-p.photoflashsquaresize p.windowrect(3) p.windowrect(4)];
drawphotoflashrect = @() Screen('FillRect',p.window,p.white,p.photoflashrect);

%fixation fuction
p.fixationlocation = .5*p.windowrect([3 3 3 3;4 4 4 4])+.5*[-p.fixationsize p.fixationsize 0 0;0 0 -p.fixationsize p.fixationsize];
drawfixation = @() Screen('DrawLines', p.window, p.fixationlocation,2,p.white);
drawfixationred = @() Screen('DrawLines', p.window, p.fixationlocation,2,[200 0 0]);
drawfixationgreen = @() Screen('DrawLines', p.window, p.fixationlocation,2,[0 200 0]);


%calibrate
i = 1;keycode=[];

fprintf('\nSTARTING CALIBRATION\n')

%start playback
s=0;
while isempty(keycode) || ~keycode(KbName('space'))
    i=i+1;
    DrawFormattedText(p.window, 'CALIBRATION\n\n<SPACE> to continue or <ESCAPE> to abort', 'center', 100, p.white);
    if mod(i,2)
        drawphotoflashrect()
        if p.isMEGexperiment
            io64(p.ioObj,p.address,p.triggerstimulus-128);
        end
    elseif p.isMEGexperiment
        io64(p.ioObj,p.address,0);
    end
    if mod(i,2)
        s=1+mod(s,p.nstimuli); %next stim
        DrawFormattedText(p.window, p.stimuli{s}, 'center', 'center');
        Screen('Flip',p.window);
    else
        drawfixation();
        Screen('Flip',p.window);
    end
    [~, keycode] = KbWait([],0,GetSecs()+0.5);
    if keycode(KbName('escape'));eval(abort);end
end 

%% disable mouse cursor
if p.isrealexperiment
    HideCursor();
end
Priority(MaxPriority(p.window));

%% make sure no triggers are active
if p.isMEGexperiment
    triggerstatus=1;
    while triggerstatus
        io64(p.ioObj,p.address,0);
        triggerstatus=io64(p.ioObj,p.address);
    end
else
    triggerstatus=0;
end

%% ready to go
% Wait for key release
while KbCheck();
    WaitSecs(0.01);
end
keycode=[];
DrawFormattedText(p.window, 'Calibration Complete\n\nStart the MEG acquisition\n\n<SPACE> to start the experiment or <ESCAPE> to abort', 'center', 'center', p.white);
Screen('Flip',p.window);
while isempty(keycode) || ~keycode(KbName('space'))
    [~, keycode, ~] = KbWait();
    if keycode(KbName('escape'));eval(abort);end
end
    
%% we're go - record start time 
p.time_experiment_start = Screen('Flip', p.window); %record experiment start time
    

%% start experiment loop
nevents = size(eventlist,1);
nblocks = eventlist.blocknumber(end);
currentblock = 0;
for eventnr=1:nevents
    if eventlist.blocknumber(eventnr)>currentblock
        currentblock = eventlist.blocknumber(eventnr);
        %save temp data
        writetable(eventlist,p.csvdatafilename)
        save(p.datafilename,'p');

        % Wait for key release
        while KbCheck()
            WaitSecs(0.01);
        end
        %instructions
        if currentblock == 1
            DrawFormattedText(p.window, sprintf('Block %i\n\nPress the button for furniture words\n\n<Press any button> to start',currentblock), 'center', 'center', p.white);
        else
            idx = eventlist.blocknumber==(currentblock-1);
            hits = sum(eventlist.istarget(idx) & eventlist.response(idx));
            hitrate = 100*mean(eventlist.response(idx & eventlist.istarget==1));
            misses = sum(eventlist.istarget(idx) & ~eventlist.response(idx));
            fa = sum(~eventlist.istarget(idx) & eventlist.response(idx));
            DrawFormattedText(p.window, sprintf('Hits: %i (%.2f%%)\nMisses: %i\nFalse Alarms: %i\n\n\nBlock %i\n\n\nPress the button for furniture words\n\n<Press any button> to start',hits, hitrate, misses, fa, currentblock), 'center', 'center', p.white);
        end
        Screen('Flip', p.window);
        %wait for any keypress to start blockloop
        [~, keycode] = KbWait();
        if keycode(KbName('escape'));eval(abort);end
        drawfixation();
        time_blockstart = Screen('Flip', p.window);
        Screen('Flip', p.window, time_blockstart + p.blockstartduration - p.halfrefresh);
    end
    %wait for stimon trigger
    while triggerstatus~=p.triggerstimulus-128
        triggerstart = Screen('Flip', p.window);
        if p.isMEGexperiment
            io64(p.ioObj,p.address,p.triggerstimulus-128);
            triggerstatus=io64(p.ioObj,p.address);
        else
            triggerstatus=p.triggerstimulus-128;
        end
    end
    %prepare stim
    stim = eventlist.stim{eventnr};
    istarget = eventlist.istarget(eventnr);
    DrawFormattedText(p.window, stim, 'center', 'center');
    drawphotoflashrect();
    time_stimon = Screen('Flip',p.window, triggerstart + p.prestimblank - p.halfrefresh);
    if p.isMEGexperiment
        io64(p.ioObj,p.address,p.triggerstimulus-128);
        triggerstatus=io64(p.ioObj,p.address);
    else
        triggerstatus=p.triggerstimulus-128;
    end
    drawfixation();
    time_stimoff = Screen('Flip', p.window, time_stimon + p.stimulusduration - p.halfrefresh);
    %wait for off trigger
    while triggerstatus
        if p.isMEGexperiment
            io64(p.ioObj,p.address,0);
            triggerstatus=io64(p.ioObj,p.address);
        else
            triggerstatus=0;
        end
    end
    responsetime = min(p.responseduration) + rand()*diff(p.responseduration);
    %check for reponse
    [response,RT]=deal(0);
    while GetSecs()-time_stimon < responsetime
        [keydown, secs, keycode] = KbCheck();
        if keydown
            if ~response && any(keycode(p.responsekeys))
                response = 1;
                RT = secs-time_stimon;
                if istarget
                    drawfixationgreen();
                else
                    drawfixationred();
                end
                Screen('Flip', p.window);
            elseif keycode(KbName('escape')) %emergency exit
                eval(abort);
            end
        end
    end
    correct = response==istarget;
    %update the eventlist with results
    eventlist.response(eventnr) = response;
    eventlist.RT(eventnr) = RT;
    eventlist.correct(eventnr) = correct;
    eventlist.time_stimon(eventnr) = time_stimon-p.time_experiment_start;
    eventlist.time_stimoff(eventnr) = time_stimoff-p.time_experiment_start;
    eventlist.stimdur(eventnr) = time_stimoff-time_stimon;
    %prep next trial
    Screen('Flip', p.window);
end
% finish up 
writetable(eventlist,p.csvdatafilename)
save(p.datafilename,'p','eventlist');

hits = sum(eventlist.istarget & eventlist.response);
hitrate = 100*mean(eventlist.response(eventlist.istarget==1));
misses = sum(eventlist.istarget & ~eventlist.response);
fa = sum(~eventlist.istarget & eventlist.response);
DrawFormattedText(p.window, sprintf('Hits: %i (%.2f%%)\nMisses: %i\nFalse Alarms: %i\n\n\nExperiment complete!\n\nRelax and wait for experimenter...\n\nExperimenter: press <space> to exit',hits, hitrate, misses, fa), 'center', 'center', p.white);

Screen('Flip', p.window);
while ~keycode(KbName('space'))
    [~, keycode, ~] = KbWait();
end
Priority(0);ListenChar(0);ShowCursor();
Screen('CloseAll');

