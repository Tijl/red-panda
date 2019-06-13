function preprocessing_redpanda(subjectid)

    %% artemis stuff
    if ~ismac
        addpath('../fieldtrip')
        addpath('../CoSMoMVPA/mvpa')
    end

    %% check ft
    if isempty(which('ft_defaults'))
        addpath('~/fieldtrip')
    end
    ft_defaults;

    %% check cosmo
    if isempty(which('cosmo_wtf'))
        addpath('~/cosmomvpa/mvpa')
    end

    %% set data files
    prefix = sprintf('data/sub-%02i/meg/sub-%02i_task-words',subjectid,subjectid);
    fn_events = sprintf('%s_events.csv',prefix);
    fn_events_bids = sprintf('%s_events.tsv',prefix);
    fn_meg = sprintf('%s_meg.con',prefix);
    fn_out = sprintf('data/derivatives/cosmomvpa/sub-%02i_task-words_cosmomvpa.mat',subjectid);

    %% load behavioural file (events)
    eventlist = readtable(fn_events);

    %% load triggers
    fprintf('Reading triggers from %s\n',fn_meg);tic
    triggerdata = ft_read_data(fn_meg,'chanindx',[166,179]);
    trl = find(diff(triggerdata(1,:))>1);
    trl(diff(trl)<100)=[]; %remove possible double points
    fprintf('Finished in %.2fs\n',toc);
    assert(size(eventlist,1)==length(trl),'Number of triggers should match number of events')

    %% load and preproc MEG data 
    fprintf('Loading %s\n',fn_meg);tic
    hdr = ft_read_header(fn_meg);
    cfg = [];
    cfg.trl = [trl-100;trl+1100;repmat(-100,1,length(trl))]';
    cfg.channel = ft_channelselection('meg',hdr.label);
    cfg.dataset = fn_meg;
    ft = ft_preprocessing(cfg);
    cfg=[];
    cfg.resamplefs = 200;
    ft = ft_resampledata(cfg, ft);
    cfg=[];
    cfg.keeptrials='yes';
    cfg.covariance='no';
    ft = ft_timelockanalysis(cfg,ft);
    fprintf('Finished in %.2fs\n',toc);
    
    %% add trigger onsets to eventlist, and write a new (bids)-eventlist
    neweventlist = table();
    neweventlist.onset = trl';
    neweventlist.duration = round(eventlist.stimdur*1000);
    neweventlist = [neweventlist eventlist];
    writetable(neweventlist,fn_events_bids,'FileType','text','Delimiter','\t');
    
    %% convert to cosmo format
    fprintf('Conversion to cosmo\n');tic
    ds = cosmo_meeg_dataset(ft);
    ds.sa = table2struct(eventlist,'toscalar',1);
    cosmo_check_dataset(ds,'meeg');
    fprintf('Finished in %.2fs\n',toc);
    
    %% save
    fprintf('Saving %s\n',fn_out);
    save(fn_out, 'ds', '-v7.3')
    fprintf('Finished in %.2fs\n',toc);
    