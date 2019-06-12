function run_make_rdm(varargin)
    
    if ismac
        if isempty(which('cosmo_wtf'))
            addpath('~/CoSMoMVPA/mvpa')
        end
        nproc = 1;
    else %on HPC
        addpath('../CoSMoMVPA/mvpa');
        % start cluster, give it a unique directory
        % starting a pool can fail when 2 procs are requesting simultaneous
        % thus try again after a second until success
        pool=[];
        while isempty(pool) 
            try
                pc = parcluster('local');
                pc.JobStorageLocation=tempdir;
                pool=parpool(pc);
            catch err
                disp(err)
                delete(gcp('nocreate'));
                pause(1)
            end
        end
        nproc=cosmo_parallel_get_nproc_available();
    end
    
    opt = struct();
    opt = cosmo_structjoin(opt,varargin);
    subjectid = opt.subjectid;
    
    %% load data
    fn_ds = sprintf('data/derivatives/cosmomvpa/sub-%02i_task-words_cosmomvpa.mat',subjectid);
    fn_out = sprintf('results/sub-%02i_RDM.mat',subjectid);
    
    fprintf('Loading %s\n',fn_ds);tic
    load(fn_ds,'ds');
    fprintf('Finished in %.2fs\n',toc);

    %%
    ds = cosmo_slice(ds,~ds.sa.istarget);
    ds.sa.targets = ds.sa.stimnumber;
    ds.sa.chunks = ds.sa.blocknumber;
    nh = cosmo_interval_neighborhood(ds,'time','radius',0);
    ma = struct();
    ma.classifier = @cosmo_classify_lda;
    ma.check_partitions = false;
    ma.output = 'fold_accuracy';
    ma.nproc = nproc;
    
    %% set up partitions
    uc = unique(ds.sa.chunks);
    ut = unique(ds.sa.targets);
    words = {};
    for t=1:length(ut)
        words(t,1) = ds.sa.stim(find(ds.sa.targets==t,1));
    end
    
    combs = nchoosek(ut,2);
    ma.partitions = struct();
    ma.partitions.train_indices = {};
    ma.partitions.test_indices = {};
    sa=struct('target1',[],'target2',[],'leftoutchunk',[]);
    for i=1:size(combs,1) %pairwise comparison to test
        % find the epochs of this pair
        idx1 = ismember(ds.sa.targets,combs(i,:));
        for j=1:length(uc) % chunk to leave out
            idx2 = ds.sa.chunks==uc(j);
            % store targets in results
            sa.target1(end+1,1) = combs(i,1);
            sa.target2(end+1,1) = combs(i,2);
            % store left out chunk in result
            sa.leftoutchunk(end+1,1) = uc(j);
            % set partitions
            ma.partitions.train_indices{1,end+1} = find(idx1 & ~idx2);
            ma.partitions.test_indices{1,end+1} = find(idx1 & idx2);
        end
    end
    
    %% decode
    res = cosmo_searchlight(ds,nh,@cosmo_crossvalidation_measure,ma);
    % merge fold information into result (targets & left out chunk)
    res.sa = cosmo_structjoin(res.sa,sa);
    res.sa.word1 = words(res.sa.target1);
    res.sa.word2 = words(res.sa.target2);
    
    %% save
    fprintf('Saving %s\n',fn_out);
    save(fn_out, 'res', '-v7.3')
    fprintf('Finished in %.2fs\n',toc);
    