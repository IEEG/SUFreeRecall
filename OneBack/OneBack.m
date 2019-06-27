% Run this script first to read clustered spike data from the 1-back
% experiment and produce the xllog array, which contains properties of
% individual units, and is used by subsequent scripts.
% N.B. You need to modify line 390 on openNSx in order to get this to work
% on Unix-based OS.
% Change the home directory before use. This should be the directory under
% which all of the subjects's subdirectories are located.
% Uncomment the section "save the image" in order to produce rasters for
% each unit. Create a subdirectory under home called "Imgs3", or change the
% code as appropriate.
% Runtime: 21 min on authors' hardware, without saving figures.

% Produces data for Figs. 1, S3, S5-S11.

% Copyright Simon Khuvis, 2019. See github.com/IEEG for licensing and use
% information.

t = cputime;

Home = '/Users/simonkhuvis/Documents/MATLAB/SingleUnits';
subjectList = {'S1/VisualLocalizer/S1',...
    'S2/ProcBatch1/S2',...
    'S2/ProcBatch2/S2',...
    'S3/VisLoc/S3',...
    'S4/VisLoc/S4',...
    'S5/VisLoc/S5',...
    'S6/VisLoc/S6',...
    'S7/VisLoc/S7',...
    'S8/VisLoc/S8'};
subj = [1 2 2 3 4 5 6 7 8];
[b,a] = butter(2,[300 3000]/15000);
font = 'Helvetica';
figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);

svgsave = [4 4 1;5 3 1;6 1 1;1 14 2;7 7 1;7 5 3;9 4 3;4 4 1;9 2 3;4 5 2];
%% Count Units
% Number of units in VTC, subtracting number of visually-responsive and
% number of category-selective for Bonferroni-Holm correction. These can be
% obtained by starting at 

% noUnits = 124;
% noVisSel = 124;

% running the code repeatedly, and subtracting as appropriate, until a
% stable value is obtained.

noUnits = 124-63; % All suspected Units
noVisSel = 63-46; % Visually-selective units
%%
for subNo = 1:9
    clearvars -except xllog Home subjectList noUnits noVisSel subNo b a svgsave subj t
    chanID = 0; % For keeping track of which channel is loaded.
    while true
        cd(Home);
        cd(subjectList{subNo}); % Navigate to Subject Directory

        pt = strsplit(subjectList{subNo},'/'); pt = pt{end};

        fls = what();
        batchfile = fls.mat{cell2mat(regexp(fls.mat,['Batch\d*'  pt  '.mat']))};
        if isempty(batchfile)
            disp(['Missing Batch File: ' pt]);
        end
        batchNo = batchfile(6);

        exper = ['VisLocTest' num2str(batchNo)];
        experName = ['1-back test ' num2str(batchNo)];
        subject = ['Subject ' num2str(subj(subNo))];

        load(['./' batchfile])

        alpha = .05; % Significance cutoff

        whichexp = cellfun(@length,regexp(fileList,[exper '\.[-\w]+\.ns6']))==1;
        if sum(whichexp) ~= 1
            error(['Too few/many files with given name in ' batchfile '.']);
        end

        startInds = cumsum(transitions)+1;
        startInd = startInds(whichexp); endInd = startInds(find(whichexp)+1)-1;

        % For unit spreadsheet:
        class = 0;
        burst = false;

    %%

        if exist('unit','var')
            unit = unit + 1;
        else
            unit = 1;
        end

        if ~exist('chanNo','var')
            chanNo = 1;
        end

        if ~exist('stimon','var')
            load(['./' exper '.stimon.mat'])
        end

        if ~exist('codes','var')
            load(['./' exper '.codes.mat'])
        end

        clf
        ref = h5read(['chan' num2str(chanNo) '/sort_neg_simple/sort_cat.h5'],'/groups');
        if(all(ref(2,:)<unit))
            chanNo = chanNo+1;
            if chanNo > chans
                break;
            end
            unit = 1;

            while ~exist(['chan' num2str(chanNo)],'dir') && chanNo <= chans
                chanNo = chanNo + 1;
            end
            
            if chanNo > chans
                break;
            end
            
            ref = h5read(['chan' num2str(chanNo) '/sort_neg_simple/sort_cat.h5'],'/groups');
        end
        refhold = zeros(1,max(ref(1,:))+1);
        refhold(ref(1,:)+1)=ref(2,:);
        ref = refhold;

        loc = loc{(chanNo > 8) + 1};

        spikes = h5read(['chan' num2str(chanNo) '/data_chan' num2str(chanNo) '.h5'],'/neg/times').'*30;

        classhold = zeros(size(spikes));
        index = h5read(['chan' num2str(chanNo) '/sort_neg_simple/sort_cat.h5'],'/index')+1;
        classes = h5read(['chan' num2str(chanNo) '/sort_neg_simple/sort_cat.h5'],'/classes');
        classes = ref(classes+1);
        classhold(index) = classes;
        classes = classhold;

        spikes = spikes(classes==unit);
        spikes = spikes((spikes>=startInd)&(spikes<=endInd)) - startInd + 1;
        if isempty(spikes)
            disp(['No spikes: channel ' num2str(chanNo) ', unit ' num2str(unit) '.'])
            continue;
        end

        %spikes = cluster_class((cluster_class(:,1) == 1),2).'*30;
        numcode = strncmpi('Face',codes,4);
        subplot(11,3,[1 23])

        % Order in which trials are displayed.
        ord = [unique(codes(strncmp(codes,'face',4)));unique(codes(~strncmp(codes,'face',4)))];
        
        %% Plot the spike raster
        hold on
        set(gca,'YDir','reverse')
        rowNo = 0;
        runOrd = zeros(size(stimon));
        rastcode = zeros(size(numcode));
        stimCount = zeros(size(ord));
        
        cats = cellfun(@(x) [upper(x(1)) x(2:end)], unique(cellfun(@(x) x(isstrprop(x,'alpha')), codes,'Uni',0),'stable'),'Uni',0);
        colorcode = zeros(5,3);
        colorcode(strcmpi(cats,'face'),:) = [94.5 23.5 36.1]/100;
        colorcode(~strcmpi(cats,'face'),:) = [33.7 41.6 99.2; 58 62.7 98;16.1 17.3 25.9;7.8 9.8 22]/100;
        
        for stimIt = 1:length(ord)
            stimID = find(strcmp(codes,ord{stimIt}));
            if(stimIt<length(stimCount))
                stimCount(stimIt+1) = length(stimID)+stimCount(stimIt);
            end
            stimCount(stimIt) = length(stimID)/2+stimCount(stimIt);
            for instIt = 1:length(stimID)
                rowNo = rowNo + 1; % Which row to plot on the raster
                runNo = stimID(instIt); % Which run number
                runOrd(rowNo) = runNo; % Keep track of which trial
                rast = (spikes((spikes>(stimon(runNo)-1.5e4)) & (spikes<(stimon(runNo)+2.9e4)))-...
                    stimon(runNo))/3e4;
                
                if ~isempty(rast) 
                    plot([1;1]*rast,[rowNo-.85 rowNo-.15],'Color',...
                        colorcode(strncmpi(ord{stimIt},cats,4),:),'LineWidth',3)
                end
                rastcode(rowNo) = find(strncmpi(ord{stimIt},cats,4),1);
            end
        end
        plot([0 0],[0 length(stimon)],'k--')
        plot([0.5 0.5],[0 length(stimon)],'k--')
        axis([-.5 2.9/3 0 length(stimon)])
        set(gca,'YTick',stimCount)
        set(gca,'YTickLabel',cellfun(@(x) x(end-1:end),ord,'Uni',0))
        title([subject ': ' experName ', channel ' num2str(chanNo) ' (' loc '), unit ' num2str(unit)])
        ylabel('Stimulus')
        set(gca,'box','off')
        set(gca,'TickLength',[.002 .002])
        
        %% Is it an ON or OFF cell?
        
        addme = -1.5e4:0.1e4:2.9e4;

        intervals = bsxfun(@plus,repmat(stimon,length(addme),1),addme.');
        spikecount = @(val) sum(spikes<=val);
        map = diff(arrayfun(spikecount,intervals)).';
        
        timecourse = {30*map(numcode,:);30*map(~numcode,:)};

        subplot(11,3,[1 23])
        
        onmap = downsample(diff(...
            arrayfun(spikecount,reshape([stimon+0.3e4;stimon+1.5e4],1,2*length(stimon)))),2).';

        offmap = downsample(diff(...
            arrayfun(spikecount,reshape([stimon+1.8e4;stimon+3.0e4],1,2*length(stimon)))),2).';

        premap = downsample(diff(...
            arrayfun(spikecount,reshape([stimon-1.2e4;stimon],1,2*length(stimon)))),2).';
        
        if max(sum(offmap>0),sum(onmap>0)) < 8
            disp(['Insufficient spikes: channel ' num2str(chanNo) ', unit ' num2str(unit) '.'])
            continue;
        end
            
        
        [pOn,statsOn,stOn]  = friedman(abs(diff(map(:,16:29),1,2)),1,'off'); % Responsiveness of firing rate to image presentation
        [pOff,statsOff,stOff] = friedman(abs(diff(map(:,31:44),1,2)),1,'off'); % Responsiveness of firing rate to image offset
        
        if (pOn < .05/noUnits/2) % Is it an ON cell?
            fill([.1 .5 .5 .1],[0 0 length(stimon)*[1 1]],[.92 .92 .92],'EdgeColor','None');
            set(gca,'children',flipud(get(gca,'children')))
            ti = 'Firing rate, 100 to 500 ms post-stimulus';
            tt = 'ON Type ';
            instmap = onmap;
            oncell = true;
            statsFin = statsOn;
            
            class = 1;
            [p,~,stat] = kruskalwallis(onmap,cellfun(@(x) x(isstrprop(x,'alpha')), codes,'Uni',0),'off');
            [pAmp,~,statAmp] = kruskalwallis(abs(log((onmap+.1)./(premap+.1))),cellfun(@(x) x(isstrprop(x,'alpha')), codes,'Uni',0),'off');
        elseif (pOff < .05/noUnits/2) % Is it an OFF cell?
            fill([.6 1 1 .6],[0 0 length(stimon)*[1 1]],[.92 .92 .92],'EdgeColor','None');
            set(gca,'children',flipud(get(gca,'children')))
            ti = 'Firing rate, 600 to 1000 ms post-offset';
            tt = 'OFF Type ';
            instmap = offmap;
            oncell = false;
            statsFin = statsOff;
            
            class = 1;
            [p,~,stat] = kruskalwallis(offmap,cellfun(@(x) x(isstrprop(x,'alpha')), codes,'Uni',0),'off');
            [pAmp,~,statAmp] = kruskalwallis(abs(log((offmap+.1)./(onmap+.1))),cellfun(@(x) x(isstrprop(x,'alpha')), codes,'Uni',0),'off');
        else
            fill([.1 .5 .5 .1],[0 0 length(stimon)*[1 1]],[.92 .92 .92],'EdgeColor','None');
            set(gca,'children',flipud(get(gca,'children')))
            ti = 'Firing rate, 100 to 500 ms post-stimulus';
            tt = 'Non-responsive or artifactual';
            instmap = onmap;
            oncell = true;
            statsFin = statsOn;
            
            class = 0;
            catsel = 0;
        end
        
        responsiveness = statsFin(2,5);
        supersig = false;
        pval = nan;
        if class % If it's visually responsive...
            
            if p < alpha/noVisSel %Check if it's category-selective
                
                mcMat = multcompare(stat,'Alpha',alpha,'Display','off','CType','lsd');
                ampMat = multcompare(statAmp,'Alpha',alpha,'Display','off','CType','lsd');

                % Does each class have a greater or lower firing rate than
                % every other?
                gtan = arrayfun(@(x) all(([-mcMat(mcMat(:,2)==x,4);mcMat(mcMat(:,1)==x,4)]*[1 -1])>0)*[1;-1],1:5);
                % Is it significant?
                sigtan = arrayfun(@(x) all([mcMat(mcMat(:,2)==x,6);mcMat(mcMat(:,1)==x,6)]<alpha), 1:5);
                
                % Does any class have a greater magnitude change in firing
                % rate than every other?
                gtanAmp = arrayfun(@(x) all(([-ampMat(ampMat(:,2)==x,4);ampMat(ampMat(:,1)==x,4)])>0),1:5);
                
                gsigtan = gtan.*sigtan;
                
                if sum(abs(gsigtan))>1
                    cds = unique(cellfun(@(x) x(isstrprop(x,'alpha')), codes,'Uni',0),'stable');
                    inds1 = strcmp(cellfun(@(x) x(isstrprop(x,'alpha')), codes,'Uni',0),cds(find(abs(gsigtan),1,'first')));
                    inds2 = strcmp(cellfun(@(x) x(isstrprop(x,'alpha')), codes,'Uni',0),cds(find(abs(gsigtan),1,'last')));
                    if oncell
                        pgt = ranksum(abs(log((onmap(inds1)+.1)./(premap(inds1)+.1))),abs(log((onmap(inds2)+.1)./(premap(inds2)+.1))));
                        isgt = (median(abs(log((onmap(inds1)+.1)./(premap(inds1)+.1))))<median(abs(log((onmap(inds2)+.1)./(premap(inds2)+.1))))) + 1;
                    else
                        pgt = ranksum(abs(log((offmap(inds1)+.1)./(onmap(inds1)+.1))),abs(log((offmap(inds2)+.1)./(onmap(inds2)+.1))));
                        isgt = (median(abs(log((offmap(inds1)+.1)./(onmap(inds1)+.1))))<median(abs(log((offmap(inds2)+.1)./(onmap(inds2)+.1))))) + 1;
                    end
                    if pgt < .05
                        gsigtmp = find(abs(gsigtan));
                        gsigtmp = gsigtmp(isgt);
                        gsigtan = 0*gsigtan; gsigtan(gsigtmp) = gtan(gsigtmp);
                    else
                        gsigtan = 0*gsigtan;
                    end
                end
                
                pval = max([mcMat(mcMat(:,2)==find(gsigtan),6);mcMat(mcMat(:,1)==find(gsigtan),6)]);

                ei = {'Inhibited' 'Excited'};
                
                if any(gsigtan)
                    tt = [tt cats{gsigtan~=0} '-Selective ' ei{(sum(gsigtan)>0)+1}];
                    class = 3;
                    catsel = find(gsigtan);
                    supersigtan = arrayfun(@(x) all([mcMat(mcMat(:,2)==x,6);mcMat(mcMat(:,1)==x,6)]<.001), 1:5);
                    supersigtanAmp = arrayfun(@(x) all([mcMat(mcMat(:,2)==x,6);mcMat(mcMat(:,1)==x,6)]<.001), 1:5);
                    supersig = any(gsigtan&supersigtan&supersigtanAmp);
                else
                    tt = [tt 'Category-Selective Unclassified'];
                    class = 2;
                    catsel = 0;
                end
            else
                tt = [tt ' Visually-Responsive Unclassified'];
                class = 1;
                catsel = 0;
            end
        end
        %% Plot the average spike frequency per trial.

        axinst = subplot(11,3,[3 24]);
        instlim = max(3*instmap/(1.2))+1;
        presno = (1:size(instmap,1))-.5;
        for classIt = 1:5
            barh(presno(rastcode==classIt),3*instmap(runOrd(rastcode==classIt))/(1.2),1,...
                'facecolor',colorcode(classIt,:),'EdgeColor','none','BaseValue',-.05*instlim);
            hold on
        end
        set(gca,'Ydir','reverse')
        set(gca,'YTick',stimCount)
        set(gca,'YTickLabel',cellfun(@(x) x(end-1:end),ord,'Uni',0))
        set(gca,'box','off')
        set(gca,'TickLength',[.002 .002])
        ylim([0 size(instmap,1)])
        xlim(instlim*[-.05 1.1])
        title(ti)
        
        %% Calculate the response latency
        
        if class>0
            if catsel
                catofint = stat.gnames{catsel};
            else
                catofint = 'face';
            end
            
            % Intervals of .02
            addme2 = -1.5e4:0.06e4:2.9e4;

            intervals2 = bsxfun(@plus,repmat(stimon,length(addme2),1),addme2.');
            map2 = diff(arrayfun(spikecount,intervals2)).';
            
            profile = mean(map2(strcmp(cellfun(@(x) x(isstrprop(x,'alpha')),...
                codes,'Uni',0),catofint),:));
            if oncell
                latency = find(conv(profile(26:end)-mean(profile(11:25))>...
                    norminv(.95)*std(profile(11:25)),[1 1 1],'valid')==3,1)*.02+.01;
            else
                latency = find(conv(profile(26:end)-mean(profile(11:25))<...
                    -norminv(.95)*std(profile(11:25)),[1 1 1],'valid')==3,1)*.02+.01;
            end
            
            if isempty(latency)
                latency = NaN;
            end
        else
            latency = NaN;
        end

        %% Peri-stimulus time histograms
        
        if catsel
            plotcat = catsel;
            othercol = [24.7 30.2 71.8]/100;
        else
            plotcat = find(strncmp('Face',cats,4),1);
            othercol = mean([24.7 30.2 71.8;94.5 23.5 36.1]/100);
        end
        plotcode = strncmpi(cats{plotcat},codes,length(cats{plotcat}));
        
        ftcourse = subplot(11,3,[25 26 28 29]);
        meancourse = @(inst) 30*sum(map(inst,:))/(sum(inst));
        semcourse = @(inst) 30*std(map(inst,:))/(sum(inst))^.5;
        
        xvals = (addme(2:end)-(addme(2)-addme(1))/2)/3e4;
        f1 = fill([xvals fliplr(xvals)],[meancourse(plotcode) + semcourse(plotcode)...
            fliplr(meancourse(plotcode) - semcourse(plotcode))],colorcode(plotcat,:),'EdgeColor','none');
        set(f1,'FaceAlpha',.1)
        hold on
        plot(xvals,meancourse(plotcode),'Color',colorcode(plotcat,:))
        xlim([addme(1) addme(end)]/3e4)
        
        
        ylabel('Firing rate [Hz]')

        f2 = fill([xvals fliplr(xvals)],[meancourse(~plotcode) + semcourse(~plotcode)...
            fliplr(meancourse(~plotcode) - semcourse(~plotcode))],othercol,'EdgeColor','none');
        set(f2,'FaceAlpha',.1)
        plot(xvals,meancourse(~plotcode),'Color',othercol)
        xlabel('Time [s]')

        yzoom = zeros(2,2);
        yzoom(1,:) = get(ftcourse,'YLim');
        
        ylim([0 max(yzoom(:,2))]);
        plot([0 0],[0 max(yzoom(:,2))],'k--')
        plot([0.5 0.5],[0 max(yzoom(:,2))],'k--')
        
        colbreak = @(mat) mat2cell(mat,size(mat,1),ones(1,size(mat,2)));
        
        sigtimes = cellfun(@ranksum,colbreak(map(numcode,:)),colbreak(map(~numcode,:)))<.05/length(xvals);
        
        set(gca,'box','off')
        set(gca,'TickLength',[.002 .002])

        %%

        instlim = max(3*instmap/(1.2))+1;
        
        if catsel
            plotcat = catsel;
        else
            plotcat = find(strncmp('Face',cats,4),1);
        end

        axbox = subplot(11,3,[27 30]);
        set(ftcourse,'box','off')
        hold on
        niceax = get(axinst,'Position');
        oldax = get(axbox,'Position');
        boxplot(axbox,instmap,strncmpi(codes,cats{plotcat},4)+1,'orientation','horizontal','colors',[othercol;colorcode(plotcat,:)],'symbol','+')
        set(axbox,'Position',[niceax(1) oldax(2) niceax(3) oldax(4)]);
        
        set(ftcourse,'box','off')

        scale = max(instmap);
        
        if supersig
            plot([1.4 1.6]*scale,[1.5 1.5],'k*','MarkerSize',10,'LineWidth',1)
            alphahold = .001;
            plot([(max(instmap(strncmpi(codes,cats{plotcat},4))) + .2*scale) 1.5*scale*[1 1]],[2 2 1.8],'k','LineWidth',1)
            plot([(max(instmap(~strncmpi(codes,cats{plotcat},4))) + .2*scale) 1.5*scale*[1 1]],[1 1 1.2],'k','LineWidth',1)
            sigcrit = '\it p \rm < 0.001';
        elseif(class == 3) 
            plot(1.5*scale,1.5,'k*','MarkerSize',10,'LineWidth',1)
            plot([(max(instmap(strncmpi(codes,cats{plotcat},4))) + .2*scale) 1.5*scale*[1 1]],[2 2 1.8],'k','LineWidth',1)
            plot([(max(instmap(~strncmpi(codes,cats{plotcat},4))) + .2*scale) 1.5*scale*[1 1]],[1 1 1.2],'k','LineWidth',1)
            sigcrit = '\it p \rm < 0.05';
        else
            sigcrit = 'n.s.';
        end
        text(1.1*instlim,.5,sigcrit,...
            'HorizontalAlignment','right','VerticalAlignment','bottom','Interpreter','tex')
        
        xlim(instlim*[-.05 1.1])
        xlabel('Unit firing rate [Hz]')
        set(gca,'YTick',[1 2])
        set(gca,'YTickLabel',{'Other',cats{plotcat}})
        ylabel('Stimulus category')
        set(gca,'box','off')
        
         %% Plot the ISI and spike waveform
        
        %%%%% FOR SUBJECT 4 %%%%%
        if subNo == 4
            spikes = spikes + 7.5e6;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%
        
        isi1 = diff(spikes((spikes>=stimon(1)-1.5e4)&(spikes<=stimon(end)+3e4))/30);
        
        yl = 50;
        
        xinds1 = linspace(0,min(2.5*median(isi1),yl),30);
        subplot(11,6,[61 62])
        bar(xinds1+(xinds1(2)-xinds1(1))/2,histc(isi1,xinds1),1,'k')
        xlabel('Time [ms]')
        ylabel('Count [spikes]')
        title('Inter-spike interval')
        set(gca,'box','off')
        set(gca,'TickLength',[.002 .002])
        xlim([min(xinds1) max(xinds1)])
        
        if chanID ~= chanNo
            openNSx(['/' fileList{whichexp}],'read',['c:' num2str(chanNo)])
            d = filtfilt(b,a,double(NS6.Data));
            chanID = chanNo;
        end
        
        subplot(11,3,33)
        
        smooth = d(bsxfun(@plus,round(spikes((spikes>=stimon(1)-1.5e4)&(spikes<=stimon(end)+3e4))'),-30:60));
        
        plot((-30:60)/30,smooth','Color',[.6 .6 .6])
        hold on
        plot((-30:60)/30,mean(smooth),'Color','k')
        title('Mean spike waveform')
        xlabel('Time [ms]')
        ylabel('Potential [\muV]')
        axis([-1 2 2*prctile(min(smooth,[],2),20) 2*prctile(max(smooth,[],2),80)])
        set(gca,'box','off')
        set(gca,'TickLength',[.002 .002])
        %%
        % Calculate the log ISI:
        isi = log10(diff(spikes((spikes>=stimon(1)-1.5e4)&(spikes<=stimon(end)+3e4))/30));
        xinds = (min(isi)-.1):.1:prctile(isi,98);
        isiH = conv(histc(isi,xinds),[.2 .2 .2 .2 .2],'same');
        % Correct for edge effects.
        isiH(1:2) = isiH(1:2).*[5/3 5/4]; isiH(end-1:end) = isiH(end-1:end).*[5/4 5/3];

        [ypk,xpk] = findpeaks(isiH,'MinPeakDistance',2);
        if (length(xpk)<2) || (~any(xinds(xpk)<2))
            burst = false;
        else
            for xIt = (length(xpk)-1)
                gmin = min(isiH(xpk(xIt):xpk(xIt+1)));
                voidVar = 1-gmin/prod(ypk(xIt:xIt+1));
                if voidVar>.7
                    burst = true;
                end
            end
        end
        %% Log the data
        cats = ['None';cats]; %#ok<AGROW>
        codeNums = 100*cellfun(@(x) find(cellfun(@(y) all(strcmp(y,x(isletter(x)))),{'face','body','house','patterns','tool'})), codes)...
            +cellfun(@(x) str2double(x(end-1:end)),codes);

        batchname = strsplit(batchfile,'.'); batchname = batchname{1};
        if ~exist('xllog','var')
            xllog = struct('pt',pt,'batchname',batchname,'loc',loc,'chanNo',chanNo,'unit',unit,'tt',tt,'oncell',oncell,...
                'class',class,'catsel',cats{catsel+1},'burst',burst,'response',[codeNums onmap offmap],'VisResponsiveness',responsiveness,'latency',latency,'timecourse',{timecourse},'p',pval);
        else
            alreadyin = find(arrayfun(@(x) isequal(x.batchname,batchname) && (x.chanNo==chanNo) && (x.unit == unit),xllog),1,'last');
            if ~isempty(alreadyin)
                xllog(alreadyin) = struct('pt',pt,'batchname',batchname,'loc',loc,'chanNo',chanNo,'unit',unit,'tt',tt,'oncell',oncell,...
                    'class',class,'catsel',cats{catsel+1},'burst',burst,'response',[codeNums onmap offmap],'VisResponsiveness',responsiveness,'latency',latency,'timecourse',{timecourse},'p',pval);
            else
                xllog = [xllog; struct('pt',pt,'batchname',batchname,'loc',loc,'chanNo',chanNo,'unit',unit,'tt',tt,'oncell',oncell,...
                    'class',class,'catsel',cats{catsel+1},'burst',burst,'response',[codeNums onmap offmap],'VisResponsiveness',responsiveness,'latency',latency,'timecourse',{timecourse},'p',pval)...
                    ]; %#ok<AGROW>
            end
        end
        
        %% Save the image
        cd([Home '/Imgs3']);

%         if class
%             set(gcf,'PaperUnits','inches','PaperPosition',[0 0 1200 1600]/100)
%         %     set(gcf,'PaperUnits','inches','PaperPosition',[0 0 1442 804]/70)
%             print([pt exper '_Channel' num2str(chanNo) 'Unit' num2str(unit) '-' num2str(noUnits)],'-dpng','-r100')
%         end
%         
%         if ismember([subNo chanNo unit],svgsave,'rows')

%             plot2svg([pt exper '_Channel' num2str(chanNo) 'Unit' num2str(unit) '-' num2str(noUnits) '.svg'])
%         end
        
        cd(Home);
    end
end

OneBackTime = cputime-t;
disp(['Runtime: ' num2str(OneBackTime)])