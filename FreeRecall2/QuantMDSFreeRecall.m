% Run this script to produce Figs. 4B,C. Load output from
% FreeRecallPreprocessing.m and RecallSpikes.mat before running this code.
% Runtime: ~18 hrs on authors' hardware.

% You can also visualize the output without rerunning the bootstrap
% analysis. Look for lines in the code labeled "$$$" to comment or skip.

% Copyright Simon Khuvis, 2019. See github.com/IEEG for licensing and use
% information.

t = cputime;

%% Extract spike rates from recall period.

global indx
indx = @(vec,ind) vec(ind,:);

% Names of all the batches
batches = indx(unique({xllog.batchname},'stable')',[1 3 4 5 6])';

% List of images, sorted within presentation blocks.
catshelp = @(inp) cell2mat(cellfun(@sort,mat2cell(inp,14+0*(1:length(inp)/14)',1),'Uni',0));
cats = @(pt) catshelp(unique(cell2mat(arrayfun(@(x) x.response(:,1),xllog(strcmp({xllog.batchname},batches{pt})),'Uni',0)),'stable'));
% If you only want one presentation block.
runcensor = @(cat,r) cat(14*(r-1)+1:14*r,:);

expand = @(x) [x{:}];

% Extract spikes from four (or two) recall periods.
spikes = arrayfun(@(x,z) expand(arrayfun(@(y) y.spikes(:,[xllog(strcmp({xllog.batchname},z)).class]>0)',...
    expand(struct2cell(x)),'Uni',0)),RecallSpikes([1 3 4 5 6]),batches,'Uni',0);

% Exclude "prompts".
eventshelp = @(x) indx(x',strncmpi({x.ref},'Face',4)|strncmpi({x.ref},'Place',5))';
% Image IDs from recall events.
events = arrayfun(@(x) expand(arrayfun(@(y) cellfun(eventshelp,y.events,'Uni',0),...
    expand(struct2cell(x)),'Uni',0)),RecallSpikes([1 3 4 5 6]),'Uni',0);
% Spike rates from the 2s prior to vocalization.
sigs = cellfun(@(ev,sp) cellfun(@(e,s) arrayfun(@(x) histc(s,x.start+3e4*[-2 0])*[1;0]/2,e'),...
    repmat(ev,size(sp,1),1),sp,'Uni',0),events,spikes,'Uni',0);

% Get the responses out of the xllog variables.
extractLog = @(xlog) arrayfun(@(x) x.response(:,[1 3-x.oncell]),xlog([xlog.class]>0),'Uni',0);

justlog = @(extracted) cellfun(@(x) [x(:,1) log10(x(:,2)+1)],extracted,'Uni',0);

counts = arrayfun(@(pt) min(cell2mat(cellfun(@(x) histc(indx(x',1),sort(cats(pt))),...
    extractLog(xllog(strcmp({xllog.batchname},batches{pt}))),'Uni',0)))',1:5,'Uni',0);
minct = min(cellfun(@min,counts));%3;

picksome = @(inp,num)indx(inp,randperm(length(inp),num));

logxl = cellfun(@(x) justlog(extractLog(xllog(strcmp({xllog.batchname},x)))),batches,'Uni',0);
logxl = cellfun(@(pt) cellfun(@(ut) [ut(:,1) ut(:,2)-mean(indx(cell2mat(pt')',2:2:2*length(pt)))'],pt,'Uni',0),logxl,'Uni',0);

extract = @(catinp,cld,pt) cell2mat(cellfun(@(x) x(x(:,1) == catinp,2)',cld{pt},'Uni',0))';
cmdsubset = @(inp) cmdscale(pdist(inp,'euclidean'),2);
normalize2 = @(inp) bsxfun(@rdivide,(bsxfun(@minus,inp,mean(inp,2))),std(inp,0,2)+1);
xform = @(ex1,ex2) mat2cell(cmdsubset(normalize2([ex1;ex2])),...
    [size(ex1,1) size(ex2,1)],2);
w = @(inp,n1,n2) squeeze(sum(cell2mat(permute(cellfun(@(x,y) cov(x(y,:)),inp,{n1';n2'},'Uni',0),[3 2 1])),3))^-1*...
    (mean(inp{1}(n1,:))-mean(inp{2}(n2,:)))';
normz = @(inp) inp./norm(inp);
cent = @(inp,n1,n2) mean([mean(inp{1}(n1,:));mean(inp{2}(n2,:))])';
evl = @(inp,n1,n2) cellfun(@(x,c,wx,n,z) (z*dot(x(setdiff(1:length(x),n),:)-c',normz(wx)'))>0,...
    inp,repmat({cent(inp,n1,n2)},2,1),repmat({w(inp,n1,n2)},2,1),{n1;n2},{1;-1},'Uni',0);
alltogether = @(r,c,cld,pt,rn) mean(cell2mat(evl(xform(extract(indx(runcensor(cats(pt),rn),r),cld,pt),...
    extract(indx(runcensor(cats(pt),rn),c),cld,pt)),1:minct-1,1:minct-1)));

%% Fill in confusion matrix for each subject
% Don't run this section when running the code from saved output. $$$

reps = 100;
superreps = 100;
MAT = cell(1,1,superreps);
rng(1)
parfor sr = 1:superreps
    tic
    clIn = arrayfun(@(z) cellfun(@(pt,y) cellfun(@(x) cell2mat(arrayfun(@(ct) picksome(x(x(:,1) == ct,:),minct),...
        cats(pt),'Uni',0)),y,'Uni',0),num2cell(1:5),logxl,'Uni',0),1:reps,'Uni',0);
    MAT{sr} = arrayfun(@(pt) arrayfun(@(rn) [nan(1,14,reps);cell2mat(arrayfun(@(row) ...
        [cell2mat(arrayfun(@(col) permute(cellfun(@(x) alltogether(row,col,x,pt,rn),clIn),[1 3 2]),1:row-1,'Uni',0)), ...
        nan(1,14-row+1,reps)], (2:14)','Uni',0))],1:sum(indx(histc(cats(pt),[100 110 210 300]),[1 3])>[0 0]),'Uni',0),1:5,'Uni',0);
    toc
end

MAT2 = mat2cell(reshape(expand(MAT),[],superreps)',superreps,[1 1 1 1 1]);
MAT2 = cellfun(@(x) mat2cell(permute(reshape(expand(x),[],superreps),[3 1 2]),1,ones(1,size(x{1},2)),superreps),MAT2,'Uni',0);
MAT2 = cellfun(@(x) cellfun(@cell2mat,x,'Uni',0) ,MAT2,'Uni',0);
acc = cellfun(@(x) cellfun(@(y) cellfun(@(z) ...
    nanmean(reshape(z,1,1,[])),mat2cell(y,[7 7],[7 7],superreps*reps),'Uni',1),x,'Uni',0),MAT2,'Uni',0);
MAT = MAT2;

%% Surrogate confusion matrix
% Don't run this section when running the code from saved output. $$$

% Scramble stimuli in all extracted responses (to build surrogate
% distributions) within categories.
scramLog = @(extracted,pt) cellfun(@(z) cell2mat(cellfun(@(y) [y(randperm(size(y,1)),1) y(:,2)], ...
    arrayfun(@(x) z(floor(z(:,1)/10)==x,:),unique(floor(cats(pt)/10))','Uni',0),'Uni',0).'),extracted,'Uni',0);

% Inter-category scramble
scramLogInter = @(extracted,pt) cellfun(@(z) cell2mat(cellfun(@(y) [y(randperm(size(y,1)),1) y(:,2)], ...
    arrayfun(@(x) z(mod(floor(z(:,1)/10),10)==x,:),unique(mod(floor(cats(pt)/10),10))','Uni',0),'Uni',0).'),extracted,'Uni',0);

SCRlogxl = @(dummy,pt) justlog(scramLog(extractLog(xllog(strcmp({xllog.batchname},batches{pt}))),pt));
SCRlogxlInter = @(dummy,pt) justlog(scramLogInter(extractLog(xllog(strcmp({xllog.batchname},batches{pt}))),pt));
nrmz = @(inp) cellfun(@(ut) [ut(:,1) ut(:,2)-mean(indx(cell2mat(inp')',2:2:2*length(inp)))'],inp,'Uni',0);

SCRextract = @(catinp,scl) cell2mat(cellfun(@(x) picksome(x(x(:,1) == catinp,2),minct)',scl,'Uni',0))';
SCRalltogether = @(r,c,scl,pt,rn) mean(cell2mat(evl(xform(SCRextract(indx(runcensor(cats(pt),rn),r),scl),...
    SCRextract(indx(runcensor(cats(pt),rn),c),scl)),1:minct-1,1:minct-1)));

rng(1)
reps = 100;
superreps = 100;
scrMAT = cell(1,1,superreps);
parfor sr = 1:superreps
    tic
    scrvv = arrayfun(@(pt) arrayfun(@(rn) arrayfun(@(x) nrmz(SCRlogxl(x,pt)), 1:reps,'Uni',0), ...
        1:sum(indx(histc(cats(pt),[100 110 210 300]),[1 3])>[0 0]),'Uni',0),1:5,'Uni',0);
    scrMAT{sr} = arrayfun(@(pt) arrayfun(@(rn) [nan(1,14,reps);cell2mat(arrayfun(@(row) [nan(1,7*(row>7),reps) ...
        cell2mat(arrayfun(@(col) permute(cellfun(@(scrv) SCRalltogether(row,col,scrv,pt,rn),...
        scrvv{pt}{rn}),[1 3 2]),7*(row>7)+1:row-1,'Uni',0)), nan(1,14-row+1,reps)], ...
        (2:14)','Uni',0))], 1:sum(indx(histc(cats(pt),[100 110 210 300]),[1 3])>[0 0]),'Uni',0),1:5,'Uni',0);
    
    scrvv = arrayfun(@(pt) arrayfun(@(rn) arrayfun(@(x) nrmz(SCRlogxlInter(x,pt)), 1:reps,'Uni',0), ...
        1:sum(indx(histc(cats(pt),[100 110 210 300]),[1 3])>[0 0]),'Uni',0),1:5,'Uni',0);
    for ptx = 1:5
        for rnx = 1:sum(indx(histc(cats(ptx),[100 110 210 300]),[1 3])>[0 0])
            scrMAT{sr}{ptx}{rnx}(8:14,1:7,:) = cell2mat(arrayfun(@(row) cell2mat(arrayfun(@(col) ...
                permute(cellfun(@(scrv) SCRalltogether(row,col,scrv,ptx,rnx),scrvv{ptx}{rnx}),[1 3 2]),1:7,'Uni',0)),(8:14)','Uni',0));
        end
    end
    toc
end

MAT2 = mat2cell(reshape(expand(scrMAT),[],superreps)',superreps,[1 1 1 1 1]);
MAT2 = cellfun(@(x) mat2cell(permute(reshape(expand(x),[],superreps),[3 1 2]),1,ones(1,size(x{1},2)),superreps),MAT2,'Uni',0);
MAT2 = cellfun(@(x) cellfun(@cell2mat,x,'Uni',0) ,MAT2,'Uni',0);
acc = cellfun(@(x) cellfun(@(y) cellfun(@(z) ...
    nanmean(reshape(z,1,1,[])),mat2cell(y,[7 7],[7 7],superreps*reps),'Uni',1),x,'Uni',0),MAT2,'Uni',0);
scrMAT = MAT2;
%% Calculate Classifier Accuracy
rfmt = @(inp) [inp(1,1,:);inp(2,2,:);inp(2,1,:)];

acc = cellfun(@(x) cell2mat(cellfun(@(y) rfmt(cellfun(@(z) ...
    nanmean(reshape(z,1,1,[])),mat2cell(y,[7 7],[7 7],size(MAT{1}{1},3)),'Uni',1)),x,'Uni',0)),MAT,'Uni',0);
acc = [acc(1) {[acc{2:3}]} acc(4:5)];
acc = cell2mat(cellfun(@(x) mean(x,2),acc,'Uni',0));

accdistMake = @(inp) cellfun(@(pt) cellfun(@(rn) cellfun(@(x) nanmean(reshape(x,1,[])),mat2cell(rn, ...
    [7 7],[7 7],ones(1,size(MAT{1}{1},3)))),pt,'Uni',0),inp,'Uni',0);
scrDist = accdistMake(scrMAT); accDist = accdistMake(MAT);
scrDist = [scrDist(1) {[scrDist{2:3}]} scrDist(4:5)];
accDist = [accDist(1) {[accDist{2:3}]} accDist(4:5)];
comb = @(inp) cellfun(@(pt) mean(cell2mat(cellfun(@(rn) rfmt(rn),pt,'Uni',0)),2),inp,'Uni',0);
scrDist = comb(scrDist); accDist = comb(accDist);

sig = cell2mat(cellfun(@(ac,sr) median(cell2mat(cellfun(@(x) mean(bsxfun(@gt,x,sr),3),...
    mat2cell(ac,3,1,ones(1,1,size(ac,3))),'Uni',0)),3),accDist,scrDist,'Uni',0));

% Bonferroni-Holm correction
compold = sig>(1-.05/12);
comp = nan;
while compold ~= comp
    comp = compold;
    comp = sig>(1-.05/sum(sum(~compold)));
end

subplot(5,1,[2 3])
b = bar([2 1 3 4],100*acc','BaseValue',.5,'EdgeColor','none');
b(1).FaceColor = [94.5 23.5 36.1]/100;
b(2).FaceColor = [24.7 30.2 71.8]/100;
b(3).FaceColor = [55 75 35]/100;
set(gca,'XTickLabel',[2 3 6 8])
ylim([50 100])
ylabel('Cross-validation accuracy (%)','FontSize',10)
hold on
plhold = bsxfun(@plus,[-.23 0 .23]',[2 1 3 4]);
plot(plhold(comp),100*mean([1 max(max(acc))]),'k*','MarkerSize',5)
set(gca,'box','off')
set(gca,'TickLength',[.01 .01])


%% Recall Exemplar Decoding
L2 = @(inp) log10(inp+1);

extr = @(sg,en,run) cellfun(@(x) x(en), sg(:,run)); % Take one event from one run at a time.
ref = @(inp) 100*find(strcmpi(inp.ref(isstrprop(inp.ref,'alpha')),{'Face','Place'}))+...
    str2double(inp.ref(isstrprop(inp.ref,'digit')))/10; % Extract id of stimulus.

% Put first block and second block together (1st and 3rd, 2nd and 4th, when
% struct2mat is used).
recmerge = @(inp) cellfun(@(x) [inp{x}],...
    mat2cell(reshape(1:length(inp),length(inp)/2,2),ones(length(inp)/2,1),2),'Uni',0);
% Remove recall events that are from the other block (even if they occured
% in the wrong block).
reccensor = @(rec,pt,r) rec(ismember([rec.event],runcensor(cats(pt),r)));
% Create struct array for each recall event and neuronal activity.
recall = cellfun(@(ev,sg,pt) recmerge(arrayfun(@(run) reccensor(arrayfun(@(en) struct(...
    'event',ref(extr(ev,en,run)),'sig',L2(extr(sg,en,run))),(1:length(ev{run}))),pt,mod(run-1,length(ev)/2)+1),...
    1:length(ev),'Uni',0)),events,sigs,num2cell(1:5),'Uni',0);

recall = cellfun(@(pt) cellfun(@(rn) arrayfun(@(ev) struct('event',ev.event,'sig',ev.sig-mean(ev.sig)),rn),pt,'Uni',0),recall,'Uni',0);

% Spike counts: each row is an exemplar; each column is a unit.
SCpres = cellfun(@(bn,bno) cell2mat(cellfun(@(res) accumarray(arrayfun(@(x) find(x==cats(bno),1),res(:,1)),res(:,2))./ ...
    arrayfun(@(x) sum(x==res(:,1)),cats(bno)),{xllog(strcmp({xllog.batchname},bn)).response},'Uni',0)),batches,num2cell(1:5),'Uni',0);

% Limit cats to the type of ex:
typecensor = @(inp,ex) inp(ismember(floor(inp/100),floor(ex/100)));

extract = @(catinp,pt) cell2mat(cellfun(@(x) x(x(:,1) == catinp,2)',logxl{pt},'Uni',0))';
normalize2 = @(inp) bsxfun(@rdivide,(bsxfun(@minus,inp,mean(inp,2))),std(inp,0,2)+1);
cmdsubset = @(inp) cmdscale(pdist(inp,'euclidean'),2);
xform = @(ex1,ex2,n1,n2) mat2cell(cmdsubset(normalize2([ex1;ex2])),...
    [size(ex1,1) size(ex2,1)],2);
w = @(inp,n1,n2) squeeze(sum(cell2mat(permute(cellfun(@(x,y) cov(x(y,:)),inp,{n1';n2'},'Uni',0),[3 2 1])),3))^-1*...
    (mean(inp{1}(n1,:))-mean(inp{2}(n2,:)))';
normz = @(inp) inp./norm(inp);
cent = @(inp,n1,n2) mean([mean(inp{1}(n1,:));mean(inp{2}(n2,:))])';
evl = @(inp,n1,n2) cellfun(@(x,c,wx,n,z) (z*dot(x(setdiff(1:length(x),n),:)-c',normz(wx)'))>0,...
    inp(1),{cent(inp,n1,n2)},{w(inp,n1,n2)},{n1},{1});

alltogether = @(rec,p,n,pt) arrayfun(@(px,nx) evl(xform([extract(px,pt);...
    rec.sig'],extract(nx,pt),1:4,1:4),1:4,1:4),p,n);

global allevlHelp
allevlHelp = @(pt,run,ind,ex) indx(typecensor(runcensor(cats(pt),run),ex),ind);
allevl = @(rec,pt,run) [nan(1,7);cell2mat(arrayfun(@(n) [alltogether(rec,allevlHelp(pt,run,1:(n-1),rec.event)', ...
    allevlHelp(pt,run,repmat(n,1,n-1),rec.event)',pt) nan(1,8-n)],(2:7)','Uni',0))];

llHelp2 = @(pt,run,inp,rec) struct('true',rec.event,'decoded',allevlHelp(pt,run,pickmax(inp,1:7),rec.event),...
    'zval',zval(pt,run,inp,rec));
lleachtype = @(pt,run) arrayfun(@(x) llHelp2(pt,run,allevl(x,pt,run),x),recall{pt}{run});

rng(1)

allperms = arrayfun(@(pt) reshape(expand(arrayfun(@(rn) arrayfun(@(ct) ...
    perms(indx([recall{pt}{rn}.event]',floor([recall{pt}{rn}.event]/100)==ct)),...
    (1:2)','Uni',0),(1:length(recall{pt})),'Uni',0)),[2 length(recall{pt})]), 1:5,'Uni',0);
allperms = cellfun(@(pt) cellfun(@(x) indx(x,randperm(size(x,1),min([1e4 size(x,1)]))),...
    pt,'Uni',0),allperms,'Uni',0);

allpermStruct = arrayfun(@(pt) reshape(expand(arrayfun(@(rn) arrayfun(@(ct) arrayfun(@(prm)...
    struct('event',num2cell(allperms{pt}{ct,rn}(prm,:)),'sig',{recall{pt}{rn}(floor([recall{pt}{rn}.event]/100)==ct).sig}),...
    (1:size(allperms{pt}{ct,rn},1))','Uni',0),(1:2)','Uni',0),(1:length(recall{pt}))','Uni',0)),[2 length(recall{pt})]),1:5,'Uni',0);

llDist = @(pt,rn) arrayfun(@(ct) cellfun(@(y) arrayfun(@(x) llHelp2(pt,rn,allevl(x,pt,rn),x),y),...
    allpermStruct{pt}{ct,rn},'Uni',0),1:2,'Uni',0);

% $$$
%%%%%%% Comment these lines when running the code from saved output %%%%%%%
decodedDistIn = arrayfun(@(pt) arrayfun(@(rn) llDist(pt,rn),1:length(recall{pt}),'Uni',0),1:5,'Uni',0);
decodedDistIn = arrayfun(@(ct) arrayfun(@(pt) arrayfun(@(rn) sum(arrayfun(@(x) x.zval.rank,cell2mat(decodedDistIn{pt}{rn}{ct})),2),...
    1:length(recall{pt}),'Uni',0),1:5,'Uni',0),(1:2)','Uni',0);

decodedIn = arrayfun(@(pt) arrayfun(@(run) lleachtype(pt,run),1:length(recall{pt}),'Uni',0),1:5,'Uni',0);
decodedIn = arrayfun(@(ct) cellfun(@(pt) cellfun(@(rn) sum(arrayfun(@(x) x.zval.rank,rn(floor([rn.true]/100)==ct))),pt,'Uni',0),...
    decodedIn,'Uni',0),(1:2)','Uni',0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

combSubj3 = @(inp) cellfun(@(x) [x(1) {{x{2}{1},x{3}{1}}} x(4:5)],inp,'Uni',0);
combRuns = @(inp) cellfun(@(x) cellfun(@(y) y{1}+y{2},x,'Uni',0),inp,'Uni',0);
combRunsDist = @(inp) cellfun(@(x) cellfun(@(y) y{1}(randi(size(y{1},1),[1e4 1]))+y{2}(randi(size(y{2},1),[1e4 1])),...
    x,'Uni',0),inp,'Uni',0);

% Comment the following line when running the code from saved output. $$$
decodedDistIn = combRunsDist(combSubj3(decodedDistIn)); decodedIn = combRuns(combSubj3(decodedIn));
zcomb1 = cell2mat(cellfun(@(dec,dist) mean(bsxfun(@gt,cell2mat(dec),cell2mat(dist))),decodedIn,decodedDistIn,'Uni',0));

typecounts = cell2mat(arrayfun(@(ct) cellfun(@(pt) ... % Recall events per category
    sum(cellfun(@(x) sum(floor([x.event]/100)==ct),pt)),recall),(1:2)','Uni',0));
typecounts = [typecounts(:,1) sum(typecounts(:,2:3),2) typecounts(:,4:5)];
rankcomb1 = cell2mat(cellfun(@cell2mat,decodedIn,'Uni',0))./typecounts;

% Fisher's Method
totalp1 = chi2cdf(-2*sum(sum(log(zcomb1(comp(1:2,:))))),2*sum(sum(comp(1:2,:))),'Upper');

% Bonferroni-Holm correction
compold = (zcomb1>(1-.05/sum(sum(comp(1:2,:)))))&comp(1:2,:);
comp1 = nan;
while compold ~= comp1
    comp1 = compold;
    comp1 = zcomb1>(1-.05/sum(sum(~compold&comp(1:2,:))));
end
thres1 = cell2mat(cellfun(@(y,z) mean((cellfun(@(x) prctile(x,95),y)./typecounts(z,:) - 1)/6),decodedDistIn,{1;2},'Uni',0));

%% Inter-category classification

global allevlHelp2
allevlHelp2 = @(pt,run,ind,ex) indx(runcensor(cats(pt),run),ind);
allevl = @(rec,pt,run) [nan(1,14);cell2mat(arrayfun(@(n) [alltogether(rec,allevlHelp2(pt,run,1:(n-1),rec.event)', ...
    allevlHelp2(pt,run,repmat(n,1,n-1),rec.event)',pt) nan(1,15-n)],(2:14)','Uni',0))];

llHelp2 = @(pt,run,inp,rec) struct('true',rec.event,'decoded',allevlHelp2(pt,run,pickmax(inp,1:14),rec.event),...
    'zval',zvalAll(pt,run,inp,rec));
lleachtype = @(pt,run) arrayfun(@(x) llHelp2(pt,run,allevl(x,pt,run),x),recall{pt}{run});

rng(1)
allperms = arrayfun(@(pt) arrayfun(@(rn) cell2mat(arrayfun(@(it)...
    indx([recall{pt}{rn}.event]',randperm(length(recall{pt}{rn})))',...
    (1:2e3)','Uni',0)),(1:length(recall{pt})),'Uni',0), 1:5,'Uni',0);

allpermStruct = arrayfun(@(pt) arrayfun(@(rn) arrayfun(@(prm)...
    struct('event',num2cell(allperms{pt}{rn}(prm,:)),'sig',{recall{pt}{rn}.sig}),...
    (1:size(allperms{pt}{rn},1))','Uni',0),(1:length(recall{pt})),'Uni',0),1:5,'Uni',0);

llDist = @(pt,rn) cellfun(@(y) arrayfun(@(x) llHelp2(pt,rn,allevl(x,pt,rn),x),y),...
    allpermStruct{pt}{rn},'Uni',0);

combSubj3 = @(inp) [inp(1) {{inp{2}{1},inp{3}{1}}} inp(4:5)];

% $$$
%%%%%%% Comment these lines when running the code from saved output %%%%%%%
decodedDistOut = arrayfun(@(pt) arrayfun(@(rn) llDist(pt,rn),1:length(recall{pt}),'Uni',0),1:5,'Uni',0);

decodedDistOut = arrayfun(@(pt) arrayfun(@(rn) sum(arrayfun(@(x) x.zval.rank,cell2mat(decodedDistOut{pt}{rn})),2),...
    1:length(recall{pt}),'Uni',0),1:5,'Uni',0);

decodedOut = arrayfun(@(pt) arrayfun(@(run) lleachtype(pt,run),1:length(recall{pt}),'Uni',0),1:5,'Uni',0);
decodedOut = cellfun(@(pt) cellfun(@(rn) sum(arrayfun(@(x) x.zval.rank,rn)),pt,'Uni',0),...
    decodedOut,'Uni',0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

combRuns = @(inp) cellfun(@(y) y{1}+y{2},inp,'Uni',0);
combRunsDist = @(inp) cellfun(@(y) y{1}(randi(size(y{1},1),[1e4 1]))+y{2}(randi(size(y{2},1),[1e4 1])),...
    inp,'Uni',0);

% Comment the following line when running the code from saved output. $$$
decodedDistOut = combRunsDist(combSubj3(decodedDistOut)); decodedOut = combRuns(combSubj3(decodedOut));

zcomb2 = cell2mat(cellfun(@(dec,dist) mean(bsxfun(@gt,dec,dist)),decodedOut,decodedDistOut,'Uni',0));

typecounts = cellfun(@(pt) ... % Recall events per category
    sum(cellfun(@(x) length(x),pt)),recall);
typecounts = [typecounts(:,1) sum(typecounts(:,2:3),2) typecounts(:,4:5)];
rankcomb2 = cell2mat(decodedOut)./typecounts;

% Fisher's Method
totalp2 = chi2cdf(-2*sum(log(zcomb2(comp(3,:)))),2*sum(comp(3,:)),'Upper');

thres2 = mean(cellfun(@(x) prctile(x,95),decodedDistOut)./typecounts/49);

% Bonferroni-Holm correction
compold = (zcomb2>(1-.05/sum(comp(3,:))))&comp(3,:);
comp2 = nan;
while compold ~= comp2
    comp2 = compold;
    comp2 = zcomb2>(1-.05/sum(~compold&comp(3,:)));
end
%% Plot recall decoding

subplot(5,1,[4 5])
yyaxis left
hx1 = bar([2 1 3 4]',[rankcomb1;4*ones(1,4)]','EdgeColor','none','BaseValue',4);
hx1(1).FaceColor = [94.5 23.5 36.1]/100;
hx1(2).FaceColor = [24.7 30.2 71.8]/100;
set(gca,'XTickLabel',[2 3 6 8])
ylim([1 7])
ylabel('Correct exemplar rank','FontSize',10)
hold on
plhold = bsxfun(@plus,[-.23 0]',[2 1 3 4]);
plot(plhold(comp1),mean([7 6*max([(max(max(rankcomb2))-7.5)/13 (max(max(rankcomb1))-4)/6])+4]),'k*')
set(gca,'box','off')
set(gca,'TickLength',[.01 .01])
xlabel('Subject','FontSize',10)
plhold = bsxfun(@plus,[.23]',[2 1 3 4]);
plot(plhold(comp2),mean([7 6*max([(max(max(rankcomb2))-7.5)/13 (max(max(rankcomb1))-4)/6])+4]),'k*')

yyaxis right
hx2 = bar([2 1 3 4]',[7.5*ones(2,4);rankcomb2]','EdgeColor','none','BaseValue',7.5,'FaceColor',[55 75 35]/100);
set(gca,'XTickLabel',[2 3 6 8])
ylim([1 14])
ylabel('Correct exemplar rank','FontSize',10)
set(gca,'box','off')
set(gca,'TickLength',[.01 .01])
xlabel('Subject','FontSize',10)

%% Plot unit counts

subplot(5,1,1)
cellcounts = cell2mat(cellfun(@(x) [sum(strcmpi(x,'Face')) sum(strcmpi(x,'Place'))]',...
    cellfun(@(pt){xllog(strcmp({xllog.pt},pt)).catsel},unique({xllog.pt}),'Uni',0),'Uni',0));
b = bar([2 1 3 4],cellcounts','LineWidth',2);
set(b(1),'FaceColor',[90 90 90]/100);
set(b(2),'FaceColor',[30 30 30]/100);
set(gca,'XTick',1:4);set(gca,'XTickLabel',[2 3 6 8])
ylabel('Units','FontSize',10)

set(gca,'box','off')
set(gca,'TickLength',[.01 .01])

FreeRecallMDSTime = cputime-t;
disp(['Runtime: ' num2str(FreeRecallMDSTime)])

function out = zval(pt,rn,inp,rec)
    global allevlHelp indx
    [p,~,sts] = ranksum(indx(llHelp(inp),allevlHelp(pt,rn,1:7,rec.event)==rec.event),...
        indx(llHelp(inp),allevlHelp(pt,rn,1:7,rec.event)~=rec.event),'tail','right');
    out = struct('p',1-p,'sums',sts.ranksum,'rank',...
        indx(tiedrank(llHelp(inp)),allevlHelp(pt,rn,1:7,rec.event)==rec.event));
end

function out = zvalAll(pt,rn,inp,rec)
    global allevlHelp2 indx
    [p,~,sts] = ranksum(indx(llHelp(inp),allevlHelp2(pt,rn,1:14,rec.event)==rec.event),...
        indx(llHelp(inp),allevlHelp2(pt,rn,1:14,rec.event)~=rec.event),'tail','right');
    out = struct('p',1-p,'sums',sts.ranksum,'rank',...
        indx(tiedrank(llHelp(inp)),allevlHelp2(pt,rn,1:14,rec.event)==rec.event));
end

function out = llHelp(inp)
    out = nansum(inp==0,2)+nansum(inp==1,1)';
end

% Tied sums of pair ranks are tested again just against each other.
function indout = pickmax(inp,inds)
    indout = inds(max(llHelp(inp(inds,inds))) == llHelp(inp(inds,inds)));
    if isequal(indout,inds)
        indout = pickmax(inp,indout);
    end
end