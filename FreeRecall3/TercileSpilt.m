% Run this script to produce Figs. 4D,E. Load output from
% FreeRecallPreprocessing.m and RecallSpikes.mat before running this code.
% Runtime: ~3 min on authors' hardware.

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

%% Recall Exemplar Decoding

extr = @(sg,en,run) cellfun(@(x) x(en), sg(:,run)); % Take one event from one run at a time.
ref = @(inp) 100*find(strcmpi(inp.ref(isstrprop(inp.ref,'alpha')),{'Face','Place'}))+...
    str2double(inp.ref(isstrprop(inp.ref,'digit')))/10; % Extract id of stimulus.

bhhelp = @(inp) min(inp(2,inp(1,:)<inp(3,:),:));
bonfholm = @(inp,alpha) inp<(alpha/bhhelp([sortrows([inp;1:length(inp)]','MissingPlacement','last')';...
    alpha./(sum(~isnan(inp)):-1:(sum(~isnan(inp))-length(inp)+1))]));

% Put first block and second block together (1st and 3rd, 2nd and 4th, when
% struct2mat is used).
recmerge = @(inp) cellfun(@(x) [inp{x}],...
    mat2cell(reshape(1:length(inp),length(inp)/2,2),ones(length(inp)/2,1),2),'Uni',0);
% Remove recall events that are from the other block (even if they occured
% in the wrong block).
reccensor = @(rec,pt,r) rec(ismember([rec.event],runcensor(cats(pt),r)));
% Create struct array for each recall event and neuronal activity.
recall = cellfun(@(ev,sg,pt) recmerge(arrayfun(@(run) reccensor(arrayfun(@(en) struct(...
    'event',ref(extr(ev,en,run)),'sig',extr(sg,en,run)),(1:length(ev{run}))),pt,mod(run-1,length(ev)/2)+1),...
    1:length(ev),'Uni',0)),events,sigs,num2cell(1:5),'Uni',0);

recall = cellfun(@(pt) cellfun(@(rn) arrayfun(@(ev) struct('event',ev.event,'sig',ev.sig),rn),pt,'Uni',0),recall,'Uni',0);

% Spike counts: each row is an exemplar; each column is a unit.
itts = 10000;
SCrec = cellfun(@(rc,bno) cell2mat(arrayfun(@(ut) accumarray(arrayfun(@(x) find(x.event==cats(bno),1), cell2mat(rc')'), ...
    cell2mat(cellfun(@(y) arrayfun(@(z) z.sig(ut),y),rc','Uni',0))',[length(cats(bno)) 1])./...
    arrayfun(@(w) sum(w==cell2mat(cellfun(@(v) [v.event],rc','Uni',0))),cats(bno)),...
    1:length(rc{1}(1).sig),'Uni',0)),recall,num2cell(1:5),'Uni',0);

newhist = @(inp) histc(inp,unique(inp,'stable'));
indxScram = @(inp,rows) inp(accumarray(cell2mat(arrayfun(@(x) [[indx(find(rows),randperm(sum(rows)));find(~rows)] x*ones(size(inp,1),1)],...
    (1:size(inp,2))','Uni',0)),repmat([find(rows);find(~rows)],size(inp,2),1))+repmat(0:size(inp,1):size(inp,1)*(size(inp,2)-1),size(inp,1),1));

rng(1)
recsamp = arrayfun(@(it) cellfun(@(pt,ptNo) cell2mat(indxScram(mat2cell(pt,ones(size(pt,1),1),newhist([xllog(strcmp({xllog.batchname},batches{ptNo})&([xllog.class]>0)).chanNo])),~isnan(pt(:,1)))),...
    SCrec,num2cell(1:length(SCrec)),'Uni',0),1:itts,'Uni',0);

SCpres = cellfun(@(bn,bno) cell2mat(cellfun(@(res) accumarray(arrayfun(@(x) find(x==cats(bno),1),res(:,1)),res(:,2))./ ...
    arrayfun(@(x) sum(x==res(:,1)),cats(bno)),{xllog(strcmp({xllog.batchname},bn)&([xllog.class]>0)).response},'Uni',0)),...
    batches,num2cell(1:5),'Uni',0);

presSplit = cellfun(@(p,r) ((p>prctile(p(~isnan(r(:,1)),:),200/3))&(~isnan(r))) - ...
    ((p<prctile(p(~isnan(r(:,1)),:),100/3))&(~isnan(r))),SCpres,SCrec,'Uni',0);

recSplit = cellfun(@(p,r) (r>prctile(r,200/3)) - (r<prctile(r,100/3)),SCpres,SCrec,'Uni',0);

recSplitDist = cellfun(@(it) cellfun(@(p,r) (r>prctile(r,200/3)) - (r<prctile(r,100/3)),SCpres,it,'Uni',0),recsamp,'Uni',0);

reorg = @(inp) [inp(1) [inp{2:3}] inp(4:5)];
match = cellfun(@(r,p) (sum((r==1)&(p==1)) + sum((r==-1)&(p==-1)))./sum(abs(r)==1),recSplit,presSplit,'Uni',0);
match = reorg(match);
matchDist = cellfun(@(rsd) cellfun(@(r,p) (sum((r==1)&(p==1))+ sum((r==-1)&(p==-1)))./sum(abs(r)==1),rsd,presSplit,'Uni',0),recSplitDist,'Uni',0);
matchDist = cellfun(@cell2mat,mat2cell(reshape([matchDist{:}],length(recSplit),[])',itts,ones(1,length(recSplit))),'Uni',0);
matchDist = reorg(matchDist);

noFaceUnits = cellfun(@(bn) sum(strcmp({xllog(ismember({xllog.batchname},bn)&([xllog.class]>0)).catsel},'Face')),...
    {batches(1) batches(2:3) batches(4) batches(5)});
noUnits = cellfun(@(bn) sum(ismember({xllog.batchname},bn)&([xllog.class]>0)),...
    {batches(1) batches(2:3) batches(4) batches(5)});
%%
subplot(5,1,1)
title('Number of Units Recorded per Subject')
plotord = 1+sum(noFaceUnits>0)-tiedrank(noFaceUnits(noFaceUnits>0));
bar(plotord, [noUnits(noFaceUnits>0);noFaceUnits(noFaceUnits>0)]')
legend({'Visually-Responsive Units','Face-Selective Units'})
ylabel('Units')
xlabel('Subject')
set(gca,'XTickLabel',num2cell(accumarray(plotord',indx(cellfun(@(x) str2num(x(2)), unique({xllog.pt},'stable'))',noFaceUnits>0))))
set(gca,'TickLength',[0 .01])

prHelp = @(r,p,mr,mp,scr) nanmean(sum((r==mr)&(p==mp)&(~isnan(scr)))./sum((r==mr)&~isnan(scr)));
fprHelp = @(r,p,mr,mp,scr,bn) nanmean(sum((r==mr)&(p==mp)&(~isnan(scr))&...
    repmat(strcmp({xllog(ismember({xllog.batchname},bn)&([xllog.class]>0)).catsel},'Face'),size(r,1),1))...
    ./sum((r==mr)&~isnan(scr)&repmat(strcmp({xllog(ismember({xllog.batchname},bn)&([xllog.class]>0)).catsel},'Face'),size(r,1),1)));
prop = cellfun(@(r,p,scr) arrayfun(@(mr,mp) prHelp(r,p,mr,mp,scr),meshgrid(-1:1),meshgrid(-1:1)'),reorg(recSplit),reorg(presSplit),reorg(SCrec),'Uni',0);

fprop = cellfun(@(r,p,scr,bn) arrayfun(@(mr,mp) fprHelp(r,p,mr,mp,scr,bn),...
    meshgrid(-1:1),meshgrid(-1:1)'),reorg(recSplit),reorg(presSplit),reorg(SCrec),...
    {batches(1) batches(2:3) batches(4) batches(5)},'Uni',0);

a2 = subplot(5,1,[2 3]);
subsel = @(inp) cellfun(@(x) x(:,[1 3]),inp,'Uni',0); % For just showing the top and bottom terciles.
bar(reshape(bsxfun(@plus,([1.5 2.5])',(5*plotord-5)),2*length(plotord),1),cell2mat(subsel(fprop(noFaceUnits>0)))','stacked');
set(a2,'XTick',sort(reshape(bsxfun(@plus,([1.5 2.5])',(5*plotord-5)),2*length(plotord),1)))
set(a2,'XTickLabel',{'Bottom','Top'})
title('Mean Distribution of Item Recall Response by Persentation Response Tercile (Face Units Only)')
hold on

% Face units, all exemplars
allexF = cellfun(@(m,md,bn) mean(nanmean(m(strcmp({xllog(ismember({xllog.batchname},bn)&([xllog.class]>0)).catsel},'Face')))<= ...
    nanmean(md(:,strcmp({xllog(ismember({xllog.batchname},bn)&([xllog.class]>0)).catsel},'Face')),2)), match,matchDist,{batches(1) batches(2:3) batches(4) batches(5)});
allexiF = bonfholm(allexF.*(noFaceUnits./noFaceUnits),.05);
plot(a2,-3+5*plotord(allexiF(noFaceUnits>0)),1.1,'k*')
text(a2,-3+5*plotord,repmat(1.2,1,sum(noFaceUnits>0)),...
    cellfun(@(x) num2str(x,1),num2cell(allexF(noFaceUnits>0)),'Uni',0),'HorizontalAlignment','center','VerticalAlignment','bottom')
ylim(a2,[0 1.4])
set(a2,'TickLength',[0 .01])

% Restrict to just face trials

presSplitF = cellfun(@(p,r,pt) ((p>prctile(p(~isnan(r(:,1))&(cats(pt)<200),:),200/3))&(~isnan(r)&repmat(cats(pt)<200,1,size(r,2)))) - ...
    ((p<prctile(p(~isnan(r(:,1))&(cats(pt)<200),:),100/3))&(~isnan(r)&repmat(cats(pt)<200,1,size(r,2)))),SCpres,SCrec,num2cell(1:5),'Uni',0);

recSplitF = cellfun(@(p,r,pt) ((r>prctile(r(cats(pt)<200,:),200/3))&(~isnan(r)&repmat(cats(pt)<200,1,size(r,2)))) - ...
    ((r<prctile(r(cats(pt)<200,:),100/3))&(~isnan(r)&repmat(cats(pt)<200,1,size(r,2)))),SCpres,SCrec,num2cell(1:5),'Uni',0);

recSplitDistF = cellfun(@(it) cellfun(@(p,r,pt) ((r>prctile(r(cats(pt)<200,:),200/3))&(~isnan(r)&repmat(cats(pt)<200,1,size(r,2)))) - ...
    ((r<prctile(r(cats(pt)<200,:),100/3))&(~isnan(r)&repmat(cats(pt)<200,1,size(r,2)))),SCpres,it,num2cell(1:5),'Uni',0),recsamp,'Uni',0);

matchF = cellfun(@(r,p) (sum((r==1)&(p==1)) + sum((r==-1)&(p==-1)))./sum(abs(r)==1),recSplitF,presSplitF,'Uni',0);
matchF = reorg(matchF);
matchDistF = cellfun(@(rsd) cellfun(@(r,p) (sum((r==1)&(p==1)) + sum((r==-1)&(p==-1)))./sum(abs(r)==1),rsd,presSplitF,'Uni',0),recSplitDistF,'Uni',0);
matchDistF = cellfun(@cell2mat,mat2cell(reshape([matchDistF{:}],length(recSplit),[])',itts,ones(1,length(recSplit))),'Uni',0);
matchDistF = reorg(matchDistF);

propF = cellfun(@(r,p,scr) arrayfun(@(mr,mp) prHelp(r,p,mr,mp,scr),meshgrid(-1:1),meshgrid(-1:1)'),reorg(recSplit),reorg(presSplit),...
    reorg(cellfun(@(s,pt) bsxfun(@plus,s,(1.0*(cats(pt)<200)./(1.0*(cats(pt)<200)))),SCrec,num2cell(1:5),'Uni',0)),'Uni',0);

fpropF = cellfun(@(r,p,scr,bn) arrayfun(@(mr,mp) fprHelp(r,p,mr,mp,scr,bn),...
    meshgrid(-1:1),meshgrid(-1:1)'),reorg(recSplit),reorg(presSplit),...
    reorg(cellfun(@(s,pt) bsxfun(@plus,s,(1.0*(cats(pt)<200)./(1.0*(cats(pt)<200)))),SCrec,num2cell(1:5),'Uni',0)),...
    {batches(1) batches(2:3) batches(4) batches(5)},'Uni',0);

a4 = subplot(5,1,[4 5]);
bar(reshape(bsxfun(@plus,([1.5 2.5])',(5*plotord-5)),2*length(plotord),1),cell2mat(subsel(fpropF(noFaceUnits>0)))','stacked');
set(a4,'XTickLabel',{'Bottom','Top'})
set(a4,'XTick',sort(reshape(bsxfun(@plus,([1.5 2.5])',(5*plotord-5)),2*length(plotord),1)))
title('Mean Distribution of Face Recall Response by Persentation Response Tercile (Face Units Only)')
legend({'Bottom Tercile Presentation','Middle Tercile Presentation','Top Tercile Presentation'})
hold on

% Face units
fexF = cellfun(@(m,md,bn) mean(nanmean(m(strcmp({xllog(ismember({xllog.batchname},bn)&([xllog.class]>0)).catsel},'Face')))<=...
    nanmean(md(:,strcmp({xllog(ismember({xllog.batchname},bn)&([xllog.class]>0)).catsel},'Face')),2)), matchF,matchDistF,...
    {batches(1) batches(2:3) batches(4) batches(5)});
fexiF = bonfholm(fexF.*(noFaceUnits./noFaceUnits).*...
    ((1.0*allexiF)./(1.0*allexiF)),.05);
plot(a4,-3+5*plotord(fexiF(noFaceUnits>0)),1.1,'k*')
text(a4,-3+5*plotord(allexiF(noFaceUnits>0)),repmat(1.2,1,sum(allexiF(noFaceUnits>0))),cellfun(@(x) num2str(x,1),num2cell(fexF((noFaceUnits>0)&allexiF)),'Uni',0),...
    'HorizontalAlignment','center','VerticalAlignment','bottom')
text(a4,-3+5*plotord(~allexiF(noFaceUnits>0)),repmat(1.2,1,sum(~allexiF&(noFaceUnits>0))),...
    cellfun(@(x) num2str(x,1),num2cell(fexF((noFaceUnits>0)&~allexiF)),'Uni',0),...
    'HorizontalAlignment','center','VerticalAlignment','bottom','Color',.5*[1 1 1])
ylim(a4,[0 1.4])
set(a4,'TickLength',[0 .01])

TercileSplitTime = cputime-t;
disp(['Runtime: ' num2str(TercileSplitTime)])