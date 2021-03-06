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

indxScram = @(inp,rows) inp(accumarray([find(rows);find(~rows)],[indx(find(rows),randperm(sum(rows)));find(~rows)]),:);

rng(1)

recsamp = arrayfun(@(it) cellfun(@(pt,ptNo) indxScram(pt,~isnan(pt(:,1))),...
    SCrec,num2cell(1:length(SCrec)),'Uni',0),1:itts,'Uni',0);
recsampF = arrayfun(@(it) cellfun(@(pt,ptNo) indxScram(pt,(cats(ptNo)<200)&~isnan(pt(:,1))),...
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

% Face-selective units
fsel = reorg(cellfun(@(bn) strcmp({xllog(ismember({xllog.batchname},bn)&([xllog.class]>0)).catsel},'Face'),batches,'Uni',0));

a2 = subplot(5,1,[2 3]);
subsel = @(inp) cellfun(@(x) x(:,[1 3]),inp,'Uni',0); % For just showing the top and bottom terciles.
bar(reshape(bsxfun(@plus,([1.5 2.5])',(5*plotord-5)),2*length(plotord),1),cell2mat(subsel(fprop(noFaceUnits>0)))','stacked');
set(a2,'XTick',sort(reshape(bsxfun(@plus,([1.5 2.5])',(5*plotord-5)),2*length(plotord),1)))
set(a2,'XTickLabel',{'Bottom','Top'})
title('Mean Distribution of Item Recall Response by Persentation Response Tercile (Face Units Only)')
hold on

% Face units, all exemplars
allexF = cellfun(@(m,md,nx) mean(nanmean(m(fsel{nx}))<= ...
    nanmean(md(:,fsel{nx}),2)), match,matchDist,num2cell(1:4));
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
    ((r<prctile(r(cats(pt)<200,:),100/3))&(~isnan(r)&repmat(cats(pt)<200,1,size(r,2)))),SCpres,it,num2cell(1:5),'Uni',0),recsampF,'Uni',0);

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
fexF = cellfun(@(m,md,nx) mean(nanmean(m(fsel{nx}))<=...
    nanmean(md(:,fsel{nx}),2)), matchF,matchDistF,...
    num2cell(1:4));
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

%% Plot surrogate distributions

figure

ds = @(x) linspace(min(x),max(x),30);
b1 = {batches(1) batches(2:3) batches(4) batches(5)};
bx = {1 [2 3] 4 5};
n1 = [3 2 6 8];
ord = [1 3 NaN 2];
% fsel = reorg(cellfun(@(bn) strcmp({xllog(ismember({xllog.batchname},bn)&([xllog.class]>0)).catsel},'Face'),batches,'Uni',0));
for p1 = [1 2 4]
    d1 = nanmean(matchDistF{p1}(:,fsel{p1}),2);
    subplot(2,sum(~isnan(ord)),3+ord(p1))
    bar(ds(d1),histc(d1,ds(d1)),1)
    title(['Subject ' num2str(n1(p1)) ': '...
        num2str(prod(arrayfun(@(x) factorial(sum(~isnan(SCrec{x}(:,1))&(cats(x)<200))),bx{p1})),5) ' unique permuatations']);
    xlabel('Concordance, face trials')
    ylabel('Count')
    ylim([0 2000])
    hold on
    r1 = nanmean(matchF{p1}(:,fsel{p1}),2);
    plot(r1*[1 1],get(gca,'YLim'),'r')
    text(r1,max(get(gca,'YLim')),num2str(fexF(p1),2),'HorizontalAlignment','left','VerticalAlignment','top')
    hold off
    set(gca,'box','off')

    d2 = nanmean(matchDist{p1}(:,fsel{p1}),2);
    subplot(2,sum(~isnan(ord)),ord(p1))
    bar(ds(d2),histc(d2,ds(d2)),1)
    title(['Subject ' num2str(n1(p1)) ': '...
        num2str(prod(arrayfun(@(x) factorial(sum(~isnan(SCrec{x}(:,1)))),bx{p1})),5) ' unique permuatations']);
    xlabel('Concordance, all trials')
    ylabel('Count')
    ylim([0 2000])
    hold on
    r2 = nanmean(match{p1}(fsel{p1}),2);
    plot(r2*[1 1],get(gca,'YLim'),'r')
    text(r2,max(get(gca,'YLim')),num2str(allexF(p1),2),'HorizontalAlignment','left','VerticalAlignment','top')
    hold off
    set(gca,'box','off')
end

TercileSplitTime = cputime-t;
disp(['Runtime: ' num2str(TercileSplitTime)])