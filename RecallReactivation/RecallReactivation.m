% Run this script to produce Figs. 3B-F. Load output from
% FreeRecallPreprocessing.m and RecallSpikes.mat before running this code.
% Runtime: 2 sec on authors' hardware.

% Copyright Simon Khuvis, 2019. See github.com/IEEG for licensing and use
% information.

t = cputime;

% Extract the spike counts around the times of interest
extractFun = @(y,z) -.25*diff(arrayfun(@(a) sum(y>a),bsxfun(@plus,[-2;2]*3e4,arrayfun(@(x) x.start, z))),1,1); % 2 sec before and after speech onset
PSTHFun = @(y,z) -10*diff(arrayfun(@(a) sum(y>a),bsxfun(@plus,(-8.5:.1:8.5)'*3e4,arrayfun(@(x) x.start, z))),1,1); % Histogram around speech onset

batchNames = unique(arrayfun(@(x) x.batchname, xllog, 'Uni', 0),'stable');
stimCat = {'Face','Place'};
Units.Face = cell(1,5);
Units.Place = cell(1,5);
faceH = cell(1,5); placeH = cell(1,5);
facecum = []; placecum = []; overallHold = [];

indx = @(vec,inp) vec(inp,:);

for batchNo = [1 3 4 5 6]
    % Which indicies of xllog are associated with this batch?
    xlinds = arrayfun(@(x) strcmp(x.batchname,batchNames{batchNo}),xllog);

    % Which recall events are of faces? Places?
    faceEvents = structfun(@(z) cellfun(@(y) strncmp(arrayfun(@(x) x.ref, y,'Uni',0),'Face',4),z.events,'Uni',0),RecallSpikes(batchNo),'Uni',0);
    placeEvents = structfun(@(z) cellfun(@(y) strncmp(arrayfun(@(x) x.ref, y,'Uni',0),'Place',5),z.events,'Uni',0),RecallSpikes(batchNo),'Uni',0);

    for stimType = 1
        unitInd = arrayfun(@(x) x.class>0,xllog(xlinds)); % Visually-responsive cells.
        spikeCount = structfun(@(x) cellfun(extractFun,x.spikes(:,unitInd),... % Spike counts from those cells in the intervals of interest
            repmat(x.events',1,sum(unitInd)),'Uni',0), RecallSpikes(batchNo),'Uni',0);
        PTSH = structfun(@(x) cellfun(PSTHFun,x.spikes(:,unitInd),...
            repmat(x.events',1,sum(unitInd)),'Uni',0), RecallSpikes(batchNo),'Uni',0);

        faceRates.Face = cellfun(@(x,y) x(y)' ,spikeCount.Face,repmat(faceEvents.Face',1,sum(unitInd)),'Uni',0);
        faceRates.Place = cellfun(@(x,y) x(y)' ,spikeCount.Place,repmat(faceEvents.Place',1,sum(unitInd)),'Uni',0);
        placeRates.Face = cellfun(@(x,y) x(y)' ,spikeCount.Face,repmat(placeEvents.Face',1,sum(unitInd)),'Uni',0);
        placeRates.Place = cellfun(@(x,y) x(y)' ,spikeCount.Place,repmat(placeEvents.Place',1,sum(unitInd)),'Uni',0);
        
        facePST.Face = cellfun(@(x,y) x(:,y)' ,PTSH.Face,repmat(faceEvents.Face',1,sum(unitInd)),'Uni',0);
        facePST.Place = cellfun(@(x,y) x(:,y)' ,PTSH.Place,repmat(faceEvents.Place',1,sum(unitInd)),'Uni',0);
        placePST.Face = cellfun(@(x,y) x(:,y)' ,PTSH.Face,repmat(placeEvents.Face',1,sum(unitInd)),'Uni',0);
        placePST.Place = cellfun(@(x,y) x(:,y)' ,PTSH.Place,repmat(placeEvents.Place',1,sum(unitInd)),'Uni',0);

        FaceRecalls = cell2mat(struct2cell(structfun(@cell2mat,faceRates,'Uni',0)));
        PlaceRecalls = cell2mat(struct2cell(structfun(@cell2mat,placeRates,'Uni',0)));
        
        facePSTH = cellfun(@(x,y) [x;y],facePST.Face,facePST.Place,'Uni',0);
        if size(facePSTH,1)==2
            facePSTH = cellfun(@(x,y) [x;y],facePSTH(1,:),facePSTH(2,:),'Uni',0);
        end
        placePSTH = cellfun(@(x,y) [x;y],placePST.Face,placePST.Place,'Uni',0);
        if size(placePSTH,1)==2
            placePSTH = cellfun(@(x,y) [x;y],placePSTH(1,:),placePSTH(2,:),'Uni',0);
        end
        overallHold = [overallHold cell2mat(cellfun(@(x,y) mean(mean([x(:,1:65);y(:,1:65)].^2)).^.5, ...
            facePSTH, placePSTH,'Uni',0))];
        faceH{batchNo} = cell2mat(permute(facePSTH,[1 3 2])); placeH{batchNo} = cell2mat(permute(placePSTH,[1 3 2]));
        facePSTH = cell2mat(cellfun(@mean,facePSTH,'Uni',0)');
        placePSTH = cell2mat(cellfun(@mean,placePSTH,'Uni',0)');
        facecum = [facecum;facePSTH];
        placecum = [placecum;placePSTH];

        % Ranksum spike rates from place vs. face recall
        dp = @(inp1,inp2) 2^.5*norminv(sum(indx(tiedrank([inp1;inp2]),1:length(inp1)))/(length(inp1)*length(inp2))-...
                (1+length(inp1))/(2*length(inp2)),0,1);
        Units.(stimCat{stimType}){batchNo} = cellfun(@(x,y) dp(x,y),...
            mat2cell(FaceRecalls,size(FaceRecalls,1),ones(1,sum(unitInd))),...
            mat2cell(PlaceRecalls,size(PlaceRecalls,1),ones(1,sum(unitInd))));
    end
end

rArrange = @(inp) [inp(1) inp(3:6)];
faceH = rArrange(faceH); placeH = rArrange(placeH);

presPref = structfun(@(x) (cell2mat(x)),Units,'Uni',0);

%%
indx = @(vec,inp) vec(inp,:);
sem = @(inp) std(inp,0,1)./size(inp,1)^.5;
revflip = @(inp) max([sem(inp)+mean(inp) fliplr(-sem(inp)+mean(inp))],0);
font = 'Helvetica';
set(gcf, 'Position', [0 0 1500 1000])
subplot(2,3,[1 2 4 5])

responsive = normcdf([xllog(arrayfun(@(y) (y.class>0)&&any(strcmp(y.batchname,batchNames([1 3 4 5 6]))),xllog)).VisResponsiveness])<.001;

selectivity = [xllog(arrayfun(@(y) (y.class>0)&&any(strcmp(y.batchname,batchNames([1 3 4 5 6]))),xllog)).dprime]';

[facethres,thresInd] = max(selectivity([xllog(arrayfun(@(y) (y.class>0)&&any(strcmp(y.batchname,batchNames([1 3 4 5 6]))),xllog)).class]==1));
facethres = mean([indx(selectivity,thresInd+1) facethres]);

% Just stringent response criterion.
enoughDat = ~isnan(presPref.Face)&responsive;
hold on
plot(selectivity(enoughDat),...
    presPref.Face(enoughDat)','k.','MarkerSize',20)

% Less strict responsiveness criterion.
enoughDat = ~isnan(presPref.Face);
FaceSel = selectivity'>facethres;

plot(selectivity(enoughDat),...
    presPref.Face(enoughDat)','.','MarkerSize',20,'Color',[.7 .7 .7])

recalldat = presPref.Face(enoughDat & FaceSel)';
fill([max(get(gca,'XLim')) facethres*[1 1] max(get(gca,'XLim'))],...
    [[1 1]*max(get(gca,'YLim')) [1 1]*min(get(gca,'YLim'))],[0 0 0],'EdgeColor','none')
alpha(.1)
set(gca,'children',flipud(get(gca,'children')))

hold off
set(gca,'box','off')
set(gca,'TickLength',[.01 .01])

pft = polyfit(selectivity(enoughDat),presPref.Face(enoughDat)',1);
xdat = get(gca,'XLim');
line(xdat,polyval(pft,xdat),'Color',[0 0 0])
% ylabel('Recall face preference','FontSize',10)
ylabel('\sl D''\rm sensitivity for face images during recall','FontSize',10)
xlabel('\sl D''\rm sensitivity for face images during presentation','FontSize',10)
[rho,p] = corr(selectivity(enoughDat),presPref.Face(enoughDat)','Type','Spearman');
text(xdat(2),get(gca,'YLim')*[0;1],{[' \rho = ' num2str(rho,2)]...
    ['\sl p\rm = ' num2str(p,1)]},'HorizontalAlignment','right','VerticalAlignment','top','FontSize',10)

ftcourse = subplot(2,3,3);

faceunits = cellfun(@(bn,fs) cellfun(@(z) mean(cell2mat(z),3),mat2cell(permute(indx([xllog(arrayfun(@(y) ...
    (y.class>0)&&any(strcmp(y.batchname,bn)),xllog)).timecourse]',fs),[3 2 1]),1,[1 1],sum(fs)),'Uni',0),...
    batchNames([1 3 4 5 6]),mat2cell(FaceSel,1,cellfun('size',faceH,3))','Uni',0);
meancourse = [mean(cell2mat(cellfun(@(x) x{1},faceunits,'Uni',0)));mean(cell2mat(cellfun(@(x) x{2},faceunits,'Uni',0)))];
sem = @(inp) std(inp)./size(inp,1).^.5;
stdcourse = [sem(cell2mat(cellfun(@(x) x{1},faceunits,'Uni',0)));sem(cell2mat(cellfun(@(x) x{2},faceunits,'Uni',0)))];

addme = -1.5e4:0.1e4:6.75e4;

xvals = (addme(2:end)-(addme(2)-addme(1))/2)/3e4;
f1 = fill([xvals fliplr(xvals)],[meancourse(1,:) + stdcourse(1,:)...
    fliplr(meancourse(1,:) - stdcourse(1,:))],[94.5 23.5 36.1]/100,'EdgeColor','none');
set(f1,'FaceAlpha',.1)
hold on
plot(xvals,meancourse(1,:),'Color',[94.5 23.5 36.1]/100)
xlim([addme(1) addme(end)]/3e4)


ylabel('Mean unit firing rate (Hz)','FontSize',10)

f2 = fill([xvals fliplr(xvals)],[meancourse(2,:) + stdcourse(2,:)...
    fliplr(meancourse(2,:) - stdcourse(2,:))],[24.7 30.2 71.8]/100,'EdgeColor','none');
set(f2,'FaceAlpha',.1)
plot(xvals,meancourse(2,:),'Color',[24.7 30.2 71.8]/100)
xlabel('Time after image presentation (s)','FontSize',10)

yzoom = zeros(2,2);
yzoom(1,:) = get(ftcourse,'YLim');

ylim([0 max(yzoom(:,2))]);
plot([0 0],[0 max(yzoom(:,2))],'k--')
plot([1.5174 1.5174],[0 max(yzoom(:,2))],'k--')

set(gca,'box','off')
set(gca,'TickLength',[.01 .01])



subplot(2,3,6)
hold on

Brates = arrayfun(@(a) structfun(@(x) cellfun(@(y,z) length(y)/z,x.spikes,repmat(x.duration',1,size(x.spikes,2))),a,'Uni',0),RecallSpikes([1 3 4 5 6]));
tempRate = arrayfun(@(x) structfun(@(y) mean(y,1),x,'Uni',0),Brates);
meanRate(1,:) = [tempRate.Face]; meanRate(2,:) = [tempRate.Place];
meanRate = repmat(mean(meanRate),2,1); % Place and face recalls can occur during any interval, so I average the baseline rates together to avoid bias.
BLInds = arrayfun(@(x) x.class>0, xllog(arrayfun(@(y) any(strcmp(y.batchname,batchNames([1 3 4 5 6]))),xllog)));

Bvar = arrayfun(@(a) structfun(@(x) cellfun(@(y,z) var(diff(y)/3e4)/(length(y)-1),x.spikes,repmat(x.duration',1,size(x.spikes,2))),a,'Uni',0),RecallSpikes([1 3 4 5 6]));
tempVar = arrayfun(@(x) structfun(@(y) mean(y,1)/2,x,'Uni',0),Bvar);
avgVar = repmat(mean([[tempVar.Face];[tempVar.Place]]),2,1);

face2 = cell2mat(cellfun(@(x,ind,bl) mean(bsxfun(@rdivide,x(:,:,ind),permute(bl(ind),[3 2 1])),3),faceH,mat2cell(FaceSel,1,cellfun('size',faceH,3)),mat2cell(meanRate(1,BLInds)',cellfun('size',faceH,3),1)','Uni',0)');
face2 = face2(any(~isnan(face2),2),:);
place2 = cell2mat(cellfun(@(x,ind,bl) mean(bsxfun(@rdivide,x(:,:,ind),permute(bl(ind),[3 2 1])),3),placeH,mat2cell(FaceSel,1,cellfun('size',placeH,3)),mat2cell(meanRate(1,BLInds)',cellfun('size',placeH,3),1)','Uni',0)');
place2 = place2(any(~isnan(place2),2),:);

plot((-8:.1:8)+.05,mean(conv2(face2,ones(1,10)/10,'valid')),'Color',[94.5 23.5 36.1]/100)
fill([-8:.1:8 8:-.1:-8]+.05,revflip(conv2(face2,ones(1,10)/10,'valid')),...
    ([94.5 23.5 36.1])/100,'EdgeColor','none')
alpha(.2)
plot((-8:.1:8)+.05,mean(conv2(place2,ones(1,10)/10,'valid')),'Color',[24.7 30.2 71.8]/100)
fill([-8:.1:8 8:-.1:-8]+.05,revflip(conv2(place2,ones(1,10)/10,'valid')),...
    ([24.7 30.2 71.8])/100,'EdgeColor','none')
alpha(.2)

[~,ptsigv,~,stts1] = ttest2((conv2(face2,ones(1,10)/10,'valid')),...
    (conv2(place2,ones(1,10)/10,'valid')));
ptsigv = ptsigv(41:5:121);
ptsig = ptsigv<=0.05; % Uncorrected

[~,fsel,~,stts2] = ttest2(mean(indx((conv2(face2,ones(1,10)/10,'valid'))',65:105)),...
    mean(indx((conv2(place2,ones(1,10)/10,'valid'))',65:105)));
disp(['d.f. = ' num2str(stts2.df) '; t statistic = ' num2str(stts2.tstat) '.'])

hold on

text(4,get(gca,'YLim')*[0;1],['\sl p\rm = ' num2str(fsel,1)],'HorizontalAlignment','right','VerticalAlignment','top','FontSize',10)

top = max(get(gca,'YLim'));
fill([-2 -2 2 2],[0 [1 1]*top 0],[0 0 0],'EdgeColor','none')
plot(bsxfun(@plus,indx((-4:.5:4)',ptsig),[-.25 .25]),[.95 .95]*top,'k','LineWidth',3)
alpha(.1)
set(gca,'children',flipud(get(gca,'children')))
xlim([-4 4])
ylim([0 top])
set(gca,'box','off')
set(gca,'TickLength',[.01 .01])
hold off
ylabel('Baseline-corrected firing rate','FontSize',10)
xlabel('Time after utterance onset (s)','FontSize',10)

%% Baseline shift?

Brates = arrayfun(@(a) structfun(@(x) cellfun(@(y,z) length(y)/z,x.spikes,repmat(x.duration',1,size(x.spikes,2))),a,'Uni',0),RecallSpikes([1 3 4 5 6]));
tempRate = arrayfun(@(x) structfun(@(y) mean(y,1),x,'Uni',0),Brates);
meanRate(1,:) = [tempRate.Face]; meanRate(2,:) = [tempRate.Place];
BLInds = arrayfun(@(x) x.class>0, xllog(arrayfun(@(y) any(strcmp(y.batchname,batchNames([1 3 4 5 6]))),xllog)));
BLFacePref = -diff(meanRate(:,BLInds))./sum(meanRate(:,BLInds));

font = 'Helvetica';
figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
set(gcf, 'Position', [0 0 700 500])

plot(selectivity,BLFacePref','k.','MarkerSize',20)
pft = polyfit(selectivity,BLFacePref',1);
xdat = get(gca,'XLim');
line(xdat,polyval(pft,xdat),'Color',[0 0 0])
ylabel('Preference for face recall periods','FontSize',10)
xlabel('\sl D''\rm sensitivity for face images during presentation','FontSize',10)
[rho,p] = corr(selectivity,BLFacePref','Type','Spearman');
text(xdat(2),get(gca,'YLim')*[0;1],{[' \rho = ' num2str(rho,2)]...
    ['\sl p\rm = ' num2str(p,1)]},'HorizontalAlignment','right','VerticalAlignment','top')
title('Units by preference for face recall periods','FontSize',10)

set(gca,'box','off')
set(gca,'TickLength',[.01 .01])

%% Delta firing rates presentation-recall

figure
enoughDat = ~isnan(presPref.Face)&responsive;

ratio = @(inp) inp(:,1)./inp(:,2);
pickFace = @(inp) arrayfun(@(x) mean(-diff(x.response(x.response(:,1)<200,[2 4]),1,2)/.4),inp);
fpres = pickFace(xllog(BLInds));
frec = mean(indx(facecum',65:105))'-meanRate(1,BLInds)';
[rho,p] = corr(fpres,frec,'Type','Spearman');
plot(fpres,frec,'.','MarkerSize',20,'Color',[.7 .7 .7])
hold on
plot(fpres(enoughDat),frec(enoughDat),'k.','MarkerSize',20)
pft = polyfit(fpres,frec,1);
xdat = get(gca,'XLim');
line(xdat,polyval(pft,xdat),'Color',[0 0 0])
ylabel('\Delta firing rate: face recall (Hz)','FontSize',10)
xlabel('\Delta firing rate: face presentation (Hz)','FontSize',10)
text(xdat(2),get(gca,'YLim')*[0;1],{[' \rho = ' num2str(rho,2)]...
    ['\sl p\rm = ' num2str(p,1)]},'HorizontalAlignment','right','VerticalAlignment','top')

fInc = (mean(indx(facecum',65:85))'-meanRate(1,BLInds)')./(norminv(1-.05/sum(BLInds)).*(avgVar(1,BLInds)'+var(indx(facecum',65:85))'/size(facecum,1)).^.5);
pInc = (mean(indx(placecum',65:85))'-meanRate(1,BLInds)')./(norminv(1-.05/sum(BLInds)).*(avgVar(1,BLInds)'+var(indx(placecum',65:85))'/size(placecum,1)).^.5);

fIncf = indx((mean(indx(facecum',65:85))'-meanRate(1,BLInds)')./(norminv(1-.05/sum(FaceSel)).*(avgVar(1,BLInds)'+var(mean(indx(facecum',65:85)))'/size(facecum,1)).^.5),FaceSel);
pIncf = indx((mean(indx(placecum',65:85))'-meanRate(1,BLInds)')./(norminv(1-.05/sum(FaceSel)).*(avgVar(1,BLInds)'+var(mean(indx(placecum',65:85)))'/size(placecum,1)).^.5),FaceSel);

set(gca,'box','off')
set(gca,'TickLength',[.01 .01])

OneBackTime = cputime-t;
disp(['Runtime: ' num2str(OneBackTime)])