% Run this script to produce Fig. 2A. Load output from OneBack.m before
% running this code.
% Runtime: ~12 hrs on authors' hardware.

% Copyright Simon Khuvis, 2019. See github.com/IEEG for licensing and use
% information.

t = cputime;

% load('XLLog-FriedmanStrictTimecouse.mat')

indx = @(vec,ind) vec(ind,:);
% Get the responses out of the xllog variables.
extractLog = @(xlog) arrayfun(@(x) x.response(:,[1 3-x.oncell]),xlog([xlog.class]>0),'Uni',0);

justlog = @(extracted) cellfun(@(x) [x(:,1) log10(x(:,2)+1)],extracted,'Uni',0);

cats = unique(cell2mat(arrayfun(@(x) x.response(:,1),xllog,'Uni',0)));
counts = min(cell2mat(cellfun(@(x) histc(indx(x',1),cats),extractLog(xllog),'Uni',0)))';
minct = min(counts); % Min num of reps per exemplar

picksome = @(inp,num)indx(inp,randperm(length(inp),num)); % Randomly pick num from inp

logxl = justlog(extractLog(xllog));

extract = @(catinp,cld) cell2mat(cellfun(@(x) x(x(:,1) == catinp,2)',cld,'Uni',0))';
cmdsubset = @(inp) cmdscale(pdist(inp,'euclidean'),2);
xform = @(ex1,ex2) mat2cell(cmdsubset([ex1;ex2]),...
    [size(ex1,1) size(ex2,1)],2);
w = @(inp,n1,n2) squeeze(sum(cell2mat(permute(cellfun(@(x,y) cov(x(y,:)),inp,{n1';n2'},'Uni',0),[3 2 1])),3))^-1*...
    (mean(inp{1}(n1,:))-mean(inp{2}(n2,:)))';
normz = @(inp) inp./norm(inp);
cent = @(inp,n1,n2) mean([mean(inp{1}(n1,:));mean(inp{2}(n2,:))])';
evl = @(inp,n1,n2) cellfun(@(x,c,wx,n,z) (z*dot(x(setdiff(1:length(x),n),:)-c',normz(wx)'))>0,...
    inp,repmat({cent(inp,n1,n2)},2,1),repmat({w(inp,n1,n2)},2,1),{n1;n2},{1;-1},'Uni',0);
alltogether = @(r,c,cld) mean(cell2mat(evl(xform(extract(cats(r),cld),extract(cats(c),cld)),1:minct-1,1:minct-1)));

%% Fill in confusion matrix
reps = 100;
superreps = 100;
MAT = cell(1,1,superreps);
parfor sr = 1:superreps
    tic
    clIn = arrayfun(@(y) cellfun(@(x) cell2mat(arrayfun(@(ct) picksome(x(x(:,1) == ct,:),minct),...
        cats,'Uni',0)),logxl,'Uni',0),1:reps,'Uni',0);
    MAT{sr} = [nan(1,length(cats),reps);cell2mat(arrayfun(@(row) [cell2mat(arrayfun(@(col) ...
        permute(cellfun(@(x) alltogether(row,col,x),clIn),[1 3 2]),...
        1:row-1,'Uni',0)), nan(1,length(cats)-row+1,reps)], (2:length(cats))','Uni',0))];
    toc
end

% Surrogate confusion matrix
% Scramble stimuli in all extracted responses (to build surrogate
% distributions) within categories.
scramLog = @(extracted) cellfun(@(z) cell2mat(cellfun(@(y) [y(randperm(size(y,1)),1) y(:,2)], ...
    arrayfun(@(x) z(floor(z(:,1)/100)==x,:),unique(floor(cats/100))','Uni',0),'Uni',0).'),extracted,'Uni',0);

scramLogIntr = @(extracted,c1,c2) cellfun(@(z) [indx(z(ismember(floor(z(:,1)/100),[c1,c2]),1), ...
    randperm(sum(ismember(floor(z(:,1)/100),[c1,c2])))) z(ismember(floor(z(:,1)/100),[c1,c2]),2)],extracted,'Uni',0);

SCRlogxl = @(dummy) justlog(scramLog(extractLog(xllog)));
SCRlogxlIntr = @(dummy,c1,c2) justlog(scramLogIntr(extractLog(xllog),c1,c2));

SCRalltogether = @(r,c,scl) mean(cell2mat(evl(xform(extract(cats(r),scl),...
    extract(cats(c),scl)),1:minct-1,1:minct-1)));

rng(1)
reps = 100;
superreps = 100;
scrMAT = cell(1,1,superreps);
parfor sr = 1:superreps
    tic
    scramvalsP = arrayfun(@(x) SCRlogxl(x), 1:reps,'Uni',0);
    scramvals = cellfun(@(y) cellfun(@(x) cell2mat(arrayfun(@(ct) picksome(x(x(:,1) == ct,:),minct),...
        cats,'Uni',0)),y,'Uni',0),scramvalsP,'Uni',0);
    scrMAT{sr} = [nan(1,length(cats),reps);cell2mat(arrayfun(@(row) ...
        [nan(1,find(floor(cats/100)==floor(cats(row)/100),1)-1,reps) ...
        cell2mat(arrayfun(@(col) permute(cellfun(@(x) SCRalltogether(row,col,x),scramvals),[1 3 2]),...
        find(floor(cats/100)==floor(cats(row)/100),1):row-1,'Uni',0)), ...
        nan(1,length(cats)-row+1,reps)], (2:length(cats))','Uni',0))];
    for cl1 = 2:5
        for cl2 = 1:cl1-1
            rsup = find(floor(cats/100)==cl1,1)-1;
            csup = find(floor(cats/100)==cl2,1)-1;
            scramvalsP = arrayfun(@(x) SCRlogxlIntr(x,cl1,cl2), 1:reps,'Uni',0);
            scramvals = cellfun(@(y) cellfun(@(x) cell2mat(arrayfun(@(ct) picksome(x(x(:,1) == ct,:),minct),...
                cats(ismember(floor(cats/100),[cl1;cl2])),'Uni',0)),y,'Uni',0),scramvalsP,'Uni',0);
            scrMAT{sr}(floor(cats/100)==cl1,floor(cats/100)==cl2,:) = cell2mat(arrayfun(@(row) ...
                cell2mat(arrayfun(@(col) permute(cellfun(@(x) SCRalltogether(row+rsup,col+csup,x),...
                scramvals),[1 3 2]),1:sum(floor(cats/100)==cl2),'Uni',0)), ...
                (1:sum(floor(cats/100)==cl1))','Uni',0));
        end
    end
    toc
end

%% Calculate Classifier Accuracy
% load('output.mat')

acc = cellfun(@(x) nanmean(reshape(x,1,[])),mat2cell(nanmean(cell2mat(MAT),3),...
    histc(floor(cats/100),1:5),histc(floor(cats/100),1:5)));

accdistMake = @(inp) cellfun(@(x) nanmean(reshape(x,1,[])),mat2cell(cell2mat(inp), ...
    histc(floor(cats/100),1:5),histc(floor(cats/100),1:5),ones(1,sum(cellfun('size',MAT,3)))));
scrDist = accdistMake(scrMAT); accDist = accdistMake(MAT);

gtDist = cell2mat(arrayfun(@(it) mean(bsxfun(@gt,accDist(:,:,it),scrDist),3),permute(1:size(accDist,3),[1 3 2]),'Uni',0));
sig = diag(median(gtDist,3));

diag3 = @(x) cell2mat(arrayfun(@(it) diag(x(:,:,it)),permute(1:size(x,3),[1 3 2]),'Uni',0));
difs = @(x) x-permute(x,[2 1 3]);
scrdifsDist = difs(repmat(diag3(scrDist),1,5)); difsDist = difs(repmat(diag3(accDist),1,5));
dDist = cell2mat(arrayfun(@(it) mean(bsxfun(@gt,difsDist(:,:,it),scrdifsDist),3),permute(1:size(difsDist,3),[1 3 2]),'Uni',0));
difsig = triu(median(dDist,3));

%% Plot, calculate significance.

font = 'Helvetica';
figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
set(gcf, 'Position', [0 0 700 600])

mMAT = cell2mat(cellfun(@(x) [nan(1,12);nan(10,1) x nan(10,1); nan(1,12)], ...
    mat2cell(mean(cell2mat(MAT),3),repmat(10,5,1),repmat(10,5,1)),'Uni',0));
imagesc(mMAT,'alphadata',~isnan(mMAT),[.5 1])

% Bonferroni-Holm correction
compold = median(gtDist,3)>(1-.05/15);
comp5 = nan;
while compold ~= comp5
    comp5 = compold;
    comp5 = median(gtDist,3)>(1-.05/(sum(sum(~comp5))-10));
end

compold = median(gtDist,3)>(1-.001/15);
comp1 = nan;
while compold ~= comp1
    comp1 = compold;
    comp1 = median(gtDist,3)>(1-.001/(sum(sum(~comp1))-10));
end

hold on

for row = 1:5
    for col = 1:5
        if comp1(row,col)
            rectangle('Position',[12*(col-1)+1,12*(row-1)+1,11,11],'LineWidth',5);
        elseif comp5(row,col)
            rectangle('Position',[12*(col-1)+1,12*(row-1)+1,11,11],'LineWidth',1);
        end
    end
end
set(gca,'XTick',6:12:66)
set(gca,'YTick',6:12:66)
set(gca,'XTickLabel',{'Face','Body','House','Pattern','Tool'},'FontSize',10)
set(gca,'YTickLabel',{'Face','Body','House','Pattern','Tool'},'FontSize',10)
xlabel('Stimulus','FontSize',10)
ylabel('Stimulus','FontSize',10)

c = colorbar;
c.Label.String = 'Classifier Accuracy';

OneBackMDSTime = cputime-t;
disp(['Runtime: ' num2str(OneBackMDSTime)])