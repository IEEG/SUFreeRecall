% Load xllog, the output of OneBack.m, and run this script to produce a 2-D
% visualization of the responses to the first three face and place
% exemplars across the population of units.
% Runtime: 4 sec on authors' hardware.

% Produces data for Figs. 2B-D.

% Copyright Simon Khuvis, 2019. See github.com/IEEG for licensing and use
% information.

tic
indx = @(vec,ind) vec(ind);

% Code numbers of stimuli
cats = unique(cell2mat(arrayfun(@(x) x.response(:,1),xllog,'Uni',0)));
% Input the test response and a training array with column 1 corresponding
% to code numbers and column 2 to responses. Performs Naive Bayes.

% Get the responses out of the xllog variables.
extractLog = @(xlog) arrayfun(@(x) x.response(:,[1 3-x.oncell]),xlog([xlog.class]>0),'Uni',0);

justlog = @(extracted) cellfun(@(x) [x(:,1) log10(x(:,2)+1)],extracted,'Uni',0);

expand = @(inp) cat(1,inp{:});
%% Plot Exemplars on Lower-Dimensional Space
rng(1)
font = 'Helvetica';
figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);

plotwhich = [101:103 301:303];
stimResp = cellfun(@(a) cell2mat(cellfun(@(z) z(1:min(min(cellfun(@length,a)),6)),a,'Uni',0)),...
    mat2cell(expand(cellfun(@(y) arrayfun(@(x) (y(y(:,1)==x,2).'-mean(y(ismember(y(:,1),plotwhich),2)))...
    /std(y(ismember(y(:,1),plotwhich),2)),cats,'Uni',0).',...
    justlog(extractLog(xllog)),'Uni',0)),63,ones(1,50)),'Uni',0);

subplot(3,2,[3 5])

plotwhich = 101:103;
distmat = cmdscale(squareform(pdist(cell2mat(stimResp(ismember(cats,plotwhich))).','euclidean')),2);
% distmat = pca(cell2mat(stimResp(ismember(cats,plotwhich))));
szs = [0 cumsum(cellfun('size',stimResp(ismember(cats,plotwhich)),2))];
colcodes = [100 23.9 37.3;99.6 52.9 61.2;30.2 7.8 11.8;35.3 19.6 22.4]/100;
markers = {'o','s','^','p'};
% scatter(distmat(:,1),distmat(:,2),50*reshape(meshgrid([1 2 6 1 2 6],1:4),[],1),[repmat([94.5 23.5 36.1],12,1);repmat([24.7 30.2 71.8],12,1)]/100)
hold on
for it = 1:3
    scatter(distmat(szs(it)+1:szs(it+1),1),...
        distmat(szs(it)+1:szs(it+1),2),150,colcodes(it,:),'filled',markers{it})
end
title('Multi-Dimensional Scaling: First Three Face Stimuli')
xlabel('Dimension 1 [AU]')
ylabel('Dimension 2 [AU]')
legend({'Face 01','Face 02','Face 03'})
legend('boxoff')

set(gca,'box','off')
set(gca,'TickLength',[.002 .002])

subplot(3,2,[4 6])

plotwhich = 301:303;
distmat = cmdscale(squareform(pdist(cell2mat(stimResp(ismember(cats,plotwhich))).','euclidean')),2);
% distmat = pca(cell2mat(stimResp(ismember(cats,plotwhich))));
szs = [0 cumsum(cellfun('size',stimResp(ismember(cats,plotwhich)),2))];
colcodes = [33.7 41.6 99.2;58 62.7 98;7.8 9.8 22;16.1 17.3 25.9]/100;
markers = {'o','s','^','p'};
% scatter(distmat(:,1),distmat(:,2),50*reshape(meshgrid([1 2 6 1 2 6],1:4),[],1),[repmat([94.5 23.5 36.1],12,1);repmat([24.7 30.2 71.8],12,1)]/100)
hold on
for it = 1:3
    scatter(distmat(szs(it)+1:szs(it+1),1),...
        distmat(szs(it)+1:szs(it+1),2),150,colcodes(it,:),'filled',markers{it})
end
title('Multi-Dimensional Scaling: First Three House Stimuli')
xlabel('Dimension 1 [AU]')
ylabel('Dimension 2 [AU]')
legend({'House 01','House 02','House 03'})
legend('boxoff')

set(gca,'box','off')
set(gca,'TickLength',[.002 .002])

subplot(3,2,[1 2])

plotwhich = [101:103 301:303];
distmat = cmdscale(squareform(pdist(cell2mat(stimResp(ismember(cats,plotwhich))).','euclidean')),2);
% distmat = pca(cell2mat(stimResp(ismember(cats,plotwhich))));
szs = [0 cumsum(cellfun('size',stimResp(ismember(cats,plotwhich)),2))];
colcodes = [100 23.9 37.3;99.6 52.9 61.2;30.2 7.8 11.8;35.3 19.6 22.4;...
    33.7 41.6 99.2;58 62.7 98;7.8 9.8 22;16.1 17.3 25.9]/100;
markers = {'o','s','^','o','s','^'};
% scatter(distmat(:,1),distmat(:,2),50*reshape(meshgrid([1 2 6 1 2 6],1:4),[],1),[repmat([94.5 23.5 36.1],12,1);repmat([24.7 30.2 71.8],12,1)]/100)
hold on
for it = 1:6
    scatter(distmat(szs(it)+1:szs(it+1),1),...
        distmat(szs(it)+1:szs(it+1),2),150,colcodes(it,:),'filled',markers{it})
end
title('Multi-Dimensional Scaling: First Three Face and House Stimuli')
xlabel('Dimension 1 [AU]')
ylabel('Dimension 2 [AU]')

set(gca,'box','off')
set(gca,'TickLength',[.002 .002])
toc