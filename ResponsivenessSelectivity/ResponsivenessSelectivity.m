% Plots Responsiveness-Selectivity Relationship
indx = @(vec,ind) vec(ind,:);

extractAll = @(xlog) arrayfun(@(x) x.response(:,[1 3-x.oncell]),xlog,'Uni',0);
batches = unique({xllog.batchname},'stable');
exCount = cell2mat(cellfun(@(y) histc(y(:,1),unique(y(:,1))),...
    {xllog(cellfun(@(x) find(strcmp({xllog.batchname},x),1),batches)).response},'Uni',0));
plotInds = ~any(bsxfun(@minus,exCount,median(exCount,2)));

% Only plot units from those subjects who took the stardard version of the
% task.
logSameSize = xllog(~strcmp({xllog.batchname},batches(~plotInds)));
[~,anovatab] = cellfun(@(x) kruskalwallis(x(:,2),floor(x(:,1)/100),'off'),...
    extractAll(logSameSize),'Uni',0);
selectivity2 = cell2mat(cellfun(@(y) y(2,5),anovatab));

facesel = 2^.5*norminv(cellfun(@(x) sum(indx(tiedrank(x(:,2)),x(:,1)<200))/(sum(x(:,1)>200)*sum(x(:,1)<200))-...
    (1+sum(x(:,1)<200))/(2*sum(x(:,1)>200)),extractAll(logSameSize)),0,1);

font = 'Helvetica';
figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
subplot(2,1,1)
loglog([logSameSize.VisResponsiveness]',selectivity2,'.k','MarkerSize',20)

pft = polyfit(log10([logSameSize.VisResponsiveness])',log10(selectivity2),1);
xdat = get(gca,'XLim');
line(xdat,10.^polyval(pft,log10(xdat)),'Color',[0 0 0])

trueslope = diff(10.^polyval(pft,log10(xdat)))/diff(xdat);

[rho,p] = corr(log10([logSameSize.VisResponsiveness])',log10(selectivity2),'Type','Spearman');

text(xdat(1),get(gca,'YLim')*[0;1],{['   \rho = ' num2str(rho,2)],...
    ['  \sl p\rm = ' num2str(p,2)]},'HorizontalAlignment','left','VerticalAlignment','top')

ylabel('Category-selectivity index')
set(gca,'box','off')

subplot(2,1,2)
semilogx([logSameSize.VisResponsiveness],facesel,'.k','MarkerSize',20)

pft = polyfit(log10([logSameSize.VisResponsiveness])',facesel,1);
xdat = get(gca,'XLim');
line(xdat,polyval(pft,log10(xdat)),'Color',[0 0 0])

trueslope = diff(10.^polyval(pft,log10(xdat)))/diff(xdat);

[rho,p] = corr(log10([logSameSize.VisResponsiveness])',facesel,'Type','Spearman');

text(xdat(1),get(gca,'YLim')*[0;1],{['   \rho = ' num2str(rho,2)],...
    ['  \sl p\rm = ' num2str(p,2)]},'HorizontalAlignment','left','VerticalAlignment','top')

ylabel('D'' sensitivity for face images')
xlabel('Responsiveness index')

set(gca,'box','off')
set(gcf, 'Position', [0 0 1500 1000])
