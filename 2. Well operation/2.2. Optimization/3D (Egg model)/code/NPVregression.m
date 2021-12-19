%%
figure;

plotregression2(-abd, -abc);
set(gca, 'FontWeight', 'Normal');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
axis square;
fig = gcf;
fig.PaperPosition = [5 5 8 8];
print('-r300', '-dpng', 'regressionplot150.png');