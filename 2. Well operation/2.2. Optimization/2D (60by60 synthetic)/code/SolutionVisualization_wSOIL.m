function [REP]=SolutionVisualization_wSOIL(REP, perm, filename, tstep)
global N Nvar nx ny maxgen ...
    direc_fig


% load PERM
%     PERMX.ref = GetPermx([permfile '.DATA']);
PERMX = perm;


% Repository sorting by mean(NPV)

% figure
fig = figure();
set(fig, 'Units', 'centimeters', 'Position', [3 3 8.3 8.3]);
set(fig, 'Units', 'centimeters', 'Position', [3 3 6 5]);
standardcoord = 0.05;
lastblank     = 0.1;
blank_space = 0.07;
ncol = 1;   nrow = 1;
sub_w = (1-((ncol-1)*blank_space+standardcoord+lastblank))/ncol;
sub_h = (1-((nrow-1)*blank_space+standardcoord+lastblank))/nrow;
mymap = repmat([0:0.25:1]',1,3);

%     h = imagesc(log(reshape(mean(PERMX,2),ny,nx)'));
h = imagesc(reshape(mean(PERMX,2),ny,nx)');
scatterSolution(REP(1,1:2*N),REP(1,2*N+1:end));
set(gca, 'XTick', [], 'YTick', []);
set(gca, 'XLim', [0 60], 'YLim', [0 60]);
set(gca, 'XTick', [0, 20, 40, 60], 'YTick', [0, 20, 40, 60]);
set(gca, 'FontName', 'Arial', 'FontSize', 9);
%     set(gcf, 'color', [0.8 0.8 0.8]);
%     title(['NPV: ' num2str(-REP.pos_fit(:,1), '%10.4e')], 'FontSize', 20);
xlabel('X Grid', 'FontName', 'Arial', 'FontSize', 9);
ylabel('Y Grid', 'FontName', 'Arial', 'FontSize', 9);
axis fill square;
% title([int2str(tstep) ' days'], 'FontName', 'Arial', 'FontSize', 9);
c = colorbar('eastoutside');
caxis([0.25 0.75]);
title(c, 'So', 'FontName', 'Arial', 'FontSize', 8);
set(c, 'YLim', [0.25 0.75], 'YTick', [0.25 0.75]);
tightfig(fig);
set(gcf, 'color', 'w');
%     colormap jet
%     colormap(mymap)
%     colormap(flipud(colormap))
%         colorbar
%     shading flat


saveas(fig, [filename '.fig']); saveas(fig, [filename '.png']);
print('-r300','-dpng',[filename '.png'])
% close(fig);
% fclose all;

end


function scatterSolution(pos,type)
global N

pos    = reshape(pos, 2, N);
grid.x = pos(1,:);
grid.y = pos(2,:);

colorIndex  = {'r', 'b', [1 1 0.2]};

hold on
for i = 1:N
    if type(1,i) >= 1/3
        s=scatter(grid.x(1,i), grid.y(1,i),40,'filled');
        %              s.MarkerEdgeColor = 'w';
        s.MarkerFaceColor = 'r';
        s.LineWidth = 2;
    elseif type(1,i) <= -1/3
        s=scatter(grid.x(1,i), grid.y(1,i),40,'filled');
        %              s.MarkerEdgeColor = 'w';
        s.MarkerFaceColor = 'k';
        s.LineWidth = 2;
    end
end
end