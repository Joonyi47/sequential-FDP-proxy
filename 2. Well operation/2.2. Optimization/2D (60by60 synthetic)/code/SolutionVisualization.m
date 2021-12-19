function [REP]=SolutionVisualization(REP, perm, filename)
    global N Nvar nx ny maxgen ...
           direc_fig

    gen         = REP.gen;
    
    % load PERM
%     PERMX.ref = GetPermx([permfile '.DATA']);
    PERMX = perm;
        

    % Repository sorting by mean(NPV)
    [~, Index_sorted] = sort(REP.pos_fit(:,1),'ascend');
    REP.pos_sorted = REP.pos(Index_sorted,:);
    REP.pos_fit_sorted = [REP.pos_fit(Index_sorted,1)];

    % figure
    fig = figure();
    set(fig, 'Units', 'centimeters', 'Position', [3 3 8.3 8.3]);
    standardcoord = 0.05;
    lastblank     = 0.1;
    blank_space = 0.07;
    ncol = 1;   nrow = 1;
    sub_w = (1-((ncol-1)*blank_space+standardcoord+lastblank))/ncol;
    sub_h = (1-((nrow-1)*blank_space+standardcoord+lastblank))/nrow;
    mymap = repmat([0:0.25:1]',1,3);

    h = imagesc(log(reshape(mean(PERMX,2),ny,nx)'));
    scatterSolution(REP.pos_sorted(1,1:2*N),REP.pos_sorted(1,2*N+1:Nvar));
    set(gca, 'XTick', [], 'YTick', []);
    set(gca, 'XTick', [20, 40, 60], 'YTick', [20, 40, 60]);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 10);
%     set(gcf, 'color', [0.8 0.8 0.8]);
%     title(['NPV: ' num2str(-REP.pos_fit(:,1), '%10.4e')], 'FontSize', 20);
    xlabel('X Grid', 'FontName', 'Times New Roman', 'FontSize', 10);
    ylabel('Y Grid', 'FontName', 'Times New Roman', 'FontSize', 10);
    axis fill tight square;
    c = colorbar('eastoutside');
    title(c, 'ln(md)', 'FontName', 'Times New Roman', 'FontSize', 10);
    tightfig(fig);
%     colormap jet
%     colormap(mymap)
%     colormap(flipud(colormap))
    %         colorbar
%     shading flat
    
    if gen < maxgen
        saveas(fig, [direc_fig '/Solution #' int2str(gen) '.png']);
        close(fig);
    else
        saveas(fig, [filename '.fig']); saveas(fig, [filename '.png']);
        print('-r300','-dpng',[filename '.png'])
        close(fig);
    end
    fclose all;
    
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
             s=scatter(grid.x(1,i), grid.y(1,i),80,'filled');
%              s.MarkerEdgeColor = 'w';
             s.MarkerFaceColor = 'r';
             s.LineWidth = 2;
         elseif type(1,i) <= -1/3
             s=scatter(grid.x(1,i), grid.y(1,i),80,'filled');
%              s.MarkerEdgeColor = 'w';
             s.MarkerFaceColor = 'k';
             s.LineWidth = 2;
         end
     end
end