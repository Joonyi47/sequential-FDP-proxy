function [REP]=SolutionVisualization_wSOIL(REP, perm, filename)
    global N Nvar nx ny maxgen ...
           direc_fig ...
           ACTIVE

    gen         = REP.gen;
    
    % load PERM
%     PERMX.ref = GetPermx([permfile '.DATA']);
           
    active = reshape(ACTIVE, nx*ny,2);
    active = active(:,1).*active(:,2);
    active = reshape(active, nx, ny)';
    active(active == 0 ) = NaN;
    
    PERMX = zeros(nx*ny*2,10);
    PERMX(ACTIVE == 1,:) = perm;

    % Repository sorting by mean(NPV)
    [~, Index_sorted] = sort(REP.pos_fit(:,1),'ascend');
    REP.pos_sorted = REP.pos(Index_sorted,:);
    REP.pos_fit_sorted = [REP.pos_fit(Index_sorted,1)];

    % figure
    fig=figure();
    set(fig, 'Units', 'centimeters', 'Position', [3 3 8.3 8.3]);
    set(fig, 'Units', 'centimeters', 'Position', [3 3 6 5]);
    standardcoord = 0.05;
    lastblank     = 0.1;
    blank_space = 0.07;
    ncol = 1;   nrow = 1;
    sub_w = (1-((ncol-1)*blank_space+standardcoord+lastblank))/ncol;
    sub_h = (1-((nrow-1)*blank_space+standardcoord+lastblank))/nrow;
    mymap = repmat([0:0.25:1]',1,3);

    imagesc(reshape(mean(PERMX(1:nx*ny,:),2),ny,nx)'.*active);
    colormap( [1 1 1; parula(256)] );
    scatterSolution(REP.pos_sorted(1,1:2*N),REP.pos_sorted(1,2*N+1:Nvar));
    set(gca, 'XTick', [], 'YTick', []);
    set(gca, 'XLim', [0 60], 'YLim', [0 60]);
    set(gca, 'XTick', [0, 20, 40, 60], 'YTick', [0, 20, 40, 60]);
    set(gca, 'FontName', 'Arial', 'FontSize', 9);
%     set(gcf, 'color', [0.8 0.8 0.8]);
%     title(['NPV: ' num2str(-REP.pos_fit(:,1), '%10.4e')], 'FontSize', 20);
    xlabel('X Grid', 'FontName', 'Arial', 'FontSize', 9);
    ylabel('Y Grid', 'FontName', 'Arial', 'FontSize', 9);
    axis fill square;
    c = colorbar('eastoutside');
    caxis([0.1 0.9]);
    title(c, 'So', 'FontName', 'Arial', 'FontSize', 8);
    set(c, 'YLim', [0.1 0.9], 'YTick', [0.1 0.9]);
    tightfig(fig);
%     colormap jet
%     colormap(mymap)
%     colormap(flipud(colormap))
%             colorbar
    shading flat
    
    if gen < maxgen
        saveas(fig, [direc_fig '/Solution #' int2str(gen) '.png']);
        close(fig);
    else
        saveas(fig, [filename '.fig']); 
        print('-r300','-dpng',[filename '.png']);
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
             s=scatter(grid.x(1,i), grid.y(1,i),40,'filled');
%              s.MarkerEdgeColor = 'r';
             s.MarkerFaceColor = 'r';
             s.LineWidth = 2;
         elseif type(1,i) <= -1/3
             s=scatter(grid.x(1,i), grid.y(1,i),40,'filled');
%              s.MarkerEdgeColor = 'k';
             s.MarkerFaceColor = 'k';
             s.LineWidth = 2;
         end
     end
end