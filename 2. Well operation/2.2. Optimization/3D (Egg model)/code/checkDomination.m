function dom_vector = checkDomination(fitness, violation)
    Np = size(fitness,1);
    dom_vector = zeros(Np,1);
    all_perm = nchoosek(1:Np,2);    % Possible permutations
    all_perm = [all_perm; [all_perm(:,2) all_perm(:,1)]];
    x = [fitness(all_perm(:,1),:), violation(all_perm(:,1),:)];
    y = [fitness(all_perm(:,2),:), violation(all_perm(:,2),:)];
    d = dominates(x,y);
    dominated_particles = unique(all_perm(d==1,2));
    dom_vector(dominated_particles) = 1;
end