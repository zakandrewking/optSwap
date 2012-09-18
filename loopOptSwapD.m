function loopOptSwapD

    % setup Cleaner
    % cleaner = onCleanup(@() cleanup);
    % global run status
    % status = 'starting';
    % run = '';
    
    sets = [
        10,0;
        5,0;
        3,1;
        3,2;
        3,3;
        3,4;
        3,5;
        3,10;
        3,20;
        5,1;
        5,2;
        5,3;
        5,4;
        5,5;
        5,10;
        5,20;];

    for i=1:size(sets,1)
        % for j=1:length(swaps)
        opt.kos = sets(i,1);
        opt.swaps = sets(i,2);
        runOptSwapD(opt);
        % end
    end

end