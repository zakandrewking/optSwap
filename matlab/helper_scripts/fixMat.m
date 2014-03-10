function ret = fixMat(aCell)
    aCell(:,end+1) = aCell(size(aCell,1),1);
    for i=1:size(aCell,1)
        t = aCell(i,1);
        x = t{1};
        f = x.f;
        aCell{i,end} = f;
    end
    ret = aCell;
end