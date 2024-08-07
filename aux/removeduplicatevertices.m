function XYsimplified = removeduplicatevertices(XY)
    isduplicated = [false; all(XY(1:end - 1, :) == XY(2:end, :), 2)];
    XYsimplified = XY(~isduplicated, :);
end
