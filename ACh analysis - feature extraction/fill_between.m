
function fill_between(x, y1, y2, alpha)
    if any(isnan(y1)) || any(isnan(y2)), hold on; return; end
    xv = [x; flipud(x)];
    yv = [y1; flipud(y2)];
    p = patch(xv, yv, [0 0 0], 'EdgeColor','none');
    if verLessThan('matlab','9.5'), set(p,'FaceAlpha',alpha); else, p.FaceAlpha = alpha; end
end

