function add_state_patch(t, mask, color, yl, ~)
    dd = diff([false; mask; false]);
    on  = find(dd==1);
    off = find(dd==-1)-1;
    for i=1:numel(on)
        x0 = t(on(i)); x1 = t(off(i));
        patch([x0 x1 x1 x0], [yl(1) yl(1) yl(2) yl(2)], color, ...
              'EdgeColor','none', 'FaceAlpha', 0.12);
    end
end
