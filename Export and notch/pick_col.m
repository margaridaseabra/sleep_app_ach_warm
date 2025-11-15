function c = pick_col(COL, st)
switch upper(st)
    case 'WK',   c = COL.WK;
    case 'NREM', c = COL.NREM;
    case 'REM',  c = COL.REM;
    case 'MA',   c = COL.MA;
    otherwise,   c = [0.5 0.5 0.5];
end
end
