function str = p_to_star(p)
if p < 0.001; str = '***';
elseif p < 0.01; str = '**';
elseif p < 0.05; str = '*';
else; str = 'ns';
end
end