function info = parse_info_from_filename(fname)
% Expect: <date>-<condition>-<mouse>-<APP_or_WT>_scores_1Hz.csv
% Returns only standardized fields used by META (no "raw" field).
info = struct('ok',false,'date','','condition','','mouse','','genotype','');
name = erase(fname, '_scores_1Hz.csv');
parts = split(name, '-');
if numel(parts) < 4, return; end
info.date      = string(parts{1});
info.condition = string(parts{2});
info.mouse     = string(parts{3});
geno = upper(string(parts{4}));
if contains(geno,'APP'), info.genotype = "APP";
elseif contains(geno,'WT'), info.genotype = "WT";
else, info.genotype = geno;
end
info.ok = true;
end
