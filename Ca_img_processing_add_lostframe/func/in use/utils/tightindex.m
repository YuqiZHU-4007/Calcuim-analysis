function ind = tightindex(vec)

[~, ind] = histc(vec, unique(vec));

