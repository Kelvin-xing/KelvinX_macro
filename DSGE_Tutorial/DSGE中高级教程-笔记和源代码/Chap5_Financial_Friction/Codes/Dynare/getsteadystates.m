Ne = M_.orig_endo_nbr;
for indexvar = 1:Ne
    varname = deblank(M_.endo_names(indexvar,:));
        eval([  varname '= oo_.steady_state(' int2str(indexvar) ');']);
end

