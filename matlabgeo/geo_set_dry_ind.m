% geo_set_dry_ind sets dry_ind = NaN in cells where the
% depth is below cutoff and to 1 elsewhere.

cutoff= 1.0e-4;
dry_ind=ones(size(X));
dry_ind(find((h2./cutoff)<1))=NaN;
