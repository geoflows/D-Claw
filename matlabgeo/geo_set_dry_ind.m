% geo_set_dry_ind sets dry_ind = NaN in cells where the
% depth is below cutoff and to 1 elsewhere.

if ~exist('dry_cutoff');
	dry_cutoff = 1.0e-4;
end

dry_ind=ones(size(X));
dry_ind(find((h2./dry_cutoff)<1))=NaN;

if ~exist('alpha_cutoff');
	alpha_cutoff = 0.5;
end

alpha_ind =ones(size(X));
alpha_ind(find((h2./alpha_cutoff)<1.))=(h2(find((h2./alpha_cutoff)<1.))./alpha_cutoff).^2;%.^2;