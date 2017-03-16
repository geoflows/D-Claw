
% set_covered_ind sets the array covered_ind to NaN for cells
% that are covered by a finer grid that will be plotted later.
% covered_ind = 1 in cells not covered.

level = amrdata(ng).level;
covered_ind=ones(size(X));
covered_indbig=covered_ind;
for og=1:ngrids
    levelog=amrdata(og).level;
    if ((levelog>level)&(PlotData(levelog)==1))
        xlowog=amrdata(og).xlow;
        ylowog=amrdata(og).ylow;
        xhighog=xlowog + amrdata(og).mx * amrdata(og).dx;
        yhighog=ylowog + amrdata(og).my * amrdata(og).dy;
        
        erasei=find(xedge>=xlowog-amrdata(og).dx & xedge<=xhighog+amrdata(og).dx);
        erasej=find(yedge>=ylowog-amrdata(og).dy& yedge<=yhighog+amrdata(og).dy);
        covered_ind([erasei(2:end-1)],[erasej(2:end-1)])=0.0;
    end
end
