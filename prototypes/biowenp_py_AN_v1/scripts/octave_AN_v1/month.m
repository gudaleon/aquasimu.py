% Helper for the Octave equivalence run (MATLAB's month() is not in core Octave). Not part of BioWEN-P.
function m = month(d)
  v = datevec(d);
  m = v(:,2);
end
