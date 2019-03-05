function [] = draw_nz_stable(psize)
%--------------------------------------------------------------------------
% [] = draw_nz_stable(psize)
% Scatter plot of the stable and long lived isotopes.
% Inputs>  psize: point size
% Outputs: none
%--------------------------------------------------------------------------
  file_form='%f %f %*f %*q %*f';
  [zstable,astable]=textread('../../../initial_abundance/anders_grevesse_89.txt',file_form,286,'headerlines',5)
  nstable=astable - zstable
  scatter(nstable,zstable,psize,'kx')
end
