function [hi, eij, fi, fij, meta] = read_eij(paramfile)
%READ_EIJ reads a binary file of parameters for a pairwise maximum entropy
%  model (undirected graphical model) estimated by plm. The model describes
%  the distribution of sequences of length L drawn from an alphabet with q
%  characters. The outputs are:
%   
%   Object      Description                         Dimensions
%   hi          sitewise fields     hi(i,Ai)        L x q
%   eij         pairwise couplings  eij(i,j,Ai,Aj)  L x L x q x q
%   fi          sitewise marginals  fi(i,Ai)        L x q
%   fij         pairwise marginals  fij(i,j,Ai,Aj)  L x L x q x q
%   meta        metadata
%
%   Note that both the eij and fij arrays are output in dense form, but 
%   will be symmetric under (i,j,ai,aj) <-> (j,i,aj,ai)
%

PRECISION = 'single';

f_eij = fopen(paramfile, 'r');
%
L = fread(f_eij, 1, 'int');
q = fread(f_eij, 1, 'int');
target_seq = char(fread(f_eij, L, 'char'))';
offset_map = fread(f_eij, L, 'int');

fi = fread(f_eij, [q L], PRECISION)';
hi = fread(f_eij, [q L], PRECISION)';

%
% Collect metadata in structure
%
meta.N = L;
meta.target_seq = target_seq;
meta.offset_map = offset_map;

fij = zeros(L, L, q, q);
eij = zeros(L, L, q, q);
for i=1:(L-1)
    for j=(i+1):L
        ij = fread(f_eij, [2 1], 'int');
        fij(i,j,:,:) = fread(f_eij, [q q], PRECISION)';
        fij(j,i,:,:) = squeeze(fij(i,j,:,:))';
        eij(i,j,:,:) = squeeze(fread(f_eij, [q q], PRECISION))';
        eij(j,i,:,:) = squeeze(eij(i,j,:,:))';
    end
end

fclose(f_eij);
end