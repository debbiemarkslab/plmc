function params = read_params(paramfile)
%READ_PARAMS reads a binary file of parameters for a pairwise maximum entropy
%  model (undirected graphical model) estimated by plmc. The model describes
%  the distribution of sequences of length L drawn from an alphabet with q
%  characters. The major parameter outputs are:
%   
%   Object      Description                         Type Dimensions
%   hi          sitewise fields     hi(i,Ai)        L x q
%   Jij         pairwise couplings  Jij(i,j,Ai,Aj)  L x L x q x q
%   fi          sitewise marginals  fi(i,Ai)        L x q
%   fij         pairwise marginals  fij(i,j,Ai,Aj)  L x L x q x q
%
%   Note that both the eij and fij arrays are output in dense form, but 
%   will be symmetric under (i,j,ai,aj) <-> (j,i,aj,ai)
%
PRECISION = 'single';
f_params = fopen(paramfile, 'r');
% Alignment dimensions
% 1: Number of sites
L = fread(f_params, 1, 'int');
% 2: Number of codes in alphabet
q = fread(f_params, 1, 'int');
% 3: Number of valid sequences in alignment
params.numSeqs = fread(f_params, 1, 'int');
% 4: Number of invalid sequences in alignment
params.numInvalidSeqs = fread(f_params, 1, 'int');
% 5: Number of iterations
params.numIter = fread(f_params, 1, 'int');
% 6: Number of iterations
params.theta = fread(f_params, 1, PRECISION);
% 7: Number of invalid sequences in alignment
params.lambda_h = fread(f_params, 1, PRECISION);
% 8: Number of invalid sequences in alignment
params.lambda_J = fread(f_params, 1, PRECISION);
% 9: Number of invalid sequences in alignment
params.lambda_group = fread(f_params, 1, PRECISION);
% 10: Number of invalid sequences in alignment
params.nEff = fread(f_params, 1, PRECISION)';
% 11: Number of invalid sequences in alignment
params.alphabet = char(fread(f_params, q, 'char'));
% 12: Number of invalid sequences in alignment
params.weights = fread(f_params, params.numSeqs + params.numInvalidSeqs, PRECISION);
% 13: Number of invalid sequences in alignment
params.target_seq = char(fread(f_params, L, 'char'))';
% 14: Number of invalid sequences in alignment
params.offset_map = fread(f_params, L, 'int');
% 15: Single-site marginals
params.fi = fread(f_params, [q L], PRECISION)';
% 16: Single-site biases
params.hi = fread(f_params, [q L], PRECISION)';
params.Jij = zeros(L, L, q, q);
% 17: Pairwise marginals
block = fread(f_params, L*(L-1)*q*q/2, PRECISION);
params.fij = zeros(L, L, q, q);
offset = 0;
for i=1:(L-1)
    for j=(i+1):L
        params.fij(j,i,:,:) = reshape(block(offset + [1:q*q]), [q q]);
        params.fij(i,j,:,:) = reshape(block(offset + [1:q*q]), [q q])';
        offset = offset + q*q;
    end
end
% 18: Couplings Jij
block = fread(f_params, L*(L-1)*q*q/2, PRECISION);
params.Jij = zeros(L, L, q, q);
offset = 0;
for i=1:(L-1)
    for j=(i+1):L
        params.Jij(j,i,:,:) = reshape(block(offset + [1:q*q]), [q q]);
        params.Jij(i,j,:,:) = reshape(block(offset + [1:q*q]), [q q])';
        offset = offset + q*q;
    end
end
fclose(f_params);
end