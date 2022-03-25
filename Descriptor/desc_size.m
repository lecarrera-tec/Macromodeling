function [nrows, ncols] = desc_size(DF)

% Usage [nrows, ncols] = desc_size(DF)
% Size of the descriptor system.

nrows = size(DF.C, 1);
ncols = size(DF.B, 2);

endfunction;
