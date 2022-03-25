clear -a;
%ignore_function_time_stamp('all');
more off;

addpath('Data', 'Model', 'MOR', 'Iparams', 'Sparams', 'State', 'VFit');


disp('Reading arguments ...');
arg_list = argv ();
if nargin ~= 1
  disp('octave <me> prefix');
  return;
end

%prefix = arg_list{1};
%file = strcat('~/Desktop/Pruebas/', prefix, 'x', prefix, 'g_H.oct');

file = arg_list{1};

disp('Loading Data File ...');
load (file);

n = size(freqs,2);
P = size(Hs,1)
H = round(0.5 * P)
assert(P == size(Hs,2));
assert(n == size(Hs,3));

A = zeros(n, 1 + H * H);
A(:,1) = freqs';
k = 2;
for i = 1:H
  for j = H+1:P
    hi = reshape(Hs(i,j,:), [], 1);
    A(:,k) = 20 * log(abs(hi));
    k = k + 1;
  end
end

disp('Saving file ...');
file = strcat(prefix, 'x', prefix, 'g_full.data');
save('-ascii', file, 'A');
disp('');
