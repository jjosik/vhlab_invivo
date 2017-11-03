function [bs,bsimg] = rc_bootstrap_set(avgstim, numsims, alpha)
%RC_BOOTSTRAP - Bootstrap for reverse correlation responses
%  [BS,BSIMG] = RC_BOOTSTRAP(AVG_STIM, NUMSIMS)
%
%  Returns the percentage of bootstrap simulations that are positive for each
%  position in the reverse correlation "average_stimlus".  The AVG_STIM is 
%  assumed to be a 3-d matrix, with the first 2 dimensions corresponding to the
%  number of grid pixels and the reconstructed times, respectively, and the third
%  dimension corresponding to individual averages of short bouts of of the reverse correlation
%  stimulus.
%
%  BSIMG is a hot/cold RGB image with the same data.

bs = [];

N = size(avgstim,3);

for i=1:numsims,
	r = 1+floor(N*rand(1,N));
	bs_frame = mean(avgstim(:,:,r),3);
	if i==1,
		bs = double(bs_frame>0);
	else,
		bs = bs + double(bs_frame>0);
	end;
end;

bs = 1000 * bs / numsims;

p_bonf = max(1,1000*alpha/100); %max(0,1000* 0.05/(size(avgstim,1)*size(avgstim,2))),

map = gray(1000);
map(1:p_bonf,:) = repmat([0 0 1],length(1:p_bonf),1);
map(1000-p_bonf+1:1000,:) = repmat([1 0 0],length(1:p_bonf),1);


bsimg = ind2rgb(round(bs),map);
