% script to generate source time function

% gaussian (see fun_gauss in ./code/srcF/mod_src.F90 and SeisSource.conf)
t0=3.0; % force_stf_timefactor in SeisSource.conf
a=1.25; % force_stf_freqfactor in ...

dtrec=0.25; %time interval of observed time series

tsrc=[0:dtrec:2*t0];
nt=length(tsrc);
for i=1:nt
fsrc(i)=exp(-(tsrc(i)-t0)^2/(a^2))/(sqrt(pi)*a);
end

%plot(tsrc,f);
