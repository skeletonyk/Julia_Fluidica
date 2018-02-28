using BenchmarkTools, Compat
using FFTW
FFTW.set_num_threads(Sys.CPU_CORES * 2)
n = 361
a = rand(n,n);
p = plan_rfft(a, flags = FFTW.MEASURE);
ip = plan_irfft(rfft(a), n, flags = FFTW.MEASURE);

tic()
for i = 1:10000
    ip*(p*a)
end
toc()
#pr = FFTW.plan_r2r(a, FFTW.REDFT10,1);
#@btime dst(dst($a)')'
#btime r2r($a, FFTW.REDFT10);
#@btime $pr * $a;
