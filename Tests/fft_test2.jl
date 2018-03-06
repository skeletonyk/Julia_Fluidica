using BenchmarkTools, Compat
using FFTW
FFTW.set_num_threads(1)#Sys.CPU_CORES * 2)
n = 605
a = rand(n,n);
p = plan_rfft(a, flags = FFTW.MEASURE);
ip = plan_irfft(rfft(a), n, flags = FFTW.MEASURE);

tic()
for i = 1:100
    @time ip*(p*a)
end
toc()
#pr = FFTW.plan_r2r(a, FFTW.REDFT10,1);
#@btime dst(dst($a)')'
#btime r2r($a, FFTW.REDFT10);
#@btime $pr * $a;
