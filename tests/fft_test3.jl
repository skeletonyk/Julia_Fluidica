using BenchmarkTools, Compat
using FFTW
using
FFTW.set_num_threads(1)#Sys.CPU_CORES * 2)
n = 605
a = rand(n,n);
p = plan_rfft(a, flags = FFTW.MEASURE);
ip = plan_irfft(rfft(a), n, flags = FFTW.MEASURE);
