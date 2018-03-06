
n = 100000;
a = rand(n)
@time tmp1 = a[2:n] - a[1:n-1]
@time diff(a)
