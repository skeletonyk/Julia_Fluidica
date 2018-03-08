push!(LOAD_PATH, pwd())

using Cylinder
println("======================= First run:")
job_run(3)

Profile.init(delay=0.02)
Profile.clear()
#clear_malloc_data()
#parms, job, solver,soln =
job_run(90)
#@btime job_run(40)
#@btime job_run(60)
#@btime job_run(80)

r = Profile.retrieve()
f = open("profile.bin", "w")
serialize(f, r)
close(f)

#using ProfileView
#f = open("profile.bin")
#r = deserialize(f);
#ProfileView.view()#r[1], lidict=r[2])
# Julia  30 - 1.75s / 60 -  7.66s /
# Matlab 30 - 2.22s / 60 - 10.00s /

# Desktop
# Julia 30 - 1.01 s / 60 - 4.58 s / 90 - 16.92 s
# Mat6ab
