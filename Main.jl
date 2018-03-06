push!(LOAD_PATH, pwd())
using Cylinder
println("======================= First run:")
job_run(3)

Profile.init(delay=0.02)
Profile.clear()
#clear_malloc_data()
parms, job, solver = job_run(50)

r = Profile.retrieve()
f = open("profile.bin", "w")
serialize(f, r)
close(f)

#using ProfileView
#f = open("profile.bin")
#r = deserialize(f);
#ProfileView.view()#r[1], lidict=r[2])
