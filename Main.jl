push!(LOAD_PATH, pwd())

using Cylinder
job_run(100)
#parms, job, solver = @fastmath job_run(100)
