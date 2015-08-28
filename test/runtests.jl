using Julimaps

filenames = filter(x -> x[end-2:end] == "jl", readdir("."))
for i in filenames
	include(i)
end