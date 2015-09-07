using Julimaps
Pkg.build("Tk")

filenames = filter(x -> x[end-1:end] == "jl", readdir("test"))
for fn in filenames
	if fn != "runtests.jl"
		println(fn)
		include(fn)
	end
end