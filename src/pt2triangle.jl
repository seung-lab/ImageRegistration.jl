# Author: Thomas Macrina
# Email: tmacrina@princeton.edu
# Date: 150807
#
# Detect which triangle in a mesh a point belongs to

using GeometricalPredicates
using Base.Test

function verts2triangle(verts)
# Use 3x2 list of coordinates to build a triangle
	Primitive(Point(verts[1,:]...),
			Point(verts[2,:]...),
			Point(verts[3,:]...))
end

function pt2triangle(pt, vertices, triangles)
# Returns index of first triangle to contain point & barycentric coordinates
	for i = 1:size(triangles,1)	
		verts = vertices[vec(triangles[i,:]), :]
		triangle = verts2triangle(verts)
		if intriangle(triangle, Point(pt...)) > 0
			w = xy2barycentric(pt, verts)
			return i, tuple(w...)
		end
	end
	return 0, (0, 0, 0)
end

function xy2barycentric(pt, verts)
# Convert Cartesian coordinate to triple of barycentric coefficients
	@assert size(verts, 1) == 3
	R = hcat(verts, ones(Int64, 3, 1))
	r = hcat(pt, 1)
	return r * R^-1
end

function pt2triangle2(pt, vertices, triangles)
	# Use randomized edge hopping to find containing triangle
	# for a, b, c in triangles: # need to loop through neighbors, not this
	# 	triangle = Primitive(Point(a...),
	# 							Point(b...),
	# 							Point(c...))
	# 	o = intriangle(triangle, pt)
	# 	if o == -1
	# 		# opposite side a
	# 	elseif o == -2
	# 		# opposite side b
	# 	elseif o == -3
	# 		# opposite side c
	# 	else
	# 		# contained in triangle
	# 	end
	# end
end