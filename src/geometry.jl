"""
Returns true if point, pt, is contained by the polygon defined by points in poly

Inputs:
  pt: 2-element array
  poly: Nx2 array as the vertices defining the polygon
          the last point is a repeat of the first
          i.e. poly[1,:] == poly[end,:]

Outputs:
  Boolean whether point is contained within the polygon

Pulled from:
https://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html

Based on Jordan curve theorem.

Draw a semi-infinite ray from the point (x infinite, y fixed) & count how many
edges that ray crosses (or actually, how many edges it is not parallel to). If
the ray crosses an odd number of edges when its slope is larger than the slope
of the edge, its contained. Only consider edges for which the point falls within
its y-range.

The semi-infinite ray method is inconsistent on boundaries, so include a test
for the border at the very end.
"""
function pt_in_poly(pt, poly)
  is_contained = false
  N = size(poly,1)
  # only count the distinct vertices
  if poly[1,:] == poly[N,:]
    N = N-1
  end
  j = N
  for i = 1:N
    i_vert = poly[i,:][:]
    j_vert = poly[j,:][:]
    # is point is on same side of edges endpoints in y?
    if (pt[2] > i_vert[2]) != (pt[2] > j_vert[2])
      # is the line from pt->i_vert parallel to j_vert->i_vert?
      if pt[1] < (j_vert[1]-i_vert[1])*(pt[2]-i_vert[2])/(j_vert[2]-i_vert[2]) + i_vert[1]
        is_contained = !is_contained
      end
    end
    j = i
  end

  if !is_contained
    j = N
    for i = 1:N
      i_vert = poly[i,:][:]
      j_vert = poly[j,:][:]
      on_line = pt_on_line_segment(pt, (i_vert, j_vert))
      if on_line
        return true
      end
    end
  end

  return is_contained
end

"""
If point exists along the line segment, the triangle it defines should be flat.
"""
function pt_on_line_segment(pt, line)
  epsilon = 1e-4
  a, b = line
  return -epsilon < norm(a-b) - (norm(a-pt) + norm(b-pt)) < epsilon
end

"""
Returns true if polyA contains polyB
"""
function poly_contains_poly(polyA, polyB)
  for i in 1:size(polyB, 1)
    pt = polyB[i,:][:]
    if pt_in_poly(pt, polyA)
      return true
    end
  end
  return false
end

"""
Returns true if polyA & polyB intersect
"""
function poly_intersects(polyA, polyB)
  return poly_contains_poly(polyA, polyB) || poly_contains_poly(polyB, polyA)
end

"""
Sutherland-Hodgman

* Assumes clip polygon is convex

https://en.wikipedia.org/wiki/Sutherland%E2%80%93Hodgman_algorithm
"""
function clip_polygon(subject, clip)
    output = [[subject[i,:]...] for i in 1:size(subject,1)]
    clip_edges = hcat(clip[1:end-1,:], clip[2:end,:])
    edge1 = line_to_vector(clip_edges[1,:])
    edge2 = line_to_vector(clip_edges[2,:])
    is_counterclockwise = cross2(edge1, edge2) > 0
    for i in 1:size(clip_edges,1)
        edge = clip_edges[i,:]
        input = output
        output = []
        start_pt = input[end]
        for j in 1:size(input,1)
            end_pt = input[j]
            if pt_on_left(end_pt, edge) == is_counterclockwise
                if pt_on_left(start_pt, edge) != is_counterclockwise
                    mid_pt = find_line_intersection([start_pt..., end_pt...], edge)
                    push!(output, mid_pt)
                end
                push!(output, end_pt)
            else
                if pt_on_left(start_pt, edge) == is_counterclockwise
                    mid_pt = find_line_intersection([start_pt..., end_pt...], edge)
                    push!(output, mid_pt)
                end                
            end
            start_pt = end_pt
        end
    end
    if length(output) > 1
      push!(output, output[1])
    end
    return hcat(output...)'
end

"""
Determine if point is to the left of an edge

"""
function pt_on_left(point, line)
    pt_vector = line_to_vector([line[1:2], point...])
    line_vector = line_to_vector(line)
    return cross2(line_vector, pt_vector) > 0
end

"""
Find intersecting point of two line segments

Inputs:
    * a: first line segment (array of two points)
    * b: second line segment (array of two points)
    
Outputs:
    * if the line segments intersect: point
    * if the line segments don't intersect: nothing

http://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
"""
function find_line_segment_intersection(a, b)
    p = a[1:2]
    r = line_to_vector(a)
    q = b[1:2]
    s = line_to_vector(b)
    if cross2(r, s) == 0
        if cross2(q-p, r) != 0
            # parallel, non-intersecting
            return nothing
        else
            # collinear
            t0 = (q-p)' * r/(r'*r)
            t1 = t0 + s' * r/(r'*r)
            if (s'*r)[1] < 0
                return p + t1[1]*r
            else
                return p + t0[1]*r
            end
        end
    else
        t = cross2(q-p, s / cross2(r, s))
        u = cross2(q-p, r / cross2(r, s))
        if (0 <= t <=1) & (0 <= u <= 1)
            # intersect @ one point
            return p + t*r
        else
            # do not intersect
            return nothing
        end
    end
end

"""
Find intersection of two lines - no protection
"""
function find_line_intersection(a, b)
    p = a[1:2]
    r = line_to_vector(a)
    q = b[1:2]
    s = line_to_vector(b)
    if cross2(r, s) == 0
      # collinear
      t0 = (q-p)' * r/(r'*r)
      t1 = t0 + s' * r/(r'*r)
      if (s'*r)[1] < 0
          return p + t1[1]*r
      else
          return p + t0[1]*r
      end
    else
      t = cross2(q-p, s / cross2(r, s))
      # u = cross2(q-p, r / cross2(r, s))
      return p + t*r
    end
end

"""
Convert 4-element line segment to 2-element vector
"""
function line_to_vector(a)
    return [a[3]-a[1], a[4]-a[2]]
end

"""
Given two 2-element vectors, return their determinant
"""
function cross2(a, b)
    return a[1]*b[2] - a[2]*b[1]
end

"""
Given an Nx2 list of convex polygon vertices, calculate polygon area
"""
function poly_area(a)
  N = size(a,1)
  if a[1,:] == a[N,:]
    N -= 1
  end
  area = 0
  for i = 1:N
    area += abs(cross2(a[i,:], a[N,:]))
  end
  return area / 2
end
