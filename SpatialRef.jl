using Base.Test

type SpatialRef
	x::Int
	y::Int
	w::Int
	h::Int
end

SpatialRef() = SpatialRef(0,0,0,0)

function +(srA::SpatialRef, srB::SpatialRef)
	x = min(srA.x, srB.x)
	y = min(srA.y, srB.y)
	w = max(srA.w+srA.x, srB.w+srB.x) - x
	h = max(srA.h+srA.y, srB.h+srB.y) - y
	return SpatialRef(x,y,w,h)
end

function test()
	A = SpatialRef(10,10,30,20)
	B = SpatialRef(5,20,20,10)
	C = A+B
	D = SpatialRef(5,10,35,20)
	@test C.x == D.x
	@test C.y == D.y
	@test C.w == D.w
	@test C.h == D.h

	@test_throws InexactError SpatialRef(0,0,-1,-1)
end