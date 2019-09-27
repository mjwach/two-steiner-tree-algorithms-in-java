package com.github.mjwach.steiner_tree_algorithms;

public class PlaneVector implements EuclideanVector
{
	private final float x;
	private final float y;

	public PlaneVector(float x, float y)
	{
		this.x = x;
		this.y = y;
	}

	public PlaneVector(double x, double y)
	{
		this((float) x, (float) y);
	}
	
	public float getX()
	{
		return x;
	}
	
	public float getY()
	{
		return y;
	}

	public PlaneVector minus(PlaneVector planeVector)
	{
		return new PlaneVector(x - planeVector.x, y - planeVector.y);
	}
	
	public PlaneVector minus(float x, float y)
	{
		return new PlaneVector(this.x - x, this.y - y);
	}
	
	public PlaneVector averagedWith(PlaneVector peer)
	{
		return new PlaneVector((x + peer.x) * .5f, (y + peer.y) * .5f);
	}

	public float dot(PlaneVector peer)
	{
		return x * peer.x + y * peer.y;
	}

	public float getDistanceFrom(PlaneVector peer)
	{
		float dx = x - peer.x;
		float dy = y - peer.y;
		return (float) Math.sqrt(dx * dx + dy * dy);
	}

	public float getDistanceSquaredFrom(PlaneVector peer)
	{
		return getDistanceSquaredFrom(peer.x, peer.y);
	}
	
	@Override
	public float getDistanceSquaredFrom(SpaceVector v)
	{
		throw new UnsupportedOperationException();
	}
	
	@Override
	public float getDistanceSquaredFrom_(EuclideanVector peer)
	{
		return peer.getDistanceSquaredFrom(this);
	}
	
	public float getDistanceSquaredFrom(float peerX, float peerY)
	{
		float dx = x - peerX;
		float dy = y - peerY;
		return dx * dx + dy * dy;
	}
	
	public double getDistanceFromPrecisely(PlaneVector peer)
	{
		double dx = x - peer.x;
		double dy = y - peer.y;
		return Math.sqrt(dx * dx + dy * dy);
	}
	
	public PlaneVector interpolatedToward(PlaneVector peer, float proportionateExtent)
	{
		return new PlaneVector(x + proportionateExtent * (peer.x - x), y + proportionateExtent * (peer.y - y));
	}

	public PlaneVector interpolatedToward(PlaneVector peer, double proportionateExtent)
	{
		return new PlaneVector(x + proportionateExtent * (peer.x - x), y + proportionateExtent * (peer.y - y));
	}

	public float crossZ(PlaneVector peer)
	{
		return x * peer.y - peer.x * y;
	}
	
	public static float crossZ(float x0, float y0, float x1, float y1)
	{
		return x0 * y1 - y0 * x1;
	}
	
	public static float crossZ(DCELGraph.Vertex v1, DCELGraph.Vertex v2, DCELGraph.Vertex v3) // == (va - vc) x (vb - vc) where (a, b, c) is one of {(1, 2, 3), (2, 3, 1), (3, 1, 2)}
	{
		return (v1.x - v3.x) * (v2.y - v3.y) - (v2.x - v3.x) * (v1.y - v3.y);
	}
	
	public static void getSteinerPoint(DCELGraph.Vertex v0, DCELGraph.Vertex v1, DCELGraph.Vertex v2, DCELGraph.Vertex repository)
	// this is the point from which segments stretching to the three vertices have a minimal cumulative length.
	// angles between these segments are bound to be 2pi/3 so the point must be on an arc of that length between any two of the
	// three vertices.  a solution, then, is to find centers and radii for two such arcs and then find the arcs' intersection.  but
	// actually that's only true when the three given points make a triangle whose angles are all less than 2pi/3; if they don't then
	// the Steiner point will be the one of the three given points that lies at the vertex of the largest angle in their triangle.
	// (in that case, one of the three segments whose lengths are to be minimized will have zero length).
	{
		// sort the points so that v0 ends up being the one that's furthest from the others; this should tend to increase the
		// accuracy of the ensuing circle-finding and intersection-finding operations for very thin triangles (maybe)
		float d02_2 = v0.getDistanceSquaredFrom(v2);
		if (d02_2 == 0) {repository.setLocation(v0); return;} // also this is a nice place to check for coincident input points,
		float d10_2 = v1.getDistanceSquaredFrom(v0);		  // which could produce dangerous NaN values in the calculations below
		if (d10_2 == 0) {repository.setLocation(v1); return;}
		float d21_2 = v2.getDistanceSquaredFrom(v1);
		if (d21_2 == 0) {repository.setLocation(v2); return;}
		
		if (d02_2 < d10_2)
		{
			DCELGraph.Vertex buffer = v2; v2 = v1; v1 = buffer;
			float buffer_ = d02_2; d02_2 = d10_2; d10_2 = buffer_;
		}
		
		if (d02_2 < d21_2)
		{
			DCELGraph.Vertex buffer = v0; v0 = v1; v1 = buffer;
			float buffer_ = d02_2; d02_2 = d21_2; d21_2 = buffer_;
		}
		
		if (d10_2 < d21_2)
		{
			DCELGraph.Vertex buffer = v0; v0 = v2; v2 = buffer;
			float buffer_ = d10_2; d10_2 = d21_2; d21_2 = buffer_;
		}
		
		// check for a very obtuse angle; the sorting will have placed it at v1
		float d10 = (float) Math.sqrt(d10_2);
		float d21 = (float) Math.sqrt(d21_2);
		
		float x10 = v1.x - v0.x, y10 = v1.y - v0.y;
		float x21 = v2.x - v1.x, y21 = v2.y - v1.y;
		
		if (x10 * x21 + y10 * y21 > d10 * d21 * .5f)
		{
			repository.setLocation(v1);
			return;
		}
		
		// find two circles through v0 on which the Steiner point must lie
		float d02 = (float) Math.sqrt(d02_2);

		float midpointX01 = .5f * (v0.x + v1.x), midpointY01 = .5f * (v0.y + v1.y);
		float midpointX02 = .5f * (v0.x + v2.x), midpointY02 = .5f * (v0.y + v2.y);
		
		float unitX01 = (v1.x - v0.x) / d10, unitY01 = (v1.y - v0.y) / d10; 
		float unitX02 = (v2.x - v0.x) / d02, unitY02 = (v2.y - v0.y) / d02;
		
		float centerOffset01Length = d10 * Numbers.oneOver2Sqrt3_f, centerOffset02Length = d02 * Numbers.oneOver2Sqrt3_f;
		float centerOffsetX01 = unitY01 * centerOffset01Length, centerOffsetY01 = -unitX01 * centerOffset01Length;
		float centerOffsetX02 = unitY02 * centerOffset02Length, centerOffsetY02 = -unitX02 * centerOffset02Length;
		
		if (centerOffsetX01 * x21 + centerOffsetY01 * y21 > 0)
		{ // these are supposed to point out from the triangle, not into it
			centerOffsetX01 = -centerOffsetX01;
			centerOffsetY01 = -centerOffsetY01;
		}
		
		if (centerOffsetX02 * x21 + centerOffsetY02 * y21 < 0)
		{
			centerOffsetX02 = -centerOffsetX02;
			centerOffsetY02 = -centerOffsetY02;
		}
		
		float center01X = midpointX01 + centerOffsetX01, center01Y = midpointY01 + centerOffsetY01;
		float center02X = midpointX02 + centerOffsetX02, center02Y = midpointY02 + centerOffsetY02;
		float radius01 = centerOffset01Length * 2;
		float radius02 = centerOffset02Length * 2;
		
		// find the circles' intersections and report the one that isn't v0
		float r0 = radius01, r1 = radius02, cx1, cy1;
		
		// cx0 = 0; // this calculation will use (center01X, center01Y) as its origin; this seems likely to increase precision, and
		// cy0 = 0; // having cx0 = cy0 = 0 simplifies the arithmetic (the calculation without this shift is shown in a comment below)
		cx1 = center02X - center01X;
		cy1 = center02Y - center01Y;
		
		boolean xySwitch = Math.abs(cx1) > Math.abs(cy1);

		if (xySwitch) // ...moreover, the "/ cy1" steps below seem to lead to dangerous imprecision when cy1 is much smaller than cx1
		{			// (as in cx1/cy1, which blows up for some ordinary cy1 values such as 0, but won't blow up for ANY ordinary
			float buffer = cx1; cx1 = cy1; cy1 = buffer; // cx1 value, given a not-too-small cy1), so here x and y are switched then
		}			// to make cx1 the smaller one, which changes nothing else about the circles and only requires a switch back later
		
		float C = ((r1 - r0) * (r0 + r1) - cx1*cx1 - cy1*cy1) * .5f;
		float D = -C / cy1;
		float E = cx1 / cy1;
		float A = 1 + E*E;
		float B = -2*E*D;
		float C_ = D*D - r0*r0;
		float discriminant = (float) Math.sqrt(B*B - 4*A*C_);
		float s0 = (-B + discriminant) / (2*A);
		float s1 = (-B - discriminant) / (2*A);
		float w0 = D - E*s0; // (s0, w0) and (s1, w1) are now the circles' two intersection points, in some coordinate plane;
		float w1 = D - E*s1; // now they just need to be translated back into the original coordinate system
		
		float x0, y0, x1, y1;
		
		if (xySwitch)
			{x0 = w0; y0 = s0; x1 = w1; y1 = s1;}
		else
			{x0 = s0; y0 = w0; x1 = s1; y1 = w1;}
		
		x0 += center01X; y0 += center01Y;
		x1 += center01X; y1 += center01Y; // this finishes the translation
		
		float dx0 = x0 - v0.x, dy0 = y0 - v0.y;
		float dx1 = x1 - v0.x, dy1 = y1 - v0.y;
		
		if (dx0*dx0 + dy0*dy0 > dx1*dx1 + dy1*dy1)
			repository.setLocation(x0, y0);
		else
			repository.setLocation(x1, y1);
		
		/*
			scratchpad:
		
			(x - cx0)^2 + (y - cy0)^2 = r0^2
			(x - cx1)^2 + (y - cy1)^2 = r1^2
			x^2 - 2x*cx0 + cx0^2 + y^2 -2y*cy0 + cy0^2 = r0^2
			x^2 - 2x*cx1 + cx1^2 + y^2 -2y*cy1 + cy1^2 = r1^2
			
			2x*cx1 - 2x*cx0 + cx0^2 - cx1^2 + 2y*cy1 - 2y*cy0 + cy0^2 - cy1^2 = r0^2 - r1^2
			2x(cx1 - cx0) + 2y(cy1 - cy0) + (cx0 - cx1)(cx0 + cx1) + (cy0 - cy1)(cy0 + cy1) = (r0 - r1)(r0 + r1)
			-2x(cx01) + -2y(cy01) + (cx01)(cx0 + cx1) + (cy01)(cy0 + cy1) = (r0 - r1)(r0 + r1)
			-2x(cx01) + -2y(cy01) = (r0 - r1)(r0 + r1) - (cx01)(cx0 + cx1) - (cy01)(cy0 + cy1)
			x(cx01) + y(cy01) = (-(r0 - r1)(r0 + r1) + (cx01)(cx0 + cx1) + (cy01)(cy0 + cy1)) / 2
				= C
			y = (C - cx01*x)/cy01
			D = C/cy01
			E = cx01 / cy01
			y = D - E*x
			
			x^2 - 2x*cx0 + cx0^2 + y^2 -2y*cy0 + cy0^2 = r0^2
			x^2 - 2cx0*x + cx0^2 + D^2 - 2EDx + E^2*x^2 - 2D*cy0 + 2E*cy0*x + cy0^2 = r0^2
			(1 + E^2)x^2 + 2(-cx0 - ED + E*cy0)x + cx0^2 + D^2 - 2D*cy0 + cy0^2 - r0^2 = 0
			Ax^2 + Bx + C_ = 0
			x = (-B +/- sqrt(B^2 - 4AC_)) / 2A
		 */
	}
	
	// octants of the plane and the rays between them; the octants are numbered in
	// counter-clockwise order starting from the x > 0, x > y octant
	private enum OctantOrBorder
	{ // the order of these matters
		border01234567, // this one is the origin
		border70, octant0, border01, octant1, border12, octant2, border23, octant3,
		border34, octant4, border45, octant5, border56, octant6, border67, octant7
	}
	
	private static OctantOrBorder getOctantOrBorder(float x, float y)
	{
		if (x > 0)
			if (y > 0)
				if (x > y)
					return OctantOrBorder.octant0;
				else if (x < y)
					return OctantOrBorder.octant1;
				else
					return OctantOrBorder.border01;
			else if (y < 0)
				if (x > -y)
					return OctantOrBorder.octant7;
				else if (x < -y)
					return OctantOrBorder.octant6;
				else
					return OctantOrBorder.border67;
			else
				return OctantOrBorder.border70;
		else if (x < 0)
			if (y > 0)
				if (x > -y)
					return OctantOrBorder.octant2;
				else if (x < -y)
					return OctantOrBorder.octant3;
				else
					return OctantOrBorder.border23;
			else if (y < 0)
				if (x > y)
					return OctantOrBorder.octant5;
				else if (x < y)
					return OctantOrBorder.octant4;
				else
					return OctantOrBorder.border45;
			else
				return OctantOrBorder.border34;
		else
			if (y > 0)
				return OctantOrBorder.border12;
			else if (y < 0)
				return OctantOrBorder.border56;
			else
				return OctantOrBorder.border01234567;
	}

	public static int compareVectorsByAngle(float x0, float y0, float x1, float y1)
	{
		OctantOrBorder octant0 = getOctantOrBorder(x0, y0);
		OctantOrBorder octant1 = getOctantOrBorder(x1, y1);
		
		int octantComparison = octant0.compareTo(octant1);
		
		if (octantComparison != 0)
			return octantComparison;
		
		switch (octant0)
		{
		case border01234567: return 0;
		case border70: case border01: return x0 < x1 ? -1 : x0 > x1 ? 1 : 0;
		case border12: case border23: return y0 < y1 ? -1 : y0 > y1 ? 1 : 0;
		case border34: case border45: return x0 > x1 ? -1 : x0 < x1 ? 1 : 0;
		case border56: case border67: return y0 > y1 ? -1 : y0 < y1 ? 1 : 0;
		case octant0:
		case octant1:
		case octant2:
		case octant3:
		case octant4:
		case octant5:
		case octant6:
		case octant7:
		{
			float z = crossZ(x0, y0, x1, y1);

			if (z > 0)
				return -1;
			else if (z < 0)
				return 1;
			
			float d2_0 = x0 * x0 + y0 * y0;
			float d2_1 = x1 * x1 + y1 * y1;
			
			return d2_0 < d2_1 ? -1 : d2_0 > d2_1 ? 1 : 0;
		}
		default:
			throw new RuntimeException("" + octant0);
		}
	}

	@Override
	public String toString()
	{
		return "<" + x  + ", " + y  + ">";
	}
}
