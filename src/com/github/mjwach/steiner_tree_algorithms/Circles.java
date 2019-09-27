package com.github.mjwach.steiner_tree_algorithms;

public class Circles
{
	public static SimpleCircle getSmallestCircleEnclosingTwoCircles(PlaneVector center0, float radius0, PlaneVector center1, float radius1)
	{
		double distanceBetweenCenters = center0.getDistanceFromPrecisely(center1);
		double radius = (distanceBetweenCenters + radius0 + radius1) * .5;
		double centerDistanceFromCenter0 = radius - radius0;
		PlaneVector center = center0.interpolatedToward(center1, centerDistanceFromCenter0 / distanceBetweenCenters);
		return new SimpleCircle(center, (float) radius);
	}
	
	public static SimpleCircle getSmallestCircleEnclosingAndTouchingThreeCircles(PlaneVector center0, float radius0, PlaneVector center1, float radius1, PlaneVector center2, float radius2)
	{
		/*
			look up the Problem of Apollonius for guidance
			
			vv scratchpad
			
			(xc - x1)^2 + (yc - y1)^2 = (rc - r1)^2
			(xc - x2)^2 + (yc - y2)^2 = (rc - r2)^2
			(xc - x3)^2 + (yc - y3)^2 = (rc - r3)^2
			
			xc^2 - 2xc*x1 + x1^2 + yc^2 - 2yc*y1 + y1^2 = rc^2 - 2rc*r1 + r1^2
			xc^2 - 2xc*x2 + x2^2 + yc^2 - 2yc*y2 + y2^2 = rc^2 - 2rc*r2 + r2^2
			
			-2xc*x1 + 2xc*x2 + x1^2 - x2^2 - 2yc*y1 + 2yc*y2 + y1^2 - y2^2 =
				-2rc*r1 + 2rc*r2 + r1^2 - r2^2
			-2(x1 - x2)xc - 2(y1 - y2)yc + x1^2 - x2^2 + y1^2 - y2^2 =
				-2(r1 - r2)rc + r1^2 - r2^2
			-2(x1 - x2)xc - 2(y1 - y2)yc = -2(r1 - r2)rc + r1^2 - r2^2 - x1^2 + x2^2 - y1^2 + y2^2
			
			2(x1 - x2)xc + 2(y1 - y2)yc = 2(r1 - r2)rc - r1^2 + r2^2 + x1^2 - x2^2 + y1^2 - y2^2
			2(x1 - x3)xc + 2(y1 - y3)yc = 2(r1 - r3)rc - r1^2 + r3^2 + x1^2 - x3^2 + y1^2 - y3^2
			
			a2*xc + b2*yc = c2
			a3*xc + b3*yc = c3
			
			a2 b2   xc   c2
			a3 b3 * yc = c3
			
			xc   (a2 b2 ^-1   c1
			yc =  a3 b3)    * c3
			
			(a b ^-1                    d -b
			 c d)    = 1 / (ad - bc) * -c  a
			 
			                           b3 -b2   A B
			1 / (a2 * b3 - a3 * b2) * -a3  a2 = C D
			 
			ei = -r1^2 + ri^2 + x1^2 - xi^2 + y1^2 - yi^2
			 
			xc = A * c1 + B * c3 = 2((r1 - r2)A + (r1 - r3)B)rc + A*e2 + B*e3 = E*rc + F
			yc = C * c1 + D * c3                                              = G*rc + H
			
			(E*rc + F - x1)^2 + (G*rc + H - y1)^2 = (rc - r1)^2
			J = F - x1
			K = H - y1
			(E*rc + J)^2 + (G*rc + K)^2 = (rc - r1)^2
			E^2*rc^2 + 2E*J*rc + J^2 + G^2*rc^2 + 2*G*K*rc + K^2 = rc^2 - 2*rc*r1 + r1^2
			(E^2 + G^2 - 1)rc^2 + 2(E*J + G*K + r1)rc + J^2 + K^2 - r1^2 = 0
			
			L = E^2 + G^2 - 1
			M = 2(E*J + G*K + r1)
			N = J^2 + K^2 - r1^2
			rc = (-M + Math.sqrt(M*M - 4*L*N)) / (2L)
			
			solve linearly for xc, yc
			put back into a topmost equation up there to get r
			put r into the xc, yc solutions (from the linear deal) to get xc, yc
		 */
		double x1 = center0.getX();
		double y1 = center0.getY();
		double r1 = radius0;
		double x2 = center1.getX();
		double y2 = center1.getY();
		double r2 = radius1;
		double x3 = center2.getX();
		double y3 = center2.getY();
		double r3 = radius2;
		
		double a2 = 2 * (x1 - x2);
		double b2 = 2 * (y1 - y2);
		double a3 = 2 * (x1 - x3);
		double b3 = 2 * (y1 - y3);
		double x12 = x1 * x1;
		double y12 = y1 * y1;
		double r12 = r1 * r1;
		double e2 = -r12 + r2*r2 + x12 - x2*x2 + y12 - y2*y2;
		double e3 = -r12 + r3*r3 + x12 - x3*x3 + y12 - y3*y3;
		double detInv = a2*b3 - a3*b2;
		double A = b3 / detInv;
		double B = -b2 / detInv;
		double C = -a3 / detInv;
		double D = a2 / detInv;
		double dr2 = r1 - r2;
		double dr3 = r1 - r3;
		double E = 2 * (dr2*A + dr3*B);
		double G = 2 * (dr2*C + dr3*D);
		double F = A*e2 + B*e3;
		double H = C*e2 + D*e3;
		double J = F - x1;
		double K = H - y1;
		double L = E*E + G*G - 1;
		double M = 2 * (E*J + G*K + r1);
		double N = J*J + K*K - r12;

		double discriminant = Math.sqrt(M*M - 4*L*N);
		double numerator1 = -M + discriminant;
		double numerator2 = -M - discriminant;
		double denominator = 2 * L;
		final double rc;
		if (denominator > 0 && numerator1 > 0 && numerator2 > 0 ||
			denominator < 0 && numerator1 < 0 && numerator2 < 0)
		{
			double rMax = Math.max(r1, Math.max(r2, r3));
			double rc1 = numerator1 / denominator;
			double rc2 = numerator2 / denominator;

			if (rc1 >= rMax)
				if (rc2 >= rMax)
					rc = Math.min(rc1, rc2);
				else
					rc = rc1;
			else
				rc = rc2;
		}
		else
			rc = (denominator > 0 == numerator1 > numerator2 ? numerator1 : numerator2) / denominator;
		/* aiming for the most positive of the two potential values of rc--if both are negative then this is harmless; and if one is
		positive then that one is surely right; and if both are positive (which does happen sometimes--note that the three equations
		that all this is based on do NOT capture the requirement that rc be greater than xi; their squaring of rc - ri discards that
		information, and that is why this choice arises) then the larger one is surely right--ah, but wait, there's at least one more
		case where both solutions may be positive:  when the three circles are all more or less in line, yet the middle one is fat
		and offset a bit and it's possible for a circle outside the three of them to touch all three in one of two ways...  well it's
		tricky to imagine without a picture, but here is an example in x, y, r,	x... format:
		-400.7035, -513.3982, 299.41135, -95.95776, -544.1875, 44.02644, -1176.4719, -952.30676, 59.30562, -946.177, 34.18043, 1072.3174, -643.36163, -750.94653, 629.1763
		--so anyway, the real strategy here is to pick the least solution that is no less than any of the input radii */
		
		double xc = E*rc + F;
		double yc = G*rc + H;
		
		return new SimpleCircle((float) xc, (float) yc, (float) rc); // well, that was monstrous
	}

	public static PlaneVector getCenterOfCircleThrough(PlaneVector point0, PlaneVector point1, PlaneVector point2)
	{
		/*
		basic idea:  points make a triangle; triangle's sides have perpendicular bisectors; circle's center is at the
				intersection of any two of these
		
		m1 = (x0 + x1) / 2
		m2 = (x0 + x2) / 2
		
		v1 = x1 - x0
		v2 = x2 - x0
		
		(h - m1) * v1 = 0
		(hx - m1x) * v1x + (hy - m1y) * v1y = 0
		hx v1x + hy v1y - m1x v1x - m1y v1y = 0
		hx v1x + hy v1y = m1x v1x + m1y v1y
		
		hx v1x + hy v1y = m1 * v1
		hx v2x + hy v2y = m2 * v2
		
		v1x v1y   hx   m1 * v1
		v2x v2y * hy = m2 * v2
		
		d = 1 / (v1x * v2y - v2x * v1y)
		hx        v2y -v1y   m1 * v1
		hy = d * -v2x  v1x * m2 * v2
		 */
		
        PlaneVector v1 = point1.minus(point0);
        PlaneVector v2 = point2.minus(point0);
        PlaneVector m1 = point1.averagedWith(point0);
        PlaneVector m2 = point2.averagedWith(point0);
		float v1m1 = v1.dot(m1);
		float v2m2 = v2.dot(m2);
		float d = 1 / (v1.crossZ(v2));

        return new PlaneVector((v2.getY() * v1m1 - v1.getY() * v2m2) * d, (-v2.getX() * v1m1 + v1.getX() * v2m2) * d);
	}

	public static boolean circleIntersectsCircle(PlaneVector location0, float radius0, boolean circle0HasBody, PlaneVector location1, float radius1, boolean circle1HasBody)
	{
		float centerDistanceSquared = location0.getDistanceSquaredFrom(location1);
		float radialSum = radius0 + radius1;
		boolean areClose = centerDistanceSquared <= radialSum * radialSum;
		
		if (areClose)
		{
			float radialDifference = radius0 - radius1;
		
			if (circle0HasBody)
				if (circle1HasBody)
					return true;
				else
					return radius0 >= radius1 || centerDistanceSquared >= radialDifference * radialDifference;
			else
				if (circle1HasBody)
					return radius0 <= radius1 || centerDistanceSquared >= radialDifference * radialDifference;
				else
					return centerDistanceSquared >= radialDifference * radialDifference;
		}
		else
			return false;
	}

	public static boolean circleIntersectsAxisAlignedBox(PlaneVector circleLocation, float circleRadius, boolean circleHasBody, float boxMinX, float boxMinY, float boxMaxX, float boxMaxY)
	{
		float x = circleLocation.getX();
		float y = circleLocation.getY();
		
		if (circleHasBody && x >= boxMinX && x <= boxMaxX && y >= boxMinY && y <= boxMaxY)
			return true;
		
		float halfWidth = (boxMaxX - boxMinX) * .5f;
		float halfHeight = (boxMaxY - boxMinY) * .5f;

		float cX = x - (boxMinX + halfWidth);
		float cY = y - (boxMinY + halfHeight);
		float r2 = circleRadius * circleRadius;
		
		if (circleHasBody && cX * cX + cY * cY <= r2)
			return true;
		
		return circlePerimeterCenteredOnCIntersectsAxisAlignedBoxCenteredOnOrigin(cX, cY, circleRadius, r2, halfWidth, halfHeight, true);
	}

	public static boolean circlePerimeterCenteredOnCIntersectsAxisAlignedBoxCenteredOnOrigin(float cX, float cY, float r, float r2, float hW, float hH, boolean boxHasBody)
	{
		if (cX >= hW)
			if (cY > hH)
				return d2(cX, cY, hW, hH) <= r2 && pointExtendsToCircle(-hW, -hH, cX, cY, r2);
			else if (cY < -hH)
				return d2(cX, cY, hW, -hH) <= r2 && pointExtendsToCircle(-hW, hH, cX, cY, r2);
			else
				return cX - hW <= r && (pointExtendsToCircle(-hW, hH, cX, cY, r2) || pointExtendsToCircle(-hW, -hH, cX, cY, r2));
		else if (cX <= -hW)
			if (cY > hH)
				return d2(cX, cY, -hW, hH) <= r2 && pointExtendsToCircle(hW, -hH, cX, cY, r2);
			else if (cY < -hH)
				return d2(cX, cY, -hW, -hH) <= r2 && pointExtendsToCircle(hW, hH, cX, cY, r2);
			else
				return -hW - cX <= r && (pointExtendsToCircle(hW, hH, cX, cY, r2) || pointExtendsToCircle(hW, -hH, cX, cY, r2));
		else if (cY > hH)
			return cY - hH <= r && (pointExtendsToCircle(-hW, -hH, cX, cY, r2) || pointExtendsToCircle(hW, -hH, cX, cY, r2));
		else if (cY < -hH)
			return -hH - cY <= r && (pointExtendsToCircle(-hW, hH, cX, cY, r2) || pointExtendsToCircle(hW, hH, cX, cY, r2));
		else
			if (boxHasBody)
				return true;
			else
			{
				boolean circleExtendsToBox = hH - cY <= r || cY + hH <= r || hW - cX <= r || cX + hW <= r;
	
				if (circleExtendsToBox)
				{
					boolean boxExtendsToCircle = pointExtendsToCircle(hW, hH, cX, cY, r2) || pointExtendsToCircle(-hW, hH, cX, cY, r2) || pointExtendsToCircle(-hW, -hH, cX, cY, r2) || pointExtendsToCircle(hW, -hH, cX, cY, r2);
					return boxExtendsToCircle;
				}
				else
					return false;
			}
	}

	private static float d2(float x0, float y0, float x1, float y1)
	{
		float dx = x0 - x1;
		float dy = y0 - y1;
		return dx * dx + dy * dy;
	}

	private static boolean pointExtendsToCircle(float x, float y, float circleCenterX, float circleCenterY, float circleRadiusSquared)
	{
		return d2(circleCenterX, circleCenterY, x, y) >= circleRadiusSquared;
	}

	public static float getDistanceBetweenCircles(PlaneVector center0, float radius0, boolean circle0HasBody, PlaneVector center1, float radius1, boolean circle1HasBody)
	{
		float d2 = center0.getDistanceSquaredFrom(center1);
		float radiusSum = radius0 + radius1;
		
		if (d2 > radiusSum * radiusSum)
			return sqrt(d2) - radiusSum;
		
		if (circle0HasBody)
			if (circle1HasBody)
				return 0;
			else
				return getDistance_circlesAreClose_AHasBody(d2, radius0, radius1);
		else
			if (circle1HasBody)
				return getDistance_circlesAreClose_AHasBody(d2, radius1, radius0);
			else
				return getDistanceBetweenClosePerimeters(d2, Math.abs(radius0 - radius1));
	}

	private static float getDistance_circlesAreClose_AHasBody(float d2, float radiusA, float radiusB)
	{
		float dr = radiusB - radiusA;
		
		if (dr <= 0)
			return 0;
		else
			return getDistanceBetweenClosePerimeters(d2, dr);
	}

	private static float getDistanceBetweenClosePerimeters(float d2, float dr)
	{
		if (d2 >= dr * dr)
			return 0;
		else
			return dr - sqrt(d2);
	}

	private static float sqrt(float number)
	{
		return (float) Math.sqrt(number);
	}

	public static boolean checkContainmentAndOrder(DCELGraph.Vertex v0, DCELGraph.Vertex v1, DCELGraph.Vertex v2, DCELGraph.Vertex v3)
	// true iff (the four points do not lie on the same circle) && ((v3 is strictly within the circle through v0, v1, and v2) ^
	// !(v0, v1, and v2 lie in counter-clockwise order on the circle))--but also, if the four points do not share a circle, then
	// transposing any two adjacent arguments (where v0 and v3 are considered adjacent) will simply reverse the result, meaning that
	// this test can be interpreted in a variety of ways
	//
	// the basic formula is taken from "Primitives for the Manipulation of General Subdivisions and the Computation of Voronoi
	// Diagrams" by Leonidas Guibas and Jorge Stolfi
	{
		// find the two points that are closest to one another by some reasonable, quick-to-compute measure
		float s01 = Math.abs(v0.x - v1.x) + Math.abs(v0.y - v1.y);
		float s02 = Math.abs(v0.x - v2.x) + Math.abs(v0.y - v2.y);
		float s03 = Math.abs(v0.x - v3.x) + Math.abs(v0.y - v3.y);
		float s12 = Math.abs(v1.x - v2.x) + Math.abs(v1.y - v2.y);
		float s13 = Math.abs(v1.x - v3.x) + Math.abs(v1.y - v3.y);
		float s23 = Math.abs(v2.x - v3.x) + Math.abs(v2.y - v3.y);
		
		float smallestDistance = s01; int smallestDistanceIndex = 0;
		if (smallestDistance > s02) {smallestDistance = s02; smallestDistanceIndex = 1;}
		if (smallestDistance > s03) {smallestDistance = s03; smallestDistanceIndex = 2;}
		if (smallestDistance > s12) {smallestDistance = s12; smallestDistanceIndex = 3;}
		if (smallestDistance > s13) {smallestDistance = s13; smallestDistanceIndex = 4;}
		if (smallestDistance > s23) {smallestDistance = s23; smallestDistanceIndex = 5;}
		
		// reorder the points so that v0 is one of the two close-together ones, making sure not to flip the result in the process
		float v0x, v0y, v1x, v1y, v2x, v2y, v3x, v3y;
		
		switch (smallestDistanceIndex)
		{
		case 0:
		case 1:
		case 2:
			v0x = v0.x; v0y = v0.y; v1x = v1.x; v1y = v1.y; v2x = v2.x; v2y = v2.y; v3x = v3.x; v3y = v3.y; break;
		case 3:
		case 4:
			v0x = v1.x; v0y = v1.y; v1x = v0.x; v1y = v0.y; v2x = v3.x; v2y = v3.y; v3x = v2.x; v3y = v2.y; break;
		case 5:
			v0x = v2.x; v0y = v2.y; v1x = v0.x; v1y = v0.y; v2x = v1.x; v2y = v1.y; v3x = v3.x; v3y = v3.y; break;
		default:
			throw new RuntimeException("this is apparently impossible; value was " + smallestDistanceIndex);
		}
		
		// shift the points so that v0 is the new origin (this could have been done without the previous reordering; but the
		// reordering will perhaps conserve precision by ensuring that any very small distances will end up in the densest
		// part of the floating-point value range...  but this hasn't been tested
		v1x -= v0x; v1y -= v0y;
		v2x -= v0x; v2y -= v0y;
		v3x -= v0x; v3y -= v0y;
		
		// apply the formula below, with simplifications afforded by the knowledge that v0 is (0, 0)
		float d1 = v1x * v1x + v1y * v1y; // TODO - somewhere around here, maybe a little earlier or later, check for especially big
		float d2 = v2x * v2x + v2y * v2y; // (or small??) coordinate values (probably after the subtraction of v0) and use something
		float d3 = v3x * v3x + v3y * v3y; // like BigDecimal with the original formula (below) to avoid infinite addends just below
		
		return
			(v2x*v1y - v1x*v2y) * d3 +
			(v1x*v3y - v3x*v1y) * d2 +
			(v3x*v2y - v2x*v3y) * d1 > 0;
		
		/*
		this is the basic formula; for many values it gives the correct result without any of the above acrobatics, but thanks
		to floating-point imprecision it's inaccurate for too many other values (and also, it doesn't even seem to be faster
		than the above)
		
		float d0 = v0.x * v0.x + v0.y * v0.y;
		float d1 = v1.x * v1.x + v1.y * v1.y;
		float d2 = v2.x * v2.x + v2.y * v2.y;
		float d3 = v3.x * v3.x + v3.y * v3.y;
		
		return
			(v0.x*v1.y - v1.x*v0.y) * (d2 - d3) +
			(v2.x*v0.y - v0.x*v2.y) * (d1 - d3) +
			(v0.x*v3.y - v3.x*v0.y) * (d1 - d2) +
			(v1.x*v2.y - v2.x*v1.y) * (d0 - d3) +
			(v3.x*v1.y - v1.x*v3.y) * (d0 - d2) +
			(v2.x*v3.y - v3.x*v2.y) * (d0 - d1) > 0;
		 */
	}
}