package com.github.mjwach.steiner_tree_algorithms;

public class SpaceVector implements EuclideanVector
{
	public static final SpaceVector zero = new SpaceVector();
	
	private final float x, y, z;

	public SpaceVector()
	{
		this(0, 0, 0);
	}
	
	public SpaceVector(float x, float y, float z)
	{
		this.x = x;
		this.y = y;
		this.z = z;
	}
	
	public SpaceVector(double x, double y, double z)
	{
		this.x = (float) x;
		this.y = (float) y;
		this.z = (float) z;
	}
	
	public SpaceVector(PlaneVector xy, float z)
	{
		this(xy.getX(), xy.getY(), z);
	}

	public float getX()
	{
		return x;
	}
	
	public float getY()
	{
		return y;
	}
	
	public float getZ()
	{
		return z;
	}
	
	public PlaneVector getXY()
	{
		return new PlaneVector(x, y);
	}
	
	public PlaneVector getXZ()
	{
		return new PlaneVector(x, z);
	}
	
	public PlaneVector getYZ()
	{
		return new PlaneVector(y, z);
	}

	public SpaceVector negated()
	{
		return new SpaceVector(-x, -y, -z);
	}
	
	public SpaceVector plus(SpaceVector peer)
	{
		return new SpaceVector(this.x + peer.x, this.y + peer.y, this.z + peer.z);
	}

	public SpaceVector plus(float x, float y, float z)
	{
		return new SpaceVector(this.x + x, this.y + y, this.z + z);
	}
	
	public SpaceVector plus(PlaneVector xy, float z)
	{
		return new SpaceVector(x + xy.getX(), y + xy.getY(), this.z + z);
	}

	public SpaceVector minus(SpaceVector peer)
	{
		return new SpaceVector(x - peer.x, y - peer.y, z - peer.z);
	}
	
	public SpaceVector minus(PlaneVector xy, float z)
	{
		return new SpaceVector(x - xy.getX(), y - xy.getY(), this.z - z);
	}

	public SpaceVector times(float multiplicand)
	{
		return new SpaceVector(x * multiplicand, y * multiplicand, z * multiplicand);
	}

	public SpaceVector averagedWith(SpaceVector peer)
	{
		return new SpaceVector((x + peer.x) * .5f, (y + peer.y) * .5f, (z + peer.z) * .5f);
	}

	public SpaceVector interpolatedToward(SpaceVector v, float t)
	{
		float ix = x + (v.x - x) * t;
		float iy = y + (v.y - y) * t;
		float iz = z + (v.z - z) * t;
		return new SpaceVector(ix, iy, iz);
	}

	public SpaceVector cross(SpaceVector peer)
	{
		return new SpaceVector(y * peer.z - peer.y * z, z * peer.x - peer.z * x, x * peer.y - peer.x * y);
	}
	
	public float getLengthSquared()
	{
		return x*x + y*y + z*z;
	}
	
	public float getLength()
	{
		return (float) Math.sqrt(getLengthSquared());
	}
	
	public SpaceVector scaled(float xScale, float yScale, float zScale)
	{
		return new SpaceVector(xScale * x, yScale * y, zScale * z);
	}

	public SpaceVector scaledTo(float length)
	{
		float factor = length / getLength();
		return new SpaceVector(factor * x, factor * y, factor * z);
	}
	
	@Override
	public float getDistanceSquaredFrom_(EuclideanVector peer)
	{
		return peer.getDistanceSquaredFrom(this);
	}
	
	@Override
	public float getDistanceSquaredFrom(SpaceVector v)
	{
		return getDistanceSquaredFrom(v.x, v.y, v.z);
	}
	
	public float getDistanceSquaredFrom(float x, float y, float z)
	{
		float dx = this.x - x;
		float dy = this.y - y;
		float dz = this.z - z;
		
		return dx*dx + dy*dy + dz*dz;
	}

	@Override
	public float getDistanceSquaredFrom(PlaneVector v)
	{
		throw new UnsupportedOperationException();
	}

	public float getDistanceFrom(SpaceVector v)
	{
		return (float) Math.sqrt(getDistanceSquaredFrom(v));
	}

	public float[] asArray()
	{
		return new float[]{x, y, z};
	}
	
	@Override
	public String toString()
	{
		StringBuilder result = new StringBuilder();
		result.append('<');
		result.append(x);
		result.append(", ");
		result.append(y);
		result.append(", ");
		result.append(z);
		result.append('>');
		return result.toString();
	}
	
	public static SpaceVector getIntersectionOfLineWithPlane(float lineBaseX, float lineBaseY, float lineBaseZ, float slopeX, float slopeY, float slopeZ, SpaceVector planeAnchor, SpaceVector planeNormal)
	{
		float lx = lineBaseX;
		float ly = lineBaseY;
		float lz = lineBaseZ;

		float sx = slopeX;
		float sy = slopeY;
		float sz = slopeZ;
		
		float ax = planeAnchor.getX();
		float ay = planeAnchor.getY();
		float az = planeAnchor.getZ();
		
		float nx = planeNormal.getX();
		float ny = planeNormal.getY();
		float nz = planeNormal.getZ();
		
		/*

p = l + st

(p-a)*n = 0
(l - a + st)*n = 0
(l-a)*n + (n*s)t = 0
t = (a-l)*n / s*n

		 */
		
		float t = ((ax-lx)*nx + (ay-ly)*ny + (az-lz)*nz) / (sx*nx + sy*ny + sz*nz);
		
		return new SpaceVector(lx + sx*t, ly + sy*t, lz + sz*t);
	}

	public float sumComponents()
	{
		return x + y + z;
	}

	public SpaceVector getPiecewiseMinimumWith(SpaceVector v)
	{
		return new SpaceVector(Math.min(x, v.x), Math.min(y, v.y), Math.min(z, v.z));
	}

	public SpaceVector getPiecewiseMaximumWith(SpaceVector v)
	{
		return new SpaceVector(Math.max(x, v.x), Math.max(y, v.y), Math.max(z, v.z));
	}

	public SpaceVector normalized()
	{
		float l = getLength();
		return new SpaceVector(x/l, y/l, z/l);
	}
}
