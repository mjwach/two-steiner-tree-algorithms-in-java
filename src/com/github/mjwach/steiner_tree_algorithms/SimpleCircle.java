package com.github.mjwach.steiner_tree_algorithms;

public class SimpleCircle
{
	public final PlaneVector center;
	public final float radius;

	public SimpleCircle(float centerX, float centerY, float radius)
	{
		this(new PlaneVector(centerX, centerY), radius);
	}
	
	public SimpleCircle(PlaneVector center, float radius)
	{
		this.center = center;
		this.radius = radius;
	}
	
	@Override
	public int hashCode()
	{
		return center.hashCode() * (Float.floatToIntBits(radius) ^ -2118117872);
	}
	
	@Override
	public String toString()
	{
		return "(" + center + ", " + radius + ")";
	}

	public SimpleCircle chooseLargerCircle(PlaneVector peerCenter, float peerRadius)
	{
		if (radius < peerRadius)
			return new SimpleCircle(peerCenter, peerRadius);
		else
			return this;
	}

	public SimpleCircle chooseLargerCircle(SimpleCircle peer)
	{
		if (radius < peer.radius)
			return peer;
		else
			return this;
	}
	
	public SimpleCircle chooseSmallerCircle(SimpleCircle peer)
	{
		if (radius > peer.radius)
			return peer;
		else
			return this;
	}
}
