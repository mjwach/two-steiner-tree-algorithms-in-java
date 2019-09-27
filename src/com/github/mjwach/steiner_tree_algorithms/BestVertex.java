package com.github.mjwach.steiner_tree_algorithms;

public abstract class BestVertex
{
	public GraphVertex value;
	private float associatedWeight;
	
	public BestVertex()
	{
		this(null, Float.MAX_VALUE);
	}
	
	public BestVertex(GraphVertex initialValue, float initialAssociatedWeight)
	{
		value = initialValue;
		associatedWeight = initialAssociatedWeight;
	}

	public abstract void update(GraphVertex v) throws ShortestPathsSearchInterruption;
	
	protected void update(GraphVertex v, float weight)
	{
		if (weight < associatedWeight)
		{
			value = v;
			associatedWeight = weight;
		}
	}

	@SuppressWarnings("unchecked")
	public <Vertex extends GraphVertex> Vertex getValue(Vertex typeRepresentative)
	{
		return (Vertex) value;
	}
}
