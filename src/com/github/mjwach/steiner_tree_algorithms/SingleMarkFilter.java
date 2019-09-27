package com.github.mjwach.steiner_tree_algorithms;

public class SingleMarkFilter implements MarkFilter
{
	private final int markIndex;

	public SingleMarkFilter(int markIndex)
	{
		this.markIndex = markIndex;
	}

	@Override
	public boolean passes(GraphVertex vertex)
	{
		return vertex.readMark(markIndex);
	}

	@Override
	public boolean passes(GraphEdge edge)
	{
		return edge.readMark(markIndex);
	}
}
