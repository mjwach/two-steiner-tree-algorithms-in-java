package com.github.mjwach.steiner_tree_algorithms;

public class OpenMarkFilter implements MarkFilter
{
	@Override
	public boolean passes(GraphVertex vertex)
	{
		return true;
	}

	@Override
	public boolean passes(GraphEdge edge)
	{
		return true;
	}
}
