package com.github.mjwach.steiner_tree_algorithms;

public interface MarkFilter
{
	boolean passes(GraphVertex vertex);
	boolean passes(GraphEdge edge);
}
