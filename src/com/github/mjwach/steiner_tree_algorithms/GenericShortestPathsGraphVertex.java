package com.github.mjwach.steiner_tree_algorithms;

import java.util.Collection;

public interface GenericShortestPathsGraphVertex extends VirtualGraphVertex
{
	Collection<? extends GenericShortestPathsGraphVertex> getNeighbors();
}
