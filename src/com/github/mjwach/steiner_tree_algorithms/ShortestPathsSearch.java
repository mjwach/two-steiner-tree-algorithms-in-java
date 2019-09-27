package com.github.mjwach.steiner_tree_algorithms;

public interface ShortestPathsSearch
{
	void advance() throws ShortestPathsSearchInterruption;
	boolean isFinished();
}
