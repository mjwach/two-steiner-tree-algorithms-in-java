package com.github.mjwach.steiner_tree_algorithms;

@SuppressWarnings("serial")
public class UncheckedShortestPathsSearchInterruption extends RuntimeException
{
	private ShortestPathsSearchInterruption interruption;

	public UncheckedShortestPathsSearchInterruption(ShortestPathsSearchInterruption interruption)
	{
		this.interruption = interruption;
	}

	public void rethrow() throws ShortestPathsSearchInterruption
	{
		throw interruption;
	}
}
