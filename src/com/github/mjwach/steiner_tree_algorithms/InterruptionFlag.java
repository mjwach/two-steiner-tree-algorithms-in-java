package com.github.mjwach.steiner_tree_algorithms;

public class InterruptionFlag implements InterruptionSignal
{
	private boolean isActive;

	public void activate()
	{
		isActive = true;
	}
	
	@Override
	public boolean isActive()
	{
		return isActive;
	}
}
