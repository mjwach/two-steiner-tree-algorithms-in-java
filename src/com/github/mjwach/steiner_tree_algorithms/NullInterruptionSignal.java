package com.github.mjwach.steiner_tree_algorithms;

public class NullInterruptionSignal implements InterruptionSignal
{
	public static final NullInterruptionSignal instance = new NullInterruptionSignal();

	@Override
	public boolean isActive()
	{
		return false;
	}
}
