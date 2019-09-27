package com.github.mjwach.steiner_tree_algorithms;

public interface Action2Throwing<Participant1, Participant2, ThrownException extends Throwable>
{
	void perform(Participant1 participant1, Participant2 participant2) throws ThrownException;
}
