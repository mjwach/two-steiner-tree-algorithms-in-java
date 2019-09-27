package com.github.mjwach.steiner_tree_algorithms;

public interface GraphEdge
{
	void setMarkFully(int markIndex);
	void clearMarkFully(int markIndex);
	boolean readMark(int markIndex);
	GraphVertex getEnd0();
	GraphVertex getEnd1();
	boolean isNotNegativelyOriented();
}
