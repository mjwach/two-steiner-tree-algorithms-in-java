package com.github.mjwach.steiner_tree_algorithms;

public interface GraphEdge
{
	void setMarkInBothDirections(int markIndex);
	void clearMarkInBothDirections(int markIndex);
	boolean readMark(int markIndex);
	GraphVertex getEnd0();
	GraphVertex getEnd1();
	GraphVertex getEnd(GraphVertex otherEnd);
	boolean isNotNegativelyOriented();
}
