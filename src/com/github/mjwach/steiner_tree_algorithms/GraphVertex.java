package com.github.mjwach.steiner_tree_algorithms;

public interface GraphVertex
{
	int getIndex();
	void setMark(int markIndex);
	boolean readMark(int markIndex);
	void clearMarkFromVerticesInConnectedSubgraph(int markIndex);
	void clearMarkFromEdgesInConnectedSubgraph(int markIndex);
	void clearMarkFromVerticesAndEdgesInConnectedSubgraph(int markIndex);
	GraphEdge findEdgeTo(GraphVertex peer);
	void removeEdge(GraphEdge edge);
	void removeAllEdges();
	GraphEdge addEdge(GraphVertex farEnd);
	void doForEachEdge(Action1<GraphEdge> action);
}
