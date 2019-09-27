package com.github.mjwach.steiner_tree_algorithms;

public interface Graph
{
	GraphVertex getVertex(int index);
	void clearAllMarks();
	void markMinimalSpanningTree_Kruskal(int markIndex);
	ShortestPathsTable findShortestPathsForAllVertexPairs(InterruptionSignal interruptionSignal) throws ShortestPathsSearchInterruption;
	float[] getVerticesAsExesMadeFromLineSegments_xyxy(float exHalfWidth, Boolean requiredMark, int markIndex);
	float[] getEdgesAsLineSegments_xyxy(Boolean requiredMark, int markIndex);
	void addDelaunayTriangulationEdges();
	float getWeightOfMarkedTree(int markIndex);
}
