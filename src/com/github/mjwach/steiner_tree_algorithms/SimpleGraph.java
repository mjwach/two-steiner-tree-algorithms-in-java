package com.github.mjwach.steiner_tree_algorithms;

public class SimpleGraph
{
	private float[] x;
	private float[] y;
	private int[] directedEdges;
	private int[] edgePartitionsBySource;
	
	public void dofangle()
	{
		x = new float[]{-70, -30, 30, 60, 32, -40, -70};
		y = new float[]{40, 50, 50, 0, -30, -34, -2};
		directedEdges = new int[]{1, 2,   0, 2, 6,   0, 1,   5,   5, 6,   3, 4,   1, 4};
		edgePartitionsBySource = new int[]{0, 2, 5, 7, 8, 10, 12, 14};
		
		for (int i = 0; i < 7; i++)
		{
			x[i] += 200;
			y[i] += 200;
		}
	}

	public float[] getEdgesAsLineSegments_xyxy()
	{
		float[] result = new float[directedEdges.length * 2];
		int nextIndex = 0;
		
		for (int sourceIndex = 0; sourceIndex < edgePartitionsBySource.length - 1; sourceIndex++)
			for (int ei0 = edgePartitionsBySource[sourceIndex],
					 ei1 = edgePartitionsBySource[sourceIndex + 1],
					edgeIndex = ei0; edgeIndex < ei1; edgeIndex++)
			{
				int destinationIndex = directedEdges[edgeIndex];
				
				if (sourceIndex < destinationIndex)
				{
					result[nextIndex++] = x[sourceIndex];
					result[nextIndex++] = y[sourceIndex];
					result[nextIndex++] = x[destinationIndex];
					result[nextIndex++] = y[destinationIndex];
				}
			}
		
		return result;
	}

	public boolean didFangle()
	{
		return directedEdges != null;
	}
}
