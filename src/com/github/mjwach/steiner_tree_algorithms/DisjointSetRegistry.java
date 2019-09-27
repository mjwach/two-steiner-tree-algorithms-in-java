package com.github.mjwach.steiner_tree_algorithms;

public class DisjointSetRegistry
{
	private int[] parents;
	private int[] ranks;

	public DisjointSetRegistry(int elementCount)
	{
		parents = new int[elementCount];
		ranks = new int[elementCount];
		
		for (int i = 0; i < elementCount; ++i)
			parents[i] = i;
	}
	
	public int find(int elementIndex)
	{
		if (parents[elementIndex] != elementIndex)
			parents[elementIndex] = find(parents[elementIndex]);
		
		return parents[elementIndex];
	}
	
	public void union(int element1Index, int element2Index)
	{
		int root1 = find(element1Index);
		int root2 = find(element2Index);
		
		if (root1 == root2)
			return;
		
		int r1 = ranks[root1];
		int r2 = ranks[root2];
		
		if (r1 < r2)
			parents[root1] = root2;
		else if (r2 < r1)
			parents[root2] = root1;
		else
		{
			parents[root1] = root2;
			++ranks[root2];
		}
	}
}
