package com.github.mjwach.steiner_tree_algorithms;

import java.util.ArrayList;
import java.util.Arrays;

import gnu.trove.list.array.TIntArrayList;

public class ShortestPathsTable
{
	private final int vertexCount;
	private final float[] pathLengths;
	private final int[] penultimateSteps;
	private final int[] sortedPeers;
	private ShortestPathsSearch[] searches;

	private static final int nullIndex = Integer.MIN_VALUE;
	private static final int yetUnknownIndex = -1;
	
	public ShortestPathsTable(int vertexCount)
	{
		this(null, vertexCount);
	}
	
	public ShortestPathsTable(ShortestPathsSearch[] searches)
	{
		this(searches, searches.length);
	}

	private ShortestPathsTable(ShortestPathsSearch[] searches, int vertexCount)
	{
		this.vertexCount = vertexCount;
		int cellCount = vertexCount * vertexCount;

		pathLengths = new float[cellCount];
		penultimateSteps = new int[cellCount];
		sortedPeers = new int[cellCount];
		this.searches = searches;

		if (searches != null)
		{
			Arrays.fill(penultimateSteps, yetUnknownIndex);
			Arrays.fill(sortedPeers, yetUnknownIndex);
		}
	}

	private void ensurePenultimateStepValueIsPresent(int sourceIndex, int cellIndex) throws ShortestPathsSearchInterruption
	{
		while (penultimateSteps[cellIndex] == yetUnknownIndex)
			searches[sourceIndex].advance();
	}

	private void ensureSortedPeersValueIsPresent(int sourceIndex, int cellIndex) throws ShortestPathsSearchInterruption
	{
		while (sortedPeers[cellIndex] == yetUnknownIndex)
			searches[sourceIndex].advance();
	}

	public GraphVertex getNextVertexOnShortestPath(Graph graph, GraphVertex source, GraphVertex destination) throws ShortestPathsSearchInterruption
	{
		int index = getNextVertexOnShortestPath(source.getIndex(), destination.getIndex());
		return index >= 0 ? graph.getVertex(index) : null;
	}
	
	public int getNextVertexOnShortestPath(int sourceIndex, int destinationIndex) throws ShortestPathsSearchInterruption
	{
		// next step from source to destination is the next-to-last step from destination to source, so:
		int cellIndex = getCellIndex(destinationIndex, sourceIndex);
		ensurePenultimateStepValueIsPresent(destinationIndex, cellIndex);
		return penultimateSteps[cellIndex];
	}

	private int getCellIndex(int si, int di)
	{
		return si * vertexCount + di;
	}

	public void setCellValue(GraphVertex source, GraphVertex destination, float pathLength, GraphVertex penultimateStep, int indexInSortingOrder)
	{
		int cellIndex = getCellIndex(source.getIndex(), destination.getIndex());
		int sortingIndex = getCellIndex(source.getIndex(), indexInSortingOrder);
		
		pathLengths[cellIndex] = pathLength;
		penultimateSteps[cellIndex] = penultimateStep != null ? penultimateStep.getIndex() : nullIndex;
		sortedPeers[sortingIndex] = destination.getIndex();
	}

	public float getPathLength(GraphVertex source, GraphVertex destination) throws ShortestPathsSearchInterruption
	{
		return getPathLength(source.getIndex(), destination.getIndex());
	}

	public float getPathLength(int sourceIndex, int destinationIndex) throws ShortestPathsSearchInterruption
	{
		int cellIndex = getCellIndex(sourceIndex, destinationIndex);
		ensurePenultimateStepValueIsPresent(sourceIndex, cellIndex);
		return pathLengths[cellIndex];
	}

	public int[] getVerticesSortedByDistanceFrom(int index)
	{
		if (searches != null)
			throw new UnsupportedOperationException();
		
		int rowIndex = index * vertexCount;
		return Arrays.copyOfRange(sortedPeers, rowIndex, rowIndex + vertexCount);
	}
	
	// returns a list of all the paths stretching from a to b, where a and b are elements of endpointsToInclude.  if a path is included in the list, then the reverse
	// of it will not be included.  the returned list is formatted ababababab... and is sorted from shortest path to longest path.
	public TIntArrayList sortShortestPaths_abab(TIntArrayList endpointsToInclude) throws ShortestPathsSearchInterruption
	{
		if (endpointsToInclude.isEmpty())
			return new TIntArrayList();
		
		boolean[] inclusionFlags = new boolean[vertexCount];
		for (int i = 0; i < endpointsToInclude.size(); ++i)
		{
			int endpoint = endpointsToInclude.getQuick(i);
			if (inclusionFlags[endpoint])
				throw new IllegalArgumentException("endpoints are to be unique");
			else
				inclusionFlags[endpoint] = true;
		}
		
		ArrayList<TIntArrayList> listsToMerge = new ArrayList<>();
		
		for (int i = 0; i < endpointsToInclude.size(); ++i)
		{
			int source = endpointsToInclude.get(i);
			TIntArrayList list = new TIntArrayList();

			for (int start = getCellIndex(source, 0), end = start + vertexCount, j = start + 1; j < end; ++j)
			{
				ensureSortedPeersValueIsPresent(source, j);
				int destination = sortedPeers[j];
				
				int pathCellIndex = start + destination;
				if (penultimateSteps[pathCellIndex] == nullIndex)
					break;
				
				if (inclusionFlags[destination] && source < destination)
				{
					list.add(source);
					list.add(destination);
				}
			}
			
			listsToMerge.add(list);
		}

		while (listsToMerge.size() > 1)
		{
			ArrayList<TIntArrayList> mergedLists = new ArrayList<>();
			
			for (int i = 0; i < listsToMerge.size(); i += 2)
			{
				TIntArrayList list1 = listsToMerge.get(i);
				TIntArrayList list2 = null;
				
				int nextI = i + 1;
				if (nextI < listsToMerge.size())
					list2 = listsToMerge.get(nextI);
				
				TIntArrayList mergedList;
				
				if (list2 == null)
					mergedList = list1;
				else
				{
					mergedList = new TIntArrayList(list1.size() + list2.size());
					
					for (int i1 = 0, i2 = 0; i1 < list1.size() || i2 < list2.size();)
					{
						if (i1 == list1.size())
						{
							mergedList.add(list2.get(i2++));
							mergedList.add(list2.get(i2++));
						}
						else if (i2 == list2.size())
						{
							mergedList.add(list1.get(i1++));
							mergedList.add(list1.get(i1++));
						}
						else
						{
							int s1 = list1.get(i1);
							int d1 = list1.get(i1 + 1);
							float l1 = getPathLength(s1, d1);
							
							int s2 = list2.get(i2);
							int d2 = list2.get(i2 + 1);
							float l2 = getPathLength(s2, d2);
							
							if (l1 <= l2)
							{
								mergedList.add(s1);
								mergedList.add(d1);
								i1 += 2;
							}
							else
							{
								mergedList.add(s2);
								mergedList.add(d2);
								i2 += 2;
							}
						}
					}
				}
				
				mergedLists.add(mergedList);
			}
			
			listsToMerge = mergedLists;
		}
		
		return listsToMerge.get(0);
	}
}
