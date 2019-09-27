package com.github.mjwach.steiner_tree_algorithms;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import gnu.trove.list.array.TIntArrayList;

public class UnfurlingPath
{
	private final TIntArrayList prefix;
	private final int destination;
	private final ShortestPathsTable paths;
	
	public UnfurlingPath(int source, int destination, ShortestPathsTable paths)
	{
		prefix = new TIntArrayList();
		prefix.add(source);
		
		this.destination = destination;

		this.paths = paths;
	}

	private int getStep(int index) throws ShortestPathsSearchInterruption
	{
		while (index >= prefix.size() && prefix.getQuick(prefix.size() - 1) != destination)
			prefix.add(paths.getNextVertexOnShortestPath(prefix.getQuick(prefix.size() - 1), destination));
		
		return index < prefix.size() ? prefix.getQuick(index) : Integer.MIN_VALUE;
	}

	public float getWeight(UnfurlingPath alreadyCountedNeighbor) throws ShortestPathsSearchInterruption
	{
		for (int i = 0; ; ++i)
		{
			int step = getStep(i);
			if (step < 0)
				return 0;
			else if (alreadyCountedNeighbor == null || step != alreadyCountedNeighbor.getStep(i))
			{
				int previousStep = getStep(Math.max(0, i - 1));
				return paths.getPathLength(previousStep, destination);
			}
		}
	}
	
	public static float getCombinedWeight(ArrayList<UnfurlingPath> outgoingPaths) throws ShortestPathsSearchInterruption
	{
		try
		{
			Collections.sort(outgoingPaths, new Comparator<UnfurlingPath>()
			{
				@Override
				public int compare(UnfurlingPath path0, UnfurlingPath path1)
				{
					try
					{
						for (int i = 0; ; ++i)
						{
							int thisStep = path0.getStep(i);
							int peerStep = path1.getStep(i);
							
							if (thisStep < peerStep)
								return -1;
							else if (thisStep > peerStep)
								return 1;
							
							if (thisStep < 0)
								return 0;
						}
					}
					catch (ShortestPathsSearchInterruption interruption)
					{
						throw new UncheckedShortestPathsSearchInterruption(interruption);
					}
				}
			});
		}
		catch (UncheckedShortestPathsSearchInterruption interruption)
		{
			interruption.rethrow();  // oh Java you wily beast
		}
		
		float weight = outgoingPaths.get(0).getWeight(null);
		
		for (int i = 1; i < outgoingPaths.size(); ++i)
		{
			UnfurlingPath thisPath = outgoingPaths.get(i);
			UnfurlingPath lastPath = outgoingPaths.get(i - 1);
			weight += thisPath.getWeight(lastPath);
		}
		
		return weight;
	}
}
