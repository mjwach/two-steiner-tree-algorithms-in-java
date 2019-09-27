package com.github.mjwach.steiner_tree_algorithms;

import java.util.ArrayList;
import java.util.Collection;

public class SteinerSubgraphShortestPathTrees
{
	public static BestVertex bestVertex(final Collection<? extends VirtualGraphVertex> anchors, final ShortestPathsTable paths)
	{
		return anchors.size() <= 3 ?
				new BestVertex()
				{
					@Override
					public void update(GraphVertex v) throws ShortestPathsSearchInterruption
					{
						update(v, getWeightOfStarTree_le3(v, anchors, paths));
					}
				} :
				new BestVertex()
				{
					@Override
					public void update(GraphVertex v) throws ShortestPathsSearchInterruption
					{
						update(v, getWeightOfStarTree_gt3(v, anchors, paths));
					}
				};
	}

	public static float getWeightOfStarTree(GraphVertex root, Collection<? extends VirtualGraphVertex> leaves, ShortestPathsTable paths) throws ShortestPathsSearchInterruption
	{
		if (leaves.size() <= 3)
			return getWeightOfStarTree_le3(root, leaves, paths);
		else
			return getWeightOfStarTree_gt3(root, leaves, paths);
	}
	
	public static float getWeightOfStarTree_le3(GraphVertex root, final Collection<? extends VirtualGraphVertex> leaves, final ShortestPathsTable paths) throws ShortestPathsSearchInterruption
	{
		float weight = 0;

		for (VirtualGraphVertex anchor: leaves)
			weight += paths.getPathLength(root, anchor.asGraphVertex());
		
		return weight;
	}

	public static float getWeightOfStarTree_gt3(GraphVertex root, final Collection<? extends VirtualGraphVertex> leaves, final ShortestPathsTable paths) throws ShortestPathsSearchInterruption
	{
		ArrayList<UnfurlingPath> arms = new ArrayList<>();
		
		for (VirtualGraphVertex anchor: leaves)
			arms.add(new UnfurlingPath(root.getIndex(), anchor.asGraphVertex().getIndex(), paths));
		
		return UnfurlingPath.getCombinedWeight(arms);
	}

	public static float overestimateWeightOfTree(GenericShortestPathsGraphVertex root, ShortestPathsTable paths) throws ShortestPathsSearchInterruption
			// actually it might be the exact weight?  some overlapping of paths is accounted for here, but maybe not all...?  anyway the point is that, at least, this doesn't underestimate
	{
		return overestimateTreeWeight(root, null, paths);
	}

	private static float overestimateTreeWeight(GenericShortestPathsGraphVertex vertex, GenericShortestPathsGraphVertex source, ShortestPathsTable paths) throws ShortestPathsSearchInterruption
	{
		if (vertex.getNeighbors().size() == 1 && source != null)
			return 0;
		
		float weight = getWeightOfOutgoingPaths(vertex, paths);

		if (source != null)
			weight -= paths.getPathLength(source.asGraphVertex(), vertex.asGraphVertex());
		
		for (GenericShortestPathsGraphVertex neighbor: vertex.getNeighbors())
			if (neighbor != source)
				weight += overestimateTreeWeight(neighbor, vertex, paths);
		
		return weight;
	}

	private static float getWeightOfOutgoingPaths(GenericShortestPathsGraphVertex vertex, ShortestPathsTable paths) throws ShortestPathsSearchInterruption
	{
		return SteinerSubgraphShortestPathTrees.getWeightOfStarTree(vertex.asGraphVertex(), vertex.getNeighbors(), paths);
	}
}
