package com.github.mjwach.steiner_tree_algorithms;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map.Entry;

import com.github.mjwach.steiner_tree_algorithms.EdgeArraysGraph.Edge;
import com.github.mjwach.steiner_tree_algorithms.EdgeArraysGraph.Vertex;

public class SubgraphSteinerVertexCache // this class handles some of the tasks involved in finding a Steiner subgraph...  it's called a cache because it was originally intended to cache some of the results of its searches, so as to save work; but
		// it turned out that that probably wouldn't actually help much (or at least, it looked like it wouldn't).  so now this is not really a cache; it's just a helper class that COULD cache some of the data passing through it.  renaming it might be
		// more honest.
		// --some work was done in here to actually perform caching, but what is cached is just hit counts recording the amount of repetition of queries, not actual query results.  so simply reactivating the obviously inactive chunks of caching code here
		// will not make the caching work, since full caching was not implemented.  the cache key classes included here should be useful in caching, but the actual caching code would need to be newly written, if caching were desired.  (the incomplete
		// caching code here was written for the purpose of evaluating the prospect of caching.)
{
	private final HashMap<ApproximateMidpointQuery, Integer> approximateMidpoints;
	private final HashMap<ExactSteinerPointQuery, Integer> exactSteinerPoints;
	private final HashMap<VertexToPathDistanceQuery, Integer> vertexToPathDistances;
	
	public SubgraphSteinerVertexCache()
	{
		approximateMidpoints = null;
		exactSteinerPoints = null;
		vertexToPathDistances = null;
	}
	
	private static class ApproximateMidpointQuery
	{
		public Vertex va;
		public Vertex vb;
		public int ratioDenominator;
		
		@Override
		public boolean equals(Object object)
		{
			ApproximateMidpointQuery peer = (ApproximateMidpointQuery) object;
			return peer.va == va && peer.vb == vb && peer.ratioDenominator == ratioDenominator;
		}
		
		@Override
		public int hashCode()
		{
	        int result = 486187739;
	        result = (result * 16777619) ^ va.getIndex();
	        result = (result * 16777619) ^ vb.getIndex();
	        result = (result * 16777619) ^ ratioDenominator;
	        return result;
		}
	}
	
	@SuppressWarnings("unused")
	public Vertex getApproximateMidpoint(Vertex va, Vertex vb, int ratioDenominator, ShortestPathsTable paths, Graph graph) throws ShortestPathsSearchInterruption
	{
		if (ratioDenominator == 2 && va.getIndex() > vb.getIndex())
		{
			Vertex v = va;
			va = vb;
			vb = v;
		}
		
		if (false)
		{
			ApproximateMidpointQuery query = new ApproximateMidpointQuery();
			query.va = va;
			query.vb = vb;
			query.ratioDenominator = ratioDenominator;
			Integer count = approximateMidpoints.get(query);
			if (count == null) count = 0;
			approximateMidpoints.put(query, count + 1);
		}
		
		float ratio = 1.0f / ratioDenominator;
		float fullDistance = paths.getPathLength(vb, va);
		float distanceToMidpoint = fullDistance * ratio;
		
		Vertex pv = va;
		Vertex nv = va;
		float lastStep = 0;
		while (distanceToMidpoint > 0 && nv != vb)
		{
			pv = nv;
			nv = (Vertex) paths.getNextVertexOnShortestPath(graph, pv, vb);
			lastStep = paths.getPathLength(pv, nv);
			distanceToMidpoint -= lastStep;
		}
		
		if (distanceToMidpoint < -.5f * lastStep)
			return pv;
		else
			return nv;
	}
	
	private static class ShortestPathsGraphEdge implements Comparable<ShortestPathsGraphEdge>
	{
		public int i1, i2;
		private float length;
		
		@Override
		public int compareTo(ShortestPathsGraphEdge peer)
		{
			return length < peer.length ? -1 : length > peer.length ? 1 : 0;
		}
	}
	
	private static class ShortestPathsGraphVertex implements GenericShortestPathsGraphVertex
	{
		public final Vertex v;
		public final ArrayList<ShortestPathsGraphVertex> neighbors;

		public ShortestPathsGraphVertex(Vertex vertex)
		{
			v = vertex;
			neighbors = new ArrayList<>();
		}
		
		@Override
		public GraphVertex asGraphVertex()
		{
			return v;
		}
		
		@Override
		public Collection<? extends GenericShortestPathsGraphVertex> getNeighbors()
		{
			return neighbors;
		}

		public void update(BestVertex bestVertex, float searchRadius, int markIndex, ShortestPathsTable paths) throws ShortestPathsSearchInterruption
		{
			Vertex searchRoot = v;

			ArrayList<Vertex> searchFringe = new ArrayList<>();  // doing this recursively is arguably simpler, but it has sometimes led to unnecessary stack overflow
			ArrayList<Vertex> nextSearchFringe = new ArrayList<>();
			
			searchFringe.add(searchRoot);
			searchRoot.setMark(markIndex);
			
			while (!searchFringe.isEmpty())
			{
				for (Vertex currentVertex: searchFringe)
				{
					float distance = paths.getPathLength(currentVertex, searchRoot);
					if (distance <= searchRadius)
					{
						bestVertex.update(currentVertex);
						
						for (Edge edge: currentVertex.edges)
						{
							Vertex neighbor = edge.getDestination(currentVertex);
							if (!neighbor.readMark(markIndex))
							{
								nextSearchFringe.add(neighbor);
								neighbor.setMark(markIndex);
							}
						}
					}
				}
				
				ArrayList<Vertex> newFringe = searchFringe;
				newFringe.clear();
				
				searchFringe = nextSearchFringe;
				nextSearchFringe = newFringe;
			}
		}

		public void clearMarks(int markIndex)
		{
			v.clearConnectedSubgraphMark(markIndex);
		}
	}
	
	private static class ExactSteinerPointQuery
	{
		public ArrayList<Vertex> vertices;
		
		@Override
		public boolean equals(Object object)
		{
			ExactSteinerPointQuery peer = (ExactSteinerPointQuery) object;
			
			if (vertices.size() != peer.vertices.size())
				return false;
			
			for (int i = 0; i < vertices.size(); ++i)
				if (vertices.get(i) != peer.vertices.get(i))
					return false;
			
			return true;
		}
		
		@Override
		public int hashCode()
		{
	        int result = 486187739;

	        for (Vertex vertex: vertices)
	        	result = (result * 16777619) ^ vertex.getIndex();
	        
	        return result;
		}
	}

	@SuppressWarnings("unused")
	public Vertex getExactSteinerPoint(final ArrayList<Vertex> vertices, int workingMarkIndex, final ShortestPathsTable paths) throws ShortestPathsSearchInterruption
	{
		if (vertices.size() < 3)
			throw new IllegalArgumentException();
		
		Collections.sort(vertices, new Comparator<Vertex>()
		{
			@Override
			public int compare(Vertex v1, Vertex v2)
			{
				int i1 = v1.getIndex();
				int i2 = v2.getIndex();
				return i1 < i2 ? -1 : i1 > i2 ? 1 : 0;
			}
		});
		
		if (false)
		{
			ExactSteinerPointQuery query = new ExactSteinerPointQuery();
			query.vertices = vertices;
			
			Integer count = exactSteinerPoints.get(query);
			if (count == null) count = 0;
			exactSteinerPoints.put(query, count + 1);
		}
		
		ArrayList<ShortestPathsGraphVertex> tree = new ArrayList<>();

		{
			ArrayList<ShortestPathsGraphEdge> edges = new ArrayList<>();
			for (int i = 0; i < vertices.size() - 1; ++i)
				for (int j = i + 1; j < vertices.size(); ++j)
				{
					ShortestPathsGraphEdge edge = new ShortestPathsGraphEdge();
					edge.i1 = i;
					edge.i2 = j;
					edge.length = paths.getPathLength(vertices.get(i), vertices.get(j));
					edges.add(edge);
				}
			Collections.sort(edges);
			
			DisjointSetRegistry sets = new DisjointSetRegistry(vertices.size());
			for (Vertex vertex: vertices) tree.add(new ShortestPathsGraphVertex(vertex));
			
			int edgeCount = 0;
			for (ShortestPathsGraphEdge edge: edges)
			{
				int a = edge.i1;
				int b = edge.i2;
	
				if (sets.find(a) != sets.find(b))
				{
					sets.union(a, b);
					ShortestPathsGraphVertex va = tree.get(a);
					ShortestPathsGraphVertex vb = tree.get(b);
					va.neighbors.add(vb);
					vb.neighbors.add(va);
					
					if (++edgeCount == tree.size() - 1)
						break;
				}
			}
		}
		
		float weight = SteinerSubgraphShortestPathTrees.overestimateWeightOfTree(tree.get(0), paths);
		
		float minimalFinalEdgeWeightRelatedQuantity = Float.MAX_VALUE;
		Edge lightestFinalEdge = null;
		for (int i = 1; i < vertices.size(); ++i)
		{
			Vertex v = vertices.get(i);
			for (Edge edge: v.edges)
			{
				float q = edge.getWeightRelatedQuantityForComparison();
				if (q < minimalFinalEdgeWeightRelatedQuantity)
				{
					minimalFinalEdgeWeightRelatedQuantity = q;
					lightestFinalEdge = edge;
				}
			}
		}
		float minimalFinalEdgeWeight = lightestFinalEdge.getWeight(minimalFinalEdgeWeightRelatedQuantity);

		float searchRadius = weight - 63.0f / 64 * minimalFinalEdgeWeight;
		
		BestVertex bestVertex = SteinerSubgraphShortestPathTrees.bestVertex(vertices, paths);
		
		tree.get(0).update(bestVertex, searchRadius, workingMarkIndex, paths);
		
		for (int i = 1; i < vertices.size(); ++i)
			bestVertex.update(vertices.get(i));
		
		tree.get(0).clearMarks(workingMarkIndex);
		
		return (Vertex) bestVertex.value;
	}
	
	private static class VertexToPathDistanceQuery
	{
		public Vertex vertex;
		public Vertex pathStart;
		public Vertex pathEnd;
		
		@Override
		public boolean equals(Object object)
		{
			VertexToPathDistanceQuery peer = (VertexToPathDistanceQuery) object;
			return vertex == peer.vertex && pathStart == peer.pathStart && pathEnd == peer.pathEnd;
		}
		
		@Override
		public int hashCode()
		{
	        int result = 486187739;
	        result = (result * 16777619) ^ vertex.getIndex();
	        result = (result * 16777619) ^ pathStart.getIndex();
	        result = (result * 16777619) ^ pathEnd.getIndex();
	        return result;
		}
	}

	@SuppressWarnings("unused")
	public float getDistanceFromVertexToNearestVertexOnPathThatIsNotTheStartOfThePath(Vertex vertex, Vertex pathStart, Vertex pathEnd, ShortestPathsTable paths) throws ShortestPathsSearchInterruption
	{
		if (false)
		{
			VertexToPathDistanceQuery query = new VertexToPathDistanceQuery();
			query.vertex = vertex;
			query.pathStart = pathStart;
			query.pathEnd = pathEnd;
			
			Integer count = vertexToPathDistances.get(query);
			if (count == null) count = 0;
			vertexToPathDistances.put(query, count + 1);
		}
		
		int shortenedPathEnd = pathEnd.getIndex();
		int shortenedPathStart = paths.getNextVertexOnShortestPath(pathStart.getIndex(), shortenedPathEnd);
		
		float ds = paths.getPathLength(vertex.getIndex(), shortenedPathStart);
		float de = paths.getPathLength(vertex.getIndex(), shortenedPathEnd);

		if (ds < de)
			return getDistanceFromVertexToNearestVertexOnPath(vertex.getIndex(), shortenedPathStart, shortenedPathEnd, de, paths);
		else
			return getDistanceFromVertexToNearestVertexOnPath(vertex.getIndex(), shortenedPathEnd, shortenedPathStart, ds, paths);
	}

	private float getDistanceFromVertexToNearestVertexOnPath(int vertex, int nearerPathEnd, int furtherPathEnd, float distanceToFurtherEnd, ShortestPathsTable paths) throws ShortestPathsSearchInterruption
	{
		float minimalDistance = Float.MAX_VALUE;
		
		int v = nearerPathEnd;
		float minimalFutureD;
		do
		{
			float d = paths.getPathLength(vertex, v);
			
			if (d < minimalDistance)
				minimalDistance = d;
			
			v = paths.getNextVertexOnShortestPath(v, furtherPathEnd);
			
			if (v < 0)
				break;
			
			float dToEnd = paths.getPathLength(v, furtherPathEnd);
			minimalFutureD = distanceToFurtherEnd - dToEnd;
		}
		while (minimalFutureD < minimalDistance);
		
		return minimalDistance;
	}

	@SuppressWarnings("unused")
	public void reportOnUsage()
	{
		if (false) reportOnUsage("AM", approximateMidpoints); // in casual tests, average was about 5, and it appeared that probably mainly small searches were getting the most cache hits per query
		if (false) reportOnUsage("ESP", exactSteinerPoints); // in casual tests, average was about 2, and it appeared that probably mainly small searches were getting the most cache hits per query
		if (false) reportOnUsage("VTPD", vertexToPathDistances); // in casual tests, average was about 6 with net crossing depth set to 1 (higher settings make for a substantially higher average)
	}

	private <Query> void reportOnUsage(String cacheName, HashMap<Query, Integer> cache)
	{
		System.out.println(cacheName + " cache size = " + cache.size());
		
		int total = 0;
		for (Entry<Query, Integer> entry: cache.entrySet())
		{
			Integer count = entry.getValue();
			System.out.print("" + count + " ");
			total += count;
		}
		System.out.println();
		
		System.out.println("total = " + total + "; average = " + total / (float) cache.size());
	}
}
