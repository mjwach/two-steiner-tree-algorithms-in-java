package com.github.mjwach.steiner_tree_algorithms;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Deque;
import java.util.Iterator;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.Stack;

import org.jheaps.AddressableHeap.Handle;
import org.jheaps.tree.PairingHeap;

import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;

public class EdgeArraysGraph implements Graph, Iterable<EdgeArraysGraph.Vertex>
{
	private final ArrayList<Vertex> vertices;
	
	public static interface EdgeFactory
	{
		Edge act(Vertex end0, Vertex end1);
	}
	
	public static interface ActionOnVertexOrEdge
	{
		void perform(Vertex vertex);
		void perform(Edge edge);
	}
	
	public static class Vertex implements GraphVertex, VirtualGraphVertex
	{
		public int index;
		public final ArrayList<Edge> edges;
		private int mark;
		
		private static final EdgeFactory defaultEdgeFactory = new EdgeFactory()
		{
			@Override
			public Edge act(Vertex end0, Vertex end1)
			{
				return new EuclideanEdge((EuclideanVertex) end0, (EuclideanVertex) end1);
			}
		};
		
		private Vertex(int index)
		{
			this.index = index;
			edges = new ArrayList<>();
		}
		
		@Override
		public String toString()
		{
			return Integer.toString(index);
		}
		
		@Override
		public GraphVertex asGraphVertex()
		{
			return this;
		}

		@Override
		public int getIndex()
		{
			return index;
		}

		@Override
		public void setMark(int markIndex)
		{
			mark |= 1 << markIndex;
		}

		@Override
		public boolean readMark(int markIndex)
		{
			return (mark & 1 << markIndex) != 0;
		}
		
		public void clearMark(int markIndex)
		{
			mark &= ~(1 << markIndex);
		}
		
		public void clearAllMarks()
		{
			mark = 0;
		}
		
		@Override
		public GraphEdge addEdge(GraphVertex farEnd)
		{
			return addEdge(farEnd, defaultEdgeFactory);
		}

		public GraphEdge addEdge(GraphVertex farEnd, EdgeFactory edgeFactory)
		{
			Vertex v = (Vertex) farEnd;
			Edge edge = edgeFactory.act(this, v);
			edges.add(edge);
			v.edges.add(edge);
			return edge;
		}
		
		@Override
		public void doForEachEdge(Action1<GraphEdge> action)
		{
			for (Edge edge: edges)
				action.perform(edge);
		}

		@Override
		public GraphEdge findEdgeTo(GraphVertex peer)
		{
			for (Edge edge: edges)
				if (edge.getDestination(this) == peer)
					return edge;
			
			return null;
		}
		
		public boolean hasEdgeTo(GraphVertex peer)
		{
			return findEdgeTo(peer) != null;
		}

		public int getEdgeCount()
		{
			return edges.size();
		}
		
		@Override
		public void removeEdge(GraphEdge edge)
		{
			Edge e = (Edge) edge;
			e.getDestination(this).edges.remove(e);
			edges.remove(e);
		}

		@Override
		public void removeAllEdges()
		{
			for (Edge edge: edges)
			{
				Vertex destination = edge.getDestination(this);
				if (destination != this)
					destination.edges.remove(edge);
			}
			
			edges.clear();
		}

		public void clearConnectedSubgraphMark(int markIndex)
		{
			if (readMark(markIndex))
			{
				clearMark(markIndex);
				for (Edge edge: edges)
					edge.getDestination(this).clearConnectedSubgraphMark(markIndex);
			}
		}

		public EdgeArraysGraph buildShortCutsGraph(InterestEvaluation interestEvaluation, MarkFilter markFilter, int graphwideVertexCount, Sleeve<ShortcutGraphVertex[]> oldVertexToNewVertexMap_out)
		{
			EdgeArraysGraph graph = new EdgeArraysGraph();
			ShortcutGraphVertex[] vertexCache = new ShortcutGraphVertex[graphwideVertexCount];
			buildShortCutsGraph(interestEvaluation, markFilter, null, null, graph, vertexCache);
			
			if (oldVertexToNewVertexMap_out != null) oldVertexToNewVertexMap_out.setIdentity(vertexCache);
			return graph;
		}

		private void buildShortCutsGraph(InterestEvaluation interestEvaluation, MarkFilter markFilter, ShortcutGraphVertex source, ArrayList<Edge> pathFromSource, EdgeArraysGraph graph, ShortcutGraphVertex[] vertexCache)
		{
			Edge sourceEdge = pathFromSource != null ? pathFromSource.get(pathFromSource.size() - 1) : null;
			boolean cycleIsCompletedHere = false;
			
			if (interestEvaluation.perform(this))
			{
				ShortcutGraphVertex destination = vertexCache[getIndex()];
				if (destination == null)
					destination = vertexCache[getIndex()] = new ShortcutGraphVertex(this, graph.vertices.size());
				else
					cycleIsCompletedHere = true;
				
				graph.vertices.add(destination);
				
				final ArrayList<Edge> pathFromSource_ = pathFromSource;
				if (source != null)
					source.addEdge(destination, new EdgeFactory()
					{
						@Override
						public Edge act(Vertex end0, Vertex end1)
						{
							return new ShortcutGraphEdge(end0, end1, pathFromSource_);
						}
					});
				
				source = destination;
				pathFromSource = null;
			}
			
			if (!cycleIsCompletedHere)
				for (Edge edge: edges)
					if (edge != sourceEdge && markFilter.passes(edge))
					{
						Vertex destination = edge.getDestination(this);
						
						if (markFilter.passes(destination))
						{
							ArrayList<Edge> pathThroughEdge = pathFromSource != null ? pathFromSource : new ArrayList<Edge>();
							pathThroughEdge.add(edge);
							
							destination.buildShortCutsGraph(interestEvaluation, markFilter, source, pathThroughEdge, graph, vertexCache);
						}
					}
		}

		public void asRootOfMarkedTreePruneUninterestingBranches(int markIndex, Predicate1<Vertex> isInteresting)
		{
			inMarkedTreePruneUninterestingBranchesAndReportOwnInterestingness(markIndex, null, isInteresting);
		}
		
		private boolean inMarkedTreePruneUninterestingBranchesAndReportOwnInterestingness(int markIndex, Edge source, Predicate1<Vertex> isInteresting)
		{
			boolean thisIsInteresting = isInteresting.isTrue(this);
			
			for (Edge edge: edges)
				if (edge != source && edge.readMark(markIndex))
					thisIsInteresting = edge.getDestination(this).inMarkedTreePruneUninterestingBranchesAndReportOwnInterestingness(markIndex, edge, isInteresting) || thisIsInteresting;
			
			if (!thisIsInteresting)
			{
				this.clearMark(markIndex);
				source.clearMark(markIndex);
			}
			
			return thisIsInteresting;
		}

		public float getWeightOfMarkedTree(int markIndex)
		{
			return getWeightOfMarkedTree(markIndex, null);
		}
		
		public float getWeightOfMarkedTree(int markIndex, Edge source)
		{
			float weight = 0;
			
			for (Edge edge: edges)
				if (edge.readMark(markIndex) && edge != source)
					weight += edge.getWeight() + edge.getDestination(this).getWeightOfMarkedTree(markIndex, edge);
			
			return weight;
		}

		public void doForEachVertexOrEdgeInMarkedTree(int markIndex, ActionOnVertexOrEdge action)
		{
			doForEachVertexOrEdgeInMarkedTree(markIndex, action, null);
		}

		private void doForEachVertexOrEdgeInMarkedTree(int markIndex, ActionOnVertexOrEdge action, Edge source)
		{
			if (readMark(markIndex))
			{
				action.perform(this);
				
				for (Edge edge: edges)
					if (edge != source && edge.readMark(markIndex))
					{
						action.perform(edge);
						edge.getDestination(this).doForEachVertexOrEdgeInMarkedTree(markIndex, action, edge);
					}
			}
		}
	}
	
	public static class EuclideanVertex extends Vertex
	{
		private EuclideanVector location;
		
		public EuclideanVertex(int index, float x, float y)
		{
			this(index, new PlaneVector(x, y));
		}
		
		public EuclideanVertex(int index, float x, float y, float z)
		{
			this(index, new SpaceVector(x, y, z));
		}
		
		public EuclideanVertex(int index, EuclideanVector location)
		{
			super(index);
			this.location = location;
		}

		public void setLocation(EuclideanVector newValue)
		{
			location = newValue;
		}

		public EuclideanVector getLocation()
		{
			return location;
		}

		public float getX()
		{
			return location.getX();
		}

		public float getY()
		{
			return location.getY();
		}

		public PlaneVector getLocationAsPlaneVector()
		{
			return (PlaneVector) location;
		}
		
		public float getDistanceSquaredFrom(EuclideanVertex v)
		{
			return getDistanceSquaredFrom(v.location);
		}

		public float getDistanceSquaredFrom(EuclideanVector v)
		{
			return location.getDistanceSquaredFrom_(v);
		}
	}
	
	public static class ShortcutGraphVertex extends Vertex
	{
		public final Vertex underlyingVertex;

		public ShortcutGraphVertex(Vertex underlyingVertex, int index)
		{
			super(index);
			this.underlyingVertex = underlyingVertex;
		}

		public void settleMarkedTreeIntoUnderlyingGraph(int markIndexInShortcutTree, int markIndexInUnderlyingGraph)
		{
			settleMarkedTreeIntoUnderlyingGraph(null, markIndexInShortcutTree, markIndexInUnderlyingGraph);
		}

		private void settleMarkedTreeIntoUnderlyingGraph(ShortcutGraphEdge source, int markIndexInShortcutTree, int markIndexInUnderlyingGraph)
		{
			underlyingVertex.setMark(markIndexInUnderlyingGraph);
			
			for (Edge e: edges)
				if (e != source && e.readMark(markIndexInShortcutTree))
				{
					ShortcutGraphEdge edge = (ShortcutGraphEdge) e;
					ShortcutGraphVertex destination = (ShortcutGraphVertex) edge.getDestination(this);

					for (Edge underlyingEdge: edge.underlyingPath)
						underlyingEdge.setMarkFullyInSelfAndEnds(markIndexInUnderlyingGraph);
					
					destination.settleMarkedTreeIntoUnderlyingGraph(edge, markIndexInShortcutTree, markIndexInUnderlyingGraph);
				}
		}
	}
	
	public static class ShortcutGraphEdge extends Edge
	{
		public final ArrayList<Edge> underlyingPath;
		private float weight;
		private boolean weightIsSet;

		public ShortcutGraphEdge(Vertex end0, Vertex end1, ArrayList<Edge> underlyingPath)
		{
			super(end0, end1);
			this.underlyingPath = underlyingPath;
		}
		
		@Override
		public float getWeight()
		{
			if (!weightIsSet)
			{
				weight = calculateWeight();
				weightIsSet = true;
			}
			
			return weight;
		}

		private float calculateWeight()
		{
			float w = 0;
			
			for (Edge edge: underlyingPath)
				w += edge.getWeight();
			
			return w;
		}
		
		@Override
		public float getWeightRelatedQuantityForComparison()
		{
			return getWeight();
		}
		
		@Override
		public float getWeight(float weightRelatedQuantity)
		{
			return getWeight();
		}
	}
	
	public static abstract class Edge implements GraphEdge
	{
		public Vertex end0;
		public Vertex end1;
		private int mark;
		
		public static final Comparator<Edge> comparatorByWeight = new Comparator<Edge>()
		{
			@Override
			public int compare(Edge edge1, Edge edge2)
			{
				return edge1.compareWeights(edge2);
			}
		};
		
		public Edge(Vertex end0, Vertex end1)
		{
			this.end0 = end0;
			this.end1 = end1;
		}

		public void setMarkFullyInSelfAndEnds(int markIndex)
		{
			setMarkFully(markIndex);
			end0.setMark(markIndex);
			end1.setMark(markIndex);
		}

		public void clearMarkFullyInSelfAndEnds(int markIndex)
		{
			clearMarkFully(markIndex);
			end0.clearMark(markIndex);
			end1.clearMark(markIndex);
		}

		@Override
		public void setMarkFully(int markIndex)
		{
			mark |= 1 << markIndex;
		}
		
		@Override
		public void clearMarkFully(int markIndex)
		{
			clearMark(markIndex);
		}

		@Override
		public boolean readMark(int markIndex)
		{
			return (mark & 1 << markIndex) != 0;
		}

		public void clearMark(int markIndex)
		{
			mark &= ~(1 << markIndex);
		}
		
		public void clearAllMarks()
		{
			mark = 0;
		}
		
		@Override
		public GraphVertex getEnd0()
		{
			return end0;
		}
		
		@Override
		public GraphVertex getEnd1()
		{
			return end1;
		}
		
		@Override
		public boolean isNotNegativelyOriented()
		{
			return true;
		}
		
		public boolean endsAt(GraphVertex v)
		{
			return end0 == v || end1 == v;
		}

		public boolean destinationHasGreaterIndex(Vertex vertex)
		{
			return getDestination(vertex).getIndex() > vertex.getIndex();
		}

		public Vertex getDestination(Vertex source)
		{
			return end0 == source ? end1 : end0;
		}
		
		public abstract float getWeightRelatedQuantityForComparison();		
		public abstract float getWeight();
		public abstract float getWeight(float weightRelatedQuantity);

		private int compareWeights(Edge peer)
		{
			float w1 = this.getWeightRelatedQuantityForComparison();
			float w2 = peer.getWeightRelatedQuantityForComparison();
			return w1 < w2 ? -1 : w1 > w2 ? 1 : 0;
		}
	}
	
	public static class EuclideanEdge extends Edge
	{
		public EuclideanEdge(EuclideanVertex end0, EuclideanVertex end1)
		{
			super(end0, end1);
		}
		
		public float getLengthSquared()
		{
			return ((EuclideanVertex) end0).getDistanceSquaredFrom((EuclideanVertex) end1);
		}

		@Override
		public float getWeightRelatedQuantityForComparison()
		{
			return getLengthSquared();
		}

		@Override
		public float getWeight()
		{
			return (float) Math.sqrt(getLengthSquared());
		}

		@Override
		public float getWeight(float weightRelatedQuantity)
		{
			return (float) Math.sqrt(weightRelatedQuantity);
		}
	}

	public EdgeArraysGraph(List<? extends Vertex> vertices)
	{
		this.vertices = new ArrayList<>(vertices);
	}
	
	public EdgeArraysGraph()
	{
		this(new ArrayList<Vertex>());
	}

	public int getVertexCount()
	{
		return vertices.size();
	}
	
	@Override
	public Iterator<Vertex> iterator()
	{
		return vertices.iterator();
	}

	@Override
	public GraphVertex getVertex(int index)
	{
		return vertices.get(index);
	}
	
	public void addVertex(Vertex vertex)
	{
		vertices.add(vertex);
	}

	@Override
	public void clearAllMarks()
	{
		for (Vertex vertex: vertices)
		{
			vertex.clearAllMarks();
			
			for (Edge edge: vertex.edges)
				edge.clearAllMarks();
		}
	}

	@Override
	public void markMinimalSpanningTree_Kruskal(int markIndex)
	{
		markMinimalSpanningTree_Kruskal(markIndex, false);
	}
	
	public void markMinimalSpanningTree_Kruskal(int markIndex, boolean markEdgelessVerticesToo)
	{
		ArrayList<Edge> edges = gatherAllEdges();
		
		Collections.sort(edges, Edge.comparatorByWeight);
		
		DisjointSetRegistry sets = new DisjointSetRegistry(vertices.size());
		
		int edgeCount = 0;
		for (Edge edge: edges)
		{
			int i1 = edge.end0.index;
			int i2 = edge.end1.index;
			
			if (sets.find(i1) != sets.find(i2))
			{
				sets.union(i1, i2);

				edge.setMarkFully(markIndex);
				edge.end0.setMark(markIndex);
				edge.end1.setMark(markIndex);
				
				if (++edgeCount == vertices.size() - 1)
					break;
			}
		}
		
		if (markEdgelessVerticesToo)
			for (Vertex vertex: vertices)
				vertex.setMark(markIndex);
	}

	@Override
	public ShortestPathsTable findShortestPathsForAllVertexPairs(InterruptionSignal interruptionSignal) throws ShortestPathsSearchInterruption
	{
		ShortestPathsTable table = new ShortestPathsTable(vertices.size());
		
		for (Vertex vertex: vertices)
		{
			ShortestPathsSearch search = new DijkstraShortestPathsSearch(vertices, vertex, table, interruptionSignal);
			while (!search.isFinished())
				search.advance();
		}
		
		return table;
	}
	
	public ShortestPathsTable findShortestPathsForAllVertexPairs_lazily(InterruptionSignal interruptionSignal)
	{
		DijkstraShortestPathsSearch[] searches = new DijkstraShortestPathsSearch[vertices.size()];
		ShortestPathsTable table = new ShortestPathsTable(searches);
		
		for (int i = 0; i < searches.length; ++i)
			searches[i] = new DijkstraShortestPathsSearch(vertices, vertices.get(i), table, interruptionSignal);
		
		return table;
	}
	
	private static class DijkstraShortestPathsSearch implements ShortestPathsSearch
	{
		private final Vertex root;
		private final ShortestPathsTable table;
		private final InterruptionSignal interruptionSignal;
		private final ArrayList<Vertex> vertices;
		private PairingHeap<Double, Vertex> heap;
		private ArrayList<Handle<Double, Vertex>> entries;
		private ArrayList<Vertex> precedingVertices;
		private int i;
		
		public DijkstraShortestPathsSearch(ArrayList<Vertex> vertices, Vertex root, ShortestPathsTable table, InterruptionSignal interruptionSignal)
		{
			this.vertices = vertices;
			this.root = root;
			this.table = table;
			this.interruptionSignal = interruptionSignal;
		}

		@Override
		public void advance() throws ShortestPathsSearchInterruption
		{
			if (heap == null)
			{
				heap = new PairingHeap<Double, EdgeArraysGraph.Vertex>();
				entries = new ArrayList<>();
				precedingVertices = new ArrayList<>();
				
				for (int i = 0; i < vertices.size(); ++i) entries.add(null);
			
				for (Vertex v: vertices)
				{
					entries.set(v.index, heap.insert(v == root ? 0 : Double.POSITIVE_INFINITY, v));
					precedingVertices.add(null);
					if (interruptionSignal.isActive()) throw new ShortestPathsSearchInterruption();
				}
			}
			
			if (isFinished())
				throw new UnsupportedOperationException("no more advancement is possible");
			
			Handle<Double, Vertex> entry = heap.deleteMin();
			Vertex v = entry.getValue();
			double distance = entry.getKey();
			
			table.setCellValue(root, v, (float) distance, precedingVertices.get(v.index), i);
			
			for (Edge edge: v.edges)
			{
				if (interruptionSignal.isActive()) throw new ShortestPathsSearchInterruption();
				
				double distanceToDestination = distance + edge.getWeight();
				Vertex destination = edge.getDestination(v);
				Handle<Double, Vertex> destinationEntry = entries.get(destination.getIndex());
				
				if (distanceToDestination < destinationEntry.getKey())
				{
					destinationEntry.decreaseKey(distanceToDestination);
					precedingVertices.set(destination.getIndex(), v);
				}
			}
			
			++i;
		}

		@Override
		public boolean isFinished()
		{
			return heap != null && heap.isEmpty();
		}
	}

	private ArrayList<Edge> gatherAllEdges()
	{
		ArrayList<Edge> edges = new ArrayList<>();
		
		for (Vertex vertex: vertices)
			for (Edge edge: vertex.edges)
				if (edge.destinationHasGreaterIndex(vertex))
					edges.add(edge);
		
		return edges;
	}

	@Override
	public float[] getVerticesAsExesMadeFromLineSegments_xyxy(float exHalfWidth, Boolean requiredMark, int markIndex)
	{
		TFloatArrayList result = new TFloatArrayList();
		for (int i = 0; i < vertices.size(); ++i)
		{
			EuclideanVertex vertex = (EuclideanVertex) vertices.get(i);
			if (requiredMark == null || vertex.readMark(markIndex) == requiredMark)
			{
				float x = vertex.getX();
				float y = vertex.getY();
				result.add(x - exHalfWidth);
				result.add(y - exHalfWidth);
				result.add(x + exHalfWidth);
				result.add(y + exHalfWidth);
				result.add(x - exHalfWidth);
				result.add(y + exHalfWidth);
				result.add(x + exHalfWidth);
				result.add(y - exHalfWidth);
			}
		}
		
		return result.toArray();
	}

	@Override
	public float[] getEdgesAsLineSegments_xyxy(Boolean requiredMark, int markIndex)
	{
		TFloatArrayList result = new TFloatArrayList();
		for (Vertex vertex: vertices)
			for (Edge edge: vertex.edges)
				if (edge.destinationHasGreaterIndex(vertex) && (requiredMark == null || edge.readMark(markIndex) == requiredMark))
				{
					result.add(((EuclideanVertex) vertex).getX());
					result.add(((EuclideanVertex) vertex).getY());
					result.add(((EuclideanVertex) edge.getDestination(vertex)).getX());
					result.add(((EuclideanVertex) edge.getDestination(vertex)).getY());
				}
		
		return result.toArray();
	}

	@Override
	public void addDelaunayTriangulationEdges()
	{
		// an intermediate DCEL is created here...  that's a bit silly; it'd probably be better to just do the Delaunay algorithm within this graph (but that'd require
		// a second version of the triangulation algorithm to be written) - TODO, handle that then
		DCELGraph triangulation = new DCELGraph(vertices.size());
		for (int i = 0; i < vertices.size(); ++i)
		{
			EuclideanVertex v = (EuclideanVertex) vertices.get(i);
			triangulation.vertices[i] = new DCELGraph.Vertex(v.getIndex(), v.getLocationAsPlaneVector());
		}
		triangulation.addDelaunayTriangulationEdges();
		
		for (DCELGraph.Vertex v: triangulation.vertices)
			v.doForEachEdge(new Action1<GraphEdge>()
			{
				@Override
				public void perform(GraphEdge edge)
				{
					if (edge.isNotNegativelyOriented())
						addEdgeIfNotPresent(edge.getEnd0().getIndex(), edge.getEnd1().getIndex());
				}
			});
	}

	public void addEdgeIfNotPresent(int end0Index, int end1Index)
	{
		Vertex v0 = vertices.get(end0Index);
		Vertex v1 = vertices.get(end1Index);

		if (v0.findEdgeTo(v1) == null)
			v0.addEdge(v1);
	}
	
	private static class SteinerSubgraphVertex
	{
		public Vertex v;
		public final ArrayList<SteinerSubgraphVertex> neighbors;
		public final Boolean isSteinerVertex;
		public boolean mark;

		public SteinerSubgraphVertex(Vertex v, Boolean isSteinerVertex)
		{
			this.v = v;
			this.isSteinerVertex = isSteinerVertex;
			neighbors = new ArrayList<>();
		}
		
		@Override
		public String toString()
		{
			return v.toString();
		}

		public void addNeighbor(Vertex neighbor, Boolean neighborIsSteinerVertex)
		{
			SteinerSubgraphVertex n = new SteinerSubgraphVertex(neighbor, neighborIsSteinerVertex);
			addNeighbor(n);
		}

		public void addNeighbor(SteinerSubgraphVertex neighbor)
		{
			neighbors.add(neighbor);
			neighbor.neighbors.add(this);
		}

		@SuppressWarnings("unused")
		public void doForEachEdgeInGraph(Action2<SteinerSubgraphVertex, SteinerSubgraphVertex> action)
		{
			doForEachEdgeInGraph(action, null);
		}
		
		private void doForEachEdgeInGraph(Action2<SteinerSubgraphVertex, SteinerSubgraphVertex> action, SteinerSubgraphVertex source)
		{
			for (SteinerSubgraphVertex neighbor: neighbors)
				if (neighbor != source)
				{
					action.perform(this, neighbor);
					neighbor.doForEachEdgeInGraph(action, this);
				}
		}

		public <ThrownException extends Throwable> void doForEachEdgeInGraph(Action2Throwing<SteinerSubgraphVertex, SteinerSubgraphVertex, ThrownException> action) throws ThrownException
		{
			doForEachEdgeInGraph(action, null);
		}
		
		private <ThrownException extends Throwable> void doForEachEdgeInGraph(Action2Throwing<SteinerSubgraphVertex, SteinerSubgraphVertex, ThrownException> action, SteinerSubgraphVertex source) throws ThrownException
		{
			for (SteinerSubgraphVertex neighbor: neighbors)
				if (neighbor != source)
				{
					action.perform(this, neighbor);
					neighbor.doForEachEdgeInGraph(action, this);
				}
		}

		public SteinerNetStats gatherHiLoPrerequisiteData(SteinerSubgraphVertex neighborInNet, Vertex newTerminal, int maximalNetNeighborDepth, SubgraphSteinerVertexCache cache, ShortestPathsTable paths) throws ShortestPathsSearchInterruption
		{
			SteinerNetStats stats1 = gatherHiLoDataInBranch(neighborInNet, newTerminal, true, maximalNetNeighborDepth, cache, paths);
			SteinerNetStats stats2 = neighborInNet.gatherHiLoDataInBranch(this, newTerminal, false, maximalNetNeighborDepth, cache, paths);
			return stats1.plus(stats2);
		}

		private SteinerNetStats gatherHiLoDataInBranch(SteinerSubgraphVertex neighborInNet, Vertex newTerminal, boolean includeThisEdge, int maximalNetNeighborDepth, SubgraphSteinerVertexCache cache, ShortestPathsTable paths) throws ShortestPathsSearchInterruption
		{
			SteinerNetStats stats = new SteinerNetStats(neighborInNet, includeThisEdge ? this : null, newTerminal, paths);
			
			if (neighborInNet.isSteinerVertex)
			{
				for (SteinerSubgraphVertex neighborBeyond: neighborInNet.neighbors)
					if (neighborBeyond != this)
						stats = stats.plus(neighborInNet.gatherHiLoDataInBranch(neighborBeyond, newTerminal, true, maximalNetNeighborDepth, cache, paths));
			}
			else if (maximalNetNeighborDepth > 0)
			{
				SteinerSubgraphVertex representativeOfNearestEdge = neighborInNet.getRepresentativeOfNearestEdge(newTerminal, this, cache, paths);
				if (representativeOfNearestEdge != null)
					stats = stats.plus(neighborInNet.gatherHiLoDataInBranch(representativeOfNearestEdge, newTerminal, true, maximalNetNeighborDepth - 1, cache, paths));
			}
			
			return stats;
		}

		public SteinerSubgraphVertex getRepresentativeOfNearestEdge(Vertex referencePoint, SteinerSubgraphVertex source, SubgraphSteinerVertexCache cache, ShortestPathsTable paths) throws ShortestPathsSearchInterruption
		{
			float minimalDistance = Float.MAX_VALUE;
			SteinerSubgraphVertex representativeOfNearestEdge = null;
			
			for (SteinerSubgraphVertex neighbor: neighbors)
				if (neighbor != source)
				{
					float distanceToEdgeMoreOrLess = cache.getDistanceFromVertexToNearestVertexOnPathThatIsNotTheStartOfThePath(referencePoint, v, neighbor.v, paths);
					if (distanceToEdgeMoreOrLess < minimalDistance)
					{
						minimalDistance = distanceToEdgeMoreOrLess;
						representativeOfNearestEdge = neighbor;
					}
				}
			
			return representativeOfNearestEdge;
		}

		public boolean settleIntoMainGraph(int markIndex, ArrayList<Edge> markedEdgeRepository, ShortestPathsTable paths, Graph graph) throws ShortestPathsSearchInterruption
		{
			v.setMark(markIndex);
			return settleIntoNewGraph(null, markIndex, markedEdgeRepository, paths, graph);
		}

		private boolean settleIntoNewGraph(SteinerSubgraphVertex source, int markIndex, ArrayList<Edge> markedEdgeRepository, ShortestPathsTable paths, Graph graph) throws ShortestPathsSearchInterruption
		{
			boolean cycleEncountered = false;
			
			for (SteinerSubgraphVertex neighbor: neighbors)
				if (neighbor != source)
				{
					boolean breakingNewGround = false;
					for (Vertex lastVertex = null, vertex = v; vertex != null; lastVertex = vertex, vertex = (Vertex) paths.getNextVertexOnShortestPath(graph, vertex, neighbor.v))
						if (lastVertex != null)
						{
							Edge edge = (Edge) vertex.findEdgeTo(lastVertex);
							
							boolean edgeMark = edge.readMark(markIndex);
							if (!edgeMark)
							{
								breakingNewGround = true;
								markedEdgeRepository.add(edge);
							}
							
							if (breakingNewGround && (edgeMark || vertex.readMark(markIndex)))
								cycleEncountered = true;
							
							edge.setMarkFully(markIndex);
							vertex.setMark(markIndex);
						}
					
					if (!cycleEncountered && neighbor.settleIntoNewGraph(this, markIndex, markedEdgeRepository, paths, graph))
						cycleEncountered = true;
				}
			
			return cycleEncountered;
		}

		public boolean impliesAnyDegenerateEdges()
		{
			return neighbors.size() > (isSteinerVertex ? 3 : 1);
		}
	}
	
	private static class SteinerNetStats
	{
		public final SteinerSubgraphVertex furthestVertex;
		public final float distanceFromTerminalToFurthestVertex;
		public final float netWeight;

		public SteinerNetStats(SteinerSubgraphVertex vertexInNet, SteinerSubgraphVertex sourceVertexForReference, Vertex newTerminal, ShortestPathsTable paths) throws ShortestPathsSearchInterruption
		{
			furthestVertex = vertexInNet;
			distanceFromTerminalToFurthestVertex = paths.getPathLength(vertexInNet.v, newTerminal);
			netWeight = sourceVertexForReference != null ? paths.getPathLength(sourceVertexForReference.v, vertexInNet.v) : 0;
		}

		private SteinerNetStats(SteinerNetStats furthestVertexSource, float netWeight)
		{
			furthestVertex = furthestVertexSource.furthestVertex;
			distanceFromTerminalToFurthestVertex = furthestVertexSource.distanceFromTerminalToFurthestVertex;
			this.netWeight = netWeight;
		}

		public SteinerNetStats plus(SteinerNetStats peer)
		{
			SteinerNetStats furthestVertexSource = distanceFromTerminalToFurthestVertex > peer.distanceFromTerminalToFurthestVertex ? this : peer;
			return new SteinerNetStats(furthestVertexSource, netWeight + peer.netWeight);
		}
	}
	
	public static class SteinerSubgraphOptions
	{
		public boolean initializeSteinerNetsByDeterministicallySearchingForDecentLocations = false;
		public int maximalJumpDepthInPass1OfSingleVertexJiggling = 3;
		public int maximalJumpDepthInMultipleVertexJiggling = 1;
		public int maximalMultipleVertexJumpCountPerPatch = 1724;
		public int maximalJumpDepthInPass2OfSingleVertexJiggling = 3;
		public boolean ifMultipleVertexJigglingIsOnThenForPatchesWithoutMultipleSteinerVerticesRunAtMostOneSingleVertexJigglingPass = true;
		public Random random;
		public int maximalDepthToTraverseInNumberOfImpliedDegenerateEdgesTouchingTerminalsToPassByWhenFindingBoundsOfLocalNetworksOfSteinerVerticesForAdjustment = 1;
	}

	/*
	 * This method marks an approximate Steiner tree within this graph, connecting all terminal vertices in a way that minimizes (approximately) the weight of the connecting subgraph.
	 * Terminals must be unique; if any is repeated then an exception will be thrown.  Weight is analyzed only via Edge's weight-related methods (the abstract ones), so if a custom
	 * weighting scheme is to be used, it suffices here to express it by implementing those methods in an appropriate subclass of Edge, and then supplying an EdgeFactory that creates
	 * that kind of edge.  Though the marked Steiner tree may or may not be exactly minimal, it will always at least be a tree - unless some of the specified terminals are in
	 * unconnected subgraphs, in which case a forest of Steiner trees will be produced, one for each maximal subset of terminals in which all terminals reside in the same connected
	 * subgraph.
	 * 
	 * Here the connection of a Steiner tree is done incrementally, with one terminal being added at a time to the growing Steiner tree.
	 * 
	 * The parameter markIndex selects the particular mark that will be left in the shape of the tree; this mark will also be used for scratch work during the algorithm.  It is
	 * expected that at the start of this call, no existing instances of this mark will exist in parts of the graph visited by the algorithm.  The natural way to ensure this, of
	 * course, is to select some mark that isn't being used anywhere at all in the graph.
	 * 
	 * The parameter interruptionSignal provides the caller with a way to stop the algorithm early and discard its potential result.  If the function is interrupted then it may leave
	 * the graph in a partially marked state.  The signal object is not treated specially with respect to synchronization; its methods are called naively, and if any synchronization
	 * is desired then the object must manage it internally.  This argument may not be null; if interruption is not desired, then something like NullInterruptionSignal.instance must
	 * be passed in.
	 * 
	 * The options parameter controls aspects of the algorithm's method for placing Steiner vertices:
	 * 
	 * An initial search may be toggled; this search is intended to improve the quality of the result, but it is not guaranteed to do so.  It is guaranteed to take some time, however.
	 * But it may also conceivably reduce the time taken by the two jiggling steps that follow, if they are turned on.  Still, it's wise to start with this option turned off.  Sometimes
	 * it eats quite a large amount of time, and in such cases it's unlikely to improve results enough to justify that.
	 * 
	 * Single-vertex jiggling is a greedy graph-refinement strategy that seems intuitively to be (probably) cheaper in general than either the other type of jiggling or the initial
	 * vertex-positioning step.  Setting its maximal jump depth to 0 will deactivate it entirely.  A higher jump depth may increase the quality of the result, but will also increase
	 * running time.
	 * 
	 * Single-vertex jiggling occurs in two passes.  These are identical with one another, except in that pass 1 happens after any initial positioning and before multiple-vertex
	 * jiggling, whereas pass 2 happens after all of that.  It may be good enough to activate only one of these passes, or neither of them, but there could conceivably be some advantage
	 * in running both.  (The question of what's best is complex and has been left unanalyzed.)  If multiple-vertex jiggling is deactivated and pass 2 has a lower depth than pass 1
	 * does, then it's probably true that pass 2 will be guaranteed to accomplish nothing; but in other scenarios it's not entirely clear that running both passes would waste work.
	 * 
	 * In any case, at most one single-jiggling pass will be run for chunks of the graph that locally skip the multiple-vertex jiggling pass - unless the boolean switch calling for
	 * that behavior is turned off.  (If multiple-vertex jiggling is deactivated entirely, then this switch has no effect; in that case exactly the specified single-vertex jiggling
	 * passes will run for every patch of the graph.)
	 * 
	 * Multiple-vertex jiggling is a strategy similar to single-vertex jiggling that moves multiple vertices simultaneously instead of just one.  Setting its maximal depth count to 0
	 * will deactivate it.  This type of jiggling will only happen in some graphs; in a few graphs it may be especially expensive (in running time).  The option setting a maximal jump
	 * count per patch is meant to prevent this from becoming a major problem.  That value should be set higher than the jump depth caps - say, in the hundreds.  (But it may be
	 * explicitly deactivated; setting it to Integer.MAX_VALUE will accomplish that.)  When this cap is in danger of being violated, only a random subset of all the possible
	 * multiple-vertex jumps will be considered; the randomness is controlled by the "random" parameter in the options, which therefore must be supplied when multiple-vertex jiggling
	 * is turned on and the cap is set.  All other aspects of this algorithm are deterministic, and if the maximal jump depth for multiple-vertex jiggling is set to 0, or if the
	 * maximal jump count per patch is set sufficiently high, then the entire algorithm will proceed without randomness.)
	 * 
	 * Finally, there is an option determining how large, in some cases, the subgraphs to be jiggled and positioned in the ways described above will be.  Generally, it makes sense,
	 * when incrementally splicing a new terminal into the growing Steiner tree, to find and readjust the section of the tree in the vicinity of the splicing point that ought to be
	 * "flexible", that is, the section comprising only Steiner vertices and the edges connecting them, as well as of course the edges connecting these to the nearest terminals.  And
	 * if there are further edges on the other side of some terminal vertex, then those should perhaps not be considered part of the network that's about to be readjusted, because
	 * terminals can't move and so the edges beyond them are seemingly insulated from any movements of the Steiner vertices on the side from which the flexing originates.
	 * 
	 * The trouble with this way of thinking (please stop reading here and treat this parameter as another "bigger value means better results and slower searches" parameter if you're
	 * not super interested in this stuff) is that in some cases, a movement on "this" side will bring one of the local edges so close to an edge on the other side of a terminal
	 * that it'd suddenly be appropriate to add a new Steiner vertex between the two edges, turning their conceptual V shape into a Y.  This addition would bring the entire little
	 * network of Steiner vertices beyond that terminal into the network that was already being readujusted.  In a Steiner tree algorithm that works with the entire Euclidean plane
	 * (and that isn't restricted to a graph), the condition triggering this addition can easily be detected as soon as it occurs:  when two edges form an angle of less than 120
	 * degrees, a new Steiner vertex should be added between them, and the two local networks thus joined should be readjusted.  But in an arbitrary graph, "angles" between edges
	 * aren't so easy to measure (and the concept of them may not even be defined in the first place).  In this particular algorithm, the solution taken is different:  When a network
	 * to be readjusted is being identified, it is assumed that any terminal encountered should be treated as a hard boundary if it has only one edge connected to it - but otherwise
	 * the terminal should be considered as a potential bridge to an additional network section beyond it, a section including exactly one of the other edges extending from it.  It
	 * would be potentially correct to include all the networks touching the terminal; if this were done recursively then the entire graph would end up being jiggled with each new
	 * terminal addition.  But that might take too long.  Instead, this algorithm chooses only one edge per terminal, to keep things simple, and during the ensuing jiggling, it
	 * pretends that an additional Steiner vertex has been added to bridge the local network with that edge's network, placing this new vertex initially right atop the terminal.  Now
	 * this too could end up including a large portion of the graph in the network to be readjusted after each terminal addition.  To prevent that, a maximal depth option is provided
	 * here.  Setting this depth to 0 prevents this type of "bridging" between networks from happening at all.  Setting it to 1 allows one level of neighboring networks to be included
	 * in the current one.  A setting of 2 will include the networks neighboring those as well...  and so on.  Regardless, only one edge per terminal, the one that seems "nearest" to
	 * the newly added terminal, will be bridged to in this way, so setting this option to a very large value will not necessarily cause the whole graph to be jiggled with each
	 * terminal addition.
	 * 
	 * Here is a specific example of why this option matters:  A six-vertex, eight-edge graph that looks like
	 * ______
	 * |>--<|
	 * 
	 * with the four "corner" vertices designated as terminals.  Suppose that the terminals happen to be added to the growing Steiner tree in clockwise order, starting at the
	 * lower-left corner.  Then the Steiner tree will start as just the leftmost vertical edge; and with the addition of the third terminal, that edge will be spliced to include a
	 * Steiner vertex, which will probably then be jiggled (or otherwise placed) directly upon the upper-left terminal vertex, since that makes for the shortest possible graph.  The
	 * degenerate edge between that terminal and the Steiner vertex will be removed in a normal cleanup step, and the Steiner tree will then look like this:
	 * 
	 * ______
	 * |
	 * 
	 * Then the fourth vertex will be added.  The long upper edge in the tree will be chosen for splicing; a Steiner vertex will be added amid it.  Any jiggling or other adjustment
	 * happening after that will perhaps be constrained from producing the full, correct Steiner tree, which (assuming the edge weights are such that the > and < here seem to possess
	 * approximately 120-degree angles) looks like
	 * 
	 * >--<
	 * 
	 * , by the fact that of the two Steiner vertices in that correct tree, only one now exists in the tree being adjusted, which therefore can come no closer to that ideal shape
	 * than something like this:
	 * 
	 * |\--<
	 * 
	 * This could have been prevented via the retention of the degenerate edge between the first Steiner vertex and the upper-left terminal that it settled atop.  The trouble with
	 * keeping degenerate edges like that is that sometimes they'll just add useless complexity to future iterations of the algorithm, and that, worse, ones that are needed won't
	 * always be the ones that actually arise.  When THREE edges happen to meet at a terminal, a Steiner vertex and a degenerate edge connecting it to the terminal could be added in
	 * order to separate two of those three from the terminal, replacing them with the one new degenerate edge; then a second Steiner vertex could be added, to replace that degenerate
	 * edge and the third of the three original edges with another new degenerate edge.  But there are multiple ways to do this, differing in the choice of which two edges to join at
	 * the start, and they produce different graph topologies that will all "look" the same as long as the degenerate edges remain degenerate, but that will stretch and flex
	 * differently when the graph is readjusted upon addition of a new terminal.  If only the degenerate edges that happen to actually arise during adjustment are kept, then the other
	 * possible topologies of degenerate networks will be ignored, and their potentials lost.
	 * 
	 * So in this algorithm, all degenerate edges are immediately collapsed, and each vertex with "extra" edges connected to it (a terminal with two or more neighbors, or a Steiner
	 * vertex with four or more neighbors) is considered to be a bank of potential new degenerate edges that will be (imperfectly) called into being when they seem to be needed.
	 * 
	 * In the example above, a maximal depth setting of 1 would be enough to join the two jigglable networks represented by the top and left edges, respectively, of the L-shaped
	 * three-terminal tree-in-progress, recreating the lost Steiner vertex that is needed to form the
	 * 
	 * >-
	 * 
	 * configuration on the left side of the tree.
	 * 
	 * It is left up to the user to determine how important this sort of logic actually is, and whether this clumsy maximal-depth mechanism is even remotely close to being good enough
	 * to generally solve the problem it's intended to solve.  (It certainly does not PERFECTLY solve the problem, at least, since in theory there is an input configuration for which
	 * the addition of just one new terminal during one of the algorithm's iterations could require EVERY other terminal to need a new Steiner vertex added atop it, and moreover, each
	 * such Steiner vertex addition could expend one entire level of the sort of "depth" that this option refers to - which means that in the worst case, it will not be possible to
	 * gain any speed via this option without sacrificing correctness.  No attempt will be made here to draw in ASCII a graph that produces this worst-case behavior; please believe
	 * that according to this author's intuition, which seems fairly clear on the matter, it does exist.)
	 * 
	 * 
	 * 
	 * This algorithm has been optimized substiantially, though it likely isn't nearly as fast as it could be.  It spends much of its time in a standard all-pairs-shortest-paths
	 * calculation that could easily have been parallelized; instead, this calculation has been made as "lazy" as possible:  Early in the algorithm, for each vertex in the graph, a
	 * search for the shortest paths from that vertex to each other vertex is started and then immediately paused.  Later, whenever the algorithm wants to know something about the
	 * shortest path between two vertices, it resumes the appropriate one of those searches, runs it just long enough to produce the requested information, and then pauses it again.
	 * This can save a very large amount of work!  But it is not trivially compatible with parallelization.
	 * 
	 * Another strategy that could be used is caching of the results of some searches for vertices meeting particular criteria.  Three such searches that seemed to be well suited to
	 * this have been encapsulated in the SubgraphSteinerVertexCache class; the intention was to store their results in hash maps or other indices.  But more recently, casual testing
	 * and some light thinking on the matter seemed to suggest that caching would not actually help very much here.  The queries that are repeated very frequently seem to be the ones
	 * that are the easiest to redo from scratch, for the most part.  There may be savings to find here, but they probably aren't huge, if they exist at all.
	 * 
	 * Parallelization is promising here.  The overall algorithm is incremental of course, but various pieces of it could easily be run in parallel.  The all-pairs-shortest-paths
	 * section, again, could be parallelized, though some clevernees might be needed to decide which parts of it to run first.  Patch constructions could probably be parallelized
	 * pretty straightforwardly.  Maybe the same goes for the identification and evaluation of potential edges for splicing.  Evaluation of different choices during Steiner vertex
	 * jiggling could be done in parallel.  There is lots of opportunity here.
	 * 
	 * Some of the above, of course, involves the manipulation of shared data, such as markings in the graph and in various metagraphs.  But it should be simple enough to make these
	 * thread-local.
	 * 
	 * There may be great potential for optimization via a reworking of the basic incremental strategy used here.  But that is outside the scope of this comment.
	 * 
	 * The text above concerns potential parallelization of this algorithm on a multi-threaded CPU.  A GPU might be the natural tool for use in solving this problem after all, though.
	 * But using a GPU would surely require a reworking of every part of the algorithm, after which it would be an entirely different thing.  So that too is outside the scope of this
	 * comment.  Probably somebody has already done good work on that topic, anyhow.
	 */
	public void markSteinerSubgraph(TIntArrayList terminals, int markIndex, InterruptionSignal interruptionSignal, SteinerSubgraphOptions options)
	{
		try
		{
		
		if (options.maximalJumpDepthInMultipleVertexJiggling > 0 && options.maximalMultipleVertexJumpCountPerPatch < Integer.MAX_VALUE && options.random == null)
			throw new IllegalArgumentException("options.random should not be null if multiple-vertex jiggling is to take place and yet be limited; either set the jump count limit to Integer.MAX_VALUE, " +
					"set the multi-vertex jump depth to 0, or set options.random");
		
		if (terminals.isEmpty())
			return;
		
		for (int i = 0; i < terminals.size(); ++i)
		{
			int terminal = terminals.get(i);
			if (terminal >= vertices.size() || vertices.get(terminal).getIndex() != terminal)
				throw new RuntimeException("indices are to be valid and aligned");
		}
		
		ShortestPathsTable paths = findShortestPathsForAllVertexPairs_lazily(interruptionSignal);
		if (interruptionSignal.isActive()) return; // this one is important; paths could be null or dangerously invalid if interruption has occurred!
		TIntArrayList sortedShortestPaths = paths.sortShortestPaths_abab(terminals);
		if (interruptionSignal.isActive()) return;
		
		{
			int[] indexInGraphToIndexInTerminalList = new int[vertices.size()];
			ArrayList<SteinerSubgraphVertex> nodes = new ArrayList<>();

			for (int i = 0; i < terminals.size(); ++i)
			{
				int terminal = terminals.getQuick(i);
				indexInGraphToIndexInTerminalList[terminal] = i;
				nodes.add(new SteinerSubgraphVertex(vertices.get(terminal), null));
			}
			
			DisjointSetRegistry sets = new DisjointSetRegistry(terminals.size());
			
			for (int i = 0, edgeCount = 0; i < sortedShortestPaths.size() && edgeCount < terminals.size() - 1; i += 2)
			{
				if (interruptionSignal.isActive()) return;
				int a = indexInGraphToIndexInTerminalList[sortedShortestPaths.getQuick(i)];
				int b = indexInGraphToIndexInTerminalList[sortedShortestPaths.getQuick(i + 1)];

				if (sets.find(a) != sets.find(b))
				{
					sets.union(a, b);

					nodes.get(a).addNeighbor(nodes.get(b));
					++edgeCount;
				}
			}

			for (int i = 0; i < terminals.size(); ++i)
			{
				SteinerSubgraphVertex terminal = nodes.get(i);

				if (!terminal.mark)
				{
					TIntArrayList sortedOriginalTerminals = spreadMarksThroughTreeFromOnePointPuttingShortestEdgesFirst(terminal, paths);
					markSteinerSubgraph_(sortedOriginalTerminals, markIndex, markIndex, options, paths, interruptionSignal);
				}
				if (interruptionSignal.isActive()) return;
			}
		}
		
		}
		catch (ShortestPathsSearchInterruption interruption)
		{
			return;
		}
	}
	
	private static class SteinerSubgraphEdge implements Comparable<SteinerSubgraphEdge>
	{
		public SteinerSubgraphVertex source;
		public SteinerSubgraphVertex destination;
		private final float length;

		public SteinerSubgraphEdge(SteinerSubgraphVertex source, SteinerSubgraphVertex destination, ShortestPathsTable paths) throws ShortestPathsSearchInterruption
		{
			this.source = source;
			this.destination = destination;
			length = getPathLength(source, destination, paths);
		}

		private SteinerSubgraphEdge(SteinerSubgraphVertex source, SteinerSubgraphVertex destination, float length)
		{
			this.source = source;
			this.destination = destination;
			this.length = length;
		}

		@Override
		public int compareTo(SteinerSubgraphEdge peer)
		{
			return length < peer.length ? -1 : length > peer.length ? 1 : 0;
		}

		public SteinerSubgraphEdge reversed()
		{
			return new SteinerSubgraphEdge(destination, source, length);
		}
	}

	private TIntArrayList spreadMarksThroughTreeFromOnePointPuttingShortestEdgesFirst(SteinerSubgraphVertex someVertexInTree, ShortestPathsTable paths) throws ShortestPathsSearchInterruption
	{
		TIntArrayList result = new TIntArrayList();
		
		if (someVertexInTree.neighbors.isEmpty())
			result.add(someVertexInTree.v.getIndex());
		else
		{
			SteinerSubgraphEdge shortestEdge = findShortestEdgeInTree(someVertexInTree, paths);

			PriorityQueue<SteinerSubgraphEdge> fringe = new PriorityQueue<>();
			fringe.offer(shortestEdge);
			fringe.offer(shortestEdge.reversed());
			
			while (!fringe.isEmpty())
			{
				SteinerSubgraphEdge edge = fringe.poll();
				edge.destination.mark = true;
				result.add(edge.destination.v.getIndex());
				
				for (SteinerSubgraphVertex neighbor: edge.destination.neighbors)
					if (neighbor != edge.source)
						fringe.add(new SteinerSubgraphEdge(edge.destination, neighbor, paths));
			}
		}
		
		return result;
	}

	private SteinerSubgraphEdge findShortestEdgeInTree(SteinerSubgraphVertex someVertexInTree, ShortestPathsTable paths) throws ShortestPathsSearchInterruption
	{
		return findShortestEdgeInTreeBranch(someVertexInTree, null, paths);
	}
	
	private SteinerSubgraphEdge findShortestEdgeInTreeBranch(SteinerSubgraphVertex someVertexInTree, SteinerSubgraphVertex source, ShortestPathsTable paths) throws ShortestPathsSearchInterruption
	{
		SteinerSubgraphEdge shortestEdge = source != null ? new SteinerSubgraphEdge(source, someVertexInTree, paths) : null;
		
		for (SteinerSubgraphVertex neighbor: someVertexInTree.neighbors)
			if (neighbor != source)
			{
				SteinerSubgraphEdge edge = findShortestEdgeInTreeBranch(neighbor, someVertexInTree, paths);
				if (shortestEdge == null || edge != null && edge.compareTo(shortestEdge) < 0)
					shortestEdge = edge;
			}
		
		return shortestEdge;
	}

	private static float getPathLength(SteinerSubgraphVertex source, SteinerSubgraphVertex destination, ShortestPathsTable paths) throws ShortestPathsSearchInterruption
	{
		return paths.getPathLength(source.v, destination.v);
	}
	
	private static class PotentialEdgeForSplicing
	{
		private final SteinerSubgraphVertex end1;
		private final SteinerSubgraphVertex end2;
		public final float hiBoundOnWeightChange;
		public final float loBoundOnWeightChange;
		public final float netWeight;

		public PotentialEdgeForSplicing(Vertex newTerminal, SteinerSubgraphVertex end1, SteinerSubgraphVertex end2, SteinerNetStats stats, ShortestPathsTable paths) throws ShortestPathsSearchInterruption
		{
			this.end1 = end1;
			this.end2 = end2;
			
			loBoundOnWeightChange = -stats.netWeight + stats.distanceFromTerminalToFurthestVertex;
			
			float distance_n1 = paths.getPathLength(newTerminal, end1.v);
			float distance_n2 = paths.getPathLength(newTerminal, end2.v);
			float distance_12 = paths.getPathLength(end1.v, end2.v);
			
			if (distance_12 > distance_n1)
				if (distance_12 > distance_n2)
					hiBoundOnWeightChange = distance_n1 + distance_n2 - distance_12;
				else
					hiBoundOnWeightChange = distance_n1;
			else
				if (distance_n2 < distance_n1)
					hiBoundOnWeightChange = distance_n2;
				else
					hiBoundOnWeightChange = distance_n1;
			
			netWeight = stats.netWeight;
		}
	}

	private static class SteinerSubgraphPatch implements Comparable<SteinerSubgraphPatch>
	{
		private static class PatchVertex implements GenericShortestPathsGraphVertex
		{
			public int index;
			public ArrayList<PatchVertex> neighbors;
			public SteinerSubgraphVertex v;
			public Vertex location;
			public Vertex guessedLocation;
			public boolean mark;
			
			public PatchVertex(SteinerSubgraphVertex v)
			{
				this(v, v.v);
			}

			public PatchVertex(SteinerSubgraphVertex v, Vertex vv)
			{
				this.v = v;
				location = vv;
				neighbors = new ArrayList<>();
			}
			
			@Override
			public GraphVertex asGraphVertex()
			{
				return location;
			}
			
			@Override
			public Collection<? extends GenericShortestPathsGraphVertex> getNeighbors()
			{
				return neighbors;
			}

			public boolean isTerminal()
			{
				return v != null && !v.isSteinerVertex;
			}

			public Vertex getBestGuessAtLocation()
			{
				return guessedLocation != null ? guessedLocation : location;
			}

			public PatchVertex markAndReturn()
			{
				mark = true;
				return this;
			}

			public void actOnVerticesWithinRange(Action1<? super Vertex> action, int rangeInEdgesInclusive, int workingMarkIndex)
			{
				if (rangeInEdgesInclusive == 1)
				{
					action.perform(location);
					for (Edge edge: location.edges)
						action.perform(edge.getDestination(location));
				}
				else if (rangeInEdgesInclusive > 1)
				{
					actOnVerticesWithinRange_withMarking(action, location, rangeInEdgesInclusive, workingMarkIndex);
					location.clearConnectedSubgraphMark(workingMarkIndex);
				}
			}

			private static void actOnVerticesWithinRange_withMarking(Action1<? super Vertex> action, Vertex vertex, int rangeInEdgesInclusive, int workingMarkIndex)
			{
				if (vertex.readMark(workingMarkIndex))
					return;
				
				vertex.setMark(workingMarkIndex);
				action.perform(vertex);
				
				if (rangeInEdgesInclusive > 0)
					for (Edge edge: vertex.edges)
						actOnVerticesWithinRange_withMarking(action, edge.getDestination(vertex), rangeInEdgesInclusive - 1, workingMarkIndex);
			}
		}
		
		private ArrayList<PatchVertex> vertices;
		private TIntArrayList edges_abab;
		private final PatchVertex terminal;
		private final PatchVertex splice;
		private final EdgeToBeReplaced edgeBeforeSplicing;
		private final ArrayList<RecreatedSteinerVertex> recreatedSteinerVertices;
		private final float approximateTotalWeight;
		private final float weightChange;

		public SteinerSubgraphPatch(PotentialEdgeForSplicing edgeForSplicing, SteinerSubgraphVertex newTerminal, int workingMarkIndex, SteinerSubgraphOptions options, final ShortestPathsTable paths, Graph graph, SubgraphSteinerVertexCache cache) throws ShortestPathsSearchInterruption
		{
			boolean singleJiggling1GenerallyOn = options.maximalJumpDepthInPass1OfSingleVertexJiggling > 0;
			boolean multiJigglingGenerallyOn = options.maximalJumpDepthInMultipleVertexJiggling > 0 && options.maximalMultipleVertexJumpCountPerPatch > 0;
			boolean singleJiggling2GenerallyOn = options.maximalJumpDepthInPass2OfSingleVertexJiggling > 0;
			int maximalTerminalCrossingDepth = options.maximalDepthToTraverseInNumberOfImpliedDegenerateEdgesTouchingTerminalsToPassByWhenFindingBoundsOfLocalNetworksOfSteinerVerticesForAdjustment;

			vertices = new ArrayList<>();
			edges_abab = new TIntArrayList();
			recreatedSteinerVertices = new ArrayList<>();
			
			float oldWeight = edgeForSplicing.netWeight;
			boolean initializingWholeNet = options.initializeSteinerNetsByDeterministicallySearchingForDecentLocations;
			
			edgeBeforeSplicing = new EdgeToBeReplaced(edgeForSplicing.end1, edgeForSplicing.end2);
			
			PatchVertex spliceNeighborA = addVerticesAndReturnFirst(edgeForSplicing.end1, edgeForSplicing.end2, newTerminal.v, maximalTerminalCrossingDepth, cache, paths);
			PatchVertex spliceNeighborB = addVerticesAndReturnFirst(edgeForSplicing.end2, edgeForSplicing.end1, newTerminal.v, maximalTerminalCrossingDepth, cache, paths);
			terminal = addVertex(newTerminal);
			splice = addVertex(new PatchVertex(null, initializingWholeNet ? null : cache.getApproximateMidpoint(edgeForSplicing.end1.v, edgeForSplicing.end2.v, 2, paths, graph)));

			link(spliceNeighborA, splice);
			link(spliceNeighborB, splice);
			link(terminal, splice);
			
			if (initializingWholeNet)
			{
				Deque<PatchVertex> initialPositioningQueue = new ArrayDeque<>();
				for (PatchVertex v: vertices)
					if (v.isTerminal())
					{
						v.mark = true;
						initialPositioningQueue.offerLast(v);
					}
				
				PatchVertex innermostSteinerVertex = null;
				while (!initialPositioningQueue.isEmpty())
					for (PatchVertex n: initialPositioningQueue.pollFirst().neighbors)
						if (!n.mark)
						{
							n.mark = true;
							innermostSteinerVertex = n;
							initialPositioningQueue.offerLast(n);
						}
				
				for (PatchVertex neighbor: innermostSteinerVertex.neighbors)
					setAndReturnGuessedLocation(innermostSteinerVertex, neighbor, paths, graph, cache);
				
				deterministicallyPlaceVerticesConsumingGuessedLocations(null, innermostSteinerVertex, workingMarkIndex, paths, cache);
			}
			
			approximateTotalWeight = oldWeight + paths.getPathLength(terminal.location, splice.location);
			
			if (singleJiggling1GenerallyOn || multiJigglingGenerallyOn || singleJiggling2GenerallyOn)
			{
				ArrayList<PatchVertex> steinerVertices = new ArrayList<>();
				for (PatchVertex vertex: vertices)
					if (!vertex.isTerminal())
						steinerVertices.add(vertex);
				
				if (singleJiggling1GenerallyOn)
					performSingleVertexJiggling(steinerVertices, options.maximalJumpDepthInPass1OfSingleVertexJiggling, workingMarkIndex, paths);
				
				boolean multiJigglingJudgedLocallyUnnecessary;
				if (multiJigglingGenerallyOn)
					if (steinerVertices.size() >= 2)
					{
						multiJigglingJudgedLocallyUnnecessary = false;
						performMultipleVertexJiggling(steinerVertices, options.maximalJumpDepthInMultipleVertexJiggling, options.maximalMultipleVertexJumpCountPerPatch, options.random, workingMarkIndex, paths);
					}
					else
						multiJigglingJudgedLocallyUnnecessary = true;
				else
					multiJigglingJudgedLocallyUnnecessary = false;
				
				boolean noMoreSingleJigglingAllowed = multiJigglingJudgedLocallyUnnecessary && options.ifMultipleVertexJigglingIsOnThenForPatchesWithoutMultipleSteinerVerticesRunAtMostOneSingleVertexJigglingPass;
				if (singleJiggling2GenerallyOn && !noMoreSingleJigglingAllowed)
					performSingleVertexJiggling(steinerVertices, options.maximalJumpDepthInPass2OfSingleVertexJiggling, workingMarkIndex, paths);
			}
			
			float newWeight = SteinerSubgraphShortestPathTrees.overestimateWeightOfTree(vertices.get(0), paths);
			weightChange = newWeight - oldWeight;
		}

		private void performSingleVertexJiggling(final ArrayList<PatchVertex> steinerVertices, int maximalJumpDepth, int workingMarkIndex, final ShortestPathsTable paths) throws ShortestPathsSearchInterruption
		{
			Stack<PatchVertex> stack = new Stack<>();

			for (PatchVertex vertex: steinerVertices)
				stack.push(vertex.markAndReturn());

			float distanceTraveled = 0;
			while (!stack.isEmpty() && distanceTraveled < getMaximalDistanceTraveled(steinerVertices, approximateTotalWeight))
			{
				final PatchVertex vertex = stack.pop();
				vertex.mark = false;
				
				BestVertex bestVertex = SteinerSubgraphShortestPathTrees.bestVertex(vertex.neighbors, paths);

				class BestVertexUpdate implements Action1<Vertex>
				{
					public ShortestPathsSearchInterruption interruption;

					@Override
					public void perform(Vertex v)
					{
						if (interruption == null)
							try
							{
								bestVertex.update(v);
							}
							catch (ShortestPathsSearchInterruption interruption)
							{
								this.interruption = interruption;
							}
					}
				}
				BestVertexUpdate bestVertexUpdate = new BestVertexUpdate();
				vertex.actOnVerticesWithinRange(bestVertexUpdate, maximalJumpDepth, workingMarkIndex);
				if (bestVertexUpdate.interruption != null)
					throw bestVertexUpdate.interruption;
				
				if (bestVertex.value != vertex.location)
				{
					Vertex newLocation = bestVertex.getValue(vertex.location);
					distanceTraveled += paths.getPathLength(vertex.location, newLocation);
					vertex.location = newLocation;

					stack.push(vertex.markAndReturn());
					for (PatchVertex neighbor: vertex.neighbors)
						if (!neighbor.isTerminal() && !neighbor.mark)
							stack.push(neighbor.markAndReturn());
				}
			}
		}
		
		private void performMultipleVertexJiggling(final ArrayList<PatchVertex> steinerVertices, int maximalJumpDepth, int maximalMultipleVertexJumpCount, Random random, int workingMarkIndex, final ShortestPathsTable paths) throws ShortestPathsSearchInterruption
		{
			final Vertex[] hypotheticalLocations = new Vertex[vertices.size()];
			for (int i = 0; i < vertices.size(); ++i)
				hypotheticalLocations[i] = vertices.get(i).location;
			
			@SuppressWarnings("serial")
			class AddingList<Element> extends ArrayList<Element> implements Action1<Element>
			{
				@Override
				public void perform(Element participant)
				{
					add(participant);
				}
			};
			
			final ArrayList<AddingList<Vertex>> neighborsWithinRange = new ArrayList<>();
			for (int i = 0; i < steinerVertices.size(); ++i)
				neighborsWithinRange.add(new AddingList<Vertex>());

			boolean aFurtherJumpMayBeHelpful;
			float distanceTraveled = 0;
			float initialTotalWeight = -1;
			do
			{
				for (int i = 0; i < steinerVertices.size(); ++i)
				{
					AddingList<Vertex> neighbors = neighborsWithinRange.get(i);
					neighbors.clear();
					steinerVertices.get(i).actOnVerticesWithinRange(neighbors, maximalJumpDepth, workingMarkIndex);
				}
				
				long combinationCountL = neighborsWithinRange.get(0).size();
				double combinationCountD = combinationCountL;
				for (int i = 1; i < neighborsWithinRange.size(); ++i)
				{
					combinationCountL *= neighborsWithinRange.get(i).size();
					combinationCountD *= neighborsWithinRange.get(i).size();
				}
				{
					double discrepancy = Math.abs(combinationCountD - combinationCountL);
					if (discrepancy >= .15 * combinationCountD)
						combinationCountL = Long.MAX_VALUE;
				}
				
				class BestMultipleVertexJump
				{
					public Vertex[] bestLocations = new Vertex[hypotheticalLocations.length];
					public float lowestWeight = Float.MAX_VALUE;
					private boolean bestPositionIsNonZero;
					private float zeroWeight;
					private boolean zeroWeightIsSet;
					
					public void update(int[] neighborIndices) throws ShortestPathsSearchInterruption
					{
						boolean nonZeroPosition = false;
						for (int i = 0; i < neighborIndices.length; ++i)
						{
							hypotheticalLocations[steinerVertices.get(i).index] = neighborsWithinRange.get(i).get(neighborIndices[i]);
							nonZeroPosition = nonZeroPosition || neighborIndices[i] != 0;
						}
						
						float weight = 0;
						
						for (int i = 0; i < edges_abab.size(); i += 2)
						{
							int a = edges_abab.getQuick(i);
							int b = edges_abab.getQuick(i + 1);
							weight += paths.getPathLength(hypotheticalLocations[a], hypotheticalLocations[b]);
						}
						
						if (!nonZeroPosition)
						{
							zeroWeight = weight;
							zeroWeightIsSet = true;
						}
						
						if (weight < lowestWeight)
						{
							System.arraycopy(hypotheticalLocations, 0, bestLocations, 0, hypotheticalLocations.length);
							lowestWeight = weight;
							bestPositionIsNonZero = nonZeroPosition;
						}
					}

					public boolean isNonZero()
					{
						return bestPositionIsNonZero;
					}

					public float getZeroWeight()
					{
						if (!zeroWeightIsSet)
							throw new UnsupportedOperationException();
						
						return zeroWeight;
					}
				}
				
				BestMultipleVertexJump bestJump = new BestMultipleVertexJump();
				int[] neighborIndices = new int[neighborsWithinRange.size()];
				
				if (combinationCountL > maximalMultipleVertexJumpCount)
				{
					bestJump.update(neighborIndices);
					
					int jumpCount = (int) Math.min(combinationCountL / 2, maximalMultipleVertexJumpCount);
					
					for (int ji = 0; ji < jumpCount; ++ji)
					{
						int movementCount;
						int sillyLimit = 5;
						do
						{
							movementCount = 0;
							for (int i = 0; i < neighborsWithinRange.size(); ++i)
							{
								neighborIndices[i] = random.nextInt(neighborsWithinRange.get(i).size());
								if (neighborIndices[i] > 0) ++movementCount;
							}
						}
						while (movementCount == 1 && --sillyLimit > 0);
						/*
						 * ^^ this is an imperfect way to do this...  it would be better to find a single random number in [0, total combination count) and then use the modulus or whatever to extract a set of neighbor indices from that, and to
						 * employ, in the selection of the initial number, a random number algorithm that would never repeat itself.  such algorithms exist, although some may have significant flaws.  (of course, it'd be possible to hack in this
						 * behavior by simply tracking in a hash set or a long bitfield all numbers tried thus far; there may be hash set classes or sparse bitfield classes that could actually make this efficient enough.  more investigation is
						 * needed.  finally, it is of course possible to generate all possible choices in order and then randomly shuffle the list of them and take a prefix of the result, but that may have an unacceptable startup cost, though
						 * such could be mitigated with clever caching, and the shuffling might be too expensive as well, although a partial shuffling should be easy enough to do...)
						 * 
						 * also, of course, there is the movementCount issue:  the aim here is to avoid doing any single-vertex jiggling, since that has (somewhat arbitrarily) been relegated to the passes designated for such work.  it's not clear
						 * that that's a good aim, but it is the aim.  (an alternative would be to include all multi-vertex movements regardless of the movement count; another would be to do that if single-vertex jiggling was not to be done, and
						 * to otherwise skip the single movements as is being done here.)  anyway, if this is to be done, then very likely there is a more efficient way to do it, and there is certainly a more correct way.
						 */
						
						bestJump.update(neighborIndices);
					}
				}
				else
				{
					boolean combinationIsValid;
					do
					{
						bestJump.update(neighborIndices);
						combinationIsValid = incrementCombination(neighborIndices, neighborsWithinRange, 0);
					}
					while (combinationIsValid);
				}
				
				if (initialTotalWeight < 0)
					initialTotalWeight = bestJump.getZeroWeight();

				aFurtherJumpMayBeHelpful = false;
				
				if (bestJump.isNonZero())
				{
					for (int i = 0; i < bestJump.bestLocations.length; ++i)
					{
						PatchVertex vertex = vertices.get(i);
						if (!vertex.isTerminal())
						{
							Vertex newLocation = bestJump.bestLocations[i];
							distanceTraveled += paths.getPathLength(vertex.location, newLocation);
							vertex.location = newLocation;
						}
					}
					
					aFurtherJumpMayBeHelpful = true;
				}
			}
			while (aFurtherJumpMayBeHelpful && distanceTraveled < getMaximalDistanceTraveled(steinerVertices, initialTotalWeight));
		}

		private float getMaximalDistanceTraveled(final ArrayList<PatchVertex> steinerVertices, float initialTotalWeight)
		{
			return steinerVertices.size() * 2.5f * initialTotalWeight;
		}

		private boolean incrementCombination(int[] neighborIndices, ArrayList<? extends ArrayList<Vertex>> neighborsWithinRange, int placeIndex)
		{
			if (placeIndex >= neighborIndices.length)
				return false;

			int index = ++neighborIndices[placeIndex];
			if (index == neighborsWithinRange.get(placeIndex).size())
			{
				neighborIndices[placeIndex] = 0;
				boolean stillValid = incrementCombination(neighborIndices, neighborsWithinRange, placeIndex + 1);
				if (!stillValid)
					return false;
			}
			
			return true;
		}

		private void link(PatchVertex va, PatchVertex vb)
		{
			va.neighbors.add(vb);
			vb.neighbors.add(va);
			edges_abab.add(va.index);
			edges_abab.add(vb.index);
		}

		private PatchVertex addVerticesAndReturnFirst(SteinerSubgraphVertex source, SteinerSubgraphVertex v, Vertex newTerminal, int remainingTerminalCrossingDepth, SubgraphSteinerVertexCache cache, ShortestPathsTable paths) throws ShortestPathsSearchInterruption
		{
			PatchVertex pv = addVertex(v);
			
			PatchVertex result = pv;
			SteinerSubgraphVertex splitAwayNeighbor = null;
			
			if (v.impliesAnyDegenerateEdges())
			{
				int terminalCrossingCost = v.isSteinerVertex ? 0 : 1;
				int nextTerminalCrossingDepth = remainingTerminalCrossingDepth - terminalCrossingCost;
				
				if (nextTerminalCrossingDepth >= 0)
				{
					splitAwayNeighbor = v.getRepresentativeOfNearestEdge(newTerminal, source, cache, paths);
					
					PatchVertex newSteinerVertex = addVertex(new PatchVertex(null, v.v));
					link(newSteinerVertex, pv);
					PatchVertex specialNeighbor = addVerticesAndReturnFirst(v, splitAwayNeighbor, newTerminal, nextTerminalCrossingDepth, cache, paths);
					link(newSteinerVertex, specialNeighbor);
					
					recreatedSteinerVertices.add(new RecreatedSteinerVertex(newSteinerVertex, pv, source, splitAwayNeighbor));
					result = newSteinerVertex;
				}
			}
			
			if (v.isSteinerVertex)
				for (SteinerSubgraphVertex n: v.neighbors)
					if (n != source && n != splitAwayNeighbor)
						link(pv, addVerticesAndReturnFirst(v, n, newTerminal, remainingTerminalCrossingDepth, cache, paths));
			
			return result;
		}
		
		private PatchVertex addVertex(SteinerSubgraphVertex v)
		{
			return addVertex(new PatchVertex(v));
		}

		private PatchVertex addVertex(PatchVertex vertex)
		{
			int index = vertices.size();
			vertices.add(vertex);
			vertex.index = index;
			
			return vertex;
		}
		
		private Vertex setAndReturnGuessedLocation(PatchVertex source, PatchVertex vertex, ShortestPathsTable paths, Graph graph, SubgraphSteinerVertexCache cache) throws ShortestPathsSearchInterruption
		{
			if (vertex.isTerminal())
				return vertex.location;
			else
			{
				Vertex guessedLocation = null;
				int incorporatedGuessCount = 0;
				
				for (PatchVertex neighbor: vertex.neighbors)
					if (neighbor != source)
					{
						Vertex neighboringGuess = setAndReturnGuessedLocation(vertex, neighbor, paths, graph, cache);
						int ratioDenominator = 1 + ++incorporatedGuessCount;
						guessedLocation = guessedLocation != null ? cache.getApproximateMidpoint(guessedLocation, neighboringGuess, ratioDenominator, paths, graph) : neighboringGuess;
					}
				
				vertex.guessedLocation = guessedLocation;
				return vertex.guessedLocation;
			}
		}
		
		private void deterministicallyPlaceVerticesConsumingGuessedLocations(PatchVertex source, PatchVertex vertex, int workingMarkIndex, ShortestPathsTable paths, SubgraphSteinerVertexCache cache) throws ShortestPathsSearchInterruption
		{
			if (!vertex.isTerminal())
			{
				ArrayList<Vertex> referencePoints = new ArrayList<>();
				for (PatchVertex neighbor: vertex.neighbors)
					referencePoints.add(neighbor.getBestGuessAtLocation());
				
				vertex.location = cache.getExactSteinerPoint(referencePoints, workingMarkIndex, paths);
				vertex.guessedLocation = null;
				
				for (PatchVertex neighbor: vertex.neighbors)
					if (neighbor != source)
						deterministicallyPlaceVerticesConsumingGuessedLocations(vertex, neighbor, workingMarkIndex, paths, cache);
			}
		}

		@Override
		public int compareTo(SteinerSubgraphPatch peer)
		{
			return weightChange < peer.weightChange ? -1 : weightChange > peer.weightChange ? 1 : 0;
		}
		
		private static class EdgeToBeReplaced
		{
			private final SteinerSubgraphVertex end0;
			private final SteinerSubgraphVertex end1;

			public EdgeToBeReplaced(SteinerSubgraphVertex end0, SteinerSubgraphVertex end1)
			{
				this.end0 = end0;
				this.end1 = end1;
			}

			public void remove()
			{
				end0.neighbors.remove(end1);
				end1.neighbors.remove(end0);
			}
		}
		
		private static class RecreatedSteinerVertex
		{
			private final PatchVertex vertex;
			private final PatchVertex base;
			private boolean isDeleted;
			private final EdgeToBeReplaced oldEdgeA;
			private final EdgeToBeReplaced oldEdgeB;

			public RecreatedSteinerVertex(PatchVertex vertex, PatchVertex base, SteinerSubgraphVertex neighborA, SteinerSubgraphVertex neighborB)
			{
				this.vertex = vertex;
				this.base = base;
				oldEdgeA = new EdgeToBeReplaced(base.v, neighborA);
				oldEdgeB = new EdgeToBeReplaced(base.v, neighborB);
			}

			public void initializeInSteinerSubgraph()
			{
				if (vertex.location == base.location)
					deleteInPatchGraph();
				else
					vertex.v = new SteinerSubgraphVertex(null, true);
			}

			private void deleteInPatchGraph()
			{
				base.neighbors.remove(vertex);
				
				for (PatchVertex neighbor: vertex.neighbors)
					if (neighbor != base)
					{
						replaceInList(neighbor.neighbors, vertex, base);
						base.neighbors.add(neighbor);
					}
				
				isDeleted = true;
			}

			public void linkIntoSteinerSubgraph()
			{
				if (!isDeleted)
				{
					oldEdgeA.remove();
					oldEdgeB.remove();
					
					for (PatchVertex neighbor: vertex.neighbors)
					{
						addIfNotPresent(neighbor.v.neighbors, vertex.v);
						addIfNotPresent(vertex.v.neighbors, neighbor.v);
					}
				}
			}
		}

		public void apply()
		{
			// these may be invalidated momentarily; best to prevent their further use
			vertices = null;
			edges_abab = null;
			
			splice.v = new SteinerSubgraphVertex(null, true);
			for (RecreatedSteinerVertex v: recreatedSteinerVertices)
				v.initializeInSteinerSubgraph();

			edgeBeforeSplicing.remove();
			for (PatchVertex neighbor: splice.neighbors)
			{
				addIfNotPresent(neighbor.v.neighbors, splice.v);
				addIfNotPresent(splice.v.neighbors, neighbor.v);
			}
			
			for (RecreatedSteinerVertex v: recreatedSteinerVertices)
				v.linkIntoSteinerSubgraph();

			placeAndJoinBranch(terminal, null);
		}

		private static void placeAndJoinBranch(PatchVertex v, PatchVertex source)
		{
			v.v.v = v.location;
			
			if (source != null && v.v.v == source.v.v)
				join(v, source);
			
			for (PatchVertex n: v.neighbors)
				if (n != source)
					placeAndJoinBranch(n, v);
		}

		private static void join(PatchVertex a, PatchVertex b)
		{
			if (b.isTerminal())
			{
				PatchVertex c = a;
				a = b;
				b = c;
			}
			
			SteinerSubgraphVertex av = a.v;
			SteinerSubgraphVertex bv = b.v;
			b.v = av;
			
			av.neighbors.remove(bv);
			bv.neighbors.remove(av);
			
			for (SteinerSubgraphVertex neighbor: bv.neighbors)
				replaceInList(neighbor.neighbors, bv, av);
			av.neighbors.addAll(bv.neighbors);
			bv.neighbors.clear();
		}
		
		private static <Element> void replaceInList(ArrayList<Element> list, Element oldElement, Element newElement)
		{
			list.set(list.indexOf(oldElement), newElement);
		}
		
		private static <Element> void addIfNotPresent(ArrayList<Element> list, Element element)
		{
			if (!list.contains(element))
				list.add(element);
		}
	}

	private void markSteinerSubgraph_(TIntArrayList sortedConnectedTerminals, int markIndex, int workingMarkIndex, final SteinerSubgraphOptions options, final ShortestPathsTable paths, InterruptionSignal interruptionSignal) throws ShortestPathsSearchInterruption
	{
		SteinerSubgraphVertex subgraphRoot = new SteinerSubgraphVertex(vertices.get(sortedConnectedTerminals.get(0)), false);
		
		if (sortedConnectedTerminals.size() >= 2)
		{
			Vertex terminal1 = vertices.get(sortedConnectedTerminals.get(1));
			subgraphRoot.addNeighbor(terminal1, false);
			
			final SubgraphSteinerVertexCache cache = new SubgraphSteinerVertexCache();
			
			for (int i = 2; i < sortedConnectedTerminals.size(); ++i)
			{
				final Vertex terminalI = vertices.get(sortedConnectedTerminals.get(i));
				
				final ArrayList<PotentialEdgeForSplicing> potentialEdgesForSplicing = new ArrayList<>();
				
				subgraphRoot.doForEachEdgeInGraph(new Action2Throwing<SteinerSubgraphVertex, SteinerSubgraphVertex, ShortestPathsSearchInterruption>()
				{
					@Override
					public void perform(SteinerSubgraphVertex v1, SteinerSubgraphVertex v2) throws ShortestPathsSearchInterruption
					{
						SteinerNetStats stats = v1.gatherHiLoPrerequisiteData(v2, terminalI, options.maximalDepthToTraverseInNumberOfImpliedDegenerateEdgesTouchingTerminalsToPassByWhenFindingBoundsOfLocalNetworksOfSteinerVerticesForAdjustment, cache, paths);
						potentialEdgesForSplicing.add(new PotentialEdgeForSplicing(terminalI, v1, v2, stats, paths));
					}
				});
				if (interruptionSignal.isActive()) return;
				
				PotentialEdgeForSplicing edgeWithLowestHi = potentialEdgesForSplicing.get(0);
				for (int j = 1; j < potentialEdgesForSplicing.size(); ++j)
				{
					PotentialEdgeForSplicing e = potentialEdgesForSplicing.get(j);
					if (e.hiBoundOnWeightChange < edgeWithLowestHi.hiBoundOnWeightChange)
						edgeWithLowestHi = e;
				}
				
				SteinerSubgraphVertex terminalInSubgraph = new SteinerSubgraphVertex(terminalI, false);
				SteinerSubgraphPatch leastExpensivePatch = null;
				
				for (PotentialEdgeForSplicing edge: potentialEdgesForSplicing)
					if (edge == edgeWithLowestHi || edge.loBoundOnWeightChange < edgeWithLowestHi.hiBoundOnWeightChange) // the == check is for the case in which lo and hi are somehow exactly the same for edgeWithLowestHi (using <= instead of < would
																										// presumably make this unnecessary, but that < seems potentially slightly faster in general, and doesn't seem to significantly affect correctness...)
					{
						SteinerSubgraphPatch patch = new SteinerSubgraphPatch(edge, terminalInSubgraph, workingMarkIndex, options, paths, this, cache);
						if (leastExpensivePatch == null || patch.compareTo(leastExpensivePatch) < 0)
							leastExpensivePatch = patch;
						if (interruptionSignal.isActive()) return;
					}
				
				leastExpensivePatch.apply();
				if (interruptionSignal.isActive()) return;
			}
			
			cache.reportOnUsage();
		}

		ArrayList<Edge> markedEdges = new ArrayList<>();
		boolean cycleEncountered = subgraphRoot.settleIntoMainGraph(markIndex, markedEdges, paths, this);
		
		if (cycleEncountered && !interruptionSignal.isActive())
		{
			final boolean[] isTerminal = new boolean[vertices.size()];
			for (int i = 0; i < sortedConnectedTerminals.size(); ++i)
				isTerminal[sortedConnectedTerminals.getQuick(i)] = true;

			Sleeve<ShortcutGraphVertex[]> shortcutVerticesByOriginalVertex = new Sleeve<>();
			EdgeArraysGraph effectiveSteinerGraph = getShortcutsGraph(isTerminal, markIndex, subgraphRoot.v, shortcutVerticesByOriginalVertex);
			for (Edge edge: markedEdges)
				edge.clearMarkFullyInSelfAndEnds(markIndex);
			
			int mstMarkIndex = 0;
			effectiveSteinerGraph.markMinimalSpanningTree_Kruskal(mstMarkIndex);
			ShortcutGraphVertex someTerminalNodeInShortcutGraph = shortcutVerticesByOriginalVertex.getIdentity()[sortedConnectedTerminals.get(0)];
			someTerminalNodeInShortcutGraph.asRootOfMarkedTreePruneUninterestingBranches(mstMarkIndex, new Predicate1<Vertex>()
			{
				@Override
				public boolean isTrue(Vertex v)
				{
					ShortcutGraphVertex sv = (ShortcutGraphVertex) v;
					return isTerminal[sv.underlyingVertex.getIndex()];
				}
			});

			someTerminalNodeInShortcutGraph.settleMarkedTreeIntoUnderlyingGraph(mstMarkIndex, markIndex);
		}
	}
	
	public static class InterestEvaluation
	{
		private final boolean[] interestFlags;

		public InterestEvaluation(boolean[] interestFlags)
		{
			this.interestFlags = interestFlags;
		}

		public boolean perform(Vertex vertex)
		{
			return interestFlags != null && interestFlags[vertex.getIndex()] || vertex.edges.size() >= 3;
		}
	}

	/*
	 * this produces a metagraph that has one vertex for each vertex of interest in this graph, and one edge for each stretch of edges and vertices joining two vertices of interest.  a vertex of interest is any vertex that has three or more incident edges
	 * or that is specified as such in the argument (which may be null).  note that the strings of edges that are collapsed into edges in the shortcut graph are entirely linear, in the sense that they do not branch--because a branch is a point where a
	 * vertex has three or more neighbors, and any such point ends the shortcut, and so cannot be within it.  so this graph straightforwardly represents the original graph.
	 */
	public EdgeArraysGraph getShortcutsGraph(boolean[] specialInterestByVertexIndex, int filteringMarkIndex, Vertex suggestedRoot, Sleeve<ShortcutGraphVertex[]> oldVertexToNewVertexMap_out)
	{
		Vertex root = null;
		
		InterestEvaluation interestEvaluation = new InterestEvaluation(specialInterestByVertexIndex);
		MarkFilter markFilter = filteringMarkIndex >= 0 ? new SingleMarkFilter(filteringMarkIndex) : new OpenMarkFilter();
		
		if (suggestedRoot != null && interestEvaluation.perform(suggestedRoot) && markFilter.passes(suggestedRoot))
			root = suggestedRoot;
		else
			for (Vertex vertex: vertices)
				if (interestEvaluation.perform(vertex) && markFilter.passes(vertex))
				{
					root = vertex;
					break;
				}
		
		if (root == null)
			return null;
		
		return root.buildShortCutsGraph(interestEvaluation, markFilter, vertices.size(), oldVertexToNewVertexMap_out);
	}

	@Override
	public float getWeightOfMarkedTree(int markIndex)
	{
		for (Vertex vertex: vertices)
			if (vertex.readMark(markIndex))
				return vertex.getWeightOfMarkedTree(markIndex, null);
		
		return 0;
	}

	public SpaceVector asSpatialGraphGetBoundingBoxDimensions()
	{
		if (vertices.isEmpty())
			return SpaceVector.zero;
		
		SpaceVector lo = new SpaceVector(Float.MAX_VALUE, Float.MAX_VALUE, Float.MAX_VALUE);
		SpaceVector hi = new SpaceVector(-Float.MAX_VALUE, -Float.MAX_VALUE, -Float.MAX_VALUE);
		
		for (Vertex vertex: vertices)
		{
			SpaceVector l = (SpaceVector) ((EuclideanVertex) vertex).getLocation();
			lo = lo.getPiecewiseMinimumWith(l);
			hi = hi.getPiecewiseMaximumWith(l);
		}
		
		return hi.minus(lo);
	}
}
