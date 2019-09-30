package com.github.mjwach.steiner_tree_algorithms;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Queue;

import org.jheaps.AddressableHeap.Handle;
import org.jheaps.tree.PairingHeap;

import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;

public class DCELGraph implements Graph
{
	public Vertex[] vertices;
	private final static Comparator<Vertex> horizontalComparator;
	private final static Comparator<Vertex> verticalComparator;
	
	public static class Vertex implements GraphVertex
	{
		public float x;
		public float y;
		public int index;
		public Edge edge;
		private int mark;
		protected int xIndex;
		protected int yIndex;
		
		public Vertex(int index, PlaneVector value)
		{
			this(index, value.getX(), value.getY());
		}
		
		public Vertex(int index, float x, float y)
		{
			this.index = index;
			this.x = x;
			this.y = y;
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
		
		public void clearMark(int markIndex)
		{
			mark &= ~(1 << markIndex);
		}

		public void clearAllMarks()
		{
			mark = 0;
		}

		@Override
		public boolean readMark(int markIndex)
		{
			return (mark & 1 << markIndex) != 0;
		}

		public boolean sharesLocationWith(Vertex vertex)
		{
			return x == vertex.x && y == vertex.y;
		}
		
		public PlaneVector getLocation()
		{
			return new PlaneVector(x, y);
		}

		public void setLocation(Vertex model)
		{
			this.x = model.x;
			this.y = model.y;
		}
		
		public void setLocation(PlaneVector xy)
		{
			x = xy.getX();
			y = xy.getY();
		}

		public void setLocation(float x, float y)
		{
			this.x = x;
			this.y = y;
		}

		public float getDistanceSquaredFrom(Vertex vertex)
		{
			return getDistanceSquaredFrom(vertex.x, vertex.y);
		}

		public float getDistanceSquaredFrom(float vx, float vy)
		{
			float dx = x - vx;
			float dy = y - vy;
			return dx * dx + dy * dy;
		}

		public double getDistanceFrom(Vertex vertex)
		{
			return Math.sqrt(getDistanceSquaredFrom(vertex));
		}

		public float getDistanceSquaredFrom(PlaneVector location)
		{
			return getDistanceSquaredFrom(location.getX(), location.getY());
		}
		
		@Override
		public String toString()
		{
			return "<" + x + ", " + y + ">";
		}

		public void markVerticesWithinNSteps(int n, int markIndex, TIntArrayList markedVertexIndexRepository, Queue<Vertex> buffer)
		{
			buffer.clear();
			buffer.add(this);
			
			for (int i = 0; i <= n; i++)
				for (int ringSize = buffer.size(), j = 0; j < ringSize; j++)
				{
					Vertex vertex = buffer.poll();
					vertex.setMark(markIndex);
					markedVertexIndexRepository.add(vertex.index);
				
					if (i < n)
						for (Edge nextEdge = vertex.edge, flagEdge = null; nextEdge != vertex.edge || nextEdge != flagEdge; nextEdge = nextEdge.leftNeighbor(), flagEdge = nextEdge)
							if (!nextEdge.destination.readMark(markIndex))
								buffer.offer(nextEdge.destination);
				}
		}

		public void markVerticesWithinNStepsIncreasingNAsNecessaryToSatisfyPredicateAtLeastOnce(int minimalN, int markIndex, TIntArrayList markedVertexIndexRepository,
				Queue<Vertex> buffer, Predicate1<Vertex> predicate)
		{
			buffer.clear();
			buffer.add(this);
			
			boolean predicateHasBeenSatisfied = false;
			for (int i = 0; i <= minimalN || !predicateHasBeenSatisfied; ++i)
				for (int ringSize = buffer.size(), j = 0; j < ringSize; ++j)
				{
					Vertex vertex = buffer.poll();
					vertex.setMark(markIndex);
					markedVertexIndexRepository.add(vertex.index);
					predicateHasBeenSatisfied = predicateHasBeenSatisfied || predicate.isTrue(vertex);
				
					if (i < minimalN || !predicateHasBeenSatisfied)
						for (Edge nextEdge = vertex.edge, flagEdge = null; nextEdge != vertex.edge || nextEdge != flagEdge; nextEdge = nextEdge.leftNeighbor(), flagEdge = nextEdge)
							if (!nextEdge.destination.readMark(markIndex))
								buffer.offer(nextEdge.destination);
				}
		}
		
		@Override
		public void clearMarkFromVerticesInConnectedSubgraph(int markIndex)
		{
			if (readMark(markIndex))
			{
				clearMark(markIndex);
				
				for (Edge nextEdge = edge, flagEdge = null; nextEdge != edge || nextEdge != flagEdge; nextEdge = nextEdge.leftNeighbor(), flagEdge = nextEdge)
					nextEdge.destination.clearMarkFromVerticesInConnectedSubgraph(markIndex);
			}
		}
		
		@Override
		public void clearMarkFromVerticesAndEdgesInConnectedSubgraph(int markIndex)
		{
			if (readMark(markIndex))
			{
				clearMark(markIndex);
				
				for (Edge nextEdge = edge, flagEdge = null; nextEdge != edge || nextEdge != flagEdge; nextEdge = nextEdge.leftNeighbor(), flagEdge = nextEdge)
					if (nextEdge.readMark(markIndex))
					{
						nextEdge.clearMarkInBothDirections(markIndex);
						nextEdge.destination.clearMarkFromVerticesAndEdgesInConnectedSubgraph(markIndex);
					}
			}
		}
		
		@Override
		public void clearMarkFromEdgesInConnectedSubgraph(int markIndex)
		{
			for (Edge nextEdge = edge, flagEdge = null; nextEdge != edge || nextEdge != flagEdge; nextEdge = nextEdge.leftNeighbor(), flagEdge = nextEdge)
				if (nextEdge.readMark(markIndex))
				{
					nextEdge.clearMarkInBothDirections(markIndex);
					nextEdge.destination.clearMarkFromEdgesInConnectedSubgraph(markIndex);
				}
		}
		
		@Override
		public GraphEdge addEdge(GraphVertex destination)
		{
			Vertex v0 = this;
			Vertex v1 = (Vertex) destination;
			
			Edge e0 = new Edge(v1);
			Edge e1 = new Edge(v0);
			e0.opposite = e1;
			e1.opposite = e0;
			Edge leftNeighbor0 = v0.findLeftNeighborForNewEdge(v1);
			Edge rightNeighbor0 = leftNeighbor0 != null ? leftNeighbor0.rightNeighbor() : null;
			Edge leftNeighbor1 = v1.findLeftNeighborForNewEdge(v0);
			Edge rightNeighbor1 = leftNeighbor1 != null ? leftNeighbor1.rightNeighbor() : null;
			if (v0.edge == leftNeighbor0) v0.edge = e0;
			if (v1.edge == leftNeighbor1) v1.edge = e1;
			if (leftNeighbor0 != null) leftNeighbor0.opposite.leftSuccessor = e0;
			if (rightNeighbor0 != null) rightNeighbor0.opposite.rightSuccessor = e0;
			if (leftNeighbor1 != null) leftNeighbor1.opposite.leftSuccessor = e1;
			if (rightNeighbor1 != null) rightNeighbor1.opposite.rightSuccessor = e1;
			if (leftNeighbor0 != null) e1.rightSuccessor = leftNeighbor0; else e1.rightSuccessor = e0;
			if (rightNeighbor0 != null) e1.leftSuccessor = rightNeighbor0; else e1.leftSuccessor = e0;
			if (leftNeighbor1 != null) e0.rightSuccessor = leftNeighbor1; else e0.rightSuccessor = e1;
			if (rightNeighbor1 != null) e0.leftSuccessor = rightNeighbor1; else e0.leftSuccessor = e1;
			
			return e0;
		}
		
		@Override
		public void doForEachEdge(Action1<GraphEdge> action)
		{
			for (Edge e = edge, flagEdge = null; e != flagEdge; e = e.leftNeighbor(), flagEdge = edge)
				action.perform(e);
		}

		public Edge findLeftNeighborForNewEdge(Vertex destination)
		{
			float ox = destination.x - x;
			float oy = destination.y - y;
			
			for (Edge e = edge, flagEdge = null; e != edge || e != flagEdge; e = e.leftNeighbor(), flagEdge = e)
				if (e.crossBodyVectorWith(ox, oy) < 0)
					return e;
			
			return edge;
		}

		@Override
		public GraphEdge findEdgeTo(GraphVertex v)
		{
			for (Edge e = edge, flagEdge = null; e != edge || e != flagEdge; e = e.leftNeighbor(), flagEdge = e)
				if (e.destination == v)
					return e;
			
			return null;
		}
		
		@Override
		public void removeEdge(GraphEdge edge)
		{
			throw new RuntimeException("not implemented yet");
		}
		
		@Override
		public void removeAllEdges()
		{
			throw new RuntimeException("not implemented yet");
		}
	}

	public static class Edge implements GraphEdge
	{
		public Vertex destination;
		public Edge leftSuccessor;
		public Edge rightSuccessor;
		public Edge opposite;
		public boolean isInGraph;
		private int mark;
		
		public static Comparator<Edge> comparatorByLength = new Comparator<Edge>()
		{
			@Override
			public int compare(Edge edge1, Edge edge2)
			{
				float l1 = edge1.getLengthSquared();
				float l2 = edge2.getLengthSquared();
				return l1 < l2 ? -1 : l1 > l2 ? 1 : 0;
			}
		};
		
		public Edge(Vertex destination)
		{
			this.destination = destination;
			isInGraph = true;
		}

		public void setMark(int markIndex)
		{
			mark |= 1 << markIndex;
		}
		
		public void setMark(int markIndex, boolean newValue)
		{
			if (newValue)
				setMark(markIndex);
			else
				clearMark(markIndex);
		}
		
		@Override
		public void setMarkInBothDirections(int markIndex)
		{
			setMark(markIndex);
			opposite.setMark(markIndex);
		}
		
		@Override
		public void clearMarkInBothDirections(int markIndex)
		{
			clearMark(markIndex);
			opposite.clearMark(markIndex);
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
		public boolean readMark(int markIndex)
		{
			return (mark & 1 << markIndex) != 0;
		}

		public Edge successor(boolean left)
		{
			return left ? leftSuccessor : rightSuccessor;
		}

		public void setSuccessor(boolean left, Edge newSuccessor)
		{
			if (left) leftSuccessor = newSuccessor; else rightSuccessor = newSuccessor;
		}

		public Edge leftNeighbor()
		{
			return opposite.rightSuccessor;
		}
		
		public Edge rightNeighbor()
		{
			return opposite.leftSuccessor;
		}

		public Edge neighbor(boolean left)
		{
			return left ? leftNeighbor() : rightNeighbor();
		}
		
		public Edge distinctNeighbor()
		{
			Edge neighbor = leftNeighbor();
			return neighbor != this ? neighbor : null;
		}

		public Edge defaultingToOtherEdgeWithSameDestination()
		{
			return isInGraph ? this : destination.edge.opposite;
		}

		public Edge rightPredecessor()
		{
			return opposite.leftSuccessor.opposite;
		}

		public Edge leftPredecessor()
		{
			return opposite.rightSuccessor.opposite;
		}

		public Edge times(boolean sign)
		{
			return sign ? this : opposite;
		}

		public void setMarks(int markIndex, boolean value)
		{
			setMark(markIndex, value);
			opposite.setMark(markIndex, value);
		}
		
		@Override
		public String toString()
		{
			String sourceString = opposite != null ? opposite.destination.toString() + " " + opposite.destination.index : "null";
			String destinationString = Integer.toString(destination.index) + " " + destination;
			return sourceString + " -> " + destinationString;
		}

		public float getLength()
		{
			return (float) Math.sqrt(getLengthSquared());
		}
		
		public float getLengthSquared()
		{
			Vertex base = opposite.destination;
			float bx = destination.x - base.x;
			float by = destination.y - base.y;
			
			return bx * bx + by * by;
		}

		public float crossBodyVectorWith(float x, float y)
		{
			Vertex base = opposite.destination;
			float bx = destination.x - base.x;
			float by = destination.y - base.y;
			
			return bx * y - x * by;
		}

		private boolean neighboringTriangleContainsGraph(boolean leftNeighbor, BoundingEdgeArray graphBounds)
		{ // this assumes the specified neighboring face is either an ordinary triangular face within the graph somewhere, or the triangular convex hull of the entire graph
			boolean leftBoundHit = false, rightBoundHit = false, bottomBoundHit = false, topBoundHit = false;
			
			Edge edge = this;
			do
			{
				leftBoundHit = leftBoundHit || edge.destination == graphBounds.ccwHullEdgeWithLeftmostDestination.destination;
				rightBoundHit = rightBoundHit || edge.destination == graphBounds.cwHullEdgeWithRightmostDestination.destination;
				bottomBoundHit = bottomBoundHit || edge.destination == graphBounds.cwHullEdgeWithBottommostDestination.destination;
				topBoundHit = topBoundHit || edge.destination == graphBounds.ccwHullEdgeWithTopmostDestination.destination;
				edge = edge.successor(leftNeighbor);
			}
			while (edge != this);
			
			return leftBoundHit && rightBoundHit && bottomBoundHit && topBoundHit;
		}
		
		@Override
		public GraphVertex getEnd0()
		{
			return opposite.destination;
		}
		
		@Override
		public GraphVertex getEnd1()
		{
			return destination;
		}
		
		@Override
		public GraphVertex getEnd(GraphVertex otherEnd)
		{
			return otherEnd == getEnd1() ? getEnd0() : getEnd1();
		}
		
		@Override
		public boolean isNotNegativelyOriented()
		{
			return opposite.destination.getIndex() < destination.getIndex();
		}
	}
	
	static
	{
		horizontalComparator = new Comparator<Vertex>()
		{
			@Override
			public int compare(Vertex vertex1, Vertex vertex2)
			{
				return vertex1.x < vertex2.x ? -1 : vertex1.x > vertex2.x ? 1 : vertex1.y < vertex2.y ? -1 : vertex1.y > vertex2.y ? 1 : 0;
			}
		};

		verticalComparator = new Comparator<Vertex>()
		{
			@Override
			public int compare(Vertex vertex1, Vertex vertex2)
			{
				return vertex1.y < vertex2.y ? -1 : vertex1.y > vertex2.y ? 1 : vertex1.x < vertex2.x ? -1 : vertex1.x > vertex2.x ? 1 : 0;
			}
		};
	}
	
	public DCELGraph(List<? extends Vertex> vertices)
	{
		this(vertices.size());
		
		for (int i = 0; i < vertices.size(); ++i)
			this.vertices[i] = vertices.get(i);
	}
	
	public DCELGraph(List<PlaneVector> vertices, Dummy dummy)
	{
		this(vertices.size());
		
		for (int v = 0; v < vertices.size(); v++)
			this.vertices[v] = new Vertex(v, vertices.get(v));
	}
	
	public DCELGraph(int prospectiveVertexCount)
	{
		this.vertices = new Vertex[prospectiveVertexCount];
	}
	
	public GraphVertex getVertex(int index)
	{
		return vertices[index];
	}
	
	private static class SteinerTreeVertex extends Vertex
	{
		private float indexfulD2FromCenter;
		private final boolean isOriginal;
		public int movementCount;
		public float xBackup, yBackup;

		public SteinerTreeVertex(Vertex model)
		{
			this(model.index, model.x, model.y);
		}
		
		public SteinerTreeVertex(int index, float x, float y)
		{
			super(index, x, y);
			isOriginal = true;
			xBackup = x;
			yBackup = y;
		}

		public SteinerTreeVertex(int index)
		{
			super(index, 0, 0);
			isOriginal = false;
			xBackup = yBackup = Float.NaN; // this will make it clear to any x == xBackup tests or the like that these backups are not to be trusted, so that any movement of the vertex will surely register as actual movement
			x = y = Float.NaN; // this assignment isn't intended to prevent any bugs, but it does help keep things clear
		}

		public SteinerTreeEdge edge()
		{
			return (SteinerTreeEdge) edge;
		}

		public void reportEdgesNearMarkedTriangulationVertices(ArrayList<SteinerTreeEdge> repository, DCELGraph delaunayTriangulation, int markIndex)
		{
			for (SteinerTreeEdge nextEdge = edge(), flagEdge = null; nextEdge != edge() || nextEdge != flagEdge; nextEdge = nextEdge.leftNeighbor_(), flagEdge = nextEdge)
				if (!nextEdge.readMark(markIndex) && delaunayTriangulation.vertices[nextEdge.closestOriginalVertex.index].readMark(markIndex))
				{
					repository.add(nextEdge);
					nextEdge.setMarks(markIndex, true);
					nextEdge.destination().reportEdgesNearMarkedTriangulationVertices(repository, delaunayTriangulation, markIndex);
				}
		}

		public void takeLocallyCorrectPosition(ArrayList<SteinerTreeVertex> affectedNeighborsRepository, float insignificantMovementThresholdSquared, int markIndex)
		{
			Edge edge_ = edge;
			Vertex v0 = edge_.destination;
			edge_ = edge_.leftNeighbor();
			Vertex v1 = edge_.destination;
			edge_ = edge_.leftNeighbor();
			Vertex v2 = edge_.destination;
			
			float oldX = x, oldY = y;
			PlaneVector.getSteinerPoint(v0, v1, v2, this);
			float dx = x - oldX, dy = y - oldY;

			movementCount++;
			
			for (SteinerTreeEdge nextEdge = edge(), flagEdge = null; nextEdge != edge() || nextEdge != flagEdge; nextEdge = nextEdge.leftNeighbor_(), flagEdge = nextEdge)
			{
				nextEdge.movementCount = true;
				nextEdge.opposite().movementCount = true;
			}
			
			if (!(dx * dx + dy * dy <= insignificantMovementThresholdSquared)) // dx and dy may be NaN in which case a more straightforward test involving > would fail
				for (SteinerTreeEdge nextEdge = edge(), flagEdge = null; nextEdge != edge() || nextEdge != flagEdge; nextEdge = nextEdge.leftNeighbor_(), flagEdge = nextEdge)
				{
					SteinerTreeVertex neighbor = nextEdge.destination();
				
					if (!neighbor.isOriginal && !neighbor.readMark(markIndex))
					{
						affectedNeighborsRepository.add(neighbor);
						neighbor.setMark(markIndex);
					}
				}
		}

		public void reportMovements(ArrayList<SteinerTreeVertex> vertexRepository, ArrayList<SteinerTreeEdge> edgeRepository)
		{
			if (movementCount > 0)
			{
				movementCount = 0;
				vertexRepository.add(this);
				
				for (SteinerTreeEdge nextEdge = edge(), flagEdge = null; nextEdge != edge() || nextEdge != flagEdge; nextEdge = nextEdge.leftNeighbor_(), flagEdge = nextEdge)
					if (nextEdge.movementCount)
					{
						nextEdge.movementCount = false;
						nextEdge.opposite().movementCount = false;
						edgeRepository.add(nextEdge);
						
						nextEdge.destination().reportMovements(vertexRepository, edgeRepository);
					}
			}
		}

		public int removeDegenerateSteinerNeighborsAndSortEdges(SteinerTreeVertex[] vertexRepository, ArrayList<SteinerTreeEdge> edgeBuffer, int markIndex)
		{
			return removeDegenerateSteinerNeighborsAndSortEdges(vertexRepository, edgeBuffer, null, markIndex);
		}
		
		private int removeDegenerateSteinerNeighborsAndSortEdges(SteinerTreeVertex[] vertexRepository, ArrayList<SteinerTreeEdge> edgeBuffer, SteinerTreeEdge callSource_opposite, int markIndex)
		{
			// first, use this method recursively to clean up degenerate edges beyond this vertex's degenerate edges; leave this
			// vertex's degenerate edges marked
			int degenerateVertexCount = 0;
			
			for (SteinerTreeEdge nextEdge = edge(), flagEdge = null; nextEdge != edge() || nextEdge != flagEdge; nextEdge = nextEdge.leftNeighbor_(), flagEdge = nextEdge)
				if (!nextEdge.readMark(markIndex))
				{
					SteinerTreeVertex destination = nextEdge.destination();
					
					if (!destination.isOriginal && destination.x == x && destination.y == y)
					{
						degenerateVertexCount++;
						vertexRepository[destination.index] = null; // note that this change is noticed below; thus it must happen before the recursive call immediately below
						
						nextEdge.setMarks(markIndex, true); // incidentally:  marked edges are degenerate and will be removed, so there'll be no need to unmark them
						degenerateVertexCount += destination.removeDegenerateSteinerNeighborsAndSortEdges(vertexRepository, edgeBuffer, nextEdge.opposite(), markIndex);
					}
				}

			// now gather the non-degenerate edges that will end up connected to this vertex once degenerate edges have been removed
			edgeBuffer.clear();
			ArrayList<SteinerTreeEdge> combinedEdgeList = edgeBuffer;
			
			for (SteinerTreeEdge nextEdge = edge(), flagEdge = null; nextEdge != edge() || nextEdge != flagEdge; nextEdge = nextEdge.leftNeighbor_(), flagEdge = nextEdge)
				if (!nextEdge.readMark(markIndex) || nextEdge == callSource_opposite)
					combinedEdgeList.add(nextEdge);
				else
					for (SteinerTreeEdge neighborEdge0 = nextEdge.opposite(), neighborEdge = neighborEdge0, nFlagEdge = null; neighborEdge != neighborEdge0 || neighborEdge != nFlagEdge; neighborEdge = neighborEdge.leftNeighbor_(), nFlagEdge = neighborEdge)
						if (!neighborEdge.readMark(markIndex))
							combinedEdgeList.add(neighborEdge);
			
			// finally, correctly link those edges to this vertex, sorting them first since that needs to be done (except where this
			// vertex is itself a redundant one that's been reached via a degenerate edge) and it's handy to do it now when all the
			// edges are together in a random-access list
			for (int i = 0; i < combinedEdgeList.size(); i++)
				combinedEdgeList.get(i).opposite().destination = this;
			
			if (vertexRepository[index] != null)
				Collections.sort(combinedEdgeList);
			
			edge = combinedEdgeList.get(0);

			for (int i0 = 0; i0 < combinedEdgeList.size(); i0++)
			{
				int i1 = i0 + 1;
				if (i1 >= combinedEdgeList.size()) i1 = 0;
				
				SteinerTreeEdge edge0 = combinedEdgeList.get(i0);
				SteinerTreeEdge edge1 = combinedEdgeList.get(i1);
				
				edge0.opposite.rightSuccessor = edge1;
				edge1.opposite.leftSuccessor = edge0;
			}
			
			return degenerateVertexCount;
		}
	}
	
	private static class SteinerTreeEdge extends Edge implements Comparable<SteinerTreeEdge>
	{
		public SteinerTreeVertex closestOriginalVertex;
		private int closestOriginalVertexDistance;
		private boolean closestOriginalVertexDirection; // true iff the vertex is closer to the destination end of the edge than to the source end
		public boolean movementCount;
		private double length;

		public SteinerTreeEdge(SteinerTreeVertex destination)
		{
			super(destination);
			length = 0; // it's assumed elsewhere that new edges know their length to be 0, so it's nice to set that explicitly
		}

		public SteinerTreeVertex destination()
		{
			return (SteinerTreeVertex) destination;
		}
		
		public SteinerTreeEdge opposite()
		{
			return (SteinerTreeEdge) opposite;
		}

		public SteinerTreeEdge leftSuccessor()
		{
			return (SteinerTreeEdge) leftSuccessor;
		}
		
		public SteinerTreeEdge rightSuccessor()
		{
			return (SteinerTreeEdge) rightSuccessor;
		}
		
		public SteinerTreeEdge leftNeighbor_()
		{
			return (SteinerTreeEdge) leftNeighbor();
		}

		public double getLength_()
		{
			SteinerTreeVertex v0 = destination();
			SteinerTreeVertex v1 = opposite().destination();
			
			if (v0.x == v0.xBackup && v0.y == v0.yBackup && v1.x == v1.xBackup && v1.y == v1.yBackup)
				return length;
			else
				return v0.getDistanceFrom(v1);
		}

		public void setClosestOriginalVertex(SteinerTreeVertex vertex, int distanceToVertex, boolean direction)
		{
			SteinerTreeEdge opposite_ = opposite();

			this.closestOriginalVertex = opposite_.closestOriginalVertex = vertex;
			this.closestOriginalVertexDistance = opposite_.closestOriginalVertexDistance = distanceToVertex;
			
			this.closestOriginalVertexDirection = direction;
			opposite_.closestOriginalVertexDirection = !direction;
		}

		public void updateOldClosestOriginalVertexValuesAheadOfSplicingZone(SteinerTreeEdge partnerAheadOfSplice)
		{
			if (closestOriginalVertexDirection == false &&
					(partnerAheadOfSplice.closestOriginalVertex != closestOriginalVertex || !partnerAheadOfSplice.closestOriginalVertexDirection))
				updateClosestOriginalVertexValues(closestOriginalVertex, closestOriginalVertex, closestOriginalVertexDistance + 1, partnerAheadOfSplice);
		}
		
		private void updateClosestOriginalVertexValues(SteinerTreeVertex affectedClosestVertex, SteinerTreeVertex replacingVertex, int replacingVertexDistance,
				SteinerTreeEdge otherPossibleReplacementSource)
		{
			if (closestOriginalVertex == affectedClosestVertex)
			{
				if (otherPossibleReplacementSource.closestOriginalVertex != replacingVertex && otherPossibleReplacementSource.closestOriginalVertexDistance + 1 < replacingVertexDistance)
				{
					replacingVertex = otherPossibleReplacementSource.closestOriginalVertex;
					replacingVertexDistance = otherPossibleReplacementSource.closestOriginalVertexDistance + 1;
				}

				SteinerTreeEdge opposite_ = opposite();
				SteinerTreeEdge left = leftSuccessor();
				SteinerTreeEdge right = rightSuccessor();
				
				boolean replacingDirection = closestOriginalVertexDirection;
				
				if (left != opposite_ && left.closestOriginalVertex != affectedClosestVertex && left.closestOriginalVertexDistance + 1 < replacingVertexDistance)
				{
					replacingVertex = left.closestOriginalVertex;
					replacingVertexDistance = left.closestOriginalVertexDistance + 1; // this could conceivably be a "- 1", depending on left's closest vertex direction--but it turns out that in any case where that's true, the value being calculated here will not actually get used
					replacingDirection = true;
				}
				
				if (right != opposite_ && right.closestOriginalVertex != affectedClosestVertex && right.closestOriginalVertexDistance + 1 < replacingVertexDistance)
				{
					replacingVertex = right.closestOriginalVertex;
					replacingVertexDistance = right.closestOriginalVertexDistance + 1;
					replacingDirection = true;
				}
				
				setClosestOriginalVertex(replacingVertex, replacingVertexDistance, replacingDirection);
				
				if (left != opposite_) left.updateClosestOriginalVertexValues(affectedClosestVertex, replacingVertex, replacingVertexDistance + 1, right);
				if (right != opposite_) right.updateClosestOriginalVertexValues(affectedClosestVertex, replacingVertex, replacingVertexDistance + 1, left);
			}
		}

		public void propagateNewClosestOriginalVertexRegionForward()
		{
			propagateClosestOriginalVertexRegionForward_(closestOriginalVertex, closestOriginalVertexDistance);
		}

		private void propagateClosestOriginalVertexRegionForward_(SteinerTreeVertex vertexForRegion, int vertexDistance)
		{
			int nextDistance = vertexDistance + 1;
			leftSuccessor().propagateClosestOriginalVertexRegionForward(vertexForRegion, nextDistance);
			rightSuccessor().propagateClosestOriginalVertexRegionForward(vertexForRegion, nextDistance);
		}

		private void propagateClosestOriginalVertexRegionForward(SteinerTreeVertex vertexForRegion, int vertexDistance)
		{
			if (vertexDistance < closestOriginalVertexDistance) // a <= would work here too, but then the last edge's opposite might qualify, and that'd be an extra thing to think about and maybe guard against...
			{
				setClosestOriginalVertex(vertexForRegion, vertexDistance, false);
				propagateClosestOriginalVertexRegionForward_(vertexForRegion, vertexDistance);
			}
		}

		@Override
		public int compareTo(SteinerTreeEdge peer)
		{ // this makes the most sense when the edges share sources (or destinations)
			Vertex source0 = opposite.destination;
			Vertex source1 = peer.opposite.destination;
			
			float x0 = destination.x - source0.x;
			float y0 = destination.y - source0.y;
			float x1 = peer.destination.x - source1.x;
			float y1 = peer.destination.y - source1.y;
			
			return PlaneVector.compareVectorsByAngle(x0, y0, x1, y1);
		}
	}
	
	public enum SteinerDuplicatePointPolicy
	{
		includeOnlyOneVertexPerLocation,
		includeOneVertexPerInputPointWithCoincidentVerticesArtificiallyConnectedViaFansOfZeroLengthEdges
	}
	
	public static DCELGraph makeSteinerTree(List<PlaneVector> vertices)
	{
		return makeSteinerTree(vertices, SteinerDuplicatePointPolicy.includeOnlyOneVertexPerLocation);
	}

	public static DCELGraph makeSteinerTree(List<PlaneVector> vertices, SteinerDuplicatePointPolicy duplicatePointPolicy)
	{
		return makeSteinerTreeAndDelaunayTriangulation(vertices, duplicatePointPolicy).steinerTree;
	}

	public static CalculatedSteinerTree makeSteinerTreeAndDelaunayTriangulation(List<PlaneVector> vertices, SteinerDuplicatePointPolicy duplicatePointPolicy)
	{
		return makeSteinerTreeAndDelaunayTriangulation(vertices, duplicatePointPolicy, 4, .0008f);
	}
	
	public static class CalculatedSteinerTree
	{
		public DCELGraph steinerTree;
		public DCELGraph delaunayTriangulation;
		public int[] mapFromInputVerticesToSteinerTreeVertices;
	}
	
	public static CalculatedSteinerTree makeSteinerTreeAndDelaunayTriangulation(List<PlaneVector> vertices, SteinerDuplicatePointPolicy duplicatePointPolicy,
			int steinerPointPositionSearchDistance, float steinerRelativeInsignificantMovementThreshold)
	// The basic algorithm used here is a modification of an "Incremental Optimization Heuristic" presented in "Two Heuristics for
	// the Euclidean Steiner Tree Problem" by Derek R. Dreyer and Michael L. Overton (June 11, 1997); this version greatly reduces
	// the running time of the original algorithm by first computing a Delaunay triangulation, which it then traverses, whenever a
	// vertex is to be added to the growing Steiner tree, in order to find edges in the tree that are "near" enough to the new vertex
	// to serve as potential linking sites where that vertex may join the tree - the original algorithm would consider every edge in
	// the tree to be a potential linking site, so if the maximal degree of a vertex in the triangulation is not too great, then this
	// modification may remove something like a O(n) factor from the running time.
	//
	// The numeric parameters here are a "search distance" determining the radius (in edges traversed) of the portion of the Delaunay
	// triangulation to be searched for potential linking sites, and a "relative movement threshold" that's invoked when a portion of
	// the Steiner tree is to be given approximately optimal vertex positions:  The Steiner vertices in that part of the graph are
	// repeatedly positioned in locally optimal ways, so that the whole subgraph tends to gradually shimmy into a nearly optimal
	// position; and during that process, vertices are considered to have done enough shimmying when their per-step movement is
	// smaller than the distance specified by this parameter (which distance is actually multiplied by a factor that grows with the
	// graph's bounding box's size and shrinks as vertex count grows, so that ideally the same value of the parameter will serve
	// equally well for vertex sets covering many different ranges and having many different point densities).
	//
	// In case more control is desired, two more parameters for the "shimmying" system may be set via another override of this method.
	{
		return makeSteinerTreeAndDelaunayTriangulation(vertices, duplicatePointPolicy, steinerPointPositionSearchDistance, steinerRelativeInsignificantMovementThreshold, 1.003f, 48);
	}
	
	public static CalculatedSteinerTree makeSteinerTreeAndDelaunayTriangulation(List<PlaneVector> vertices, SteinerDuplicatePointPolicy duplicatePointPolicy,
			int steinerPointPositionSearchDistance, float steinerRelativeInsignificantMovementThreshold, float insignificantMovementThresholdGrowthFactor, int hardShimmyingIterationLimit)
	{
		if (steinerPointPositionSearchDistance <= 0)
			throw new IllegalArgumentException("search distance " + steinerPointPositionSearchDistance);

		if (vertices.size() <= 2)
		{
			CalculatedSteinerTree result = new CalculatedSteinerTree();
			result.mapFromInputVerticesToSteinerTreeVertices = new int[vertices.size()];
			result.delaunayTriangulation = makeDelaunayTriangulation(vertices);
			if (vertices.size() == 0)
				result.steinerTree = makeDelaunayTriangulation(vertices);
			else
				result.steinerTree = makeOneOrTwoPointSteinerTree(vertices.get(0), vertices.size() == 2 ? vertices.get(1) : null, duplicatePointPolicy, result.mapFromInputVerticesToSteinerTreeVertices);
			
			return result;
		}

		final SteinerTreeVertex[] treeVerticesByX = new SteinerTreeVertex[vertices.size() * 2 - 2];
		
		final int[] originalPointsSorted = new int[vertices.size()];
		
		Action1<Vertex[]> originalPointOrdering = new Action1<Vertex[]>()
		{
			@Override
			public void perform(Vertex[] allSortedVerticesIncludingDuplicates)
			{
				for (int v = 0; v < allSortedVerticesIncludingDuplicates.length; v++)
					originalPointsSorted[v] = allSortedVerticesIncludingDuplicates[v].index;
			}
		};
		
		Action1<Vertex> treeVertexCreation = new Action1<DCELGraph.Vertex>()
		{
			@Override
			public void perform(Vertex vertex)
			{
				SteinerTreeVertex treeVertex = new SteinerTreeVertex(vertex);
				treeVertex.xIndex = vertex.xIndex;
				treeVertex.yIndex = vertex.yIndex;
				
				treeVerticesByX[treeVertex.xIndex] = treeVertex;
			}
		};

		int[] mapFromOldIndicesToNewIndices = new int[vertices.size()];
		DCELGraph delaunayTriangulation = makeDelaunayTriangulation(vertices, DelaunayDuplicatePointPolicy.makeOneVertexPerInputLocationButInvolveOnlyTheFirstVertexPerLocationInTheTriangulation, mapFromOldIndicesToNewIndices,
				treeVertexCreation, originalPointOrdering);
		
		Vertex[] oldTriangulationVertices = Arrays.copyOf(delaunayTriangulation.vertices, delaunayTriangulation.vertices.length);
		int vertexCount = 0;
		
		{ // temporarily removing unused vertices from the Delaunay triangulation to simplify Steiner tree development
			for (int v = 0; v < delaunayTriangulation.vertices.length; v++)
			{
				Vertex vertex = delaunayTriangulation.vertices[v];
	
				if (vertex.edge != null)
				{
					vertex.index = vertexCount;
					delaunayTriangulation.vertices[vertex.index] = vertex;
					vertexCount++;
				}
			}
			
			if (vertexCount == 0) // the test above for a vertex that is in the triangulation simply checks for an edge, which is incorrect in one case:  the case of a single-vertex "triangulation".  such a triangulation arises when all points have been
				vertexCount = 1;  // found to be duplicates of one another.  in that case, the aforementioned test is incorrect; any one (but only one) of the duplicates could be considered to be a valid member of the triangulation.  setting count to 1 now
								  // expresses that the first vertex is the one so chosen.  it is already in the right spot in its array, and its index need not be altered.
			
			for (int i = vertexCount; i < delaunayTriangulation.vertices.length; i++)
				delaunayTriangulation.vertices[i] = null;
		}
		
		int maximalSteinerVertexCount = Math.max(0, vertexCount - 2); // zero, one, or two vertices in a graph can form the skeleton of a Steiner tree with no additional vertices; beyond that, each additional point may require one more Steiner point
		int maximalSteinerTreeVertexCount = vertexCount + maximalSteinerVertexCount;
		SteinerTreeVertex[] treeVertices = new SteinerTreeVertex[maximalSteinerTreeVertexCount];
		
		for (int v = 0; v < vertexCount; v++) // realigning indices of vertices in the Steiner tree to match any changes made while duplicate
		{									// locations were being removed from the triangulation
			SteinerTreeVertex treeVertex = treeVerticesByX[v];
			int newIndex = oldTriangulationVertices[treeVertex.index].index;
			treeVertex.index = newIndex;
			treeVertices[newIndex] = treeVertex;
		}
		
		DCELGraph steinerTree;

		if (vertexCount >= 3)
		{
			float maximalY = -Float.MAX_VALUE, minimalY = Float.MAX_VALUE;
			
			for (int v = 0; v < vertexCount; v++)
			{
				Vertex vertex = treeVertices[v];
				if (vertex.y < minimalY) minimalY = vertex.y;
				if (vertex.y > maximalY) maximalY = vertex.y;
			}
			
			float graphWidth = maximalY - minimalY;
			float graphHeight = treeVerticesByX[vertexCount - 1].x - treeVerticesByX[0].x;
			
			float threshold = steinerRelativeInsignificantMovementThreshold * (graphWidth + graphHeight) / (float) Math.sqrt(vertexCount);
			float insignificantMovementThreshold = Math.max(threshold, Float.MIN_VALUE * 8);
	
			steinerTree = makeSteinerTree(treeVertices, treeVerticesByX, vertexCount, delaunayTriangulation, steinerPointPositionSearchDistance, insignificantMovementThreshold, 1.003f, 48);
		}
		else // too many duplicate vertices have been stripped out for the main algorithm to be appropriate
			steinerTree = makeOneOrTwoPointSteinerTree(treeVertices[0].getLocation(), vertexCount == 2 ? treeVertices[1].getLocation() : null, duplicatePointPolicy, null);
		
		for (int v = 0; v < oldTriangulationVertices.length; v++) oldTriangulationVertices[v].index = v;
		
		if (vertexCount != vertices.size() && duplicatePointPolicy == SteinerDuplicatePointPolicy.includeOneVertexPerInputPointWithCoincidentVerticesArtificiallyConnectedViaFansOfZeroLengthEdges)
		{
			int[] includedVerticesByCoincidentVertex = new int[vertices.size()];
			
			for (int v = 0; v < vertices.size(); v++)
			{
				int p = originalPointsSorted[v];
				Vertex includedVertex = oldTriangulationVertices[p];
				
				if (includedVertex.edge != null)
				{
					for (int i = v; 	i < originalPointsSorted.length && includedVertex.sharesLocationWith(oldTriangulationVertices[originalPointsSorted[i]]); i++)
						includedVerticesByCoincidentVertex[originalPointsSorted[i]] = p;
					for (int i = v - 1; i >= 0 							&& includedVertex.sharesLocationWith(oldTriangulationVertices[originalPointsSorted[i]]); i--)
						includedVerticesByCoincidentVertex[originalPointsSorted[i]] = p;
				}
			}
			
			Vertex[] newSteinerVertices = new Vertex[steinerTree.vertices.length + vertices.size() - vertexCount];
			
			for (int v = 0, n = vertices.size(); v < steinerTree.vertices.length; v++)
			{ // this aligns Steiner tree vertices with the original input, leaving gaps for duplicate input locations
				SteinerTreeVertex vertex = (SteinerTreeVertex) steinerTree.vertices[v];
				int newIndex;

				if (vertex.isOriginal)
				{
					newIndex = delaunayTriangulation.vertices[vertex.index].index;
					mapFromOldIndicesToNewIndices[newIndex] = newIndex;
				}
				else
					newIndex = n++;

				vertex.index = newIndex;
				newSteinerVertices[newIndex] = vertex;
			}
			
			for (int v = 0; v < newSteinerVertices.length; v++)
				if (newSteinerVertices[v] == null)
				{
					Vertex vertex = new Vertex(v, vertices.get(v));
					newSteinerVertices[v] = vertex;

					int fanRootIndex = includedVerticesByCoincidentVertex[v];
					Vertex root = newSteinerVertices[fanRootIndex];
					Edge edge0 = new Edge(root);
					Edge edge1 = new Edge(vertex);
					vertex.edge = edge0;
					edge0.opposite = edge1;
					edge1.opposite = edge0;
					edge1.leftSuccessor = edge1.rightSuccessor = edge1.opposite;
					edge0.rightSuccessor = root.edge;
					edge0.leftSuccessor = root.edge.rightNeighbor();
					edge0.rightSuccessor.opposite.leftSuccessor = edge1;
					edge0.leftSuccessor.opposite.rightSuccessor = edge1;
					root.edge = edge1;
				}
			
			steinerTree.vertices = newSteinerVertices;
		}
		
		delaunayTriangulation.vertices = oldTriangulationVertices;
		
		CalculatedSteinerTree result = new CalculatedSteinerTree();
		result.steinerTree = steinerTree;
		result.delaunayTriangulation = delaunayTriangulation;
		result.mapFromInputVerticesToSteinerTreeVertices = mapFromOldIndicesToNewIndices;
		
		return result;
	}
	
	private static DCELGraph makeOneOrTwoPointSteinerTree(PlaneVector vertex0, PlaneVector vertex1, SteinerDuplicatePointPolicy duplicatePointPolicy, int[] mapFromInputVerticesToSteinerTreeVertices)
	{
		boolean twoVertices = vertex1 != null &&
			(!vertex0.equals(vertex1) || duplicatePointPolicy == SteinerDuplicatePointPolicy.includeOneVertexPerInputPointWithCoincidentVerticesArtificiallyConnectedViaFansOfZeroLengthEdges);

		DCELGraph graph = new DCELGraph(twoVertices ? 2 : 1);
		
		graph.vertices[0] = new SteinerTreeVertex(0, vertex0.getX(), vertex0.getY());
		if (mapFromInputVerticesToSteinerTreeVertices != null) mapFromInputVerticesToSteinerTreeVertices[0] = 0;
		
		if (twoVertices)
		{
			graph.vertices[1] = new SteinerTreeVertex(1, vertex1.getX(), vertex1.getY());
			if (vertex1 != null && mapFromInputVerticesToSteinerTreeVertices != null) mapFromInputVerticesToSteinerTreeVertices[1] = 1;
		}
		else
			if (vertex1 != null && mapFromInputVerticesToSteinerTreeVertices != null) mapFromInputVerticesToSteinerTreeVertices[1] = 0;
		
		if (twoVertices)
		{
			Vertex v0 = graph.vertices[0];
			Vertex v1 = graph.vertices[1];
			
			v0.edge = new Edge(v1);
			v1.edge = new Edge(v0);
			v0.edge.opposite = v0.edge.leftSuccessor = v0.edge.rightSuccessor = v1.edge;
			v1.edge.opposite = v1.edge.leftSuccessor = v1.edge.rightSuccessor = v0.edge;
		}
		
		return graph;
	}

	private static DCELGraph makeSteinerTree(SteinerTreeVertex[] vertices, SteinerTreeVertex[] verticesByX, int vertexCount, DCELGraph delaunayTriangulation, int searchDistance, float insignificantMovementThreshold, float insignificantMovementThresholdGrowthFactor, int hardShimmyingIterationLimit)
	{
		int markIndex = 0;
		
		// find a vertex that's near the center of the graph in some sense
		float halfVertexCount = vertexCount * .5f;

		for (int v = 0; v < vertexCount; v++)
		{
			SteinerTreeVertex vertex = vertices[v];
			float x = vertex.xIndex - halfVertexCount;
			float y = vertex.yIndex - halfVertexCount;
			vertex.indexfulD2FromCenter = x * x + y * y;
		}
		
		SteinerTreeVertex nearestVertexToCenter = null;
		for (int i = 0, c = vertexCount >> 1; i < vertexCount; i++)
		{
			int v = (i & 1) == 1 ? c + i : c - i;
			SteinerTreeVertex vertex = verticesByX[v];
			
			if (nearestVertexToCenter == null || vertex.indexfulD2FromCenter < nearestVertexToCenter.indexfulD2FromCenter)
				nearestVertexToCenter = vertex;
			
			float x = v - halfVertexCount;
			
			if (nearestVertexToCenter.indexfulD2FromCenter <= 2 * x * x)
				break;
		}
		
		// build an array of vertices, more or less ordered by distance from the graph's center; these will join the tree one by one
		ArrayList<SteinerTreeVertex> verticesByDistanceFromCenter = new ArrayList<>();
		
		{
			Vertex nearestTriangulationVertexToCenter = delaunayTriangulation.vertices[nearestVertexToCenter.index];
			
			Queue<Vertex> verticesToExplore = new ArrayDeque<>(vertexCount);
			nearestTriangulationVertexToCenter.setMark(markIndex);
			verticesToExplore.add(nearestTriangulationVertexToCenter);
	
			while (!verticesToExplore.isEmpty())
			{
				Vertex vertex = verticesToExplore.poll();
				verticesByDistanceFromCenter.add(vertices[vertex.index]);
				
				int edgeCount = 0;
				Edge bestEdge = vertex.edge;
				for (Edge edge = vertex.edge, flagEdge = null; edge != flagEdge || edge != vertex.edge; edge = edge.leftNeighbor(), flagEdge = edge, edgeCount++)
					if (vertices[edge.destination.index].indexfulD2FromCenter < vertices[bestEdge.destination.index].indexfulD2FromCenter)
						bestEdge = edge;
				
				for (Edge leftEdge = bestEdge, rightEdge = bestEdge.rightNeighbor(); edgeCount > 0; leftEdge = leftEdge.leftNeighbor(), rightEdge = rightEdge.rightNeighbor())
				{
					if (!leftEdge.destination.readMark(markIndex))
					{
						leftEdge.destination.setMark(markIndex);
						verticesToExplore.offer(leftEdge.destination);
					}
					edgeCount--;

					if (!rightEdge.destination.readMark(markIndex) && edgeCount > 0)
					{
						rightEdge.destination.setMark(markIndex);
						verticesToExplore.offer(rightEdge.destination);
					}
					edgeCount--;
				}
			}
			
			for (int v = 0; v < vertexCount; v++)
				delaunayTriangulation.vertices[v].clearMark(markIndex);
		}
		
		// build a Steiner tree from the first three vertices
		int v = 0;
		double graphLength; // the sum of all edges' lengths
		
		{
			SteinerTreeVertex v0 = verticesByDistanceFromCenter.get(0);
			SteinerTreeVertex v1 = verticesByDistanceFromCenter.get(1);
			SteinerTreeVertex v2 = verticesByDistanceFromCenter.get(2);
			
			SteinerTreeVertex steinerPoint = new SteinerTreeVertex(vertexCount);
			PlaneVector.getSteinerPoint(v0, v1, v2, steinerPoint);
			steinerPoint.xBackup = steinerPoint.x; steinerPoint.yBackup = steinerPoint.y;
			vertices[vertexCount++] = steinerPoint;
			verticesByX = null; // this has probably just become invalid; it isn't needed now anyway, so it shouldn't be accessible
			
			v0.edge = new SteinerTreeEdge(steinerPoint);
			Edge e0 = new SteinerTreeEdge(v0);
			v1.edge = new SteinerTreeEdge(steinerPoint);
			Edge e1 = new SteinerTreeEdge(v1);
			v2.edge = new SteinerTreeEdge(steinerPoint);
			Edge e2 = new SteinerTreeEdge(v2);
			steinerPoint.edge = e0;
			
			v0.edge.opposite = e0;
			e0.opposite = v0.edge = e0.leftSuccessor = e0.rightSuccessor = v0.edge;
			v1.edge.opposite = e1;
			e1.opposite = v1.edge = e1.leftSuccessor = e1.rightSuccessor = v1.edge;
			v2.edge.opposite = e2;
			e2.opposite = v2.edge = e2.leftSuccessor = e2.rightSuccessor = v2.edge;
			
			v0.edge.leftSuccessor = e2;
			v0.edge.rightSuccessor = e1;
			v1.edge.leftSuccessor = e0;
			v1.edge.rightSuccessor = e2;
			v2.edge.leftSuccessor = e1;
			v2.edge.rightSuccessor = e0; // if the three vertices are in clockwise order then these assignments are wrong--but that'll be corrected later
			
			v0.edge().setClosestOriginalVertex(v0, 0, false);
			v1.edge().setClosestOriginalVertex(v1, 0, false);
			v2.edge().setClosestOriginalVertex(v2, 0, false);
			
			double distance0 = v0.getDistanceFrom(steinerPoint);
			double distance1 = v1.getDistanceFrom(steinerPoint);
			double distance2 = v2.getDistanceFrom(steinerPoint);
			v0.edge().length = v0.edge().opposite().length = distance0;
			v1.edge().length = v1.edge().opposite().length = distance1;
			v2.edge().length = v2.edge().opposite().length = distance2;
			
			graphLength = distance0 + distance1 + distance2;
			
			v += 3;
		}
		
		// add the remaining vertices to the tree
		class SteinerVertexAdditionResult
		{
			private final SteinerTreeVertex newOriginalVertex;
			private final SteinerTreeVertex newSteinerVertex;
			private final ArrayList<SteinerTreeVertex> movedVertices;
			private final ArrayList<SteinerTreeEdge> movedEdges;
			private final float[] vertexLocationsX;
			private final float[] vertexLocationsY;
			private final double[] edgeLengths;
			public final double graphLength;

			public SteinerVertexAdditionResult(SteinerTreeVertex newOriginalVertex, SteinerTreeVertex newSteinerVertex, SteinerTreeEdge splicedEdge, double oldGraphLength)
			{
				this.newOriginalVertex = newOriginalVertex;
				this.newSteinerVertex = newSteinerVertex;
				
				movedVertices = new ArrayList<SteinerTreeVertex>();
				movedEdges = new ArrayList<>();
				
				newSteinerVertex.reportMovements(movedVertices, movedEdges);
				
				vertexLocationsX = new float[movedVertices.size()];
				vertexLocationsY = new float[movedVertices.size()];
				edgeLengths = new double[movedEdges.size()];
				
				for (int i = 0; i < movedVertices.size(); i++)
				{
					Vertex movedVertex = movedVertices.get(i);
					vertexLocationsX[i] = movedVertex.x;
					vertexLocationsY[i] = movedVertex.y;
				}
				
				for (int i = 0; i < movedEdges.size(); i++)
				{
					SteinerTreeEdge movedEdge = movedEdges.get(i);
					double currentEdgeLength = movedEdge.getLength_();
					edgeLengths[i] = currentEdgeLength;
					
					oldGraphLength += currentEdgeLength - movedEdge.length;
				}
				
				graphLength = oldGraphLength;
			}

			public void undoShimmying()
			{
				for (int i = 0; i < movedVertices.size(); i++)
				{
					SteinerTreeVertex movedVertex = movedVertices.get(i);
					movedVertex.x = movedVertex.xBackup;
					movedVertex.y = movedVertex.yBackup;
				}
			}

			public SteinerTreeVertex reinstallForGood()
			{
				// splice the new Steiner vertex back in
				for (Edge edge = newSteinerVertex.edge, flagEdge = null; edge != newSteinerVertex.edge || edge != flagEdge; edge = edge.leftNeighbor(), flagEdge = edge)
				{
					edge.rightSuccessor.opposite.leftSuccessor = edge.opposite;
					edge.leftSuccessor.opposite.rightSuccessor = edge.opposite;
					edge.destination.edge = edge.opposite;
				}
				
				// restore the vertex positions that were found after the splicing
				for (int i = 0; i < movedVertices.size(); i++)
				{
					SteinerTreeVertex movedVertex = movedVertices.get(i);
					movedVertex.x = movedVertex.xBackup = vertexLocationsX[i];
					movedVertex.y = movedVertex.yBackup = vertexLocationsY[i];
				}
				
				for (int i = 0; i < movedEdges.size(); i++)
				{
					SteinerTreeEdge movedEdge = movedEdges.get(i);
					movedEdge.length = movedEdge.opposite().length = edgeLengths[i];
				}

				// update any old closest original vertex values that are outside the split but still affected by it
				for (SteinerTreeEdge edge = newSteinerVertex.edge(), flagEdge = null; edge != newSteinerVertex.edge || edge != flagEdge; edge = edge.leftNeighbor_(), flagEdge = edge)
				{
					SteinerTreeEdge leftSuccessor = edge.leftSuccessor();
					SteinerTreeEdge rightSuccessor = edge.rightSuccessor();
					
					if (leftSuccessor != edge.opposite)
						leftSuccessor.updateOldClosestOriginalVertexValuesAheadOfSplicingZone(rightSuccessor);
					
					if (rightSuccessor != leftSuccessor) // this test is probably overkill, but whatever, it's cheap
						rightSuccessor.updateOldClosestOriginalVertexValuesAheadOfSplicingZone(leftSuccessor);
				}
				
				// greedily repair any newly created gap in a formerly contiguous closest-vertex region
				SteinerTreeEdge e0 = newOriginalVertex.edge().opposite();
				SteinerTreeEdge e1 = e0.leftNeighbor_();
				SteinerTreeEdge e2 = e1.leftNeighbor_();
				
				SteinerTreeVertex e1l, e1r, e2l, e2r;
				int e1ld, e1rd, e2ld, e2rd;
				
				if (e1.destination().isOriginal)
				{
					e1l = e1r = e1.destination();
					e1ld = e1rd = -1;
				}
				else
				{
					e1l = e1.leftSuccessor().closestOriginalVertex;
					e1ld = e1.leftSuccessor().closestOriginalVertexDistance;
					e1r = e1.rightSuccessor().closestOriginalVertex;
					e1rd = e1.rightSuccessor().closestOriginalVertexDistance;
				}
				
				if (e2.destination().isOriginal)
				{
					e2l = e2r = e2.destination();
					e2ld = e2rd = -1;
				}
				else
				{
					e2l = e2.leftSuccessor().closestOriginalVertex;
					e2ld = e2.leftSuccessor().closestOriginalVertexDistance;
					e2r = e2.rightSuccessor().closestOriginalVertex;
					e2rd = e2.rightSuccessor().closestOriginalVertexDistance;
				}
				
				SteinerTreeVertex e1m = null, e2m = null;
				int e1md = Integer.MIN_VALUE, e2md = Integer.MIN_VALUE;
				if (e1l == e2l)
				{
					e1m = e1l;
					e2m = e2l;
					e1md = e1ld;
					e2md = e2ld;
				}
				else if (e1l == e2r)
				{
					e1m = e1l;
					e2m = e2r;
					e1md = e1ld;
					e2md = e2rd;
				}
				else if (e1r == e2l)
				{
					e1m = e1r;
					e2m = e2l;
					e1md = e1rd;
					e2md = e2ld;
				}
				else if (e1r == e2r)
				{
					e1m = e1r;
					e2m = e2r;
					e1md = e1rd;
					e2md = e2rd;
				}
				
				if (e1m != null)
					if (e1md < e2md)
					{
						e2.setClosestOriginalVertex(e2m, e2md - 1, false);
						e1.setClosestOriginalVertex(e2m, e2md - 2, true);
					}
					else
					{
						e1.setClosestOriginalVertex(e1m, e1md - 1, false);
						e2.setClosestOriginalVertex(e1m, e1md - 2, true);
					}

				// finish healing the newly split edge by filling in any closest-original-vertex data that's still missing from the three edges that have replaced it
				e0.setClosestOriginalVertex(newOriginalVertex, 0, true);

				if (e1.closestOriginalVertex == null)
					if (e1.destination().isOriginal)
						e1.setClosestOriginalVertex(e1.destination(), 0, true);
					else
					{
						e1.setClosestOriginalVertex(newOriginalVertex, 1, false);
						e1.propagateNewClosestOriginalVertexRegionForward();
					}
				
				if (e2.closestOriginalVertex == null)
					if (e2.destination().isOriginal)
						e2.setClosestOriginalVertex(e2.destination(), 0, true);
					else
					{
						e2.setClosestOriginalVertex(newOriginalVertex, 1, false);
						e2.propagateNewClosestOriginalVertexRegionForward();
					}

				// finally, now that everything else is consistent, allow the nascent closest-vertex region surrounding the new vertex to grow into the tree
				e0.opposite().propagateNewClosestOriginalVertexRegionForward();
				
				// return
				return newSteinerVertex;
			}
		}
		
		TIntArrayList nearbyVertexIndices = new TIntArrayList();
		Queue<Vertex> nearbyVertexBuffer = new ArrayDeque<>();
		ArrayList<SteinerTreeEdge> potentialEdgesForSplicing = new ArrayList<>();
		ArrayList<SteinerTreeVertex> verticesNeedingMovement = new ArrayList<>();
		ArrayList<SteinerTreeVertex> verticesAboutToNeedMovement = new ArrayList<>();
		
		EdgeFrame savedEdgeFrame = new EdgeFrame();
		
		for (; v < verticesByDistanceFromCenter.size(); v++)
		{
			SteinerTreeVertex vertexBeingAdded = verticesByDistanceFromCenter.get(v);

			// find the Steiner tree edges that are near the vertex to be added
			Vertex triangulationVertex = delaunayTriangulation.vertices[vertexBeingAdded.index];
			
			nearbyVertexIndices.reset();
//			triangulationVertex.markVerticesWithinNSteps(searchDistance, markIndex, nearbyVertexIndices, nearbyVertexBuffer);
			triangulationVertex.markVerticesWithinNStepsIncreasingNAsNecessaryToSatisfyPredicateAtLeastOnce(searchDistance, markIndex, nearbyVertexIndices, nearbyVertexBuffer,
					new Predicate1<Vertex>()  // here the search for nearby tree vertices is artificially forced to continue until at least one actual tree edge has been
					{						 // encountered...  normally, the manner in which the vertices are ordered will guarantee that just one step out in the Delaunay
						@Override			// triangulation will reveal some portion of the growing Steiner tree; and in fact it seems to be difficult to produce a test case showing
						public boolean isTrue(Vertex vertex)  // that it's even possible for no part of the tree to be found in one step.  but it also seems difficult to PROVE that
						{									 // a tree edge will be encountered, and floating-point numbers can be unreliable when they are very close together...
							return vertices[vertex.index].edge() != null;  // so for safety's sake, the search is explicitly required to continue until an edge is found.
						}
					});

			potentialEdgesForSplicing.clear();
			for (int i = 0; i < nearbyVertexIndices.size(); i++)
				vertices[nearbyVertexIndices.getQuick(i)].reportEdgesNearMarkedTriangulationVertices(potentialEdgesForSplicing, delaunayTriangulation, markIndex);
			
			for (int i = 0; i < potentialEdgesForSplicing.size(); i++)
				potentialEdgesForSplicing.get(i).setMarks(markIndex, false);
			
			triangulationVertex.clearMarkFromVerticesInConnectedSubgraph(markIndex);
			
			// find the best old edge for the new vertex to be connected to (via one new Steiner vertex that'll split the old edge)--
			// that is, the one that makes for the shortest graph--and remember the result of that best addition, so that it can be
			// safely undone to make room for the next, potentially better addition to be tried, until all have been tried; then,
			// redo the best one
			SteinerVertexAdditionResult bestAdditionResult = null;
			
			for (int i = 0; i < potentialEdgesForSplicing.size(); i++)
			{
				// splice in the vertex being added
				SteinerTreeEdge splicedEdge = potentialEdgesForSplicing.get(i);
				savedEdgeFrame.setEdge(splicedEdge);
				
				SteinerTreeVertex newSteinerVertex = new SteinerTreeVertex(vertexCount);
				SteinerTreeVertex v0 = vertexBeingAdded;
				SteinerTreeVertex v1 = splicedEdge.destination();
				SteinerTreeVertex v2 = splicedEdge.opposite().destination();
				Edge edge0 = new SteinerTreeEdge(v0), edge1 = new SteinerTreeEdge(v1), edge2 = new SteinerTreeEdge(v2);
				edge0.opposite = new SteinerTreeEdge(newSteinerVertex);
				edge1.opposite = new SteinerTreeEdge(newSteinerVertex);
				edge2.opposite = new SteinerTreeEdge(newSteinerVertex);
				edge0.opposite.opposite = edge0;
				edge1.opposite.opposite = edge1;
				edge2.opposite.opposite = edge2;
				newSteinerVertex.edge = edge0;
				v0.edge = edge0.opposite;
				v1.edge = edge1.opposite;
				v2.edge = edge2.opposite;
				edge0.opposite.leftSuccessor = edge2;
				edge0.opposite.rightSuccessor = edge1;
				edge1.opposite.leftSuccessor = edge0;
				edge1.opposite.rightSuccessor = edge2;
				edge2.opposite.leftSuccessor = edge1;
				edge2.opposite.rightSuccessor = edge0;
				
				edge0.leftSuccessor = edge0.rightSuccessor = edge0.opposite;
				
				if (splicedEdge.rightSuccessor == splicedEdge.opposite)
					edge1.leftSuccessor = edge1.rightSuccessor = edge1.opposite;
				else
				{
					edge1.leftSuccessor = splicedEdge.leftSuccessor;
					edge1.rightSuccessor = splicedEdge.rightSuccessor;
					edge1.leftSuccessor.opposite.rightSuccessor = edge1.opposite;
					edge1.rightSuccessor.opposite.leftSuccessor = edge1.opposite;
				}
				
				if (splicedEdge.opposite.rightSuccessor == splicedEdge)
					edge2.leftSuccessor = edge2.rightSuccessor = edge2.opposite;
				else
				{
					edge2.leftSuccessor = splicedEdge.opposite.leftSuccessor;
					edge2.rightSuccessor = splicedEdge.opposite.rightSuccessor;
					edge2.leftSuccessor.opposite.rightSuccessor = edge2.opposite;
					edge2.rightSuccessor.opposite.leftSuccessor = edge2.opposite;
				}
				
				// shift Steiner vertices around until the edges are (probably) about as short as they're gonna get
				double newGraphLength = graphLength - splicedEdge.length; // the new edges know their length as 0, so the splicing reduces the graph length by the lost edge's length and increases it by 0--which is technically wrong, but consistent, and as long as it's consistent, the length update system will correct it correctly
				verticesNeedingMovement.add(newSteinerVertex);
				
				for (int shimmyingIterationCount = 0, maximalVertexMovementCount = 0; !verticesNeedingMovement.isEmpty() && shimmyingIterationCount < hardShimmyingIterationLimit; shimmyingIterationCount++)
				{
					float insignificantMovementThresholdSquared = insignificantMovementThreshold * insignificantMovementThreshold;

					for (int vnm = 0; vnm < verticesNeedingMovement.size(); vnm++)
					{
						SteinerTreeVertex vertexToBeMoved = verticesNeedingMovement.get(vnm);
						
						if (vertexToBeMoved.movementCount == maximalVertexMovementCount && maximalVertexMovementCount > 0)
							insignificantMovementThreshold *= insignificantMovementThresholdGrowthFactor; // a new wave of movement seems to be starting; try a little less hard to get it right than for the last one
						
						vertexToBeMoved.takeLocallyCorrectPosition(verticesAboutToNeedMovement, insignificantMovementThresholdSquared, markIndex);
						
						if (vertexToBeMoved.movementCount > maximalVertexMovementCount)
							maximalVertexMovementCount = vertexToBeMoved.movementCount;
					}
					
					for (int vanm = 0; vanm < verticesAboutToNeedMovement.size(); vanm++)
						verticesAboutToNeedMovement.get(vanm).clearMark(markIndex);
					
					verticesNeedingMovement.clear();
					
					ArrayList<SteinerTreeVertex> buffer = verticesNeedingMovement;
					verticesNeedingMovement = verticesAboutToNeedMovement;
					verticesAboutToNeedMovement = buffer;
				}
				
				verticesNeedingMovement.clear(); // if the hard shimmying iteration limit was hit then this list might still have entries
				
				// record and the result of the vertex addition and test its quality
				SteinerVertexAdditionResult result = new SteinerVertexAdditionResult(vertexBeingAdded, newSteinerVertex, splicedEdge, newGraphLength);
				
				if (bestAdditionResult == null || bestAdditionResult.graphLength > result.graphLength)
					bestAdditionResult = result;
				
				// undo the addition to make room for the next one
				result.undoShimmying();
				savedEdgeFrame.restore();
				vertexBeingAdded.edge = null;
			}
			
			// reconnect the vertex in the best way that was found, and update running position and length values to reflect that this
			// is the new canonical configuration of the graph
			assert closestVerticesAreRight(vertices, vertexCount);
			SteinerTreeVertex newSteinerVertex = bestAdditionResult.reinstallForGood();
			vertices[vertexCount++] = newSteinerVertex;
			assert closestVerticesAreRight(vertices, vertexCount);  // keeping these closest-vertex records correct was especially tricky!
			
			graphLength = bestAdditionResult.graphLength; // this graph length value has been through a lot of additions and
					// subtractions, so maybe it's not as accurate as it should be...  recalculating it from scratch would presumably
					// help with that, but doing so every time would guarantee a O(n^2) running time...  it could be done, say, once
					// every n/C times for some constant C; that'd change the bound to O(Cn) = O(n); and there might be other
					// solutions...  but for now it's probably good enough to trust in double to hold up under this kind of abuse
		}
		
		// the graph is basically done, but some Steiner points may have ended up atop original vertices; these must be removed--and
		// also, edges probably aren't ordered correctly, so sort them
		int degenerateVertexCount = 0;
		
		for (int i = 0; i < vertexCount; i++)
			if (vertices[i] != null)
				degenerateVertexCount += vertices[i].removeDegenerateSteinerNeighborsAndSortEdges(vertices, new ArrayList<SteinerTreeEdge>(), markIndex);
		
		// make the final graph
		int finalVertexCount = vertexCount - degenerateVertexCount;
		
		DCELGraph graph = new DCELGraph(finalVertexCount);
		
		for (int i = 0, fi = 0; i < vertexCount; i++)
		{
			SteinerTreeVertex vertex = vertices[i];
		
			if (vertex != null)
			{
				graph.vertices[fi] = vertex;
				vertex.index = fi;
				++fi;
			}
		}
		
		return graph;
	}
	
	private static boolean closestVerticesAreRight(SteinerTreeVertex[] vs, int vc)
	{
		for (int i = 0; i < vc; ++i)
			for (SteinerTreeEdge edge0 = vs[i].edge(), edge = edge0, flagEdge = null; edge != edge0 || edge != flagEdge; edge = edge.leftNeighbor_(), flagEdge = edge)
				if (!closestVertexDataIsLocallyConsistent(edge))
					return false;
		
		return true;
	}
	
	private static boolean closestVertexDataIsLocallyConsistent(SteinerTreeEdge edge)
	{
		SteinerTreeEdge forwardEdge = edge.closestOriginalVertexDirection ? edge : edge.opposite();
		SteinerTreeEdge backwardEdge = forwardEdge.opposite();
		SteinerTreeVertex cv = forwardEdge.closestOriginalVertex;
		int distance = forwardEdge.closestOriginalVertexDistance;
		
		if (backwardEdge.closestOriginalVertex != cv || backwardEdge.closestOriginalVertexDirection || backwardEdge.closestOriginalVertexDistance != distance ||
				!forwardEdge.closestOriginalVertexDirection)
			return false;
		
		if (backwardEdge.rightSuccessor().closestOriginalVertex == cv && backwardEdge.rightSuccessor().closestOriginalVertexDistance != distance + 1)
			return false;
		
		if (backwardEdge.rightSuccessor().closestOriginalVertex != cv && backwardEdge.rightSuccessor().closestOriginalVertexDistance > distance + 1)
			return false;
		
		if (backwardEdge.rightSuccessor().closestOriginalVertex == cv && backwardEdge.rightSuccessor().closestOriginalVertexDirection)
			return false;
		
		if (backwardEdge.leftSuccessor().closestOriginalVertex == cv && backwardEdge.leftSuccessor().closestOriginalVertexDistance != distance + 1)
			return false;
		
		if (backwardEdge.leftSuccessor().closestOriginalVertex != cv && backwardEdge.leftSuccessor().closestOriginalVertexDistance > distance + 1)
			return false;
		
		if (backwardEdge.leftSuccessor().closestOriginalVertex == cv && backwardEdge.leftSuccessor().closestOriginalVertexDirection)
			return false;
		
		if (forwardEdge.destination == cv)
		{
			if (distance != 0)
				return false;
		}
		else
		{
			if (forwardEdge.leftSuccessor().closestOriginalVertex != cv && forwardEdge.rightSuccessor().closestOriginalVertex != cv)
				return false;
	
			if (forwardEdge.leftSuccessor().closestOriginalVertex == cv && forwardEdge.rightSuccessor().closestOriginalVertex == cv)
			{
				SteinerTreeEdge ls = forwardEdge.leftSuccessor();
				SteinerTreeEdge rs = forwardEdge.rightSuccessor();
				
				if (!(	ls.closestOriginalVertexDistance == distance && rs.closestOriginalVertexDistance == distance - 1 ||
						ls.closestOriginalVertexDistance == distance - 1 && rs.closestOriginalVertexDistance == distance))
					return false;
				
				SteinerTreeEdge fs, bs;
				if (ls.closestOriginalVertexDistance == distance - 1)
				{
					fs = ls;
					bs = rs;
				}
				else
				{
					fs = rs;
					bs = ls;
				}
				
				if (!fs.closestOriginalVertexDirection || bs.closestOriginalVertexDirection)
					return false;
			}
			else
			{
				if (forwardEdge.leftSuccessor().closestOriginalVertex != cv && forwardEdge.leftSuccessor().closestOriginalVertexDistance > distance)
					return false;
				
				if (forwardEdge.leftSuccessor().closestOriginalVertex == cv && forwardEdge.leftSuccessor().closestOriginalVertexDistance != distance - 1)
					return false;
				
				if (!forwardEdge.leftSuccessor().closestOriginalVertexDirection)
					return false;
	
				if (forwardEdge.rightSuccessor().closestOriginalVertex != cv && forwardEdge.rightSuccessor().closestOriginalVertexDistance > distance)
					return false;
				
				if (forwardEdge.rightSuccessor().closestOriginalVertex == cv && forwardEdge.rightSuccessor().closestOriginalVertexDistance != distance - 1)
					return false;
				
				if (!forwardEdge.rightSuccessor().closestOriginalVertexDirection)
					return false;
			}
		}
		
		return true;
	}
	
	private static class OneWayEdgeFrame
	{
		private Vertex destination;
		private Edge destinationEdge;
		private Edge rightSuccessorOpposite;
		private Edge rightSuccessorOppositeLeftSuccessor;
		private Edge leftSuccessorOpposite;
		private Edge leftSuccessorOppositeRightSuccessor;

		public void setEdge(Edge edge)
		{
			destination = edge.destination;
			destinationEdge = edge.destination.edge;
			rightSuccessorOpposite = edge.rightSuccessor.opposite;
			rightSuccessorOppositeLeftSuccessor = rightSuccessorOpposite.leftSuccessor;
			leftSuccessorOpposite = edge.leftSuccessor.opposite;
			leftSuccessorOppositeRightSuccessor = leftSuccessorOpposite.rightSuccessor;
		}

		public void restore()
		{
			destination.edge = destinationEdge;
			rightSuccessorOpposite.leftSuccessor = rightSuccessorOppositeLeftSuccessor;
			leftSuccessorOpposite.rightSuccessor = leftSuccessorOppositeRightSuccessor;
		}
	}
	
	private static class EdgeFrame
	{
		private final OneWayEdgeFrame forwardFrame;
		private final OneWayEdgeFrame backwardFrame;
		
		public EdgeFrame()
		{
			forwardFrame = new OneWayEdgeFrame();
			backwardFrame = new OneWayEdgeFrame();
		}
		
		public void setEdge(Edge edge)
		{
			forwardFrame.setEdge(edge);
			backwardFrame.setEdge(edge.opposite);
		}

		public void restore()
		{
			forwardFrame.restore();
			backwardFrame.restore();
		}
	}

	public static DCELGraph makeDelaunayTriangulation(List<PlaneVector> vertices)
	{
		return makeDelaunayTriangulation(vertices, DelaunayDuplicatePointPolicy.makeOneVertexPerInputLocationButInvolveOnlyTheFirstVertexPerLocationInTheTriangulation);
	}
	
	public static DCELGraph makeDelaunayTriangulation(List<PlaneVector> vertices, DelaunayDuplicatePointPolicy delaunayDuplicatePointPolicy)
	{
		return makeDelaunayTriangulation(vertices, delaunayDuplicatePointPolicy, null, null, null);
	}
	
	private static DCELGraph makeDelaunayTriangulation(List<PlaneVector> vertices, DelaunayDuplicatePointPolicy delaunayDuplicatePointPolicy, int[] mapFromOldIndicesToNewIndices,
			Action1<Vertex> sortedPointRegistration, Action1<Vertex[]> rawSortedVerticesRegistration)
	{
		DCELGraph graph = new DCELGraph(vertices, Dummy.instance);
		graph.addDelaunayTriangulationEdges(delaunayDuplicatePointPolicy, mapFromOldIndicesToNewIndices, sortedPointRegistration, rawSortedVerticesRegistration);
		
		return graph;
	}
	
	public enum DelaunayDuplicatePointPolicy
	{
		makeOneVertexPerInputLocationButInvolveOnlyTheFirstVertexPerLocationInTheTriangulation,
		expectAndAssumeAllInputLocationsAreDistinct
	}
	
	public void addDelaunayTriangulationEdges()
	{
		addDelaunayTriangulationEdges(DelaunayDuplicatePointPolicy.makeOneVertexPerInputLocationButInvolveOnlyTheFirstVertexPerLocationInTheTriangulation);
	}

	public void addDelaunayTriangulationEdges(DelaunayDuplicatePointPolicy delaunayDuplicatePointPolicy)
	{
		for (Vertex vertex: vertices)
			if (vertex.edge != null)
				throw new UnsupportedOperationException("this probably won't work?  maybe it will.  probably not though");
		
		addDelaunayTriangulationEdges(delaunayDuplicatePointPolicy, null, null, null);
	}
	
	private void addDelaunayTriangulationEdges(DelaunayDuplicatePointPolicy delaunayDuplicatePointPolicy, int[] mapFromOldIndicesToNewIndices, Action1<Vertex> sortedVertexRegistration, Action1<Vertex[]> rawSortedVerticesRegistration)
	{
		if (vertices.length >= 2 || sortedVertexRegistration != null)
		{
			Vertex[] sortedVertices = Arrays.copyOf(vertices, vertices.length);
			int sortedVertexCount = sortedVertices.length;
			Arrays.sort(sortedVertices, horizontalComparator);
			
			if (rawSortedVerticesRegistration != null)
				rawSortedVerticesRegistration.perform(sortedVertices);
			
			if (delaunayDuplicatePointPolicy == DelaunayDuplicatePointPolicy.makeOneVertexPerInputLocationButInvolveOnlyTheFirstVertexPerLocationInTheTriangulation)
				for (int oldV = 0, newV = 1; newV < sortedVertices.length;)
				{
					while (newV < sortedVertices.length && sortedVertices[oldV].sharesLocationWith(sortedVertices[newV]))
					{
						newV++;
						sortedVertexCount--;
					}
					
					if (newV < sortedVertices.length)
					{
						int ov = ++oldV;
						int nv = newV++;
						sortedVertices[ov] = sortedVertices[nv];
						if (mapFromOldIndicesToNewIndices != null) mapFromOldIndicesToNewIndices[ov] = nv;
					}
				}
			
			for (int v = 0; v < sortedVertexCount; v++) sortedVertices[v].xIndex = v;
			Arrays.sort(sortedVertices, 0, sortedVertexCount, verticalComparator);
			for (int v = 0; v < sortedVertexCount; v++) sortedVertices[v].yIndex = v;
			
			if (sortedVertexRegistration != null)
				for (int v = 0; v < sortedVertexCount; v++)
					sortedVertexRegistration.perform(sortedVertices[v]);
			
			if (sortedVertexCount >= 2)
				addDelaunayTriangulationEdges(sortedVertices, 0, sortedVertexCount, false, new Vertex[sortedVertexCount]);
		}
	}

	private BoundingEdgeArray addDelaunayTriangulationEdges(Vertex[] sortedVertices, int rangeEnd0, int rangeEnd1, boolean sortedByXNotY, Vertex[] temporaryStorage)
	{ // this algorithm comes from "Primitives for the Manipulation of General Subdivisions and the Computation of Voronoi Diagrams" by Leonidas Guibas and Jorge Stolfi
		int vertexCount = rangeEnd1 - rangeEnd0;
		
		switch (vertexCount)
		{
		case 0:
		case 1:
			throw new IllegalArgumentException("" + vertexCount);
		case 2:
		{
			Vertex v0 = sortedVertices[rangeEnd0 + 0];
			Vertex v1 = sortedVertices[rangeEnd0 + 1];
			Edge edge0 = new Edge(v1);
			Edge edge1 = new Edge(sortedVertices[rangeEnd0 + 0]);
			v0.edge = edge0;
			v1.edge = edge1;
			edge0.opposite = edge0.leftSuccessor = edge0.rightSuccessor = edge1;
			edge1.opposite = edge1.leftSuccessor = edge1.rightSuccessor = edge0;
			
			BoundingEdgeArray boundingEdges = new BoundingEdgeArray();
			
			if (edge0.destination.xIndex < edge1.destination.xIndex)
			{
				boundingEdges.ccwHullEdgeWithLeftmostDestination = edge0;
				boundingEdges.cwHullEdgeWithRightmostDestination = edge1;
			}
			else
			{
				boundingEdges.ccwHullEdgeWithLeftmostDestination = edge1;
				boundingEdges.cwHullEdgeWithRightmostDestination = edge0;
			}
			
			if (edge0.destination.yIndex < edge1.destination.yIndex)
			{
				boundingEdges.cwHullEdgeWithBottommostDestination = edge0;
				boundingEdges.ccwHullEdgeWithTopmostDestination = edge1;
			}
			else
			{
				boundingEdges.cwHullEdgeWithBottommostDestination = edge1;
				boundingEdges.ccwHullEdgeWithTopmostDestination = edge0;
			}
			
			return boundingEdges;
		}
		case 3:
		{
			Vertex v0 = sortedVertices[rangeEnd0 + 0];
			Vertex v1 = sortedVertices[rangeEnd0 + 1];
			Vertex v2 = sortedVertices[rangeEnd0 + 2];
			
			if (PlaneVector.crossZ(v0, v1, v2) < 0)
			{
				Vertex v_ = v0;
				v0 = v1;
				v1 = v_;
			} // now the vertices are in CCW order
			
			Edge edge01 = new Edge(v1);
			Edge edge10 = new Edge(v0);
			Edge edge02 = new Edge(v2);
			Edge edge20 = new Edge(v0);
			Edge edge12 = new Edge(v2);
			Edge edge21 = new Edge(v1);
			
			v0.edge = edge01;
			v1.edge = edge12;
			v2.edge = edge20;

			edge01.leftSuccessor = edge01.rightSuccessor = edge12;
			edge12.leftSuccessor = edge12.rightSuccessor = edge20;
			edge20.leftSuccessor = edge20.rightSuccessor = edge01;
			
			edge02.leftSuccessor = edge02.rightSuccessor = edge21;
			edge21.leftSuccessor = edge21.rightSuccessor = edge10;
			edge10.leftSuccessor = edge10.rightSuccessor = edge02;
			
			edge01.opposite = edge10;
			edge10.opposite = edge01;
			edge12.opposite = edge21;
			edge21.opposite = edge12;
			edge20.opposite = edge02;
			edge02.opposite = edge20;
			
			BoundingEdgeArray boundingEdges = new BoundingEdgeArray();
			
			for (int i = rangeEnd0; i < rangeEnd1; i++)
			{
				Vertex vertex = sortedVertices[i];
				
				switch (vertex.xIndex - rangeEnd0)
				{
				case 0: boundingEdges.ccwHullEdgeWithLeftmostDestination = vertex.edge.rightPredecessor(); break;
				case 2: boundingEdges.cwHullEdgeWithRightmostDestination = vertex.edge.opposite; break;
				}
				
				switch (vertex.yIndex - rangeEnd0)
				{
				case 0: boundingEdges.cwHullEdgeWithBottommostDestination = vertex.edge.opposite; break;
				case 2: boundingEdges.ccwHullEdgeWithTopmostDestination = vertex.edge.rightPredecessor(); break;
				}
			}

			return boundingEdges;
		}
		default:
			// build a triangulation graph in each half of the vertex list
			int vertexCount0 = vertexCount / 2;
			int rangeMiddle = rangeEnd0 + vertexCount0;
			
			resort(sortedVertices, rangeEnd0, rangeMiddle, rangeEnd0, rangeEnd1, sortedByXNotY, temporaryStorage); // if this is parallelized, it could be safely
			resort(sortedVertices, rangeMiddle, rangeEnd1, rangeEnd0, rangeEnd1, sortedByXNotY, temporaryStorage); // run in three steps:  1) these two calls; then
			
			BoundingEdgeArray bounds0 = addDelaunayTriangulationEdges(sortedVertices, rangeEnd0, rangeMiddle, !sortedByXNotY, temporaryStorage); // 2) this call; then

			BoundingEdgeArray bounds1 = addDelaunayTriangulationEdges(sortedVertices, rangeMiddle, rangeEnd1, !sortedByXNotY, temporaryStorage); // 3) this one.
			
			// walk an imaginary segment stretching between the "innermost" points of the two new graphs toward the end of the gap
			// between the graphs so that it becomes one of the two missing sections of the convex hull of the combined graph
			Edge innerEdge0, innerEdge1;
			
			if (sortedByXNotY)
			{
				innerEdge0 = bounds0.cwHullEdgeWithRightmostDestination;
				innerEdge1 = bounds1.ccwHullEdgeWithLeftmostDestination;
			}
			else
			{
				innerEdge0 = bounds1.cwHullEdgeWithBottommostDestination;
				innerEdge1 = bounds0.ccwHullEdgeWithTopmostDestination;
			}
			
			for (boolean edge0Changed = true, edge1Changed = true; edge1Changed;)
			{
				{
					while (PlaneVector.crossZ(innerEdge1.destination, innerEdge0.destination, innerEdge0.leftSuccessor.destination) > 0)
					{
						innerEdge0 = innerEdge0.leftSuccessor;
						edge0Changed = true;
					}
					
					edge1Changed = false;
				}
				
				if (edge0Changed)
				{
					while (PlaneVector.crossZ(innerEdge0.destination, innerEdge1.destination, innerEdge1.rightSuccessor.destination) < 0)
					{
						innerEdge1 = innerEdge1.rightSuccessor;
						edge1Changed = true;
					}
				
					edge0Changed = false;
				}
			}
			
			// add the missing convex hull section to the graph
			Edge newEdge0 = new Edge(innerEdge0.destination);
			Edge newEdge1 = new Edge(innerEdge1.destination);
			
			newEdge0.leftSuccessor = innerEdge0.leftSuccessor;
			newEdge0.rightSuccessor = innerEdge0.opposite;
			newEdge1.leftSuccessor = innerEdge1.opposite;
			newEdge1.rightSuccessor = innerEdge1.rightSuccessor;
			newEdge0.opposite = newEdge1;
			newEdge1.opposite = newEdge0;
			
			innerEdge0.leftSuccessor = newEdge1;
			newEdge0.leftSuccessor.opposite.rightSuccessor = newEdge1;
			innerEdge1.rightSuccessor = newEdge0;
			newEdge1.rightSuccessor.opposite.leftSuccessor = newEdge0;
			
			Edge newHullSection01 = newEdge0;
			
			// zip up the rest of the gap between the graphs
			Edge newHullSection10;
			
			for (newHullSection10 = null; newHullSection10 == null;)
			{
				Edge gapWall0 = newEdge0.rightSuccessor;
				Edge gapWall1 = newEdge1.leftSuccessor;
				float gapTriangle0Cross = PlaneVector.crossZ(newEdge1.destination, gapWall0.destination, newEdge0.destination);
				float gapTriangle1Cross = PlaneVector.crossZ(newEdge1.destination, gapWall1.destination, newEdge0.destination);
				
				if (gapTriangle0Cross <= 0) gapWall0 = null; // the gap actually ended, it turns out
				if (gapTriangle1Cross <= 0) gapWall1 = null;
				
				if (gapWall0 != null)
					while (gapWall0.leftSuccessor.leftSuccessor.leftSuccessor == gapWall0 && Circles.checkContainmentAndOrder(newEdge0.destination, newEdge1.destination, gapWall0.destination, gapWall0.leftNeighbor().destination))
					{
						Edge oldWall = gapWall0;
						gapWall0 = gapWall0.leftNeighbor();
						remove(oldWall);
					}
				
				if (gapWall1 != null)
					while (gapWall1.rightSuccessor.rightSuccessor.rightSuccessor == gapWall1 && Circles.checkContainmentAndOrder(newEdge0.destination, newEdge1.destination, gapWall1.destination, gapWall1.rightNeighbor().destination))
					{
						Edge oldWall = gapWall1;
						gapWall1 = gapWall1.rightNeighbor();
						remove(oldWall);
					}
				
				boolean nextEdgeIsOnWall0 = false;
				
				if (gapWall0 != null)
					if (gapWall1 != null)
						if (Circles.checkContainmentAndOrder(newEdge0.destination, newEdge1.destination, gapWall1.destination, gapWall0.destination))
							nextEdgeIsOnWall0 = true;
						else
							nextEdgeIsOnWall0 = false;
					else
						nextEdgeIsOnWall0 = true;
				else
					if (gapWall1 != null)
						nextEdgeIsOnWall0 = false;
					else
						newHullSection10 = newEdge0;
				
				if (newHullSection10 == null)
					if (nextEdgeIsOnWall0)
					{
						Edge lastNewEdge0 = newEdge0;
						Edge lastNewEdge1 = newEdge1;
						newEdge0 = new Edge(gapWall0.destination);
						newEdge1 = new Edge(lastNewEdge1.destination);
						newEdge0.opposite = newEdge1;
						newEdge1.opposite = newEdge0;
						
						newEdge0.leftSuccessor = gapWall0.opposite;
						newEdge0.rightSuccessor = gapWall0.rightSuccessor;
						newEdge1.rightSuccessor = lastNewEdge0;
						newEdge1.leftSuccessor = lastNewEdge1.leftSuccessor;

						gapWall0.rightSuccessor.opposite.leftSuccessor = newEdge1;
						gapWall0.rightSuccessor = newEdge1;
						lastNewEdge1.leftSuccessor.opposite.rightSuccessor = newEdge0;
						lastNewEdge1.leftSuccessor = newEdge0;
					}
					else
					{
						Edge lastNewEdge0 = newEdge0;
						Edge lastNewEdge1 = newEdge1;
						newEdge0 = new Edge(lastNewEdge0.destination);
						newEdge1 = new Edge(gapWall1.destination);
						newEdge0.opposite = newEdge1;
						newEdge1.opposite = newEdge0;
						
						newEdge0.leftSuccessor = lastNewEdge1;
						newEdge0.rightSuccessor = lastNewEdge0.rightSuccessor;
						newEdge1.leftSuccessor = gapWall1.leftSuccessor;
						newEdge1.rightSuccessor = gapWall1.opposite;
						
						gapWall1.leftSuccessor.opposite.rightSuccessor = newEdge0;
						gapWall1.leftSuccessor = newEdge0;
						lastNewEdge0.rightSuccessor.opposite.leftSuccessor = newEdge1;
						lastNewEdge0.rightSuccessor = newEdge1;
					}
			}
			
			// find the two bounds of the new graph that are at either end of the former gap
			Comparator<Vertex> minorComparator;
			Edge boundCandidate01, boundCandidate10;
			
			if (sortedByXNotY)
			{
				minorComparator = verticalComparator;
				boundCandidate01 = newHullSection01;
				boundCandidate10 = newHullSection10;
			}
			else
			{
				minorComparator = horizontalComparator;
				boundCandidate01 = newHullSection10.opposite;
				boundCandidate10 = newHullSection01.opposite;
			}
			
			while (minorComparator.compare(boundCandidate01.leftSuccessor.destination, boundCandidate01.destination) < 0)
				boundCandidate01 = boundCandidate01.leftSuccessor;
			
			while (minorComparator.compare(boundCandidate01.leftPredecessor().destination, boundCandidate01.destination) < 0)
				boundCandidate01 = boundCandidate01.leftPredecessor();
			
			while (minorComparator.compare(boundCandidate10.rightSuccessor.destination, boundCandidate10.destination) > 0)
				boundCandidate10 = boundCandidate10.rightSuccessor;
			
			while (minorComparator.compare(boundCandidate10.rightPredecessor().destination, boundCandidate10.destination) > 0)
				boundCandidate10 = boundCandidate10.rightPredecessor();
			
			// set all four of the graph's bounds, refusing to fully trust old bounds values (bounding edges' destinations are still
			// bounding vertices, but the actual edges may have been removed from the graph; and those that weren't removed may have
			// been moved to the interior of the graph by new convex hull edges that were added)
			BoundingEdgeArray boundingEdges = new BoundingEdgeArray();
			
			if (sortedByXNotY)
			{
				Edge ccwL = bounds0.ccwHullEdgeWithLeftmostDestination.defaultingToOtherEdgeWithSameDestination();
				Edge cwR = bounds1.cwHullEdgeWithRightmostDestination.defaultingToOtherEdgeWithSameDestination();
				boundingEdges.ccwHullEdgeWithLeftmostDestination = ccwL;
				boundingEdges.cwHullEdgeWithRightmostDestination = cwR;
				boundingEdges.cwHullEdgeWithBottommostDestination = boundCandidate01;
				boundingEdges.ccwHullEdgeWithTopmostDestination = boundCandidate10;
				
				while (ccwL.rightSuccessor.rightSuccessor.rightSuccessor == ccwL && !ccwL.neighboringTriangleContainsGraph(false, boundingEdges)) ccwL = ccwL.rightSuccessor.opposite;
				while (cwR.leftSuccessor.leftSuccessor.leftSuccessor == cwR && !cwR.neighboringTriangleContainsGraph(true, boundingEdges)) cwR = cwR.leftSuccessor.opposite;
				boundingEdges.ccwHullEdgeWithLeftmostDestination = ccwL;
				boundingEdges.cwHullEdgeWithRightmostDestination = cwR;
			}
			else
			{
				Edge cwB = bounds0.cwHullEdgeWithBottommostDestination.defaultingToOtherEdgeWithSameDestination();
				Edge ccwT = bounds1.ccwHullEdgeWithTopmostDestination.defaultingToOtherEdgeWithSameDestination();
				boundingEdges.cwHullEdgeWithBottommostDestination = cwB;
				boundingEdges.ccwHullEdgeWithTopmostDestination = ccwT;
				boundingEdges.ccwHullEdgeWithLeftmostDestination = boundCandidate01.leftSuccessor.opposite;
				boundingEdges.cwHullEdgeWithRightmostDestination = boundCandidate10.rightSuccessor.opposite;
				
				while (cwB.leftSuccessor.leftSuccessor.leftSuccessor == cwB && !cwB.neighboringTriangleContainsGraph(true, boundingEdges)) cwB = cwB.leftSuccessor.opposite;
				while (ccwT.rightSuccessor.rightSuccessor.rightSuccessor == ccwT && !ccwT.neighboringTriangleContainsGraph(false, boundingEdges)) ccwT = ccwT.rightSuccessor.opposite;
				boundingEdges.cwHullEdgeWithBottommostDestination = cwB;
				boundingEdges.ccwHullEdgeWithTopmostDestination = ccwT;
			}
			
			return boundingEdges;
		}
	}
	
	private void remove(Edge edge)
	{
		if (edge.destination.edge == edge.opposite) edge.destination.edge = edge.opposite.distinctNeighbor();
		if (edge.opposite.destination.edge == edge) edge.opposite.destination.edge = edge.distinctNeighbor();
		
		Edge oldRightPredecessor = edge.rightPredecessor();
		Edge oldLeftPredecessor = edge.leftPredecessor();
		
		oldRightPredecessor.rightSuccessor = oldLeftPredecessor.opposite;
		oldLeftPredecessor.leftSuccessor = oldRightPredecessor.opposite;
		edge.rightSuccessor.opposite.leftSuccessor = edge.leftSuccessor;
		edge.leftSuccessor.opposite.rightSuccessor = edge.rightSuccessor;
		
		edge.isInGraph = false;
		edge.opposite.isInGraph = false;
	}

	private void resort(Vertex[] sortedVertices, int rangeEnd0, int rangeEnd1, int lastRangeEnd0, int lastRangeEnd1, boolean sortedByXNotY, Vertex[] temporaryStorage)
	{ // this changes the vertices in the current range, in linear time, from x-order to y-order or vice-versa, and renumbers them as
		for (int v = rangeEnd0; v < rangeEnd1; v++) // needed to keep both x-indices and y-indices within that same range
		{
			Vertex vertex = sortedVertices[v];
			temporaryStorage[sortedByXNotY ? vertex.yIndex : vertex.xIndex] = vertex;
		}
		
		for (int v = rangeEnd0, i = lastRangeEnd0; i < lastRangeEnd1; i++)
		{
			Vertex vertex = temporaryStorage[i];
		
			if (vertex != null)
			{
				int newIndex = v++;
				if (sortedByXNotY) vertex.yIndex = newIndex; else vertex.xIndex = newIndex;
				sortedVertices[newIndex] = vertex;
				temporaryStorage[i] = null;
			}
		}
	}

	public void markMinimalSpanningTree_Kruskal(int markIndex)
	{
		ArrayList<Edge> edges = gatherAllEdgesAsHalves();
		
		Collections.sort(edges, Edge.comparatorByLength);
		
		DisjointSetRegistry sets = new DisjointSetRegistry(vertices.length);
		
		for (Edge edge: edges)
		{
			int i1 = edge.opposite.destination.index;
			int i2 = edge.destination.index;
			
			if (sets.find(i1) != sets.find(i2))
			{
				sets.union(i1, i2);

				edge.setMark(markIndex);
				edge.opposite.destination.setMark(markIndex);
				edge.destination.setMark(markIndex);
			}
		}
	}
	
	public void clearAllMarks()
	{
		for (Vertex vertex: vertices)
		{
			vertex.clearAllMarks();
			
			for (Edge edge0 = vertex.edge, edge = edge0, flagEdge = null; edge != edge0 || edge != flagEdge; edge = edge.leftNeighbor(), flagEdge = edge)
				edge.clearAllMarks();
		}
	}

	public ArrayList<Edge> gatherAllEdgesAsHalves()
	{
		ArrayList<Edge> edges = new ArrayList<>();
		
		for (Vertex vertex: vertices)
			for (Edge edge0 = vertex.edge, edge = edge0, flagEdge = null; edge != edge0 || edge != flagEdge; edge = edge.leftNeighbor(), flagEdge = edge)
				if (vertex.index < edge.destination.index)
					edges.add(edge);
		
		return edges;
	}

	@Override
	public ShortestPathsTable findShortestPathsForAllVertexPairs(InterruptionSignal interruptionSignal)
	{
		ShortestPathsTable table = new ShortestPathsTable(vertices.length);
		
		for (Vertex v: vertices)
			findShortestPathsFromVertex_Dijkstra(v, table);
		
		return table;
	}

	private void findShortestPathsFromVertex_Dijkstra(Vertex root, ShortestPathsTable table)
	{
		PairingHeap<Double, Vertex> heap = new PairingHeap<>();
		ArrayList<Handle<Double, Vertex>> entries = new ArrayList<>();
		ArrayList<Vertex> precedingVertices = new ArrayList<>();
		
		for (int i = 0; i < vertices.length; ++i) entries.add(null);
		
		for (Vertex v: vertices)
		{
			entries.set(v.index, heap.insert(v == root ? 0 : Double.POSITIVE_INFINITY, v));
			precedingVertices.add(null);
		}

		for (int i = 0; !heap.isEmpty(); ++i)
		{
			Handle<Double, Vertex> entry = heap.deleteMin();
			Vertex v = entry.getValue();
			double distance = entry.getKey();
			
			table.setCellValue(root, v, (float) distance, precedingVertices.get(v.index), i);
			
			for (Edge e = v.edge, flag = null; e != flag; e = e.leftNeighbor(), flag = v.edge)
			{
				double distanceToDestination = distance + e.getLength();
				Handle<Double, Vertex> destinationEntry = entries.get(e.destination.index);
				
				if (distanceToDestination < destinationEntry.getKey())
				{
					destinationEntry.decreaseKey(distanceToDestination);
					precedingVertices.set(e.destination.index, v);
				}
			}
		}
	}

	private static class BoundingEdgeArray
	{
		public Edge cwHullEdgeWithRightmostDestination;
		public Edge ccwHullEdgeWithTopmostDestination;
		public Edge ccwHullEdgeWithLeftmostDestination;
		public Edge cwHullEdgeWithBottommostDestination;
	}

	@Override
	public float[] getEdgesAsLineSegments_xyxy(Boolean requiredMark, int markIndex)
	{
		TFloatArrayList result = new TFloatArrayList();
		
		for (int i = 0; i < vertices.length; i++)
		{
			Vertex vertex = vertices[i];
		
			for (Edge edge0 = vertex.edge, edge = edge0, flagEdge = null; edge != edge0 || edge != flagEdge; edge = edge.leftNeighbor(), flagEdge = edge)
				if (i < edge.destination.index && (requiredMark == null || edge.readMark(markIndex) == requiredMark))
				{
					result.add(vertex.x);
					result.add(vertex.y);
					result.add(edge.destination.x);
					result.add(edge.destination.y);
				}
		}
		
		return result.toArray();
	}

	@Override
	public float[] getVerticesAsExesMadeFromLineSegments_xyxy(float exHalfWidth, Boolean requiredMark, int markIndex)
	{
		TFloatArrayList result = new TFloatArrayList();
		for (int i = 0; i < vertices.length; ++i)
			if (requiredMark == null || vertices[i].readMark(markIndex) == requiredMark)
			{
				float x = vertices[i].x;
				float y = vertices[i].y;
				result.add(x - exHalfWidth);
				result.add(y - exHalfWidth);
				result.add(x + exHalfWidth);
				result.add(y + exHalfWidth);
				result.add(x - exHalfWidth);
				result.add(y + exHalfWidth);
				result.add(x + exHalfWidth);
				result.add(y - exHalfWidth);
			}
		
		return result.toArray();
	}

	public PlaneVector locateVertex(int index)
	{
		return vertices[index].getLocation();
	}
	
//	@Override
//	public void markSteinerSubgraph(TIntArrayList terminals, int markIndex, InterruptionSignal interruptionSignal, SteinerSubgraphOptions options)
//	{
//		throw new UnsupportedOperationException();
//	}
	
	@Override
	public float getWeightOfMarkedTree(int markIndex)
	{
		throw new UnsupportedOperationException();
	}
}
