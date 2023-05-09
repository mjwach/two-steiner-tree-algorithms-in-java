package com.github.mjwach.steiner_tree_algorithms.examples;

import com.github.mjwach.steiner_tree_algorithms.EdgeArraysGraph;
import com.github.mjwach.steiner_tree_algorithms.EdgeArraysGraph.Edge;
import com.github.mjwach.steiner_tree_algorithms.EdgeArraysGraph.EdgeFactory;
import com.github.mjwach.steiner_tree_algorithms.EdgeArraysGraph.EuclideanVertex;
import com.github.mjwach.steiner_tree_algorithms.EdgeArraysGraph.SteinerSubgraphOptions;
import com.github.mjwach.steiner_tree_algorithms.EdgeArraysGraph.Vertex;
import com.github.mjwach.steiner_tree_algorithms.NullInterruptionSignal;
import com.github.mjwach.steiner_tree_algorithms.Numbers;
import com.github.mjwach.steiner_tree_algorithms.PlaneVector;

import gnu.trove.list.array.TIntArrayList;

public class SubgraphManhattanDistanceExample
{
	private static class CityEdge extends Edge
	{
		public CityEdge(Vertex end0, Vertex end1)
		{
			super(end0, end1);
		}

		@Override
		public float getWeightRelatedQuantityForComparison()
		{
			return getWeight();
		}

		@Override
		public float getWeight()
		{
			EuclideanVertex v0 = (EuclideanVertex) end0;
			EuclideanVertex v1 = (EuclideanVertex) end1;
			PlaneVector l0 = v0.getLocationAsPlaneVector();
			PlaneVector l1 = v1.getLocationAsPlaneVector();
			return Math.abs(l0.getX() - l1.getX()) + Math.abs(l0.getY() - l1.getY());
		}

		@Override
		public float getWeight(float weightRelatedQuantity)
		{
			return weightRelatedQuantity;
		}
	}
	
	public static void main(String[] args)
	{
		Printing.produceMathematicaCode = false;  // change this to true if you want
		
		if (Printing.produceMathematicaCode) System.out.print("{");
		
		printGraph("Graph", true, "Green", EdgeArraysGraph.defaultEdgeFactory);
		System.out.println(Printing.produceMathematicaCode ? "," : "");
		printGraph("Steiner tree of graph (Euclidean metric)", false, "Blue", EdgeArraysGraph.defaultEdgeFactory);
		System.out.println(Printing.produceMathematicaCode ? "," : "");
		printGraph("Steiner tree of graph (Manhattan metric)", false, "Magenta", new EdgeFactory()
		{			
			@Override
			public Edge act(Vertex end0, Vertex end1)
			{
				return new CityEdge(end0, end1);
			}
		});

		if (Printing.produceMathematicaCode) System.out.println("(* Paste the entire output of this program into Mathematica to see the graphs!  Note that copying from the Eclipse console can add weird formatting - paste into a plain text editor first to strip that away! *)}");
	}

	private static void printGraph(String graphName, boolean drawGraphNotSteinerTree, String color, EdgeFactory edgeFactory)
	{
		PlaneVector v0 = new PlaneVector(0, 1);
		PlaneVector v1 = v0.rotatedBy(Numbers.piTimes2_f / 3);
		PlaneVector v2 = v1.rotatedBy(Numbers.piTimes2_f / 3);
		
		PlaneVector v3 = v1.rotatedBy(Numbers.piTimes2_f / 6).times(.2f);
		PlaneVector v4 = v3.rotatedBy(Numbers.piTimes2_f / 3);
		PlaneVector v5 = v4.rotatedBy(Numbers.piTimes2_f / 3);

		PlaneVector v6 = v3.averagedWith(v4);
		PlaneVector v7 = v4.averagedWith(v5);
		PlaneVector v8 = v5.averagedWith(v3);

		EdgeArraysGraph graph = new EdgeArraysGraph();
		graph.addVertex(new EuclideanVertex(graph.getVertexCount(), v0));
		graph.addVertex(new EuclideanVertex(graph.getVertexCount(), v1));
		graph.addVertex(new EuclideanVertex(graph.getVertexCount(), v2));
		graph.addVertex(new EuclideanVertex(graph.getVertexCount(), v3));
		graph.addVertex(new EuclideanVertex(graph.getVertexCount(), v4));
		graph.addVertex(new EuclideanVertex(graph.getVertexCount(), v5));
		graph.addVertex(new EuclideanVertex(graph.getVertexCount(), v6));
		graph.addVertex(new EuclideanVertex(graph.getVertexCount(), v7));
		graph.addVertex(new EuclideanVertex(graph.getVertexCount(), v8));

		graph.getEdgeArrayVertex(0).addEdge(graph.getEdgeArrayVertex(1), edgeFactory);
		graph.getEdgeArrayVertex(1).addEdge(graph.getEdgeArrayVertex(2), edgeFactory);
		graph.getEdgeArrayVertex(2).addEdge(graph.getEdgeArrayVertex(0), edgeFactory);

		graph.getEdgeArrayVertex(3).addEdge(graph.getEdgeArrayVertex(6), edgeFactory);
		graph.getEdgeArrayVertex(6).addEdge(graph.getEdgeArrayVertex(4), edgeFactory);
		graph.getEdgeArrayVertex(4).addEdge(graph.getEdgeArrayVertex(7), edgeFactory);
		graph.getEdgeArrayVertex(7).addEdge(graph.getEdgeArrayVertex(5), edgeFactory);
		graph.getEdgeArrayVertex(5).addEdge(graph.getEdgeArrayVertex(8), edgeFactory);
		graph.getEdgeArrayVertex(8).addEdge(graph.getEdgeArrayVertex(3), edgeFactory);

		graph.getEdgeArrayVertex(0).addEdge(graph.getEdgeArrayVertex(7), edgeFactory);
		graph.getEdgeArrayVertex(1).addEdge(graph.getEdgeArrayVertex(8), edgeFactory);
		graph.getEdgeArrayVertex(2).addEdge(graph.getEdgeArrayVertex(6), edgeFactory);
		
		TIntArrayList terminals = new TIntArrayList();
		terminals.add(0);
		terminals.add(1);
		terminals.add(2);

		int steinerMarkIndex = 0;
		
		SteinerSubgraphOptions options = new SteinerSubgraphOptions();
		options.maximalMultipleVertexJumpCountPerPatch = Integer.MAX_VALUE;
		
		graph.markSteinerSubgraph(terminals, steinerMarkIndex, NullInterruptionSignal.instance, options);
		
		Printing.performOn(graph, drawGraphNotSteinerTree ? -1 : steinerMarkIndex, steinerMarkIndex + 1, graphName, color);
	}
}
