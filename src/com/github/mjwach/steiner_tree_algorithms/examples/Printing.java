package com.github.mjwach.steiner_tree_algorithms.examples;

import com.github.mjwach.steiner_tree_algorithms.Action1;
import com.github.mjwach.steiner_tree_algorithms.DCELGraph;
import com.github.mjwach.steiner_tree_algorithms.EdgeArraysGraph;
import com.github.mjwach.steiner_tree_algorithms.EdgeArraysGraph.EuclideanVertex;
import com.github.mjwach.steiner_tree_algorithms.GraphEdge;
import com.github.mjwach.steiner_tree_algorithms.GraphVertex;
import com.github.mjwach.steiner_tree_algorithms.PlaneVector;
import com.github.mjwach.steiner_tree_algorithms.SpaceVector;

public class Printing
{
	public static boolean produceMathematicaCode;
	
	public static void performOn(DCELGraph graph, String graphName, String color)
	{
		System.out.println(produceMathematicaCode ? "(* " + graphName + " edges: *)" : graphName + " edges:");
		
		if (graph.vertices.length == 0)
		{
			System.out.println(produceMathematicaCode ? "\"(empty graph)\"" : "none");
			return;
		}

		if (produceMathematicaCode)
			System.out.println("Graphics[{" + color + ",");
		int workingMarkIndex = 0;
		printUnprintedEdges(graph.vertices[0], -1, workingMarkIndex);
		graph.vertices[0].clearMarkFromEdgesInConnectedSubgraph(workingMarkIndex);
		if (produceMathematicaCode)
			System.out.println("}]");
	}

	public static void performOn(EdgeArraysGraph graph, int markIndexOfSubgraphToBePrinted, int workingMarkIndex, String graphName, String color)
	{
		System.out.println(produceMathematicaCode ? "(* " + graphName + " edges: *)" : graphName + " edges:");
		
		GraphVertex startingVertex = null;
		for (GraphVertex vertex: graph)
			if (markIndexOfSubgraphToBePrinted < 0 || vertex.readMark(markIndexOfSubgraphToBePrinted))
			{
				startingVertex = vertex;
				break;
			}
		
		if (startingVertex == null)
		{
			System.out.println(produceMathematicaCode ? "\"(empty graph)\"" : "none");
			return;
		}

		if (produceMathematicaCode)
			System.out.println((is3DEuclideanGraph(graph) ? "Graphics3D" : "Graphics") + "[{" + color + ",");
		printUnprintedEdges(startingVertex, markIndexOfSubgraphToBePrinted, workingMarkIndex);
		startingVertex.clearMarkFromEdgesInConnectedSubgraph(workingMarkIndex);
		if (produceMathematicaCode)
			System.out.println("}]");
	}

	private static boolean is3DEuclideanGraph(EdgeArraysGraph graph)
	{
		if (graph.getVertexCount() == 0)
			return false;
		
		GraphVertex vertex = graph.getVertex(0);
		if (vertex instanceof EuclideanVertex)
			return ((EuclideanVertex) vertex).getLocation() instanceof SpaceVector;
		else
			return false;
	}

	private static void printUnprintedEdges(GraphVertex vertex, int filterMarkIndex, int workingMarkIndex)
	{
		vertex.doForEachEdge(new Action1<GraphEdge>()
		{
			@Override
			public void perform(GraphEdge edge)
			{
				if ((filterMarkIndex < 0 || edge.readMark(filterMarkIndex)) && !edge.readMark(workingMarkIndex))
				{
					edge.setMarkInBothDirections(workingMarkIndex);
					printEdge(edge);
					printUnprintedEdges(edge.getEnd(vertex), filterMarkIndex, workingMarkIndex);
				}
			}
		});
	}

	private static void printEdge(GraphEdge edge)
	{
		if (produceMathematicaCode) System.out.print("Line[{");
		printVertex(edge.getEnd0());
		System.out.print(produceMathematicaCode ? ", " : " -> ");
		printVertex(edge.getEnd1());
		System.out.println(produceMathematicaCode ? "}]," : "");
	}

	private static void printVertex(GraphVertex vertex)
	{
		if (vertex instanceof DCELGraph.Vertex)
			print2DVector(((DCELGraph.Vertex) vertex).getLocation());
		else if (vertex instanceof EdgeArraysGraph.EuclideanVertex)
		{
			EuclideanVertex euclideanVertex = (EdgeArraysGraph.EuclideanVertex) vertex;
			if (euclideanVertex.getLocation() instanceof PlaneVector)
				print2DVector(euclideanVertex.getLocationAsPlaneVector());
			else
				print3DVector((SpaceVector) euclideanVertex.getLocation());
		}
		else if (vertex instanceof SubgraphArbitraryWeightsExample.CustomVertex)
		{
			SubgraphArbitraryWeightsExample.CustomVertex customVertex = (SubgraphArbitraryWeightsExample.CustomVertex) vertex;
			print2DVector(new PlaneVector(customVertex.x, customVertex.y));
		}
		else
			throw new UnsupportedOperationException();
	}

	private static void print2DVector(PlaneVector v)
	{
		if (produceMathematicaCode)
			System.out.print("{" + p(v.getX()) + ", " + p(v.getY()) + "}");
		else
			System.out.print(v);
	}

	private static void print3DVector(SpaceVector v)
	{
		if (produceMathematicaCode)
			System.out.print("{" + p(v.getX()) + ", " + p(v.getY()) + ", " + p(v.getZ()) + "}");
		else
			System.out.print(v);
	}
	
	private static String p(float n)
	{
		return String.format("%.9f", n);
	}
}
