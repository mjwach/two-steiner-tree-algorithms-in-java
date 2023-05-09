package com.github.mjwach.steiner_tree_algorithms.examples;

import com.github.mjwach.steiner_tree_algorithms.EdgeArraysGraph;
import com.github.mjwach.steiner_tree_algorithms.NullInterruptionSignal;

import gnu.trove.list.array.TIntArrayList;

import com.github.mjwach.steiner_tree_algorithms.EdgeArraysGraph.EuclideanVertex;
import com.github.mjwach.steiner_tree_algorithms.EdgeArraysGraph.SteinerSubgraphOptions;

public class Subgraph2DExample
{
	public static void main(String[] args)
	{
		Printing.produceMathematicaCode = false;  // change this to true if you want
		
		EdgeArraysGraph graph = new EdgeArraysGraph();
		int size = 6;
		for (int x = 0; x < size; ++x)
			for (int y = 0; y < size; ++y)
				graph.addVertex(new EuclideanVertex(graph.getVertexCount(), x, y));
		
		// this grid of vertices will be given edges that form a square grid with exes in the squares; here are the horizontal grid edges:
		for (int y = 0; y < size; ++y)
			for (int x = 0; x < size - 1; ++x)
				graph.getVertex(size * x + y).addEdge(graph.getVertex(size * (x + 1) + y));
		
		// and here are the vertical ones:
		for (int x = 0; x < size; ++x)
			for (int y = 0; y < size - 1; ++y)
				graph.getVertex(size * x + y).addEdge(graph.getVertex(size * x + y + 1));
		
		// and here are the diagonals:
		for (int x = 0; x < size - 1; ++x)
			for (int y = 0; y < size - 1; ++y)
			{
				graph.getVertex(size * x + y).addEdge(graph.getVertex(size * (x + 1) + y + 1));
				graph.getVertex(size * (x + 1) + y).addEdge(graph.getVertex(size * x + y + 1));
			}

		int terminal0X = 3, terminal0Y = 0;
		int terminal1X = 0, terminal1Y = 3; // just some arbitrary grid locations that make a somewhat interesting-looking tree
		int terminal2X = 5, terminal2Y = 2;
		int terminal3X = 3, terminal3Y = 5;
		
		TIntArrayList terminals = new TIntArrayList();
		terminals.add(size * terminal0X + terminal0Y);
		terminals.add(size * terminal1X + terminal1Y);
		terminals.add(size * terminal2X + terminal2Y);
		terminals.add(size * terminal3X + terminal3Y);

		int steinerMarkIndex = 0;
		
		SteinerSubgraphOptions options = new SteinerSubgraphOptions();
		options.maximalMultipleVertexJumpCountPerPatch = Integer.MAX_VALUE;  // this graph is small, so there's little need to potentially introduce randomness by limiting this
		
		graph.markSteinerSubgraph(terminals, steinerMarkIndex, NullInterruptionSignal.instance, options);
		
		if (Printing.produceMathematicaCode) System.out.print("{");
		Printing.performOn(graph, steinerMarkIndex, steinerMarkIndex + 1, "Steiner tree", "Magenta");
		if (Printing.produceMathematicaCode) System.out.println("(* Paste the entire output of this program into Mathematica to see the graph!  Note that copying from the Eclipse console can add weird formatting - paste into a plain text editor first to strip that away! *)}");
	}
}
