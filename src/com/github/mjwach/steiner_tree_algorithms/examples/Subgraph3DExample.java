package com.github.mjwach.steiner_tree_algorithms.examples;

import java.util.Random;

import com.github.mjwach.steiner_tree_algorithms.EdgeArraysGraph;
import com.github.mjwach.steiner_tree_algorithms.NullInterruptionSignal;

import gnu.trove.list.array.TIntArrayList;

import com.github.mjwach.steiner_tree_algorithms.EdgeArraysGraph.EuclideanVertex;
import com.github.mjwach.steiner_tree_algorithms.EdgeArraysGraph.SteinerSubgraphOptions;
import com.github.mjwach.steiner_tree_algorithms.GraphVertex;

public class Subgraph3DExample
{
	public static void main(String[] args)
	{
		Printing.produceMathematicaCode = false;  // change this to true if you want
		
		EdgeArraysGraph graph = new EdgeArraysGraph();
		int size = 22;
		for (int x = 0; x < size; ++x)
			for (int y = 0; y < size; ++y)
				for (int z = 0; z < size; ++z)
					graph.addVertex(new EuclideanVertex(graph.getVertexCount(), x, y, z));

		int M = size - 1;
		
		class GridGraph
		{
			GraphVertex getVertex(int x, int y, int z)
			{
				return graph.getVertex(size * (size * x + y) + z);
			}
			
			int getIndex(int x, int y, int z)
			{
				return getVertex(x, y, z).getIndex();
			}

			public void addDiagonalEdge(int x, int y, int z, int random0123)
			{
				if (x < M && y < M && z < M)
					switch (random0123)
					{
					case 0: getVertex(x, y, z).addEdge(getVertex(x + 1, y + 1, z + 1)); break;
					case 1: getVertex(x + 1, y, z).addEdge(getVertex(x, y + 1, z + 1)); break;
					case 2: getVertex(x, y + 1, z).addEdge(getVertex(x + 1, y, z + 1)); break;
					case 3: getVertex(x, y, z + 1).addEdge(getVertex(x + 1, y + 1, z));
					}
			}
		}
		
		GridGraph grid = new GridGraph();
		
		// this grid of vertices will be given edges that form a cubic grid with two of four possible diagonal edges randomly chosen to stretch across each cell's center
		Random random = new Random(618474);
		
		for (int y = 0; y < size; ++y)
			for (int z = 0; z < size; ++z)
				for (int x = 0; x < size; ++x)
				{
					GraphVertex vertex = grid.getVertex(x, y, z);
					if (x < M)
						vertex.addEdge(grid.getVertex(x + 1, y, z));
					if (y < M)
						vertex.addEdge(grid.getVertex(x, y + 1, z));
					if (z < M)
						vertex.addEdge(grid.getVertex(x, y, z + 1));

					int random0123a = random.nextInt(4);
					int random0123b = (random0123a + random.nextInt(3)) % 4;
					grid.addDiagonalEdge(x, y, z, random0123a);
					grid.addDiagonalEdge(x, y, z, random0123b);
				}

		int steinerMarkIndex = 0;
		
		int h = size / 2;
		TIntArrayList terminals = new TIntArrayList();
		terminals.add(grid.getIndex(0, 0, 0));
		terminals.add(grid.getIndex(M, 0, 0));
		terminals.add(grid.getIndex(h, M, 0));
		terminals.add(grid.getIndex(0, M, M));
		terminals.add(grid.getIndex(M, M, M));
		terminals.add(grid.getIndex(h, 0, M));
		
		SteinerSubgraphOptions options = new SteinerSubgraphOptions();
		options.maximalJumpDepthInMultipleVertexJiggling = 5;  // these values were not chosen carefully to produce an optimal combination of short running time and low tree weight...
		options.maximalMultipleVertexJumpCountPerPatch = 460;  // you must play around with all the options yourself to find good values for a given type of graph!
		options.random = random;
		
		graph.markSteinerSubgraph(terminals, steinerMarkIndex, NullInterruptionSignal.instance, options);
		
		if (Printing.produceMathematicaCode) System.out.print("{");
		Printing.performOn(graph, steinerMarkIndex, steinerMarkIndex + 1, "Steiner tree", "Magenta");
		if (Printing.produceMathematicaCode) System.out.println("(* Paste the entire output of this program into Mathematica to see the graph!  Note that copying from the Eclipse console can add weird formatting - paste into a plain text editor first to strip that away! *)}");
	}
}
