package com.github.mjwach.steiner_tree_algorithms.examples;

import java.util.Random;

import com.github.mjwach.steiner_tree_algorithms.EdgeArraysGraph;
import com.github.mjwach.steiner_tree_algorithms.NullInterruptionSignal;
import com.github.mjwach.steiner_tree_algorithms.PlaneVector;

import gnu.trove.list.array.TIntArrayList;

import com.github.mjwach.steiner_tree_algorithms.EdgeArraysGraph.Edge;
import com.github.mjwach.steiner_tree_algorithms.EdgeArraysGraph.EdgeFactory;
import com.github.mjwach.steiner_tree_algorithms.EdgeArraysGraph.SteinerSubgraphOptions;
import com.github.mjwach.steiner_tree_algorithms.EdgeArraysGraph.Vertex;

public class SubgraphArbitraryWeightsExample
{
	public static class CustomVertex extends Vertex
	{
		public final int x;
		public final int y;

		public CustomVertex(int index, int x, int y)
		{
			super(index);
			this.x = x;
			this.y = y;
		}
	}
	
	private static class CustomEdge extends Edge
	{
		private final float arbitraryWeight;

		public CustomEdge(Vertex end0, Vertex end1, float arbitraryWeight)
		{
			super(end0, end1);
			this.arbitraryWeight = arbitraryWeight;
		}

		@Override
		public float getWeightRelatedQuantityForComparison()
		{
			return arbitraryWeight;
		}

		@Override
		public float getWeight()
		{
			return arbitraryWeight;
		}

		@Override
		public float getWeight(float weightRelatedQuantity)
		{
			return weightRelatedQuantity;
		}
	}
	
	private static interface WeightFunction
	{
		float actOn(PlaneVector v, int gridLength);
	}
	
	public static void main(String[] args)
	{
		Printing.produceMathematicaCode = false;  // change this to true if you want
		
		if (Printing.produceMathematicaCode) System.out.print("{");
	
		printGraph("Steiner tree (where the main graph is a grid with trivial edge weights)", new WeightFunction()
		{
			@Override
			public float actOn(PlaneVector v, int gridLength)
			{
				return 1;
			}
		});
		
		System.out.println(Printing.produceMathematicaCode ? "," : "");
	
		printGraph("Steiner tree (where the main graph is a grid with arbitrarily sine-waving edge weights)", new WeightFunction()
		{
			@Override
			public float actOn(PlaneVector v, int gridLength)
			{
				float x = v.getX();
				float y = v.getY();
				return 2.01f + (float) (Math.sin(x / 5.5f) * Math.sin(y / 6.4f) + Math.sin(x / 4.5f) * Math.sin(y / 2.9f));
			}
		});
		
		System.out.println(Printing.produceMathematicaCode ? "," : "");
	
		printGraph("Steiner tree (where the main graph is a grid with weights that smoothly and rapidly grow with distance from graph center)", new WeightFunction()
		{
			@Override
			public float actOn(PlaneVector v, int gridLength)
			{
				float h = ((float) gridLength - 1) / 2;
				float d = v.getDistanceFrom(new PlaneVector(h, h));
				return d * d * d;
			}
		});
		
		System.out.println(Printing.produceMathematicaCode ? "," : "");
	
		printGraph("Steiner tree (where the main graph is a grid with random weights)", new WeightFunction()
		{
			private Random random = new Random(888888);
			
			@Override
			public float actOn(PlaneVector v, int gridLength)
			{
				return random.nextFloat();
			}
		});

		if (Printing.produceMathematicaCode) System.out.println("(* Paste the entire output of this program into Mathematica to see the graphs!  Note that copying from the Eclipse console can add weird formatting - paste into a plain text editor first to strip that away! *)}");
	}

	private static void printGraph(String graphName, WeightFunction weightFunction)
	{
		EdgeArraysGraph graph = new EdgeArraysGraph();
		int size = 64;
		for (int x = 0; x < size; ++x)
			for (int y = 0; y < size; ++y)
				graph.addVertex(new CustomVertex(graph.getVertexCount(), x, y));
		
		EdgeFactory edgeFactory = new EdgeFactory()
		{	
			@Override
			public Edge act(Vertex end0, Vertex end1)
			{
				CustomVertex v0 = (CustomVertex) end0;
				CustomVertex v1 = (CustomVertex) end1;
				PlaneVector center = new PlaneVector(v0.x, v0.y).averagedWith(new PlaneVector(v1.x, v1.y));
				
				return new CustomEdge(end0, end1, weightFunction.actOn(center, size));
			}
		};
		
		class GridGraph
		{
			public Vertex getVertex(int x, int y)
			{
				return graph.getEdgeArrayVertex(size * x + y);
			}

			public int getIndex(int x, int y)
			{
				return getVertex(x, y).getIndex();
			}
		}
		
		GridGraph grid = new GridGraph();
		
		// this grid of vertices will be given edges that form a square grid, with no diagonal edges
		for (int y = 0; y < size; ++y)
			for (int x = 0; x < size - 1; ++x)
				grid.getVertex(x, y).addEdge(grid.getVertex(x + 1, y), edgeFactory);
		
		for (int x = 0; x < size; ++x)
			for (int y = 0; y < size - 1; ++y)
				grid.getVertex(x, y).addEdge(grid.getVertex(x, y + 1), edgeFactory);

		int M = size - 1;
		int h = size / 2;
		
		TIntArrayList terminals = new TIntArrayList();
		terminals.add(grid.getIndex(0, 0));
		terminals.add(grid.getIndex(h, 0));
		terminals.add(grid.getIndex(M, 0));
		terminals.add(grid.getIndex(0, h));
		terminals.add(grid.getIndex(h, h));
		terminals.add(grid.getIndex(M, h));
		terminals.add(grid.getIndex(0, M));
		terminals.add(grid.getIndex(h, M));
		terminals.add(grid.getIndex(M, M));

		int steinerMarkIndex = 0;
		
		SteinerSubgraphOptions options = new SteinerSubgraphOptions();
		options.random = new Random();
		
		graph.markSteinerSubgraph(terminals, steinerMarkIndex, NullInterruptionSignal.instance, options);
		
		Printing.performOn(graph, steinerMarkIndex, steinerMarkIndex + 1, graphName, "Magenta");
	}
}
