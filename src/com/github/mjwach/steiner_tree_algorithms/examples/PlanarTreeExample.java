package com.github.mjwach.steiner_tree_algorithms.examples;

import java.util.ArrayList;

import com.github.mjwach.steiner_tree_algorithms.DCELGraph;
import com.github.mjwach.steiner_tree_algorithms.PlaneVector;
import com.github.mjwach.steiner_tree_algorithms.DCELGraph.CalculatedSteinerTree;
import com.github.mjwach.steiner_tree_algorithms.DCELGraph.SteinerDuplicatePointPolicy;

public class PlanarTreeExample
{
	public static void main(String[] arguments)
	{
		Printing.produceMathematicaCode = false;  // change this to true if you want
		
		ArrayList<PlaneVector> points = new ArrayList<>();
		points.add(new PlaneVector(1, 0));
		points.add(new PlaneVector(-1, 0));
		points.add(new PlaneVector(1.5, 2));
		points.add(new PlaneVector(-1.5, 2));
		points.add(new PlaneVector(0, 3.5));
		
		CalculatedSteinerTree graphs = DCELGraph.makeSteinerTreeAndDelaunayTriangulation(points, SteinerDuplicatePointPolicy.includeOnlyOneVertexPerLocation);
		
		if (Printing.produceMathematicaCode) System.out.println("{");
		
		Printing.performOn(graphs.steinerTree, "Steiner tree", "Red");
		if (!Printing.produceMathematicaCode) System.out.println();

		if (Printing.produceMathematicaCode) System.out.println(", ");
		
		Printing.performOn(graphs.delaunayTriangulation, "Delaunay triangulation", "Green");

		if (Printing.produceMathematicaCode) System.out.println("}\n(* Paste the entire output of this program into Mathematica to see the graphs! *)");
	}
}
