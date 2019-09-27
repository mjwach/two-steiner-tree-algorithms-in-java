package com.github.mjwach.steiner_tree_algorithms;

import java.io.Serializable;
import java.util.Comparator;

public class Pair<FirstElement, SecondElement> implements Serializable
{
	private static final long serialVersionUID = 2462276754751613724L;
	private FirstElement firstElement;
	private SecondElement secondElement;

	public static <FirstElement extends Comparable<FirstElement>, SecondElement> Comparator<Pair<FirstElement, SecondElement>> getFirstElementComparator()
	{
		return new Comparator<Pair<FirstElement, SecondElement>>()
		{
			@Override
			public int compare(Pair<FirstElement, SecondElement> pair1, Pair<FirstElement, SecondElement> pair2)
			{
				return pair1.getFirstElement().compareTo(pair2.getFirstElement());
			}
		};
	}

	public Pair()
	{
	}

	public Pair(FirstElement firstElement, SecondElement secondElement)
	{
		this.firstElement = firstElement;
		this.secondElement = secondElement;
	}

	public FirstElement getFirstElement()
	{
		return firstElement;
	}

	public SecondElement getSecondElement()
	{
		return secondElement;
	}
	
	public FirstElement get0()
	{
		return firstElement;
	}
	
	public SecondElement get1()
	{
		return secondElement;
	}

	public <Element_ extends FirstElement> Element_ setFirstElement(Element_ value)
	{
		firstElement = value;
		return value;
	}

	public <Element_ extends SecondElement> Element_ setSecondElement(Element_ value)
	{
		secondElement = value;
		return value;
	}
	
	public void set(FirstElement newFirstElement, SecondElement newSecondElement)
	{
		firstElement = newFirstElement;
		secondElement = newSecondElement;
	}
	
	@Override
	public boolean equals(Object object)
	{
		if (object instanceof Pair)
		{
			Pair<?, ?> p = (Pair<?, ?>) object;
			return elementsMatch(firstElement, p.firstElement) && elementsMatch(secondElement, p.secondElement);
		}
		else
			return false;
	}
	
	private static boolean elementsMatch(Object e1, Object e2)
	{
		return e1 == e2 || e1 != null && e1.equals(e2);
	}

	@Override
	public String toString()
	{
		return "(" + firstElement + "," + secondElement + ")";
	}
}
