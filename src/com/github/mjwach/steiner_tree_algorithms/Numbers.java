package com.github.mjwach.steiner_tree_algorithms;

import java.util.Random;

public class Numbers
{
	public static final int floatStoredMantissaSizeInBits = 23;
	public static final int floatStoredExponentSizeInBits = 8;
	public static final int floatMantissaMask =			0b00000000011111111111111111111111;
	public static final int floatMantissaImplied1Mask =	0b00000000100000000000000000000000;
	public static final int floatExponentMask =			0b01111111100000000000000000000000;
	public static final int floatExponentSignMask =		0b01000000000000000000000000000000;
	public static final int floatExponentIndicatingDenormalization = Math.getExponent(0f);
	public static final int doubleStoredMantissaSizeInBits = 52;
	public static final long greatestLongInTheMainBodyOfLongsFittingInADouble = 1L << (Numbers.doubleStoredMantissaSizeInBits + 1);
	public static final long leastLongInTheMainBodyOfLongsThatFitInADouble = -greatestLongInTheMainBodyOfLongsFittingInADouble;
	public static final float greatestFloatThatFitsInAnInt = 2147483520;
	public static final float leastFloatThatFitsInAnInt = Integer.MIN_VALUE;
	public static final float greatestFloatThatFitsInALong = 9.2233715E18f;
	public static final float leastFloatThatFitsInALong = -9.223372E18f;
	public static final double pi = Math.PI;
	public static final double piOver2 = pi * .5;
	public static final double piTimes1pt5 = pi * 1.5;
	public static final double piTimes2 = pi * 2;
	public static final float pi_f = (float) pi;
	public static final float piOver2_f = (float) piOver2;
	public static final float piTimes1pt5_f = (float) piTimes1pt5;
	public static final float piTimes2_f = (float) piTimes2;
	public static final double sqrt2 = Math.sqrt(2.0);
	public static final float sqrt2_f = (float) sqrt2;
	public static final double oneOverSqrt2 = Math.sqrt(.5);
	public static final float oneOverSqrt2_f = (float) oneOverSqrt2;
	public static final double sqrt3_over2 = Math.sqrt(3.0) / 2;
	public static final float sqrt3_over2_f = (float) sqrt3_over2;
	public static final double oneOver2Sqrt3 = .5 / Math.sqrt(3);
	public static final float oneOver2Sqrt3_f = (float) oneOver2Sqrt3;
	public static final int doubleSizeInBytes = Double.SIZE / 8;
	public static final int floatSizeInBytes = Float.SIZE / 8;
	public static final int intSizeInBytes = Integer.SIZE / 8;
	public static final int shortSizeInBytes = Short.SIZE / 8;
	public static final int subLongAddressSize = 6;
	
	public static int getMantissaBits(float f)
	{
		return Float.floatToRawIntBits(f) & floatMantissaMask;
	}
	
	public static int getMantissaBitsIncludingAnImplied1(float f)
	{
		return Float.floatToRawIntBits(f) & floatMantissaMask | floatMantissaImplied1Mask;
	}

	public static long getMantissaBitsIncludingAnImplied1(double d)
	{
		return Double.doubleToRawLongBits(d) & 0b1111111111111111111111111111111111111111111111111111L
			   |							  0b10000000000000000000000000000000000000000000000000000L;
	}

	public static boolean isInTheMainBodyOfLongsFittingInADouble(long n)
	{
		return n <= greatestLongInTheMainBodyOfLongsFittingInADouble && n >= leastLongInTheMainBodyOfLongsThatFitInADouble;
	}

	public static int bitCountCappedAt2(int n)
	{
		long lowest = Integer.lowestOneBit(n);
		return lowest == 0 ? 0 : lowest == Integer.highestOneBit(n) ? 1 : 2;
	}
	
	public static int bitCountCappedAt2(long n)
	{
		long lowest = Long.lowestOneBit(n);
		return lowest == 0 ? 0 : lowest == Long.highestOneBit(n) ? 1 : 2;
	}

	public static long roundGeometrically(double n)
	{
		long floor = (long) Math.floor(n);
		long ceiling = floor + 1;
		double mean = Math.sqrt(floor * ceiling);
		if (n < 0) mean = -mean;
		
		return n < mean ? floor : ceiling;
	}

	public static int signum(int n)
	{
		return n > 0 ? 1 : n == 0 ? 0 : -1;
	}

	private static final int[] shellSortGaps = new int[]{701, 301, 132, 57, 23, 10, 4, 1};
	
	public static void shellSort(int[] array, IntComparator comparator)
	{
		for (int g = array.length < 57 ? 4 : 0; g < shellSortGaps.length; g++)
			for (int gap = shellSortGaps[g], i = gap; i < array.length; i++)
			{
				int ai = array[i];
				
				int j;
				for (j = i; j >= gap && comparator.lessThan(ai, array[j - gap]); j -= gap)
					array[j] = array[j - gap];
				
				array[j] = ai;
			}
	}

	public static void insertionSort(int[] array, IntComparator comparator)
	{
		insertionSort(array, 0, array.length, comparator);
	}
	
	public static void insertionSort(int[] array, int i0, int i1, IntComparator comparator)
	{
		for (int i = i0; i < i1; i++)
		{
			int buffer = array[i];
			
			int j;
			for (j = i; j > 0 && comparator.lessThan(buffer, array[j - 1]); j--)
				array[j] = array[j - 1];
			
			array[j] = buffer;
		}
	}
	
	public static void insertionSort(int[] array, int[] passenger, int i0, int i1, IntComparator comparator)
	{
		for (int i = i0; i < i1; i++)
		{
			int buffer = array[i]; int buffer2 = passenger[i];
			
			int j;
			for (j = i; j > 0 && comparator.lessThan(buffer, array[j - 1]); j--)
			{
				array[j] = array[j - 1];
				passenger[j] = passenger[j - 1];
			}
			
			array[j] = buffer;
			passenger[j] = buffer2;
		}
	}

	public static void quicksort(int[] array, IntComparator comparator, Random random)
	{
		quicksort(array, 0, array.length, comparator, random);
	}
	
	public static void quicksort(int[] array, int i0, int i1, IntComparator comparator, Random random)
	{
		int count = i1 - i0;
		
		if (count > quicksortInsertionSortBoundary)
		{
			int pivotIndex = i0 + random.nextInt(count);

			int pivot = array[pivotIndex];
			array[pivotIndex] = array[i0];
			// array[0] will now be ignored for a bit; by rights it should contain the pivot but there's no need to actually set that

			int left = i0, right = i1;
			
			for (;;)
			{
				do ++left;
				while (left < i1 && comparator.lessThan(array[left], pivot));
				
				do --right;
				while (right > i0 && comparator.lessThan(pivot, array[right]));
				
				if (left > right)
					break;

				{
					int buffer = array[left];
					array[left] = array[right];
					array[right] = buffer;
				}
			}

			array[i0] = array[right];
			array[right] = pivot;
			
			quicksort(array, i0, right, comparator, random);
			quicksort(array, left, i1, comparator, random);
		}
		else
			insertionSort(array, i0, i1, comparator);
	}

	private static final int quicksortInsertionSortBoundary = 24;

	public static void quicksort(int[] array, int[] passenger, IntComparator comparator, Random random)
	{
		quicksort(array, passenger, 0, array.length, comparator, random);
	}
	
	public static void quicksort(int[] array, int[] passenger, int i0, int i1, IntComparator comparator, Random random)
	{
		int count = i1 - i0;
		
		if (count > quicksortInsertionSortBoundary)
		{
			int pivotIndex = i0 + random.nextInt(count);

			int pivot = array[pivotIndex]; int pivotPassenger = passenger[pivotIndex];
			array[pivotIndex] = array[i0]; passenger[pivotIndex] = passenger[i0];

			int left = i0, right = i1;
			
			for (;;)
			{
				do ++left;
				while (left < i1 && comparator.lessThan(array[left], pivot));
				
				do --right;
				while (right > i0 && comparator.lessThan(pivot, array[right]));
				
				if (left > right)
					break;

				{
					int buffer = array[left];
					array[left] = array[right];
					array[right] = buffer;
					
					buffer = passenger[left];
					passenger[left] = passenger[right];
					passenger[right] = buffer;
				}
			}

			array[i0] = array[right]; passenger[i0] = passenger[right];
			array[right] = pivot; passenger[right] = pivotPassenger;
			
			quicksort(array, passenger, i0, right, comparator, random);
			quicksort(array, passenger, left, i1, comparator, random);
		}
		else
			insertionSort(array, passenger, i0, i1, comparator);
	}

	public static float square(float x)
	{
		return x * x;
	}

	public static float mix(float x0, float x1, float t)
	{
		return x0 + t * (x1 - x0);
	}
}
