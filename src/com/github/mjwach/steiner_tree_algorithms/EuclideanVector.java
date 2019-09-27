package com.github.mjwach.steiner_tree_algorithms;

public interface EuclideanVector
{
	float getX();
	float getY();
	float getDistanceSquaredFrom_(EuclideanVector peer);
	float getDistanceSquaredFrom(PlaneVector v);
	float getDistanceSquaredFrom(SpaceVector v);
}
