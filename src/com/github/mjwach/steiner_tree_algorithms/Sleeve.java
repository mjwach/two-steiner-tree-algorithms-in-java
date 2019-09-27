package com.github.mjwach.steiner_tree_algorithms;

import java.io.Serializable;

public class Sleeve<Identity> implements Serializable
{
	private static final long serialVersionUID = -7302225733982077943L;
	private Identity identity;

	public Sleeve(Identity identity)
	{
		this.identity = identity;
	}
	
	public Sleeve()
	{
	}
	
	public Identity getIdentity()
	{
		return identity;
	}

	public boolean hasIdentity()
	{
		return identity != null;
	}
	
	public Identity setIdentity(Identity identity)
	{
		this.identity = identity;
		return identity;
	}
	
	@Override
	public boolean equals(Object object)
	{
		return object instanceof Sleeve && identity == ((Sleeve<?>) object).identity;
	}
	
	@Override
	public int hashCode()
	{
		return System.identityHashCode(identity);
	}
	
	@Override
	public String toString()
	{
		return identity.toString();
	}

	public Sleeve<Identity> copy()
	{
		return new Sleeve<Identity>(identity);
	}
}
