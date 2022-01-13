/*
 * SCHype - Java package for spectral clustering in hypergraphs
 * 
 * Copyright (C) 2012 Tom Michoel (The Roslin Institute, University of Edinburgh)
 * 
 */

import java.util.*;

public class Edge {
	
	public ArrayList<String> vertices; // ArrayList ensures that we output nodes edges in same order as we get them
	
	public ArrayList<String> sourceVertices; // only used in directed hypergraphs
	
	public ArrayList<String> targetVertices; // only used in directed hypergraphs
	
	public double weight;
	
	public Edge(){
		
	}
	
	/**
	 * Undirected edge
	 * @param vertices
	 */
	public Edge(ArrayList<String> vertices){
		this.vertices = vertices;
		this.weight = 1.0;
	}
	
	/**
	 * Undirected edge
	 * @param vertices
	 * @param weight
	 */
	public  Edge(ArrayList<String> vertices, double weight){
		this.vertices = vertices;
		this.weight = weight;
	}

	
	/**
	 * Directed edge
	 * @param vertices
	 */
	public Edge(ArrayList<String> sourceVertices, ArrayList<String> targetVertices){
		this.sourceVertices = sourceVertices;
		this.targetVertices = targetVertices;
		this.vertices = new ArrayList<String>();
		for (String node : this.sourceVertices)
			this.vertices.add(node);
		for (String node : this.targetVertices)
			this.vertices.add(node);
		this.weight = 1.0;
	}
	
	/**
	 * Directed edge
	 * @param vertices
	 */
	public Edge(ArrayList<String> sourceVertices, ArrayList<String> targetVertices, double weight){
		this.sourceVertices = sourceVertices;
		this.targetVertices = targetVertices;
		this.vertices = new ArrayList<String>();
		for (String node : this.sourceVertices)
			this.vertices.add(node);
		for (String node : this.targetVertices)
			this.vertices.add(node);
		this.weight = weight;
	}
	
	/**
	 * Check if an edge is a subset of a cluster
	 * @param clust
	 * @return
	 */
	public boolean belongsTo (HashSet<String> clust){
		boolean isSubSet = true;
		for (String node : this.vertices){
			isSubSet = isSubSet & clust.contains(node);
			if (!isSubSet)
				break;
		}
		return isSubSet;
	}
	
	public boolean belongsTo (HashSet<String> sourceClust, HashSet<String> targetClust){
		boolean isSubSet = true;
		for (String node : this.sourceVertices){
			isSubSet = isSubSet & sourceClust.contains(node);
			if (!isSubSet){
				return isSubSet;
			}
		}
		for (String node : this.targetVertices){
			isSubSet = isSubSet & targetClust.contains(node);
			if (!isSubSet){
				return isSubSet;
			}
		}
		return isSubSet;
	}
	
	public String toString(boolean withWeight){
		String out = "";
		for (int k=0; k<this.vertices.size(); k++)
			out += this.vertices.get(k) + "\t";
		if (withWeight)
			out += ">\t" + String.valueOf(this.weight);
		out.trim();
		return out;
	}
	
	public String toStringDirected(boolean withWeight){
		String out = "";
		for (int k=0; k<this.sourceVertices.size(); k++)
			out += this.sourceVertices.get(k) + "\t";
		out += "|\t";
		for (int k=0; k<this.targetVertices.size(); k++)
			out += this.targetVertices.get(k) + "\t";
		if (withWeight)
			out += ">\t" + String.valueOf(this.weight);
		out.trim();
		return out;
	}
}
