/*
 * SCHype - Java package for spectral clustering in hypergraphs
 * 
 * Copyright (C) 2012 Tom Michoel (The Roslin Institute, University of Edinburgh)
 * 
 */

import java.util.*;

public class Cluster {

	public HyperGraph hyperGraph;
	
	public HashSet<String> vertices;
	
	public HashSet<String> sourceVertices;
	
	public HashSet<String> targetVertices;
	
	public HashSet<Edge> edges;
	
	public double score;
	
	public double scoreBound;
	
	public HashMap<String, Double> perronVector;
	
	public HashMap<String, Double> perronVectorX;
	
	public HashMap<String, Double> perronVectorY;
	
	public Cluster(){
		
	}
	
	public Cluster(HyperGraph hg, HashSet<String> elements){
		this.hyperGraph = hg;
		this.vertices = elements;
		this.setEdges();
	}
	
	public Cluster(HyperGraph hg, HashSet<String> sourceVertices, HashSet<String> targetVertices){
		this.hyperGraph = hg;
		this.sourceVertices = sourceVertices;
		this.targetVertices = targetVertices;
		this.vertices = new HashSet<String>();
		for (String node : this.sourceVertices)
			this.vertices.add(node);
		for (String node : this.targetVertices)
			this.vertices.add(node);
		this.setEdges();
	}
	
	public Cluster(HyperGraph hg, HashSet<String> elements, double p){
		this.hyperGraph = hg;
		this.vertices = elements;
		this.setEdges();
		this.score = this.clusterScore(p);
	}
	
	public Cluster(HyperGraph hg, HashSet<String> sourceVertices, HashSet<String> targetVertices, double p, double q){
		this.hyperGraph = hg;
		this.sourceVertices = sourceVertices;
		this.targetVertices = targetVertices;
		this.setEdges();
		this.vertices = new HashSet<String>();
		for (String node : this.sourceVertices)
			this.vertices.add(node);
		for (String node : this.targetVertices)
			this.vertices.add(node);
		this.score = this.clusterScore(p,q);
	}
	
	public Cluster(HyperGraph hg, HashSet<Edge> edges, boolean dummy){
		this.hyperGraph = hg;
		this.edges = edges;
		this.setNodes();
	}
	
	/**
	 * Score of a cluster
	 */
	public double clusterScore(double p){
		for (Edge edge : this.edges)
			score += edge.weight;
		if (!this.hyperGraph.directed)
			score = score/Math.pow((double)this.vertices.size(), 1.0/p);
		else
			score = score/Math.pow((double)this.sourceVertices.size()*(double)this.targetVertices.size(), 0.5/p);
		return score;
	}

	public double clusterScore(double p, double q){	
		for (Edge edge : this.edges)
			score += edge.weight;
		if (!this.hyperGraph.directed)
			score = score/Math.pow((double)this.vertices.size(), 1.0/p);
		else
			score = score/(Math.pow((double)this.sourceVertices.size(), 0.5/p)*Math.pow((double)this.targetVertices.size(), 0.5/q));
		return score;
	}
	
	/**
	 * Create set of vertices
	 */
	public void setNodes(){
		this.vertices = new HashSet<String>();
		if (this.hyperGraph.directed){
			this.sourceVertices = new HashSet<String>();
			this.targetVertices = new HashSet<String>();
		}
		for (Edge edge : this.edges){
			for (String node : edge.vertices)
				this.vertices.add(node);
			if (this.hyperGraph.directed){
				for (String node : edge.sourceVertices)
					this.sourceVertices.add(node);
				for (String node : edge.targetVertices)
					this.targetVertices.add(node);
			}
		}
	}
	
	/**
	 * Create set of edges belong to this cluster (for undirected hypergraphs)
	 */
	public void setEdges(){
		this.edges = new HashSet<Edge>();
		HashSet<Edge> marked = new HashSet<Edge>();
		if (!this.hyperGraph.directed){
			for (Edge edge : this.hyperGraph.edges){
				boolean isSubSet = edge.belongsTo(this.vertices);
				if (isSubSet)
					marked.add(edge);
			}
		} else {
			for (Edge edge : this.hyperGraph.edges){
				boolean isSubSet = edge.belongsTo(this.sourceVertices, this.targetVertices);
				if (isSubSet)
					marked.add(edge);
			}
		}
		for (Edge edge : marked)
			this.edges.add(edge);
	}
	
	/**
	 * Convert cluster to a hypergraph
	 * @return
	 */
	public HyperGraph toHyperGraph(){
		HyperGraph hg = new HyperGraph(this.edges,this.hyperGraph.directed, this.hyperGraph.weighted);
		System.out.println("Converted cluster to hypergraph with " + hg.vertices.size() + " vertices and " + hg.edges.size() + " edges.");
		return hg;
	}
	

}
