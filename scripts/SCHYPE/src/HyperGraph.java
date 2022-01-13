/*
 * SCHype - Java package for spectral clustering in hypergraphs
 * 
 * Copyright (C) 2012 Tom Michoel (The Roslin Institute, University of Edinburgh)
 * 
 */


import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;


public class HyperGraph {

	public HashSet<String> vertices;
	
	public HashSet<Edge> edges;
	
	public HashMap<String, HashSet<Edge>> nodeEdges; // for each node, to which edges does it belong
	
	public HashMap<String, HashSet<Edge>> nodeSourceEdges; // only used for directed hypergraphs
	
	public HashMap<String, HashSet<Edge>> nodeTargetEdges; // only used for directed hypergraphs
	
	public HashMap<String, Double> perronVector; // for undirected hypergraphs
	
	public HashMap<String, Double> perronVectorX; // for directed hypergraphs
	
	public HashMap<String, Double> perronVectorY; // for directed hypergraphs
	
	public boolean directed;
	
	public boolean weighted;
	
	/**
	 * Empty constructor
	 */
	public HyperGraph(){
		this.vertices = new HashSet<String>();
		this.edges = new HashSet<Edge>();
	}
	
	/**
	 * Create from list of edges
	 */
	public HyperGraph(HashSet<Edge> edges, boolean directed, boolean weighted){
		this.edges = edges;
		this.directed = directed;
		this.setNodes();
		this.setEdgeNodes();
	}
	
	public HyperGraph(String file, boolean directed, boolean weighted){
		this.weighted = weighted;
		if (!directed)
			this.readHyperGraphUndirected(file);
		else
			this.readHyperGraphDirected(file);
	}
	
	/**
	 * Read undirected hypergraph from file
	 */
	public void readHyperGraphUndirected (String file) {
		this.edges = new HashSet<Edge>();
		try {
			Scanner fileScanner = new Scanner(new File(file)).useDelimiter("\\n");
			// walk through file, skip comment lines (starting with #)
			ArrayList<String> edge;
			String node;			
			double weight;
			while (fileScanner.hasNext()){
				Scanner line = new Scanner(fileScanner.next().toString()).useDelimiter("\\s");
				// reset weight
				weight = 1.0;
				// read first node
				node = line.next();
				node = node.trim();
				// skip if comment line
				if (!node.startsWith("#")){
					edge = new ArrayList<String>();
					edge.add(node);
					//read other nodes
					while (line.hasNext()){
						node = line.next();
						node = node.trim();
						if (!node.equals(">")) // separator between nodes and edge weight
							edge.add(node);
						else if (this.weighted) {
							node = line.next();
							weight = Double.valueOf(node);
						}
					}
					// put edge if it has more than one node
					if (edge.size()>1){
						Edge e = new Edge(edge, weight);
						this.edges.add(e);
					}
				}
			}
			this.directed = false;
			this.setNodes();
			this.setEdgeNodes();
			System.out.println("Read undirected hypergraph with " + this.vertices.size() + " nodes and " + this.edges.size() + " edges.");
		} catch (FileNotFoundException e) {
			System.out.println(e);
		}
	}
	
	/**
	 * Read directed hypergraph from file
	 */
	public void readHyperGraphDirected (String file) {
		this.edges = new HashSet<Edge>();
		try {
			Scanner fileScanner = new Scanner(new File(file)).useDelimiter("\\n");
			// walk through file, skip comment lines (starting with #)
			ArrayList<String> sourceVertices;
			ArrayList<String> targetVertices;
			String node;			
			double weight;
			while (fileScanner.hasNext()){
				Scanner line = new Scanner(fileScanner.next().toString()).useDelimiter("\\s");
				// reset weight
				weight = 1.0;
				// read first node
				node = line.next();
				node = node.trim();
				// skip if comment line
				boolean source;
				if (!node.startsWith("#")){
					source = true;
					sourceVertices = new ArrayList<String>();
					targetVertices = new ArrayList<String>();
					sourceVertices.add(node);
					//read other nodes
					while (line.hasNext()){
						node = line.next();
						node = node.trim();
						if (!node.equals(">")){ // separator between nodes and edge weight
							if (node.equals("|")){ // separator between source and target nodes
								source = false;
							} else {
								if (source){
									sourceVertices.add(node);
								} else{
									targetVertices.add(node);
								}
							}
						}else if (this.weighted){
							node = line.next();
							weight = Double.valueOf(node);
						}
					}
					Edge e = new Edge(sourceVertices, targetVertices, weight);
					this.edges.add(e);
				}
			}
			this.directed = true;
			this.setNodes();
			this.setEdgeNodes();
			System.out.println("Read directed hypergraph with " + this.vertices.size() + " nodes and " + this.edges.size() + " edges.");
		} catch (FileNotFoundException e) {
			System.out.println(e);
		}
	}
	/**
	 * Gather nodes from edge list
	 */
	public void setNodes(){
		this.vertices = new HashSet<String>();
		for (Edge edge : this.edges){
			for (String node : edge.vertices){
				this.vertices.add(node);
			}
		}
	}
	
	/**
	 * Make map from nodes to edges containing this node
	 */
	public void setEdgeNodes(){
		this.nodeEdges = new HashMap<String, HashSet<Edge>>();
		for (String node : vertices)
			this.nodeEdges.put(node, new HashSet<Edge>());
		if (this.directed){
			this.nodeSourceEdges = new HashMap<String, HashSet<Edge>>();
			this.nodeTargetEdges = new HashMap<String, HashSet<Edge>>();
			for (String node : vertices){
				this.nodeSourceEdges.put(node, new HashSet<Edge>());
				this.nodeTargetEdges.put(node, new HashSet<Edge>());
			}
		}
		for (Edge edge : this.edges){
			for (String node : edge.vertices){
				this.nodeEdges.get(node).add(edge);
			}
			if (this.directed){
				for (String node : edge.sourceVertices)
					this.nodeSourceEdges.get(node).add(edge);
				for (String node : edge.targetVertices)
					this.nodeTargetEdges.get(node).add(edge);
			}
		}
	}
	
	
	/**
	 * Remove set of edges if all its elements belong to a cluster (undirected case)
	 */
	public void removeEdges(HashSet<String> cluster){
		HashSet<Edge> marked = new HashSet<Edge>();
		for (Edge edge : this.edges){
			boolean isSubSet = true;
			for (String node : edge.vertices){
				isSubSet = isSubSet & cluster.contains(node);
			}
			if (isSubSet)
				marked.add(edge);
		}
		for (Edge edge : marked)
			this.edges.remove(edge);
		this.setNodes();
		this.setEdgeNodes();
	}
	
	/**
	 * Remove set of edges if all its elements belong to a cluster (directed case)
	 */
	public void removeEdges(HashSet<String> sourceCluster, HashSet<String> targetCluster){
		HashSet<Edge> marked = new HashSet<Edge>();
		for (Edge edge : this.edges){
			boolean isSubSet = true;
			for (String node : edge.sourceVertices){
				isSubSet = isSubSet & sourceCluster.contains(node);
			}
			for (String node : edge.targetVertices){
				isSubSet = isSubSet & targetCluster.contains(node);
			}
			if (isSubSet)
				marked.add(edge);
		}
		for (Edge edge : marked)
			this.edges.remove(edge);
		this.setNodes();
		this.setEdgeNodes();
	}
	
	/**
	 * Remove set of edges if one of its nodes belongs to a cluster (undirected case)
	 */
	public void removeNodes(HashSet<String> cluster){
		HashSet<Edge> marked = new HashSet<Edge>();
		for (Edge edge : this.edges){
			boolean isSubSet = false;
			for (String node : edge.vertices){
				isSubSet = isSubSet || cluster.contains(node);
			}
			if (isSubSet)
				marked.add(edge);
		}
		for (Edge edge : marked)
			this.edges.remove(edge);
		this.setNodes();
		this.setEdgeNodes();
	}
	
	/**
	 * Remove set of edges if one of its nodes belongs to a cluster (directed case)
	 */
	public void removeNodes(HashSet<String> sourceCluster, HashSet<String> targetCluster){
		HashSet<Edge> marked = new HashSet<Edge>();
		for (Edge edge : this.edges){
			boolean isSubSet = false;
			for (String node : edge.sourceVertices){
				isSubSet = isSubSet || sourceCluster.contains(node);
			}
			for (String node : edge.targetVertices){
				isSubSet = isSubSet || targetCluster.contains(node);
			}
			if (isSubSet)
				marked.add(edge);
		}
		for (Edge edge : marked)
			this.edges.remove(edge);
		this.setNodes();
		this.setEdgeNodes();
	}
	
	/**
	 * Create set of nodes forming an irreducible component
	 * @return
	 */
	public HashSet<String> irredComponent(){
		// start with node that has the largest number of edges
		int max = 0;
		String maxnode = "dummy";
		for (String node : this.nodeEdges.keySet()){
			if (this.nodeEdges.get(node).size()>max){
				max = this.nodeEdges.get(node).size();
				maxnode = node;
			}
		}
		HashSet<String> comp = new HashSet<String>();
		for (Edge edge : this.nodeEdges.get(maxnode))
			for (String node : edge.vertices)
				comp.add(node);
		// add neighbors until convergence
		int diff = 1;
		int size = comp.size();
		int sizenew;
		HashSet<String> added;
		while (diff>0){
			added = new HashSet<String>();
			for (String node : comp){
				for (Edge edge : this.nodeEdges.get(node))
					for (String nb : edge.vertices)
						added.add(nb);
			}
			for (String nb : added){
				comp.add(nb);
			}
			sizenew = comp.size();
			diff = sizenew - size;
			size = sizenew;
		}
		return comp;
	}
	
	/**
	 * Compute Perron vector (dominant eigenvector) for undirected hypergraph
	 * @return
	 */
	public void setPerronVector(double pnorm, double tolerance, int maxstep){
		
		HashMap<String,Double> pv = new HashMap<String,Double>();
		
		// initialize all to zero
		for (String node : vertices)
			pv.put(node, 0.0);
		
		// set irred component to 1
		HashSet<String> comp = this.irredComponent();
		
		for (String node : comp)
			pv.put(node, 1.0);
		
		normalize(pv, pnorm);
		// Power algorithm until convergence
		int step=0;
		double diff=1;
		double mu = 1.0, munew;
		while (diff>tolerance && step<maxstep){
			step++;
			// Hypergraph x vector multiplication
			pv = powermult(pv,pnorm);
			munew = norm(pv, pnorm);
			normalize(pv, munew);
			// Convergence parameter
			diff = Math.abs(1.0 - munew/mu);
			// Update mu
			mu = munew;
		}
		if (step >= maxstep)
			System.out.println("Maximum number of iterations reached, diff = " + diff);
		this.perronVector = pv;
	}
	
	/**
	 * Compute Perron vectors (dominant singular vectors) for directed hypergraph
	 * @return
	 */
	public void setPerronVectorXY(double pnorm, double qnorm, double tolerance, int maxstep){
		
		HashMap<String,Double> pvX = new HashMap<String,Double>();
		HashMap<String,Double> pvY = new HashMap<String,Double>();
		
		// initialize all to zero
		for (String node : vertices){
			pvX.put(node, 0.0);
			pvY.put(node, 0.0);
		}
		
		// set irred component to 1
		HashSet<String> comp = this.irredComponent();
		
		for (String node : comp){
			pvX.put(node, 1.0);
			pvY.put(node, 1.0);
		}
		
		normalize(pvX, pnorm);
		normalize(pvY, qnorm);
		
		// Power algorithm until convergence
		int step=0;
		double diff=1;
		double mu = 1.0, munew;
		double nu = 1.0, nunew;
		while (diff>tolerance && step<maxstep){
			step++;
			// Hypergraph x vector multiplication
			pvX = powermult(pvX,pvY,pnorm,1);
			munew = norm(pvX, pnorm);
			normalize(pvX, munew);
			
			pvY = powermult(pvX,pvY,qnorm,2);
			nunew = norm(pvY, qnorm);
			normalize(pvY, nunew);
			
			// Convergence parameter
			diff = Math.max(Math.abs(1.0 - munew/mu), Math.abs(1.0 - nunew/nu));
			// Update mu, nu
			mu = munew;
			nu = nunew;
		}
		if (step >= maxstep)
			System.out.println("Maximum number of iterations reached, diff = " + diff);
		
		this.perronVectorX = pvX;
		this.perronVectorY = pvY;
	}
	
	/**
	 * Hypergraph x vector 'multiplication'
	 * @param u
	 * @return
	 */
	public HashMap<String,Double> powermult(HashMap<String, Double> u, double pnorm){
		HashMap<String,Double> v = new HashMap<String,Double>();
		for (String node : this.vertices){
			double val = 0.0;
			for (Edge edge : this.nodeEdges.get(node)){
				// make product over nodes in edge
				double prod = 1.0;
				for (String nb : edge.vertices)
					prod = prod * u.get(nb);
				prod = Math.pow(prod, 1.0/(double)edge.vertices.size())/(double)edge.vertices.size();
				// add to sum
				val += edge.weight*prod;
			}
			v.put(node, Math.pow(val, 1.0/pnorm));
		}
		return v;
	}
	
	/**
	 * Hypergraph x vector 'multiplication'
	 * @param u
	 * @return
	 */
	public HashMap<String,Double> powermult(HashMap<String, Double> u1, HashMap<String, Double> u2, double pnorm, int dim){
		HashMap<String,Double> v = new HashMap<String,Double>();
		if (dim == 1){
			for (String node : this.vertices){
				double val = 0.0;
				for (Edge edge : this.nodeSourceEdges.get(node)){
					// make product over nodes in edge
					double prod1 = 1.0;
					for (String nb : edge.sourceVertices)
						prod1 = prod1 * u1.get(nb);
					prod1 = Math.pow(prod1, 0.5/(double)edge.sourceVertices.size());
					double prod2 = 1.0;
					for (String nb : edge.targetVertices)
						prod2 = prod2 * u2.get(nb);
					prod2 = Math.pow(prod2, 0.5/(double)edge.targetVertices.size());
					// add to sum
					val += edge.weight*prod1*prod2/(double)edge.sourceVertices.size();
				}
				v.put(node, Math.pow(val, 1.0/pnorm));
			}	
		} else if (dim == 2){
			for (String node : this.vertices){
				double val = 0.0;
				for (Edge edge : this.nodeTargetEdges.get(node)){
					// make product over nodes in edge
					double prod1 = 1.0;
					for (String nb : edge.sourceVertices)
						prod1 = prod1 * u1.get(nb);
					prod1 = Math.pow(prod1, 0.5/(double)edge.sourceVertices.size());
					double prod2 = 1.0;
					for (String nb : edge.targetVertices)
						prod2 = prod2 * u2.get(nb);
					prod2 = Math.pow(prod2, 0.5/(double)edge.targetVertices.size());
					// add to sum
					val += edge.weight*prod1*prod2/(double)edge.targetVertices.size();
				}
				v.put(node, Math.pow(val, 1.0/pnorm));
			}	
		}
		return v;
	}
	
	/**
	 * Return generalized Rayleigh-Ritz score of a vector (undirected case)
	 * @param v
	 * @param pnorm
	 * @return
	 */
	public double vectorScore(HashMap<String, Double> v, double pnorm){
		double score = 0.0;
		for (Edge edge : this.edges){
			double prod = 1.0;
			for (String node : edge.vertices){
				prod = prod * v.get(node);
			}
			prod = Math.pow(prod, 1.0/(double)edge.vertices.size());
			score += edge.weight*prod;
		}
		score = score/norm(v, pnorm);
		return score;
	}
	
	/**
	 * Return generalized Rayleigh-Ritz score of a vector pair (directed case)
	 * @param v
	 * @param pnorm
	 * @return
	 */
	public double vectorScore(HashMap<String, Double> v, HashMap<String, Double> w, double pnorm, double qnorm){
		double score = 0.0;
		for (Edge edge : this.edges){
			double prod1 = 1.0;
			for (String node : edge.sourceVertices){
				prod1 = prod1 * v.get(node);
			}
			prod1 = Math.pow(prod1, 0.5/(double)edge.sourceVertices.size());
			double prod2 = 1.0;
			for (String node : edge.targetVertices){
				prod2 = prod2 * w.get(node);
			}
			prod2 = Math.pow(prod2, 0.5/(double)edge.targetVertices.size());
			score += edge.weight*prod1*prod2;
		}
		score = score/Math.pow(norm(v, pnorm)*norm(w,qnorm),0.5);
		return score;
	}
	
	/**
	 * Norm of a vector
	 * @param vec
	 * @param pnorm
	 * @return
	 */
	public static double norm (HashMap<String, Double> vec, double pnorm){
		double nrm = 0.0;
		for (String node : vec.keySet())
			nrm += Math.pow(vec.get(node), pnorm);
		nrm = Math.pow(nrm, 1.0/pnorm);
		return nrm;
	}
	
	/**
	 * Normalize vector
	 * @param vec
	 */
	public static void normalize(HashMap<String,Double> vec, double norm){
		for (String node : vec.keySet())
			vec.put(node, vec.get(node)/norm);
	}
	
	/**
	 * Print edges
	 */
	public void printEdges(){
		if (this.directed)
			for (Edge edge : this.edges)
				System.out.println(edge.toStringDirected(this.weighted));
		else
			for (Edge edge : this.edges)
				System.out.println(edge.toString(this.weighted));
	}
	
	/**
	 * Print edges
	 */
	public void printEdges(String source){
		for (Edge edge : this.nodeEdges.get(source)){
			String line = "";
			for (String node : edge.vertices)
				line = line + "\t" + node;
			System.out.println(line);
		}
	}
	

}
