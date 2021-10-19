/*
 * SCHype - Java package for spectral clustering in hypergraphs
 * 
 * Copyright (C) 2012 Tom Michoel (The Roslin Institute, University of Edinburgh)
 * 
 */
 
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;

public class SCHype {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		// mandatory parameters
		String hgfile = null;
		String output = null;
		
		// optional parameters
		boolean dir = false;
		boolean weighted = false;
		double p = 1.0;
		double q = p;
		boolean clusteredges = true;
		double tolerance = 1E-5;
		int maxstep = 1000;
		int minclustsize = 1;
		double minclustscore = 0.0;
		
		// create options
		Options opts = new Options();
		opts.addOption("hgfile", true, "Hypergraph input file.");
		opts.addOption("output", true, "Output file name.");
		opts.addOption("dir", true, "Flag for directed hypergraphs, default false.");
		opts.addOption("weighted", true, "Flag for weighted hypergraphs, default false.");
		opts.addOption("p", true, "Edge-to-node scaling parameter, default 1.0.");
		opts.addOption("q", true, "Edge-to-node scaling parameter for target set in directed hypergpraphs," +
				", default equal to p.");
		opts.addOption("clusteredges", true, "Flag for clustering edges or nodes, default true.");
		opts.addOption("tolerance", true, "Tolerance for convergence of eigenvector calculation, default 1E-5.");
		opts.addOption("maxstep", true, "Maximum number of iterations in eigenvector calculation, default 1000, " +
				" message printed if exceeded.");
		opts.addOption("minclustsize", true, "In output, only keep clusters with at least this number of nodes, " +
				"default 1.");
		opts.addOption("minclustscore", true, "In output, only keep clusters with at least this score value, " +
				"default 0.0.");
		
		// parse command line
		CommandLineParser parser = new PosixParser();
		try {
			
			CommandLine cmd = parser.parse(opts, args);
			
			if (cmd.hasOption("hgfile"))
				hgfile = cmd.getOptionValue("hgfile");
			
			if (cmd.hasOption("output"))
				output = cmd.getOptionValue("output");
			
			if (cmd.hasOption("dir"))
				dir = Boolean.parseBoolean(cmd.getOptionValue("dir"));
			
			if (cmd.hasOption("weighted"))
				weighted = Boolean.parseBoolean(cmd.getOptionValue("weighted"));
			
			if (cmd.hasOption("p")){
				p = Double.parseDouble(cmd.getOptionValue("p"));
				q = p;
			}
			
			if (cmd.hasOption("q"))
				q = Double.parseDouble(cmd.getOptionValue("q"));
			
			if (cmd.hasOption("clusteredges"))
				clusteredges = Boolean.parseBoolean(cmd.getOptionValue("clusteredges"));
			
			if (cmd.hasOption("tolerance"))
				tolerance = Double.parseDouble(cmd.getOptionValue("tolerance"));
			
			if (cmd.hasOption("maxstep"))
				maxstep = Integer.parseInt(cmd.getOptionValue("maxstep"));
			
			if (cmd.hasOption("minclustsize"))
				minclustsize = Integer.parseInt(cmd.getOptionValue("minclustsize"));
			
			if (cmd.hasOption("minclustscore"))
				minclustscore = Double.parseDouble(cmd.getOptionValue("minclustscore"));
		}
		catch (ParseException exp) {
			System.out.println("Error while parsing command line:");
			System.out.println();
			exp.printStackTrace();
			System.exit(1);
		}
		
		// print banner
		printBanner();
		
		// we need input and output file
		if (hgfile==null)
			Die("Error: Hypergraph input file name must be provided.");
		if (output==null)
			Die("Error: Output file name must be provided.");
		
		// print parameter settings
		System.out.println("Parameters");
		System.out.println("----------");
		System.out.println("Hypergraph file name:                                   " + hgfile);
		System.out.println("Output file name (clusters - nodes):                    " + output + ".nodes.txt");
		System.out.println("Output file name (clusters - edges):                    " + output + ".edges.txt");
		System.out.println("Directed hypergraph:                                    " + dir);
		System.out.println("Weighted hypergraph:                                    " + weighted);
		if (!dir)
		System.out.println("Edge-to-node scaling parameter:                         " + p);
		else{
		System.out.println("Edge-to-node scaling parameter (source nodes):          " + p);
		System.out.println("Edge-to-node scaling parameter (target nodes):          " + q);
		}
		System.out.println("Cluster edges:                                          " + clusteredges);
		System.out.println("Tolerance (eigenvector calculation):                    " + tolerance);
		System.out.println("Maximum number of iterations (eigenvector calculation): " + maxstep);
		System.out.println("Minimum cluster size:                                   " + minclustsize);
		System.out.println("Minimum cluster score:                                  " + minclustscore);
		System.out.println("----------");
		System.out.println("");
		
		// create hypergraph object
		HyperGraph hg = new HyperGraph(hgfile, dir, weighted);
		// create clustering object
		PFClustering clust;
		if (dir)
			clust = new PFClustering(hg, p, q);
		else
			clust = new PFClustering(hg, p);
		// do clustering
		clust.recursivePfClustering();
		// some post-processing and printing the output
		clust.postProcessing();
		clust.writeClusters(output + ".nodes.txt");
		clust.writeClusterEdges(output + ".edges.txt");
	}
	
	public static void Die (String msg) {
		System.out.println(msg);
		System.exit(1);
	}
	
	public static void printBanner () {
		System.out.println("");
		System.out.println("SCHype - Spectral Clustering in Hypergraphs");
		System.out.println("-------------------------------------------");
		System.out.println("Version 1.0");
		System.out.println("Copyright (c) 2011-2012 Tom Michoel");
		System.out.println("");
	}

}
