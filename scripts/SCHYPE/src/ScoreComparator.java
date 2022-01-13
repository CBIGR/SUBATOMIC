/*
 * SCHype - Java package for spectral clustering in hypergraphs
 * 
 * Copyright (C) 2012 Tom Michoel (The Roslin Institute, University of Edinburgh)
 * 
 */

import java.util.Comparator;

public class ScoreComparator implements Comparator {

	public int compare(Object clust1, Object clust2){
		double score1 = ((Cluster)clust1).score;
		double score2 = ((Cluster)clust2).score;
		// standard sorting sorts in ascending order, we want descending order
		if (score1<score2)
			return 1;
		else if (score1>score2)
			return -1;
		else
			return 0;
	}
	
}
